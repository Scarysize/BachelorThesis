//
//  GridReducer.cpp
//  gridReader
//
//  Created by Franz Neubert on 19/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//
#include "Calculator.h"
#include "Classifier.hpp"
#include "Connectivity.hpp"
#include "GridReducer.hpp"
#include "EdgeCollapse.hpp"
#include "DynamicTest.hpp"
#include "Tetragrid.hpp"
#include "Edge.hpp"

using namespace std;

void GridReducer::run(double (*calculateCost)(Vertex *a, Vertex *b, Tetragrid *grid)) {
    buildQueue(calculateCost);
    doCollapse(calculateCost);
    // do collapses & recalculations
}

void GridReducer::buildQueue(double (*calculateCost)(Vertex *a, Vertex *b, Tetragrid *grid)) {
    // run classifier
    Classifier::classifyVertices(this->grid);
    cout << "DONE: vertex classification" << endl;
    
    vector<EdgeCollapse*> collapses;
    vector<Edge*> traversed;
    int counter = 0;
    for (auto cell : this->grid->cells) {
        if (!cell->deleted) {
            counter++;
            if (counter%500 == 0) {
                std::cout << "INFO: cells traversed: " << counter << std::endl;
            }
            for (auto edge : cell->edges) {
                Vertex *a = edge->getA();
                Vertex *b = edge->getB();
                if (a->isDeleted() || b->isDeleted()) {
                    continue;
                }
                Edge *checkEdge = Edge::isEdge(a, b, &traversed);
                
                // check if edge was already traversed
                if (checkEdge == nullptr) {
                    // mark edge as traversed
                    traversed.push_back(new Edge(a, b));
                    
                    // check if edge is either inner or boundary edge (everything else can´t be collapsed)
                    if (Classifier::isInnerEdge(a, b) ||
                        Classifier::isBoundaryEdge(a, b)) {
                        double cost = (*calculateCost)(a, b, this->grid);
                        collapses.push_back(new EdgeCollapse(a, b, cost));
                    }
                }
            }
        }
    }
    std::make_heap(collapses.begin(), collapses.end(), EdgeCollapse::CompareCost());
    this->prioq.clear();
    this->prioq = collapses;
    double avrgCost = 0;
    double min = 0;
    double max = 0;
    for (auto col : this->prioq) {
        avrgCost += col->getCost();
        if (min == 0 || col->getCost() < min) {
            min = col->getCost();
        }
        if (max == 0 || col->getCost() > max) {
            max = col->getCost();
        }
    }
    this->setLimit(max);
    cout << "INFO: average collapse cost: " << avrgCost / prioq.size() << endl;
    cout << "INFO: lowest collapse cost:  " << min << endl;
    cout << "INFO: highest collapse cost: " << max << endl;
}

void GridReducer::recalcQueue(double (*calculateCost)(Vertex *a, Vertex *b, Tetragrid *grid), EdgeCollapse *lastCol) {
    vector<EdgeCollapse*> recalced;
    for (auto collapse : this->prioq) {
        if (collapse->getA()->isDeleted()) {
//            double cost = (*calculateCost)(collapse->getB(), lastCol->getA(), this->grid);
//            recalced.push_back(new EdgeCollapse(collapse->getB(), lastCol->getA(), cost));
            continue;
        } else if (collapse->getB()->isDeleted()) {
//            double cost = (*calculateCost)(collapse->getA(), lastCol->getA(), this->grid);
//            recalced.push_back(new EdgeCollapse(collapse->getA(), lastCol->getA(), cost));
            continue;
        }
        else if(collapse->getA()->isModified() || collapse->getB()->isModified()) {
            double cost = (*calculateCost)(collapse->getA(), collapse->getB(), this->grid);
            recalced.push_back(new EdgeCollapse(collapse->getA(), collapse->getB(), cost));
        } else {
            recalced.push_back(collapse);
        }
    }
    // reset modified flag after recalculation
    for (auto vertex : this->grid->vertices) {
        if (!vertex->isDeleted() && vertex->isModified()) {
            vertex->setModified(false);
        }
    }
    make_heap(recalced.begin(), recalced.end(), EdgeCollapse::CompareCost());
    this->prioq.clear();
    this->prioq = recalced;
}


void GridReducer::doCollapse(double (*calculateCost)(Vertex *a, Vertex *b, Tetragrid *grid)){
    for (int x = 0; x < 3; x++) {
        for (int i = 0; !this->prioq.empty() /*&& i < 7*/; i++) {
            // 1. Get top most collapse from heap
            EdgeCollapse *top = this->prioq.front();
            
            if (top->getCost() > this->getLimit()) {
                cout << "INFO: limit reached" << endl;
                break;
            }
            
            /*
             TESTS/SIMLATION:
             - check if solid angle at collapse point = solid angle vertex a/b
             */
            if (!DynamicTest::testSolidAngle(top, this->grid)) {
                this->prioq.erase(this->prioq.begin());
                std::cout << "INFO: REJECTED: " << top->getA()->getId() << " -- " << top->getB()->getId() << std::endl;
                continue;
            }
            
            vector<Cell*> ncells = Connectivity::getNcells(top->getA(), top->getB());
            vector<Cell*> icells = Connectivity::getIcells(top->getA(), top->getB());
            cout << "INFO: collapsing: " << top->getA()->getId() << " -- " << top->getB()->getId() << endl;
            if (icells.size() == 0) {
                cout << "WARN: number of icells 0: " << endl;
            }
            
            
            // 2. Calculate coords of collapse point
            double collapsePoint[3];
            double coordsA[3];
            double coordsB[3];
            top->getA()->getCoords(coordsA);
            top->getB()->getCoords(coordsB);
            Calculator::calcMidPoint(coordsA, coordsB, collapsePoint);
            // 3. Update coords of vertex A to collapse coords
            top->getA()->setCoords(collapsePoint);
            //    Mark B as deleted, A as modified, updates incidents on A
            top->getB()->deleteVertex();
            top->getA()->setModified(true);
            top->getA()->incidents = ncells;
            
            // 4. Replace vertex B with vertex A in NCELLS
            for (auto ncell : ncells) {
                if (!ncell->deleted) {
                    vector<Vertex*> updated;
                    for (auto vertex : ncell->vertices) {
                        if (vertex == top->getB()) {
                            updated.push_back(top->getA());
                        } else {
                            updated.push_back(vertex);
                        }
                    }
                    ncell->vertices.clear();
                    ncell->vertices = updated;
                }
            }
            
            // 5. Delete all ICELLS
            for (auto icell : icells) {
                if (!icell->deleted) {
                    icell->deleteCell();
                }
            }
            ncells.clear();
            icells.clear();
            
            // 6. Remove the collapse from the queue
            this->prioq.erase(this->prioq.begin());
            // 7. Recalculate the queue
            recalcQueue(calculateCost, top);
        }
        buildQueue(calculateCost);
    }
}