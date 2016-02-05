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

void GridReducer::run(double (*calculateCost)(Vertex*, Vertex*, std::vector<Cell*>, std::vector<Vertex*>)) {
    buildQueue(calculateCost);
    doCollapse(calculateCost);
    // do collapses & recalculations
}

void GridReducer::buildQueue(double (*calculateCost)(Vertex *a, Vertex *b, std::vector<Cell*> cells, std::vector<Vertex*> vertices)) {
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
                Edge *checkEdge = Edge::isEdge(a, b, &traversed);
                
                // check if edge was already traversed
                if (checkEdge == nullptr) {
                    // mark edge as traversed
                    traversed.push_back(new Edge(a, b));

                    // check if edge is either inner or boundary edge (everything else can´t be collapsed)
                    if (Classifier::isInnerEdge(a, b) ||
                        Classifier::isBoundaryEdge(a, b)) {
                        double cost = (*calculateCost)(a, b, this->grid->cells, this->grid->vertices);
                        collapses.push_back(new EdgeCollapse(a, b, cost));
                    }
                }
            }
        }
    }
    std::make_heap(collapses.begin(), collapses.end(), EdgeCollapse::CompareCost());
    this->prioq = collapses;
}
//
//void GridReducer::recalcQueue(double (*calculateCost)(Vertex *, Vertex *, std::vector<Cell *>, std::vector<Vertex *>), EdgeCollapse *lastCollapse) {
//    std::vector<EdgeCollapse*> recalced;
//    for (auto collapse : this->prioq) {
//        // Collapse uses deleted vertex -> remove
//        if (collapse->getB() == lastCollapse->getB() ||
//            collapse->getA() == lastCollapse->getB()) {
//            continue;
//        }
//        // Collapse uses modified vertex -> recalc
//        else if (collapse->getA() == lastCollapse->getA() ||
//            collapse->getB() == lastCollapse->getA()) {
//            double cost = (*calculateCost)(this->vertices[collapse->getA()], this->vertices[collapse->getB()], this->cells, this->vertices);
//            collapse->setCost(cost);
//            recalced.push_back(collapse);
//        } else {
//            recalced.push_back(collapse);
//        }
//    }
//    std::make_heap(recalced.begin(), recalced.end(), EdgeCollapse::CompareCost());
//    this->prioq = recalced;
//}


void GridReducer::doCollapse(double (*calculateCost)(Vertex *a, Vertex *b, std::vector<Cell*> cells, std::vector<Vertex*> vertices)){
    for (int x = 0; x < 1; x++) {
        for (int i = 0; !this->prioq.empty(); i++) {
            // 1. Get top most collapse from heap
            EdgeCollapse *top = this->prioq.front();
            
            /*
             TESTS/SIMLATION:
                - check if solid angle at collapse point = solid angle vertex a/b
            */
//            if (!DynamicTest::testSolidAngle(top, this->grid->vertices, this->grid->cells)) {
//                this->prioq.erase(this->prioq.begin());
//                std::cout << "INFO: REJECTED" << std::endl;
//                continue;
//            }
        
            set<Cell*> ncells = Connectivity::getNcells(top->getA(), top->getB());
            set<Cell*> icells = Connectivity::getIcells(top->getA(), top->getB());
            cout << "INFO: collapsing" << top->getA()->getId() << " -- " << top->getB()->getId() << endl;
            
//            // 2. Calculate coords of collapse point
//            double collapsePoint[3];
//            double coordsA[3];
//            double coordsB[3];
//            this->vertices[top->getA()]->getCoords(coordsA);
//            this->vertices[top->getB()]->getCoords(coordsB);
//            Calculator::calcMidPoint(coordsA, coordsB, collapsePoint);
//            // 3. Update coords of vertex A to collapse coords
//            this->vertices[top->getA()]->setCoords(collapsePoint);
//            
//            // 4. Replace vertex B with vertex A in NCELLS
//            for (auto ncell : ncells) {
//                Cell *cell = this->cells[ncell];
//                for (int j = 0; j < cell->points.size(); j++) {
//                    if (cell->points[j] == top->getB()) {
//                        cell->points[j] = top->getA();
//                    }
//                }
//            }
//            
//            // 5. Delete all ICELLS
//            for (auto icell : icells) {
//                this->cells[icell]->deleteCell();
//            }
//            
//            // 6. Remove the collapse from the queue
//            this->prioq.erase(this->prioq.begin());
//            // 7. Recalculate the queue
//            recalcQueue(calculateCost, top);
        }
        buildQueue(calculateCost);
    }
}