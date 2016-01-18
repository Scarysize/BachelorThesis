//
//  GridReducer.cpp
//  gridReader
//
//  Created by Franz Neubert on 19/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "Classifier.hpp"
#include "GridReducer.hpp"
#include "EdgeCollapse.hpp"


void GridReducer::run(double (*calculateCost)(Vertex, Vertex, std::vector<Cell*>, std::vector<Vertex*>)) {
    buildQueue(calculateCost);
}

void GridReducer::buildQueue(double (*calculateCost)(Vertex a, Vertex b, std::vector<Cell*> cells, std::vector<Vertex*> vertices)) {
    // run classifier
    Classifier::classifyVertices(&this->vertices, &this->cells);
    std::cout << "DONE: vertex classification" << std::endl;
    
    std::vector<EdgeCollapse> collapses;
    std::set<std::vector<int>> edges;
    int counter = 0;
    for (auto cell : this->cells) {
        counter++;
        if (counter%250 == 0) {
            std::cout << "INFO: cells traversed: " << counter << std::endl;
        }
        if (!cell->deleted) {
            for (auto edge : cell->edges) {
                // check if edge was already traversed
                if (edges.find(edge) == edges.end()) {
                    // mark edge as traversed
                    edges.insert(edge);
                    
                    int vertexIdA = this->vertices[edge[0]]->getId();
                    int vertexIdB = this->vertices[edge[1]]->getId();
                    double cost = (*calculateCost)(*this->vertices[vertexIdA], *this->vertices[vertexIdB], cells, vertices);
                    collapses.push_back(*new EdgeCollapse(vertexIdA, vertexIdB, cost));
                }
                
            }
        }
    }
     std::make_heap(collapses.begin(), collapses.end(), EdgeCollapse::CompareCost());
     this->prioq = collapses;
};