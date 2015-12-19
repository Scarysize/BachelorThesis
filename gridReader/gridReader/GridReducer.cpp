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

struct vertexCompare {
    bool operator()(EdgeCollapse col1,  EdgeCollapse col2) {
        if ((col1.getA() == col2.getA() &&
             col1.getB() == col2.getB()) ||
            (col1.getA() == col2.getB() &&
             col1.getB() == col2.getA())) {
                return false;
            }
        return true;
    }
};

void GridReducer::run(double (*calculateCost)(Vertex, Vertex, std::vector<Cell>, std::vector<Vertex>)) {
    buildQueue(calculateCost);
}

void GridReducer::buildQueue(double (*calculateCost)(Vertex a, Vertex b, std::vector<Cell> cells, std::vector<Vertex> vertices)) {
    // run classifier
    Classifier::classifyVertices(&this->vertices, &this->cells);
    std::vector<EdgeCollapse> collapses;
    std::set<std::vector<int>> edges;
    for (auto cell : this->cells) {
        if (!cell.deleted) {
            for (auto edge : cell.edges) {
                if (edges.find(edge) == edges.end()) {
                    edges.insert(edge);
                    int vertexIdA = this->vertices[edge[0]].getId();
                    int vertexIdB = this->vertices[edge[1]].getId();
                    double cost = (*calculateCost)(this->vertices[vertexIdA], this->vertices[vertexIdB], cells, vertices);
                    collapses.push_back(*new EdgeCollapse(vertexIdA, vertexIdB, cost));
                } else {
                    std::cout << "in q: " << edge[0] << "---" << edge[1] << std::endl;
                }
                // std::cout << edge[0] << "--" << edge[1] << ": " << cost << std::endl;
            }
        }
    }
    // std::make_heap(collapses.begin(), collapses.end(), EdgeCollapse::CompareCost());
    // this->prioq = collapses;
};