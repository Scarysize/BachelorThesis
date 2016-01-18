//
//  Connectivity.cpp
//  gridReader
//
//  Created by Franz Neubert on 18/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "Connectivity.hpp"
#include "Cell.hpp"
#include <vector>

std::vector<int> Connectivity::cellsUsingVertex(int vertex, std::vector<Cell*> *cells) {
    std::vector<int> vertexCells;
    for (auto cell : *cells) {
        if (!cell->deleted) {
            std::set<int> pointSet = { cell->points[0], cell->points[1], cell->points[2], cell->points[3] };
            if (pointSet.find(vertex) != pointSet.end()) {
                vertexCells.push_back(cell->id);
            }
        }
    }
    return vertexCells;
}

std::set<int> Connectivity::getIcells(int vertexA, int vertexB, std::vector<Cell*> *cells) {
    std::vector<int> cellsUsingA = cellsUsingVertex(vertexA, cells);
    std::vector<int> cellsUsingB = cellsUsingVertex(vertexB, cells);
    std::set<int> cellsUsingAB;
    std::set_intersection(cellsUsingA.begin(), cellsUsingA.end(), cellsUsingB.begin(), cellsUsingB.end(), std::inserter(cellsUsingAB, cellsUsingAB.begin()));
    return cellsUsingAB;
}

std::set<int> Connectivity::getNcells(int vertexA, int vertexB, std::vector<Cell*> *cells) {
    std::vector<int> cellsUsingA = cellsUsingVertex(vertexA, cells);
    std::vector<int> cellsUsingB = cellsUsingVertex(vertexB, cells);
    std::set<int> diffAB;
    std::set<int> diffBA;
    std::set<int> unionDiffs;
    std::set_difference(cellsUsingA.begin(), cellsUsingA.end(), cellsUsingB.begin(), cellsUsingA.begin(), std::inserter(diffAB, diffAB.begin()));
    std::set_difference(cellsUsingB.begin(), cellsUsingB.end(), cellsUsingA.begin(), cellsUsingA.begin(), std::inserter(diffBA, diffBA.begin()));
    std::set_union(diffAB.begin(), diffAB.end(), diffBA.begin(), diffBA.end(), std::inserter(unionDiffs, unionDiffs.begin()));
    return unionDiffs;
}
