//
//  Connectivity.cpp
//  gridReader
//
//  Created by Franz Neubert on 18/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "Connectivity.hpp"
#include <vector>

#include "Vertex.hpp"
#include "Cell.hpp"

using namespace std;


set<Cell*> Connectivity::getIcells(Vertex *a, Vertex *b) {
    vector<Cell*> cellsUsingA = a->incidents;
    vector<Cell*> cellsUsingB = b->incidents;
    set<Cell*> cellsUsingAB;
    set_intersection(cellsUsingA.begin(), cellsUsingA.end(), cellsUsingB.begin(), cellsUsingB.end(), inserter(cellsUsingAB, cellsUsingAB.begin()));
    return cellsUsingAB;
}


set<Cell*> Connectivity::getNcells(Vertex *a, Vertex *b) {
    vector<Cell*> cellsUsingA = a->incidents;
    vector<Cell*> cellsUsingB = b->incidents;
    
    if (cellsUsingB.empty()) {
        set<Cell*> cellSet;
        copy(cellsUsingA.begin(), cellsUsingA.end(), inserter(cellSet, cellSet.begin()));
        return cellSet;
    }
    
    set<Cell*> diffAB;
    set<Cell*> diffBA;
    set<Cell*> unionDiffs;
    set_difference(cellsUsingA.begin(), cellsUsingA.end(), cellsUsingB.begin(), cellsUsingA.begin(), inserter(diffAB, diffAB.begin()));
    set_difference(cellsUsingB.begin(), cellsUsingB.end(), cellsUsingA.begin(), cellsUsingA.begin(), inserter(diffBA, diffBA.begin()));
    set_union(diffAB.begin(), diffAB.end(), diffBA.begin(), diffBA.end(), inserter(unionDiffs, unionDiffs.begin()));
    return unionDiffs;
}
