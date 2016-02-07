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


vector<Cell*> Connectivity::getIcells(Vertex *a, Vertex *b) {
    vector<Cell*> cellsUsingAB;
    for (auto cell : a->incidents) {
        if (find(b->incidents.begin(), b->incidents.end(), cell) != b->incidents.end()) {
            cellsUsingAB.push_back(cell);
        }
    }
    return cellsUsingAB;
}


vector<Cell*> Connectivity::getNcells(Vertex *a, Vertex *b) {

    if (b->incidents.empty()) {
        vector<Cell*> aOnly;
        copy(a->incidents.begin(), a->incidents.end(), back_inserter(aOnly));
        cout << "INFO: b.incidents empty" << endl;
        return aOnly;
    } else if (a->incidents.empty()) {
        vector<Cell*> bOnly;
        copy(b->incidents.begin(), b->incidents.end(), back_inserter(bOnly));
        cout << "INFO: a.incidents empty" << endl;
        return bOnly;
    }
//    for (auto cell : a->incidents) {
//        cout << cell->id << endl;
//    }
//    cout << " -------- " << endl;
//    for (auto cell : b->incidents) {
//        cout << cell->id << endl;
//    }
    vector<Cell*> symDiff;
    for (auto cell : a->incidents) {
        // cell is not in b, push_back
        if (find(b->incidents.begin(), b->incidents.end(), cell) == b->incidents.end()) {
            symDiff.push_back(cell);
        }
    }
    for (auto cell : b->incidents) {
        // cell is not in a, push_back
        if (find(a->incidents.begin(), a->incidents.end(), cell) == a->incidents.end()) {
            symDiff.push_back(cell);
        }
    }
    return symDiff;
}
