//
//  CostCalculations.cpp
//  gridReader
//
//  Created by Franz Neubert on 19/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "CostCalculations.hpp"
#include "Connectivity.hpp"
#include "Calculator.h"
#include "Tetragrid.hpp"
#include "Edge.hpp"
#include "Vertex.hpp"

using namespace std;

double CostCalculations::calcCombinedCost(Vertex *a, Vertex *b, Tetragrid *grid) {
    return   1 * calcEdgeLengthCost(a, b, grid) +
           200 * calcEdgeEquityCost(a, b, grid) +
             1 * calcScalarCost(a, b, grid);
}

double CostCalculations::calcEdgeLengthCost(Vertex *a, Vertex *b, Tetragrid *grid) {
    return Calculator::calcEdgeLength(a->getCoords(), b->getCoords());
}

double CostCalculations::calcEdgeEquityCost(Vertex *a, Vertex *b, Tetragrid *grid) {
    vector<Cell*> ncells = Connectivity::getNcells(a, b);
    double sum = 0;
    for (auto cell : ncells) {
        double preColSum = 0;
        for (auto edge : cell->edges) {
            preColSum += pow((edge->getLength() - cell->avrgEdgeLength()), 2);
        }
        // SIMULATE COLLAPSE
        double coordsA[3];
        double coordsB[3];
        double collapsePoint[3];
        a->getCoords(coordsA);
        b->getCoords(coordsB);
        Calculator::calcMidPoint(coordsA, coordsB, collapsePoint);
        a->setCoords(collapsePoint);
        double postColAvrgLen = cell->avrgEdgeLength();
        double postColSum = 0;
        for (auto edge : cell->edges) {
            postColSum += pow((edge->calcEdgeLength() - postColAvrgLen), 2);
        }
        sum += preColSum - postColSum;
        // REVERSE SIMULATION
        a->setCoords(coordsA);
    }
    return sum;
}

double CostCalculations::calcScalarCost(Vertex *a, Vertex *b, Tetragrid *grid) {
    return 0;
}