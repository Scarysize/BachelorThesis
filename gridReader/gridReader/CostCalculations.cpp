//
//  CostCalculations.cpp
//  gridReader
//
//  Created by Franz Neubert on 19/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "CostCalculations.hpp"
#include "Calculator.h"

double CostCalculations::calcEdgeLengthCost(Vertex *a, Vertex *b, vector<Cell*> *cells, vector<Vertex*> *vertices) {
    a->setModified(false);
    b->setModified(false);
    return Calculator::calcEdgeLength(a->getCoords(), b->getCoords());
}