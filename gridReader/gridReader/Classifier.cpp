//
//  Classifier.cpp
//  gridReader
//
//  Created by Franz Neubert on 19/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "Cell.hpp"
#include "Vertex.hpp"
#include "Classifier.hpp"
#include "Connectivity.hpp"
#include "Calculator.h"
#include "Tetragrid.hpp"
#include <vector>
#include <set>
bool Classifier::isBoundaryVertex(double angle) {
    if (angle < 4 * M_PI) {
        return true;
    }
    return false;
}

bool Classifier::isInnerVertex(double angle) {
    if (angle >= 4 * M_PI) {
        return true;
    }
    return false;
}

bool Classifier::isCornerVertex(double angle) {
    if (angle <= M_PI/2 ||  4 * M_PI - angle <= M_PI/2) {
        return true;
    }
    return false;
}

bool Classifier::isCurveCornerVertex(double angle) {
    if ((M_PI/2 < angle && angle <= (3 * M_PI) / 2) ||
        ((M_PI/2 < 4 * M_PI - angle) && (4 * M_PI - angle <= (3 * M_PI) / 2))) {
        return true;
    }
    return false;
}

bool Classifier::isBoundaryEdge(Vertex *A, Vertex *B) {
    if (A->isBoundary() && B->isBoundary()) {
        return true;
    }
    return false;
}

bool Classifier::isInnerEdge(Vertex *A, Vertex *B) {
    if (A->isInterior() && B->isInterior()) {
        return true;
    }
    return false;
}

void Classifier::classifyVertices(Tetragrid *grid) {
    for (auto vertex : grid->vertices) {
        double solidAngle = Classifier::calcSolidAngleSum(vertex, grid);
        if (isInnerVertex(solidAngle)) {
            vertex->setToInterior();
        } else if (isCurveCornerVertex(solidAngle) || isCornerVertex(solidAngle)) {
            vertex->setToCorner();
        } else if (isBoundaryVertex(solidAngle)) {
            vertex->setToBoundary();
        }
    }
}

double Classifier::calcSolidAngleSum(Vertex *seed, Tetragrid *grid) {
    double solidAngleSum = 0;
    for (auto cell : seed->incidents) {
        double pointCoords[3*3];
        double seedCoords[3];
        int i = 0;
        for (auto vertex : cell->vertices) {
            if (vertex != seed) {
                vertex->getCoords(&pointCoords[i]);
                i += 3;
            } else {
                vertex->getCoords(seedCoords);
            }
        }
        double a[3];
        double b[3];
        double c[3];
        Calculator::calcVectorBetweenPoints(seedCoords, &pointCoords[0], a);
        Calculator::calcVectorBetweenPoints(seedCoords, &pointCoords[3], b);
        Calculator::calcVectorBetweenPoints(seedCoords, &pointCoords[6], c);
        solidAngleSum += Calculator::calcSolidAngle(a, b, c, seedCoords);
    }
    return solidAngleSum;
}