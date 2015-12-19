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
#include <vector>
#include <set>

void Classifier::classifyVertices(std::vector<Vertex> *vertices, std::vector<Cell> *cells) {
    for (auto vertex : *vertices) {
        Classifier::calcSolidAngleSum(vertex.getId(), vertices, cells);
    }
}

double Classifier::calcSolidAngleSum(int seed, std::vector<Vertex> *vertices, std::vector<Cell> *cells) {
    std::vector<int> cellsUsingVertex = Connectivity::cellsUsingVertex(seed, *cells);
    double solidAngleSum = 0;
    for (auto cell : cellsUsingVertex) {
        double pointCoords[3*3];
        double seedCoords[3];
        int i = 0;
        for (auto point : cells->at(cell).points) {
            if (point != seed) {
                vertices->at(point).getCoords(&pointCoords[i]);
                i += 3;
            } else {
                vertices->at(point).getCoords(seedCoords);
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