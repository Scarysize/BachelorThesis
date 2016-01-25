//
//  DynamicTest.cpp
//  gridReader
//
//  Created by Franz Neubert on 18/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "DynamicTest.hpp"
#include "Classifier.hpp"


bool DynamicTest::testVolume(std::set<vtkIdType> ncells, vtkUnstructuredGrid *gridPre, vtkUnstructuredGrid *gridAfter) {
    for (auto ncell : ncells) {
        vtkSmartPointer<vtkPoints> ncellPointsPre = gridPre->GetCell(ncell)->GetPoints();
        vtkSmartPointer<vtkPoints> ncellPointsAfter = gridAfter->GetCell(ncell)->GetPoints();
        double ncellCoordsPre[3*4];
        double ncellCoordsAfter[3*4];
        for (int i = 0; i < 4; i++) {
            ncellPointsPre->GetPoint(i, &ncellCoordsPre[i*3]);
            ncellPointsAfter->GetPoint(i, &ncellCoordsAfter[i*3]);
        }
        double volumePre = Calculator::calcTetraVolume(&ncellCoordsPre[0], &ncellCoordsPre[3], &ncellCoordsPre[6], &ncellCoordsPre[9]);
        double volumeAfter = Calculator::calcTetraVolume(&ncellCoordsAfter[0], &ncellCoordsAfter[3], &ncellCoordsAfter[6], &ncellCoordsAfter[9]);
        // check for sign changes in volume --> signals self-intersection / folding
        if ((volumePre < 0 && volumeAfter > 0) ||
            (volumePre > 0 && volumeAfter < 0)) {
            return false;
        }
    }
    return true;
}

bool DynamicTest::testSolidAngle(EdgeCollapse *collapse, std::vector<Vertex *> vertices, std::vector<Cell *> cells) {
    if (vertices[collapse->getA()]->isInterior() && vertices[collapse->getB()]->isInterior()) {
        return true;
    }
    double angleA = Classifier::calcSolidAngleSum(collapse->getA(), &vertices, &cells);
    double angleB = Classifier::calcSolidAngleSum(collapse->getB(), &vertices, &cells);
    double collapsePoint[3];
    double coordsA[3];
    double coordsB[3];
    vertices[collapse->getA()]->getCoords(coordsA);
    vertices[collapse->getB()]->getCoords(coordsB);
    Calculator::calcMidPoint(coordsA, coordsB, collapsePoint);
    vertices[collapse->getA()]->setCoords(collapsePoint);
    double angleCollapse = Classifier::calcSolidAngleSum(collapse->getA(), &vertices, &cells);
    vertices[collapse->getA()]->setCoords(coordsA);
    // check for deviation in the solid angle to prevent dents in the boundary
    return (fabs(angleCollapse - angleA) <= 0);
}