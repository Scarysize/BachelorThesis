//
//  CostCalculator.cpp
//  gridReader
//
//  Created by Franz Neubert on 05/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "CostCalculator.hpp"
#include "Calculator.h"
#include "EdgeCollapse.hpp"

#include <math.h>
#include <list>
#include <set>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkTetra.h>


/*!
    Calculates the weighted change of the scalar value caused by collapsing an edge.
    \param pointA The vertex id of an edge
    \param pointB The vertex id of an edge
    \param weight The weight to be applied to the change
    \param tetraGrid The pointer to a tetrahedron grid containing the edge described by the vertices pointA and pointB
    \return The weighted change of the scalar value if pointA & pointB were to be collapsed
*/
double CostCalculator::calcScalarCost(vtkIdType pointA, vtkIdType pointB, double weight, vtkUnstructuredGrid *tetraGrid) {
    return weight * fabs(getPointData_AlphaWater(pointA, tetraGrid) - getPointData_AlphaWater(pointB, tetraGrid));
}

/*!
    Calculates the weighted change in volume of all tetrahedrons incident on an edge collapse.
    E[vol] = weight * ((Σ[T ε ncells](Vol(T) + Vol(T_new))) - Σ[T ε icells](Vol(T)))
    \param pointA The vertex id of an edge
    \param pointB The vertex id of an edge
    \param weight The weight to be applied to the change
    \param ncells The set of cell ids modified by the edge collapse (NCELLS)
    \param icells The set of cell ids removed by the edge collapse (ICELLS)
    \param tetraGrid The pointer to a tetrahedron grid containing the edge described by the vertices pointA and pointB
    \return the weighted change of the volume if pointA & pointB were to be collapsed
*/
//double CostCalculator::calcVolumeCost(vtkIdType pointA, vtkIdType pointB, double weight, std::set<vtkIdType> ncells, std::set<vtkIdType> icells, vtkUnstructuredGrid *tetraGrid) {
//    double coordsA[3];
//    double coordsB[3];
//    tetraGrid->GetPoints()->GetPoint(pointA, coordsA);
//    tetraGrid->GetPoints()->GetPoint(pointB, coordsB);
//    
//    // Σ[T ε icells](Vol(T))
//    double volume_icells = 0;
//    for (auto icell : icells) {
//        vtkSmartPointer<vtkPoints> icellPoints = tetraGrid->GetCell(icell)->GetPoints();
//        double icellCoords[4*3];
//        for (int i = 0; i < 4; i++) {
//            icellPoints->GetPoint(i, &icellCoords[i*3]);
//        }
//        volume_icells += Calculator::calcTetraVolume(&icellCoords[0], &icellCoords[3], &icellCoords[6], &icellCoords[9]);
//    }
//    
//    // Σ[T ε ncells](Vol(T) + Vol(T_new))
//    double volume_diff_ncells = 0;
//    vtkSmartPointer<vtkUnstructuredGrid> simulateGrid = EdgeCollapse::simulateNcells(pointA, pointB, tetraGrid);
//    for (auto ncell : ncells) {
//        vtkSmartPointer<vtkPoints> ncellPointsPre = tetraGrid->GetCell(ncell)->GetPoints();
//        vtkSmartPointer<vtkPoints> ncellPointsAfter = simulateGrid->GetCell(ncell)->GetPoints();
//        double ncellCoordsPre[3*4];
//        double ncellCoordsAfter[3*4];
//        for (int i = 0; i < 4; i++) {
//            ncellPointsPre->GetPoint(i, &ncellCoordsPre[i*3]);
//            ncellPointsAfter->GetPoint(i, &ncellCoordsAfter[i*3]);
//        }
//        volume_diff_ncells += Calculator::calcTetraVolume(&ncellCoordsPre[0], &ncellCoordsPre[3], &ncellCoordsPre[6], &ncellCoordsPre[9]) - Calculator::calcTetraVolume(&ncellCoordsAfter[0], &ncellCoordsAfter[3], &ncellCoordsAfter[6], &ncellCoordsAfter[9]);
//    }
//
//    return weight * (volume_diff_ncells + volume_icells);
//}

/*!
    Calculates the weighted edge length. This causes shorter edges to be collapsed earlier (due to lower costs).
    \param pointA The vertex id of an edge
    \param pointB The vertex id of an edge
    \param weight The weight to be applied to the change
    \param tetraGrid The pointer to a tetrahedron grid containing the edge described by the vertices pointA and pointB
    \return The length of the edge described by pointA & pointB multiplied with the weight
*/
double CostCalculator::calcEdgeLengthCost(int pointA, int pointB, double weight, double points[][3]) {
    return Calculator::calcEdgeLength(points[pointA], points[pointB]) * weight;
}

// may not be working correctly
//double CostCalculator::calcEdgeEquityCost(vtkIdType pointA, vtkIdType pointB, double weight, vtkUnstructuredGrid *tetraGrid) {
//    vtkSmartPointer<vtkUnstructuredGrid> simulateGrid = EdgeCollapse::simulateNcells(pointA, pointB, tetraGrid);
//    // std::cout <<  simulateGrid->GetNumberOfCells() << std::endl;
//    std::set<vtkIdType> ncells = EdgeCollapse::getNCells(pointA, pointB, tetraGrid);
//    double sumOverNcells = 0;
//    for (auto ncell : ncells) {
//        vtkCell *tetra = tetraGrid->GetCell(ncell);
//        vtkCell *simulateTetra = simulateGrid->GetCell(ncell);
//        double averageEdgeLength = Calculator::calcAverageEdgeLength(ncell, tetraGrid);
//        double simulateAvrgEdgeLength = Calculator::calcAverageEdgeLength(ncell, simulateGrid);
//        double sumPreCollapse = 0;
//        double sumAfterCollapse = 0;
//        for (vtkIdType edge = 0; edge < tetra->GetNumberOfEdges(); edge++) {
//            double coordsA[3];
//            double coordsB[3];
//            tetra->GetEdge((int)edge)->GetPoints()->GetPoint(0, coordsA);
//            tetra->GetEdge((int)edge)->GetPoints()->GetPoint(1, coordsB);
//            double coordsASimulate[3];
//            double coordsBSimulate[3];
//            simulateTetra->GetEdge((int)edge)->GetPoints()->GetPoint(0, coordsASimulate);
//            simulateTetra->GetEdge((int)edge)->GetPoints()->GetPoint(1, coordsASimulate);
//            sumAfterCollapse += pow((Calculator::calcEdgeLength(coordsASimulate, coordsBSimulate) - simulateAvrgEdgeLength), 2) ;
//            sumPreCollapse += pow((Calculator::calcEdgeLength(coordsA, coordsB) - averageEdgeLength), 2);
//        }
//        double diff = sumPreCollapse - sumAfterCollapse;
//        sumOverNcells += diff;
//    }
//    simulateGrid = NULL;
//    return sumOverNcells * weight;
//}


double CostCalculator::getPointData_AlphaWater(vtkIdType pointId, vtkUnstructuredGrid *tetraGrid) {
    vtkSmartPointer<vtkDataArray> scalars_AlphaWater = tetraGrid->GetPointData()->GetArray("alpha.water");
    return scalars_AlphaWater->GetTuple1(pointId);;
}

