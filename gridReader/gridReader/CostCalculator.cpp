//
//  CostCalculator.cpp
//  gridReader
//
//  Created by Franz Neubert on 05/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "CostCalculator.hpp"
#include "Calculator.h"

#include <math.h>
#include <list>
#include <set>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkTetra.h>

double CostCalculator::calcScalarCost(vtkCell *edge, double weight, vtkUnstructuredGrid *tetraGrid) {
    if (edge->GetCellType() != VTK_LINE) {
        std::cerr << "cell is not a line" << std::endl;
        return 0;
    }
    vtkIdType pointA = edge->GetPointId(0);
    vtkIdType pointB = edge->GetPointId(1);
    return weight * fabs(getPointData_AlphaWater(&pointA, tetraGrid) - getPointData_AlphaWater(&pointB, tetraGrid));
}

double CostCalculator::calcVolumeCost(int edgeId, vtkIdType tetraId,double weight, std::set<vtkIdType> *ncells, std::set<vtkIdType> *icells, vtkUnstructuredGrid *tetraGrid) {
    double A[3];
    double B[3];
    double midPoint[3];
    double volume_icells = 0;
    double volume_diff_ncells = 0;
    
    tetraGrid->GetCell(tetraId)->GetEdge(edgeId)->GetPoints()->GetPoint(0, A);
    tetraGrid->GetCell(tetraId)->GetEdge(edgeId)->GetPoints()->GetPoint(1, B);
    Calculator::calcMidPoint(A, B, midPoint);
    
    /*
        Σ[T ε icells](Vol(T))
     */
    for (std::set<vtkIdType>::iterator id = icells->begin(); id != icells->end(); id++) {
        
        vtkTetra *tetra = vtkTetra::SafeDownCast(tetraGrid->GetCell(*id));
        double coords[4][3];
        tetra->GetPoints()->GetPoint(0, coords[0]);
        tetra->GetPoints()->GetPoint(1, coords[1]);
        tetra->GetPoints()->GetPoint(2, coords[2]);
        tetra->GetPoints()->GetPoint(3, coords[3]);
        volume_icells += fabs(tetra->ComputeVolume(coords[0], coords[1], coords[2], coords[3]));
        tetra = NULL;
    }
    
    /*
        Σ[T ε ncells]( Vol(T) + Vol(T_new) )
     */
    for (std::set<vtkIdType>::iterator id = ncells->begin(); id != ncells->end(); id++) {
        vtkTetra *tetra = vtkTetra::SafeDownCast(tetraGrid->GetCell(*id));
        double coords[4][3];
        tetra->GetPoints()->GetPoint(0, coords[0]);
        tetra->GetPoints()->GetPoint(1, coords[1]);
        tetra->GetPoints()->GetPoint(2, coords[2]);
        tetra->GetPoints()->GetPoint(3, coords[3]);
        
        double volume_cell = fabs(tetra->ComputeVolume(coords[0], coords[1], coords[2], coords[3]));
        
        vtkIdList *tetraPointIds = tetra->GetPointIds();
        bool foundRelevantPoint = false;
        for (vtkIdType c = 0; c < tetraPointIds->GetNumberOfIds() && !foundRelevantPoint; c++) {
            if (tetraPointIds->GetId(c) == tetraGrid->GetCell(tetraId)->GetEdge(edgeId)->GetPointId(0)) {
                coords[(int) c][0] = midPoint[0];
                coords[(int) c][1] = midPoint[1];
                coords[(int) c][2] = midPoint[2];
                foundRelevantPoint = true;
            } else if(tetraPointIds->GetId(c) == tetraGrid->GetCell(tetraId)->GetEdge(edgeId)->GetPointId(1)) {
                coords[(int) c][0] = midPoint[0];
                coords[(int) c][1] = midPoint[1];
                coords[(int) c][2] = midPoint[2];
                foundRelevantPoint = true;
            }
        }
        double volume_newCell = fabs(tetra->ComputeVolume(coords[0], coords[1], coords[2], coords[3]));
        tetra = NULL;
        
        volume_diff_ncells += volume_cell - volume_newCell;
        // std::cout << "old volume: " << volume_cell << std::endl;
        // std::cout << "new volume: " << volume_newCell << std::endl;
    }
    // w * ( Σ[T ε icells]( Vol(T) ) + Σ[T ε ncells]( Vol(T) + Vol(T_new) ) )
    // free(midPoint);
    // free(A);
    // free(B);
    return weight * (volume_diff_ncells + volume_icells);
}


std::set<vtkIdType> CostCalculator::getIntroducedTetras(int edgeId, vtkIdType tetraId, vtkUnstructuredGrid *tetraGrid) {
    std::set<vtkIdType> icells;
    
    vtkSmartPointer<vtkIdList> tetrasWithEdge = vtkSmartPointer<vtkIdList>::New();
    tetraGrid->GetCellNeighbors(tetraId, tetraGrid->GetCell(tetraId)->GetEdge(edgeId)->GetPointIds(), tetrasWithEdge);
    
    for(vtkIdType c = 0; c < tetrasWithEdge->GetNumberOfIds(); c++) {
        icells.insert(tetrasWithEdge->GetId(c));
    }
    icells.insert(tetraId); // GetCellNeighbors will ignore the tetraId, so add it manually
    tetrasWithEdge = NULL;
    // std::cout << "tetra id: " << tetraId << std::endl;
    return icells;
}

std::set<vtkIdType> CostCalculator::getNonVanishingTetras(int edgeId, vtkIdType tetraId, vtkUnstructuredGrid *tetraGrid) {
    std::set<vtkIdType> ncells;
    
    for (int p = 0; p < 2; p++) {
        vtkSmartPointer<vtkIdList> edgePoints = vtkSmartPointer<vtkIdList>::New();
        edgePoints->InsertNextId(tetraGrid->GetCell(tetraId)->GetEdge(edgeId)->GetPointId(p));
        
        vtkSmartPointer<vtkIdList> tetrasWithPoint = vtkSmartPointer<vtkIdList>::New();
        tetraGrid->GetCellNeighbors(tetraId, edgePoints, tetrasWithPoint);
        
        for(vtkIdType c = 0; c < tetrasWithPoint->GetNumberOfIds(); c++) {
            ncells.insert(tetrasWithPoint->GetId(c));
        }
        edgePoints = NULL;
        tetrasWithPoint = NULL;
    }
    
    return ncells;
}


double CostCalculator::getPointData_AlphaWater(vtkIdType *pointId, vtkUnstructuredGrid *tetraGrid) {
    vtkSmartPointer<vtkDataArray> scalars_AlphaWater = tetraGrid->GetPointData()->GetArray("alpha.water");
    double *pointData = scalars_AlphaWater->GetTuple(*pointId);
    scalars_AlphaWater = NULL;
    return *pointData;
}

void CostCalculator::printCellNeighbours (std::set<vtkIdType> *neighbours) {
    for(std::set<vtkIdType>::iterator it1 = neighbours->begin(); it1 != neighbours->end(); it1++) {
        std::cout << " " << *it1;
    }
    std::cout << std::endl;
}
