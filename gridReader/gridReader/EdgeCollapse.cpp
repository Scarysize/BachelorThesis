//
//  EdgeCollapse.cpp
//  gridReader
//
//  Created by Franz Neubert on 05/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "EdgeCollapse.hpp"

#include <set>

#include <vtkIdTypeArray.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>

#include "Calculator.h"

EdgeCollapse::EdgeCollapse(vtkIdType pointA, vtkIdType pointB, double cost) {
    this->pointA = pointA;
    this->pointB = pointB;
    this->cost = cost;
}


void EdgeCollapse::calcCollapsePoint(vtkIdType pointA, vtkIdType pointB, vtkUnstructuredGrid *tetraGrid, double *midpoint) {
    Calculator::calcMidPoint(tetraGrid->GetPoints()->GetPoint(pointA), tetraGrid->GetPoints()->GetPoint(pointB), midpoint);
}

std::set<vtkIdType> EdgeCollapse::getIcells(vtkIdType pointA, vtkIdType pointB, vtkUnstructuredGrid *tetraGrid) {
    std::set<vtkIdType> pointAids;
    std::set<vtkIdType> pointBids;
    vtkSmartPointer<vtkIdList> pointACells = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> pointBCells = vtkSmartPointer<vtkIdList>::New();
    tetraGrid->GetPointCells(pointA, pointACells);
    tetraGrid->GetPointCells(pointB, pointBCells);
    
    for (vtkIdType point = 0; point < pointACells->GetNumberOfIds(); point++) {
        pointAids.insert(pointACells->GetId(point));
    }
    for (vtkIdType point = 0; point < pointBCells->GetNumberOfIds(); point++) {
        pointBids.insert(pointBCells->GetId(point));
    }
    
    std::set<vtkIdType> intersect;
    std::set_intersection(pointAids.begin(), pointAids.end(), pointBids.begin(), pointBids.end(), std::inserter(intersect, intersect.begin()));
    
    return intersect;
}

std::set<vtkIdType> EdgeCollapse::getNCells(vtkIdType pointA, vtkIdType pointB, vtkUnstructuredGrid *tetraGrid) {
    std::set<vtkIdType> pointAids;
    std::set<vtkIdType> pointBids;
    vtkSmartPointer<vtkIdList> pointACells = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> pointBCells = vtkSmartPointer<vtkIdList>::New();
    tetraGrid->GetPointCells(pointA, pointACells);
    tetraGrid->GetPointCells(pointB, pointBCells);
    
    for (vtkIdType point = 0; point < pointACells->GetNumberOfIds(); point++) {
        pointAids.insert(pointACells->GetId(point));
    }
    for (vtkIdType point = 0; point < pointBCells->GetNumberOfIds(); point++) {
        pointBids.insert(pointBCells->GetId(point));
    }
    
    std::set<vtkIdType> differenceAB;
    std::set_difference(pointAids.begin(), pointAids.end(), pointBids.begin(), pointBids.end(), std::inserter(differenceAB, differenceAB.begin()));
    
    std::set<vtkIdType> differenceBA;
    std::set_difference(pointBids.begin(), pointBids.end(), pointAids.begin(), pointAids.end(), std::inserter(differenceBA, differenceBA.begin()));
    
    std::set<vtkIdType> unionDifferences;

    std::set_union(differenceAB.begin(), differenceAB.end(), differenceBA.begin(), differenceBA.end(), std::inserter(unionDifferences, unionDifferences.begin()));
    
    return unionDifferences;
}


void EdgeCollapse::setCost(double cost) {
    this->cost = cost;
}

void EdgeCollapse::setNcells(std::set<vtkIdType> ncells) {
    this->ncells = ncells;
}

void EdgeCollapse::setIcells(std::set<vtkIdType> icells) {
    this->icells = icells;
}

double EdgeCollapse::getCost() {
    double cost = this->cost;
    return cost;
}

vtkIdType EdgeCollapse::getPointA() {
    return this->pointA;
}

vtkIdType EdgeCollapse::getPointB() {
    return this->pointB;
}
