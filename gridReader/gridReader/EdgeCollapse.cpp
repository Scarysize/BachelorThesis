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

#include "Calculator.h"

EdgeCollapse::EdgeCollapse(int edgeId, vtkIdType tetraId, double cost, std::set<vtkIdType> *icells, std::set<vtkIdType> *ncells) {
    this->edgeId = edgeId;
    this->tetraId = tetraId;
    this->cost = cost;
    this->ncells = ncells;
    this->icells = icells;
}

void EdgeCollapse::calcCollapsePoint(vtkUnstructuredGrid *tetraGrid, double *midpoint) {
    double points[2][3];
    tetraGrid->GetCell(this->tetraId)->GetEdge(this->edgeId)->GetPoints()->GetPoint(0, points[0]);
    tetraGrid->GetCell(this->tetraId)->GetEdge(this->edgeId)->GetPoints()->GetPoint(1, points[1]);
    
    Calculator::calcMidPoint(points[0], points[1], midpoint);
}

const double EdgeCollapse::getCost() {
    return this->cost;
}

struct EdgeCollapse::CompareCost {
    bool operator()(EdgeCollapse &col1, EdgeCollapse &col2) {
        return col1.getCost() < col2.getCost();
    }
};