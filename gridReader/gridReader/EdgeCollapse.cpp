//
//  EdgeCollapse.cpp
//  gridReader
//
//  Created by Franz Neubert on 05/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "EdgeCollapse.hpp"

#include <set>
#include <vector>

#include <vtkIdTypeArray.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>

#include "Calculator.h"
#include "Helper.hpp"

EdgeCollapse::EdgeCollapse(Vertex *a, Vertex *b, double cost) {
    this->A = a;
    this->B = b;
    this->cost = cost;
}

void EdgeCollapse::setCost(double cost) {
    this->cost = cost;
}

double EdgeCollapse::getCost() {
    double cost = this->cost;
    return cost;
}

Vertex *EdgeCollapse::getA() {
    return this->A;
}

Vertex *EdgeCollapse::getB() {
    return this->B;
}
