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

EdgeCollapse::EdgeCollapse(int pointA, int pointB, double cost) {
    this->A = pointA;
    this->B = pointB;
    this->cost = cost;
}

void EdgeCollapse::setCollapsePoint(std::vector<Vertex> vertices){
    Calculator::calcMidPoint(vertices[this->A].getCoords(), vertices[this->B].getCoords(), this->collapsePoint);
}

void EdgeCollapse::setCost(double cost) {
    this->cost = cost;
}

void EdgeCollapse::setNcells(std::set<int> ncells) {
    this->ncells = ncells;
}

void EdgeCollapse::setIcells(std::set<int> icells) {
    this->icells = icells;
}

double EdgeCollapse::getCost() {
    double cost = this->cost;
    return cost;
}

int EdgeCollapse::getA() {
    return this->A;
}

int EdgeCollapse::getB() {
    return this->B;
}
