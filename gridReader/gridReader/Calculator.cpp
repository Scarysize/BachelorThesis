//
//  Calculator.cpp
//  gridReader
//
//  Created by Franz Neubert on 02/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <vtkPoints.h>
#include "Calculator.h"



void Calculator::calcCrossProduct(double *x, double *y, double *crossProduct) {
    double c1[3];
    double c2[3];
    
    c1[0] = x[1] * y[2];
    c1[1] = x[2] * y[0];
    c1[2] = x[0] * y[1];
    
    c2[0] = x[2] * y[1];
    c2[1] = x[0] * y[2];
    c2[2] = x[1] * y[0];
    
    Calculator::subtractVectors(c2, c1, crossProduct);
}

double Calculator::calcScalarProduct(double *x, double *y) {
    return  x[0] * y[0] +
            x[1] * y[1] +
            x[2] * y[2];
}

void Calculator::subtractVectors(double *x, double *y, double *result) {
    result[0] = y[0] - x[0];
    result[1] = y[1] - x[1];
    result[2] = y[2] - x[2];
}

double Calculator::calcVectorLength(double *vector) {
    // std::cout << "vector length: " << sqrt(pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2)) << std::endl;
    return sqrt(pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2));
}

void Calculator::calcMidPoint(double *A, double *B, double *midpoint) {
    midpoint[0] = (A[0]+B[0])/2;
    midpoint[1] = (A[1]+B[1])/2;
    midpoint[2] = (A[2]+B[2])/2;
}

void Calculator::calcVectorBetweenPoints(double *x, double *y, double *result) {
    // y - x for vector(xy)
    result[0] = y[0] - x[0];
    result[1] = y[1] - x[1];
    result[2] = y[2] - x[2];
}

double Calculator::calcSolidAngle(double *a, double *b, double *c) {
    double bCrossC[3];
    calcCrossProduct(b, c, bCrossC);
    double tripleScalar = calcScalarProduct(a, bCrossC);
    double scalarAB = calcScalarProduct(a, b);
    double scalarAC = calcScalarProduct(a, c);
    double scalarBC = calcScalarProduct(b, c);
    double lengthA = calcVectorLength(a);
    double lengthB = calcVectorLength(b);
    double lengthC = calcVectorLength(c);
    double prodABC = lengthA * lengthB * lengthC;
    
    double solidAngle = fabs(tripleScalar) / (prodABC + scalarAB * lengthC + scalarAC * lengthB + scalarBC * lengthA);
    
    
    return atan(solidAngle / 2);
}






