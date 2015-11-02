//
//  Calculator.cpp
//  gridReader
//
//  Created by Franz Neubert on 02/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include <stdio.h>
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

void Calculator::subtractVectors(double *x, double *y, double *result) {
    result[0] = y[0] - x[0];
    result[1] = y[1] - x[1];
    result[2] = y[2] - x[2];
}

double Calculator::calcVectorLength(double *vector) {
    return sqrt(pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2));
}