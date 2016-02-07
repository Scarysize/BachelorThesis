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
#include <vtkUnstructuredGrid.h>
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
    return sqrt(pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2));
}

double Calculator::calcEdgeLength(double *pointA, double *pointB) {
    return sqrt(pow(pointB[0] - pointA[0], 2) + pow(pointB[1] - pointA[1], 2) + pow(pointB[2] - pointA[2], 2));
}

void Calculator::calcMidPoint(double A[3], double B[3], double *result) {
    result[0] = (A[0] + B[0]) / 2;
    result[1] = (A[1] + B[1]) / 2;
    result[2] = (A[2] + B[2]) / 2;
}

void Calculator::calcVectorBetweenPoints(double *x, double *y, double *result) {
    // y - x for vector(xy)
    result[0] = y[0] - x[0];
    result[1] = y[1] - x[1];
    result[2] = y[2] - x[2];
}

double Calculator::calcSolidAngle(double a[3], double b[3], double c[3], double seed[3]) {
    
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
    
    
    return atan(solidAngle) * 2;
}

double Calculator::calcTetraVolumeNoAbs(double a[3], double b[3], double c[3], double d[3]) {
    double diffAD[3];
    double diffBD[3];
    double diffCD[3];
    subtractVectors(d, a, diffAD);
    subtractVectors(d, b, diffBD);
    subtractVectors(d, c, diffCD);
    
    double crossBDCD[3];
    calcCrossProduct(diffBD, diffCD, crossBDCD);
    
    return calcScalarProduct(diffAD, crossBDCD)/6;
}

double Calculator::calcTetraVolume(double a[3], double b[3], double c[3], double d[3]) {
    return fabs(calcTetraVolumeNoAbs(a, b, c, d));
}



double Calculator::calcAverageEdgeLength(vtkIdType tetra, vtkUnstructuredGrid *tetraGrid) {
    double sum = 0;
    for (vtkIdType edge = 0; edge < tetraGrid->GetCell(tetra)->GetNumberOfEdges(); edge++) {
        sum += calcEdgeLength(tetraGrid->GetCell(tetra)->GetEdge((int)edge)->GetPoints()->GetPoint(0), tetraGrid->GetCell(tetra)->GetEdge((int)edge)->GetPoints()->GetPoint(1));
    }
    return (sum / tetraGrid->GetCell(tetra)->GetNumberOfEdges());
}






