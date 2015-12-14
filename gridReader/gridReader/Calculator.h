//
//  Calculator.h
//  gridReader
//
//  Created by Franz Neubert on 02/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#ifndef Calculator_h
#define Calculator_h

#include <vtkTetra.h>

class Calculator {
    
public:
    static double cellGradient(vtkTetra *cell);
    static double calcVectorLength(double vector[3]);
    static double calcEdgeLength(double pointA[3], double pointB[3]);
    static void calcCrossProduct(double x[3], double y[3], double crossProduct[3]);
    static double calcScalarProduct(double x[3], double y[3]);
    static void subtractVectors(double x[3], double y[3], double result [3]);
    static void calcMidPoint(double A[3], double B[3], double midpoint[3]);
    static double calcSolidAngle(double a[3], double b[3], double c[3], double seed[3]);
    static void calcVectorBetweenPoints(double x[3], double y[3], double result[3]);
    static double calcTetraVolume(double a[3], double b[3],double c[3],double d[3]);
};

#endif /* Calculator_h */
