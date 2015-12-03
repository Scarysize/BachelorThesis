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
    static float calcVectorLength(double vector[3]);
    static void calcCrossProduct(double x[3], double y[3], double crossProduct[3]);
    static double calcDotProduct(double x[3], double y[3]);
    static void subtractVectors(double x[3], double y[3], double result [3]);
    static void calcMidPoint(double A[3], double B[3], double midpoint[3]);
    static double calcSolidAngle(vtkCell *tetra, vtkIdType origin);
};

#endif /* Calculator_h */
