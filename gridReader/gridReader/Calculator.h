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
    static void calcCrossProduct(double x[3], double y[3], double crossProduct[3]);
    static void subtractVectors(double x[3], double y[3], double result [3]);
};

#endif /* Calculator_h */
