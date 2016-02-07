//
//  DynamicTest.hpp
//  gridReader
//
//  Created by Franz Neubert on 18/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#ifndef DynamicTest_hpp
#define DynamicTest_hpp

#include <stdio.h>
#include <set>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>

#include "Calculator.h"
#include "EdgeCollapse.hpp"
#include "Cell.hpp"
#include "Vertex.hpp"

class Tetragrid;
class DynamicTest {
public:
    static bool testVolume(std::set<vtkIdType> ncells, vtkUnstructuredGrid *gridPre, vtkUnstructuredGrid *gridAfter);
    static bool testSolidAngle(EdgeCollapse *collapse, Tetragrid *grid);
};

#endif /* DynamicTest_hpp */
