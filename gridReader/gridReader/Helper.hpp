//
//  Helper.hpp
//  gridReader
//
//  Created by Franz Neubert on 11/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#ifndef Helper_hpp
#define Helper_hpp

#include <stdio.h>
#include <list>
#include <vector>
#include <vtkCell.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "Tetragrid.hpp"

#endif /* Helper_hpp */

class Helper {
public:
    static bool edgesAreEqual(vtkCell *edgeA, vtkCell *edgeB);
    static std::list<vtkIdType> toStdList(vtkIdList *idList);
    static std::set<vtkIdType> toStdSet(vtkIdList *idList);
    static vtkIdList toVtkIdlist(std::set<vtkIdType> *idList);
    static void coords(std::vector<double> points, double coords[3]);
    static vtkSmartPointer<vtkUnstructuredGrid> makeGrid(Tetragrid *grid);
};