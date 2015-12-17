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
#include <vtkCell.h>

#endif /* Helper_hpp */

class Helper {
public:
    static bool edgesAreEqual(vtkCell *edgeA, vtkCell *edgeB);
    static std::list<vtkIdType> toStdList(vtkIdList *idList);
    static std::set<vtkIdType> toStdSet(vtkIdList *idList);
    static vtkIdList toVtkIdlist(std::set<vtkIdType> *idList);
    static std::set<vtkIdType> cellsUsingVertex(vtkIdType vertex, vtkUnstructuredGrid *tetraGrid);
};