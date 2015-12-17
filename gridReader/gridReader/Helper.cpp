//
//  Helper.cpp
//  gridReader
//
//  Created by Franz Neubert on 11/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "list"
#include "set"
#include "Helper.hpp"

bool Helper::edgesAreEqual(vtkCell *edgeA, vtkCell *edgeB) {
    if (edgeA->GetCellType() != VTK_LINE || edgeB->GetCellType() != VTK_LINE) {
        std::cerr << "one of the edges isn´t of type VTK_LINE" << std::endl;
        return false;
    } else {
        if (edgeA->GetPointIds() == edgeB->GetPointIds()) {
            std::cout << "same edge" << std::endl;
            return true;
        } else {
            std::cout << "not same" << std::endl;
        }
    }
    return false;
}

std::list<vtkIdType> Helper::toStdList(vtkIdList *idList) {
    std::list<vtkIdType> list;
    for (vtkIdType id = 0; id < idList->GetNumberOfIds(); id++) {
        list.push_back(idList->GetId(id));
    }
    
    return list;
}

std::set<vtkIdType> Helper::toStdSet(vtkIdList *idList) {
    std::set<vtkIdType> set;
    for (vtkIdType id = 0; id < idList->GetNumberOfIds(); id++) {
        set.insert(idList->GetId(id));
    }
    
    return set;
}

std::set<vtkIdType> Helper::cellsUsingVertex(vtkIdType vertex, vtkUnstructuredGrid *tetraGrid) {
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
    tetraGrid->GetPointCells(vertex, cellIds);
    std::set<vtkIdType> uniqueIds = Helper::toStdSet(cellIds);
    return uniqueIds;
}
