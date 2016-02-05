//
//  Helper.cpp
//  gridReader
//
//  Created by Franz Neubert on 11/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include <vtkSmartPointer.h>
#include <vtkTetra.h>

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

void Helper::coords(std::vector<double> points, double *coords) {
    if (points.size() == 3) {
        coords[0] = points[0];
        coords[1] = points[1];
        coords[2] = points[2];
    }
}

//vtkSmartPointer<vtkUnstructuredGrid> Helper::makeGrid(std::vector<Cell *> cells, std::vector<Vertex *> vertices) {
//    vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
//    ugrid->Allocate(cells.size(), 1);
//    
//    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//    points->SetNumberOfPoints(vertices.size());
//    int i = 0;
//    for (auto vertex : vertices) {
//        double coords[3];
//        vertex->getCoords(coords);
//        points->InsertPoint(i, coords[0], coords[1], coords[2]);
//        i++;
//    }
//    
//    for (auto cell : cells) {
//        if (!cell->deleted) {
//            vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
//            tetra->GetPointIds()->SetId(0, cell->points[0]);
//            tetra->GetPointIds()->SetId(1, cell->points[1]);
//            tetra->GetPointIds()->SetId(2, cell->points[2]);
//            tetra->GetPointIds()->SetId(3, cell->points[3]);
//            ugrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
//        }
//    }
//    ugrid->SetPoints(points);
//    return ugrid;
//}
