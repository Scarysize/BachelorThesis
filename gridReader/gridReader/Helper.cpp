//
//  Helper.cpp
//  gridReader
//
//  Created by Franz Neubert on 11/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>

#include <list>
#include <set>
#include "Helper.hpp"
#include "Vertex.hpp"
#include "Cell.hpp"


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

vtkSmartPointer<vtkUnstructuredGrid> Helper::makeGrid(Tetragrid *grid) {
    vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ugrid->Allocate(grid->cells.size(), 1);
    
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(grid->vertices.size());
    vtkSmartPointer<vtkDoubleArray> pointData = vtkSmartPointer<vtkDoubleArray>::New();
    vector<vtkDoubleArray*> scalarArrays;
    for (auto scalar : *grid->vertices.at(0)->getScalars()) {
        scalarArrays.push_back(vtkDoubleArray::New());
    }
    int i = 0;
    for (auto vertex : grid->vertices) {
        double coords[3];
        vertex->getCoords(coords);
        points->InsertPoint(i, coords[0], coords[1], coords[2]);
        int scalarCounter = 0;
        for (auto scalar : *vertex->getScalars()) {
            scalarArrays.at(scalarCounter)->InsertTuple1(i, scalar);
            scalarCounter++;
        }
        i++;
    }
    
    for (auto cell : grid->cells) {
        if (!cell->deleted) {
            vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId(0, cell->vertices.at(0)->getId());
            tetra->GetPointIds()->SetId(1, cell->vertices.at(1)->getId());
            tetra->GetPointIds()->SetId(2, cell->vertices.at(2)->getId());
            tetra->GetPointIds()->SetId(3, cell->vertices.at(3)->getId());
            ugrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }
    }
    ugrid->SetPoints(points);
    int x = 0;
    for (auto scalarArray : scalarArrays) {
        ugrid->GetPointData()->AddArray(scalarArray);
        ugrid->GetPointData()->GetArray(x)->SetName(grid->vertices.at(0)->getScalarNames()->at(x).c_str());
        x++;
    }
    return ugrid;
}
