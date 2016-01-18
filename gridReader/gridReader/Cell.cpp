//
//  Cell.cpp
//  gridReader
//
//  Created by Franz Neubert on 18/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "Cell.hpp"
#include "Helper.hpp"

#include <vtkSmartPointer.h>

Cell::Cell(int id, std::vector<int> points){
    this->id = id;
    this->points = points;
    this->deleted = false;
    std::vector<int> edge0 = { this->points[0], this->points[1] };
    std::vector<int> edge1 = { this->points[0], this->points[2] };
    std::vector<int> edge2 = { this->points[0], this->points[3] };
    std::vector<int> edge3 = { this->points[1], this->points[2] };
    std::vector<int> edge4 = { this->points[2], this->points[3] };
    std::vector<int> edge5 = { this->points[3], this->points[1] };
    this->edges = { edge0, edge1, edge2, edge3, edge4, edge5 };
}

std::vector<Cell*> Cell::cellsFromGrid(vtkUnstructuredGrid *grid) {
    std::vector<Cell*> cells;
    
    for (vtkIdType cell = 0; cell < grid->GetNumberOfCells(); cell++) {
        std::set<vtkIdType> cellPointIds = Helper::toStdSet(grid->GetCell(cell)->GetPointIds());
        std::vector<int> ids;
        for (auto id : cellPointIds) {
            ids.push_back((int) id);
        }
        cells.push_back(new Cell((int)cell, ids));
    }
    
    return cells;
}