//
//  Cell.hpp
//  gridReader
//
//  Created by Franz Neubert on 18/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#ifndef Cell_hpp
#define Cell_hpp

#include <stdio.h>
#include <set>
#include <vector>

#include <vtkUnstructuredGrid.h>

class Cell {
public:
    // properties
    int id;
    std::vector<int> points;
    std::vector<std::vector<int>> edges;
    bool deleted;
    
public:
    // constructor
    Cell(int id, std::vector<int> points);
    
    // member functions
    void deleteCell(){
        this->deleted = true;
    }
    
    // STATICS
    static std::vector<Cell*> cellsFromGrid(vtkUnstructuredGrid *grid);
};

#endif /* Cell_hpp */
