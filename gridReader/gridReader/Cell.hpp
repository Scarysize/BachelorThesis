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

using namespace std;

class Edge;
class Vertex;
class Cell {
public:
    // constructor
    Cell(int id, std::vector<Vertex*> vertices, vector<Edge*> edges);
    
    // properties
    int id;
    bool deleted;
    double volume;
    
    vector<Vertex*> vertices;
    vector<Edge*> edges;


    
    // member functions
    void deleteCell() {
        this->deleted = true;
    }
    
};

#endif /* Cell_hpp */
