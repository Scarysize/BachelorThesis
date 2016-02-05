//
//  Tetragrid.hpp
//  gridReader
//
//  Created by Franz Neubert on 05/02/16.
//  Copyright © 2016 Franz Neubert. All rights reserved.
//

#ifndef Tetragrid_hpp
#define Tetragrid_hpp

#include <stdio.h>
#include <vtkUnstructuredGrid.h>
using namespace std;

class Edge;
class Cell;
class Vertex;
class Tetragrid {
public:
    vector<Cell*> cells;
    vector<Edge*> edges;
    vector<Vertex*> vertices;
    
    static Tetragrid createGrid(vtkUnstructuredGrid *grid);
    void precalculations();
    
private:
    void precalcCells();
    void precalcEdges();
    void precalceVertices();
};

#endif /* Tetragrid_hpp */
