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
    
    Tetragrid(vector<Cell*> cells, vector<Edge*> edges, vector<Vertex*> vertices) {
        this->cells = cells;
        this->edges = edges;
        this->vertices = vertices;
    }
    
    // #riskystuff
    ~Tetragrid() {
        for (auto cell : this->cells) {
            delete cell;
        }
        for (auto edge : this->edges) {
            delete edge;
        }
        for (auto vertex : this-> vertices) {
            delete vertex;
        }
    }
    
    static Tetragrid *createGrid(vtkUnstructuredGrid *grid);
    void precalculations();
    
private:
    void precalcCells();
    void precalcEdges();
    void precalceVertices();
};

#endif /* Tetragrid_hpp */
