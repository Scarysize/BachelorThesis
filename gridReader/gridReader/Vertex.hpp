//
//  Vertex.hpp
//  gridReader
//
//  Created by Franz Neubert on 07/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#ifndef Vertex_hpp
#define Vertex_hpp

#include <stdio.h>
#include <vector>

#include <vtkUnstructuredGrid.h>

using namespace std;

class Cell;
class Vertex {
public:
    Vertex(int id, double coords[3]);
    ~Vertex() {
        for (auto cell : this->incidents) {
            delete cell;
        }
    }
    static std::vector<Vertex*> verticesFromGrid(vtkUnstructuredGrid *grid);
    
private:
    int id;
    bool boundary;
    bool interior;
    bool corner;
    
    bool deleted;
    bool modified;
    double coords[3];
    
public:
    vector<Cell*> incidents;
    int getId();
    double *getCoords();
    void getCoords(double coords[3]);
    void setCoords(double coords[3]);

    bool isBoundary();
    bool isInterior();
    bool isCorner();
    
    void setToBoundary();
    void setToInterior();
    void setToCorner();
    
    void deleteVertex();
    bool isDeleted();
    
    void setModified(bool modified);
    bool isModified();
    
    
};

#endif /* Vertex_hpp */
