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
    static std::vector<Vertex*> verticesFromGrid(vtkUnstructuredGrid *grid);
    
private:
    int id;
    bool boundary;
    bool interior;
    bool corner;
    
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
    
    void setModified(bool mod);
    bool wasModified();
    
    
};

#endif /* Vertex_hpp */
