//
//  Cell.cpp
//  gridReader
//
//  Created by Franz Neubert on 18/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//
#include "Cell.hpp"
#include "Edge.hpp"
#include <vtkSmartPointer.h>

Cell::Cell(int id, vector<Vertex*> vertices, vector<Edge*> edges){
    this->id = id;
    this->vertices = vertices;
    this->deleted = false;
    this->edges = edges;
}

double Cell::avrgEdgeLength() {
    double sum = 0;
    for (auto edge : this->edges) {
        sum += edge->calcEdgeLength();
    }
    return sum / this->edges.size();
}
