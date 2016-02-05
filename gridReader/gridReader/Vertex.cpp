//
//  Vertex.cpp
//  gridReader
//
//  Created by Franz Neubert on 07/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "Vertex.hpp"

Vertex::Vertex(int id, double coords[3]) {
    this->id = id;
    this->coords[0] = coords[0];
    this->coords[1] = coords[1];
    this->coords[2] = coords[2];
    this->boundary = false;
    this->interior = false;
    this->corner = false;
    this->modified = true;
    this->incidents = vector<Cell*>();
}

std::vector<Vertex*> Vertex::verticesFromGrid(vtkUnstructuredGrid *grid) {
    std::vector<Vertex*> vertices;
    for (vtkIdType point = 0; point < grid->GetNumberOfPoints(); point++) {
        double coords[3];
        grid->GetPoint(point, coords);
        vertices.push_back(new Vertex((int)point, coords));
    }
    return vertices;
}

void Vertex::setToBoundary() {
    this->boundary = true;
    this->corner = false;
    this->interior = false;
}

void Vertex::setToInterior() {
    this->boundary = false;
    this->corner = false;
    this->interior = true;
}
void Vertex::setToCorner() {
    this->boundary = false;
    this->corner = true;
    this->interior = false;
}
bool Vertex::isBoundary() {
    return this->boundary;
}

bool Vertex::isCorner() {
    return this->corner;
}

bool Vertex::isInterior() {
    return this->interior;
}

int Vertex::getId() {
    return this->id;
}

double *Vertex::getCoords() {
    return this->coords;
}

bool Vertex::wasModified() {
    return this->modified;
}

void Vertex::setModified(bool mod) {
    this->modified = mod;
}

void Vertex::getCoords(double *coords) {
    coords[0] = this->coords[0];
    coords[1] = this->coords[1];
    coords[2] = this->coords[2];
}

void Vertex::setCoords(double *coords) {
    this->coords[0] = coords[0];
    this->coords[1] = coords[1];
    this->coords[2] = coords[2];
}