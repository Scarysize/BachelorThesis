//
//  Vertex.cpp
//  gridReader
//
//  Created by Franz Neubert on 07/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include <vtkDataArray.h>
#include <vtkIndent.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include "Vertex.hpp"

Vertex::Vertex(int id, double coords[3]) {
    this->id = id;
    this->coords[0] = coords[0];
    this->coords[1] = coords[1];
    this->coords[2] = coords[2];
    this->boundary = false;
    this->interior = false;
    this->corner = false;
    this->deleted = false;
    this->incidents = vector<Cell*>();
}

Vertex::Vertex(int id, double coords[3], vector<double> scalars) {
    this->id = id;
    this->coords[0] = coords[0];
    this->coords[1] = coords[1];
    this->coords[2] = coords[2];
    this->boundary = false;
    this->interior = false;
    this->corner = false;
    this->deleted = false;
    this->incidents = vector<Cell*>();
    
    this->scalars = scalars;
    this->hasValues = true;
}

std::vector<Vertex*> Vertex::verticesFromGrid(vtkUnstructuredGrid *grid) {
    std::vector<Vertex*> vertices;
    vtkSmartPointer<vtkPointData> pointData = grid->GetPointData();
    if (pointData->GetNumberOfArrays() > 0) {
        vector<vtkDataArray*> scalarArrays;
        // collect scalar data arrays
        for (int i = 0; i < pointData->GetNumberOfArrays(); i++) {
            if (pointData->GetArray(i)->GetNumberOfComponents() == 1) {
                scalarArrays.push_back(pointData->GetArray(i));
                
            }
        }
        for (vtkIdType point = 0; point < grid->GetNumberOfPoints(); point++) {
            double coords[3];
            grid->GetPoint(point, coords);
            vector<double> pointScalars;
            vector<string> scalarNames;
            for (auto scalarArray : scalarArrays) {
                pointScalars.push_back(*scalarArray->GetTuple(point));
                scalarNames.push_back(string(scalarArray->GetName()));
            }
            Vertex *vertex = new Vertex((int)point, coords, pointScalars);
            vertex->setScalarNames(scalarNames);
            vertices.push_back(vertex);
        }
    } else {
        for (vtkIdType point = 0; point < grid->GetNumberOfPoints(); point++) {
            double coords[3];
            grid->GetPoint(point, coords);
            vertices.push_back(new Vertex((int)point, coords));
        }
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

bool Vertex::isDeleted() {
    return this->deleted;
}

void Vertex::deleteVertex() {
    this->deleted = true;
}

void Vertex::setModified(bool modified) {
    this->modified = modified;
}

bool Vertex::isModified() {
    return this->modified;
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

void Vertex::setHasValues(bool hasValue) {
    this->hasValues = hasValue;
}

bool Vertex::getHasValues() {
    return this->hasValues;
}

void Vertex::setScalars(vector<double> values) {
    this->scalars = values;
}

vector<double> *Vertex::getScalars() {
    return &this->scalars;
}