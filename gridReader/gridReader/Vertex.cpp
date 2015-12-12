//
//  Vertex.cpp
//  gridReader
//
//  Created by Franz Neubert on 07/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include "Vertex.hpp"


Vertex::Vertex(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

Vertex::Vertex(double coords[3]) {
    this->x = coords[0];
    this->y = coords[1];
    this->z = coords[2];
}

bool Vertex::hasSameCoords(Vertex vertex) {
    if (vertex.x == this->x && vertex.y == this->y && vertex.z == this->z) {
        return true;
    } else {
        return false;
    }
}