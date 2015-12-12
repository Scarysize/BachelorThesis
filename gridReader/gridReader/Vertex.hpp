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

class Vertex {
public:
    Vertex(double x, double y, double z);
    Vertex(double coords[3]);
    
public:
    bool hasSameCoords(Vertex vertex);
    
public:
    double x;
    double y;
    double z;
};

#endif /* Vertex_hpp */
