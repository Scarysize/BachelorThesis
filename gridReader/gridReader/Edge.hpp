//
//  Edge.hpp
//  gridReader
//
//  Created by Franz Neubert on 04/02/16.
//  Copyright © 2016 Franz Neubert. All rights reserved.
//

#ifndef Edge_hpp
#define Edge_hpp

#include <stdio.h>
#include "Cell.hpp"
#include "Vertex.hpp"

using namespace std;

class Edge {
public:
    Edge(Vertex *a, Vertex *b) {
        this->a = a;
        this->b = b;
    }
    
    static Edge *isEdge(Vertex *a, Vertex *b, vector<Edge*> *edges);
    
    Vertex *getA() {
        return this->a;
    }
    
    Vertex *getB() {
        return this->b;
    }
    
    void setLength(double length) {
        this->length = length;
    }
    
    double getLength() {
        return this->length;
    }
    
private:
    double length;
    Vertex *a;
    Vertex *b;
};

#endif /* Edge_hpp */
