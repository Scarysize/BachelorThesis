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
    
    bool operator< (const Edge &edge) const {
        bool result = true;
        if (min(a->getId(), b->getId()) < min(edge.a->getId(), edge.b->getId())) {
            result = true;
        } else if (min(edge.a->getId(), edge.b->getId()) < min(a->getId(), b->getId())) {
            result = false;
        } else {
            result = (max(a->getId(), b->getId()) < max(edge.a->getId(), edge.b->getId()));
        }
        return result;
    }
    
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
    
    double calcEdgeLength();
    
    struct CompareEdge {
        bool operator()(Edge *edge1, Edge *edge2) {
            if ((edge1->getA() == edge2->getA() ||
                edge1->getA() == edge2->getB())
                &&
                (edge1->getB() == edge2->getB() ||
                 edge1->getB() == edge2->getA())) {
                return true;
            }
            return false;
        }
    };

private:
    double length;
    Vertex *a;
    Vertex *b;
};

#endif /* Edge_hpp */
