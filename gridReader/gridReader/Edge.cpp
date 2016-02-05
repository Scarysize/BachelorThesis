//
//  Edge.cpp
//  gridReader
//
//  Created by Franz Neubert on 04/02/16.
//  Copyright © 2016 Franz Neubert. All rights reserved.
//

#include "Edge.hpp"

Edge *Edge::isEdge(Vertex *a, Vertex *b, vector<Edge *> *edges) {
    for (auto edge : *edges) {
        if ((edge->getA() == a || edge->getB() == a) &&
            (edge->getA() == b || edge->getB() == b)) {
            return edge;
        }
    }
    return nullptr;
}