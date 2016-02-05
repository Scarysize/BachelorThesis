//
//  Connectivity.hpp
//  gridReader
//
//  Created by Franz Neubert on 18/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#ifndef Connectivity_hpp
#define Connectivity_hpp

#include <stdio.h>
#include <set>
#include <vector>

class Cell;
class Vertex;
class Connectivity {
public:
    static std::vector<int> cellsUsingVertex(int vertex, std::vector<Cell*> *cells);
    static std::set<Cell*> getIcells(Vertex *a, Vertex *b);
    static std::set<Cell*> getNcells(Vertex *a, Vertex *b);
};

#endif /* Connectivity_hpp */
