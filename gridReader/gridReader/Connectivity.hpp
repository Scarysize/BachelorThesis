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
#include "Cell.hpp"

class Connectivity {
public:
    static std::vector<int> cellsUsingVertex(int vertex, std::vector<Cell*> *cells);
    static std::set<int> getIcells(int vertexA, int vertexB, std::vector<Cell*> *cells);
    static std::set<int> getNcells(int vertexA, int vertexB, std::vector<Cell*> *cells);
};

#endif /* Connectivity_hpp */
