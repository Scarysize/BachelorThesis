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

using namespace std;

class Cell;
class Vertex;
class Connectivity {
public:
    static vector<int> cellsUsingVertex(int vertex, std::vector<Cell*> *cells);
    static vector<Cell*> getIcells(Vertex *a, Vertex *b);
    static vector<Cell*> getNcells(Vertex *a, Vertex *b);
};

#endif /* Connectivity_hpp */
