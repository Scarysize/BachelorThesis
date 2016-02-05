//
//  GridReducer.hpp
//  gridReader
//
//  Created by Franz Neubert on 19/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#ifndef GridReducer_hpp
#define GridReducer_hpp

#include <stdio.h>
#include <vector>
#include "Cell.hpp"
#include "Vertex.hpp"
#include "EdgeCollapse.hpp"

class Tetragrid;
class GridReducer {
    
public:
    GridReducer(Tetragrid *grid) {
        this->grid = grid;
    }
    
    void run(double (*calculateCost)(Vertex *a, Vertex *b, std::vector<Cell*> cells, std::vector<Vertex*> vertices));
    
    
private:
    Tetragrid *grid;
    std::vector<EdgeCollapse*> prioq;
    
    void buildQueue(double (*calculateCost)(Vertex *a, Vertex *b, std::vector<Cell*> cells, std::vector<Vertex*> vertices));
    void doCollapse(double (*calculateCost)(Vertex *a, Vertex *b, std::vector<Cell*> cells, std::vector<Vertex*> vertices));
    void recalcQueue(double (*calculateCost)(Vertex *a, Vertex *b, std::vector<Cell*> cells, std::vector<Vertex*> vertices),
                     EdgeCollapse *lastCollapse);
  
};

#endif /* GridReducer_hpp */
