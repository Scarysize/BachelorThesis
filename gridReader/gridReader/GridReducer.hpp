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

using namespace std;

class Tetragrid;
class GridReducer {
    
public:
    GridReducer(Tetragrid *grid) {
        this->grid = grid;
    }
    ~GridReducer() {
        for (auto el : this->prioq) {
            delete el;
        }
        delete this->grid;
    }
    
    void run(double (*calculateCost)(Vertex *a, Vertex *b, Tetragrid *grid));
    
    Tetragrid *getGrid() {
        return this->grid;
    }

    
private:
    Tetragrid *grid;
    std::vector<EdgeCollapse*> prioq;
    
    void buildQueue(double (*calculateCost)(Vertex *a, Vertex *b, Tetragrid *grid));
    void doCollapse(double (*calculateCost)(Vertex *a, Vertex *b, Tetragrid *grid));
    void recalcQueue(double (*calculateCost)(Vertex *a, Vertex *b, Tetragrid *grid),
                     EdgeCollapse *lastCollapse);
  
};

#endif /* GridReducer_hpp */
