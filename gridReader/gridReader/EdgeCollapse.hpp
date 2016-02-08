//
//  EdgeCollapse.hpp
//  gridReader
//
//  Created by Franz Neubert on 05/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#ifndef EdgeCollapse_hpp
#define EdgeCollapse_hpp

#include <stdio.h>
#include <set>
#include <vector>

#include <vtkIdTypeArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include "Vertex.hpp"

class Vertex;
class EdgeCollapse {
    
public:
    EdgeCollapse(Vertex *a, Vertex *b, double cost);
    
    ~EdgeCollapse() {
        delete this->A;
        delete this->B;
    }
    
    struct CompareCost {
        bool operator()(EdgeCollapse *col1, EdgeCollapse *col2) {
            return col1->getCost() > col2->getCost();
        }
    };
    
    
public:
    // GETTER
    Vertex *getA();
    Vertex *getB();
    double getCost();
    
    // SETTER
    void setCost(double cost);

private:
    double cost;
    Vertex *A;
    Vertex *B;
};


#endif /* EdgeCollapse_hpp */
