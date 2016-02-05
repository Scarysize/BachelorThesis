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
    
    struct CompareCost {
        bool operator()(EdgeCollapse *col1, EdgeCollapse *col2) {
            return col1->getCost() > col2->getCost();
        }
    };
    
    
public:
    static bool compareVertices(EdgeCollapse &col1, EdgeCollapse &col2) {
        if ((col1.getA() == col2.getA() &&
             col1.getB() == col2.getB()) ||
            (col1.getA() == col2.getB() &&
             col1.getB() == col2.getA())) {
                return true;
            }
        return false;

    }
    std::set<int> ncells;
    std::set<int> icells;
    
    // GETTER
    Vertex *getA();
    Vertex *getB();
    double getCost();
    
    // SETTER
    void setCost(double cost);
    void setNcells(std::set<int> ncells);
    void setIcells(std::set<int> icells);


private:
    double cost;
    Vertex *A;
    Vertex *B;
    double collapsePoint[3];
};


#endif /* EdgeCollapse_hpp */
