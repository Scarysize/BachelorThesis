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

#include <vtkIdTypeArray.h>
#include <vtkUnstructuredGrid.h>

class EdgeCollapse {
    
public:
    EdgeCollapse(int edgeId, vtkIdType tetraId, double cost, std::set<vtkIdType> icells, std::set<vtkIdType> ncells);
    
public:
    void calcCollapsePoint(vtkUnstructuredGrid *tetraGrid, double *midpoint);
    double getCost();
    const vtkIdType getTetraId();
    const int getEdgeId();
    
    struct CompareCost;
    
    std::set<vtkIdType> ncells;
    std::set<vtkIdType> icells;

private:
    int edgeId;
    vtkIdType tetraId;
    double cost;
    int edgeType;
};


#endif /* EdgeCollapse_hpp */
