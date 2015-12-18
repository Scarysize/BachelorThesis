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
#include <vtkSmartPointer.h>

class EdgeCollapse {
    
public:
    EdgeCollapse(vtkIdType pointA, vtkIdType pointB, double cost);
    
public:
    std::set<vtkIdType> ncells;
    std::set<vtkIdType> icells;
    
    

    static void calcCollapsePoint(vtkIdType pointA, vtkIdType pointB, vtkUnstructuredGrid *tetraGrid,double *midpoint);
    static std::set<vtkIdType> getNCells(vtkIdType pointA, vtkIdType pointB, vtkUnstructuredGrid *tetraGrid);
    static std::set<vtkIdType> getIcells(vtkIdType pointA, vtkIdType pointB, vtkUnstructuredGrid *tetraGrid);
    static vtkSmartPointer<vtkUnstructuredGrid> simulateNcells(vtkIdType pointA, vtkIdType pointB, vtkUnstructuredGrid *tetraGrid);
    
    // GETTER
    vtkIdType getPointA();
    vtkIdType getPointB();
    
    // SETTER
    void setCost(double cost);
    void setNcells(std::set<vtkIdType> ncells);
    void setIcells(std::set<vtkIdType> icells);
    double getCost();


private:
    double cost;
    vtkIdType pointA;
    vtkIdType pointB;
};


#endif /* EdgeCollapse_hpp */
