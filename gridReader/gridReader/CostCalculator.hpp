//
//  CostCalculator.hpp
//  gridReader
//
//  Created by Franz Neubert on 05/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#ifndef CostCalculator_hpp
#define CostCalculator_hpp

#include <stdio.h>
#include <list>
#include <set>

#include <vtkCell.h>
#include <vtkUnstructuredGrid.h>

class CostCalculator {
    
public:
    static double calcScalarCost(vtkIdType pointA, vtkIdType pointB, double weight, vtkUnstructuredGrid *tetraGrid);
    static double calcVolumeCost(vtkIdType pointA, vtkIdType pointB, double weight, std::set<vtkIdType> ncells, std::set<vtkIdType> icells,vtkUnstructuredGrid *tetraGrid);
    static double calcEdgeLengthCost(vtkIdType pointA, vtkIdType pointB, double weight, vtkUnstructuredGrid *tetraGrid);
    static double calcEdgeEquityCost(vtkIdType pointA, vtkIdType pointB, double weight, vtkUnstructuredGrid *tetraGrid);
    /* ...other cost calculation methods */

private:
    static double getPointData_AlphaWater(vtkIdType pointId, vtkUnstructuredGrid *tetraGrid);
    /* ...other point data retrieval methods */
};

#endif /* CostCalculator_hpp */
