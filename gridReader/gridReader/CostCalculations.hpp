//
//  CostCalculations.hpp
//  gridReader
//
//  Created by Franz Neubert on 19/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#ifndef CostCalculations_hpp
#define CostCalculations_hpp

#include <stdio.h>
#include <vector>


class Tetragrid;
class Vertex;
class CostCalculations {
public:
    static double calcCombinedCost(Vertex *a, Vertex *b, Tetragrid *grid);
    static double calcEdgeLengthCost(Vertex *a, Vertex *b, Tetragrid *grid);
    static double calcEdgeEquityCost(Vertex *a, Vertex *b, Tetragrid *grid);
    static double calcScalarCost(Vertex *a, Vertex *b, Tetragrid *grid);
};

#endif /* CostCalculations_hpp */
