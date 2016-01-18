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

#include "Vertex.hpp"
#include "Cell.hpp"

class CostCalculations {
public:
    static double calcEdgeLengthCost(Vertex a, Vertex b, std::vector<Cell*> cells, std::vector<Vertex*> vertices);
};

#endif /* CostCalculations_hpp */
