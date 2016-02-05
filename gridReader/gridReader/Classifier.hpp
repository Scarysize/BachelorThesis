//
//  Classifier.hpp
//  gridReader
//
//  Created by Franz Neubert on 19/12/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#ifndef Classifier_hpp
#define Classifier_hpp

#include <stdio.h>
#include <vector>
#include <set>
#include "Cell.hpp"
#include "Vertex.hpp"

class Tetragrid;
class Classifier {
public:
    static void classifyVertices(Tetragrid *grid);
    static bool isBoundaryEdge(Vertex *A, Vertex *B);
    static bool isInnerEdge(Vertex *A, Vertex *B);
    static double calcSolidAngleSum(Vertex *seed, Tetragrid *grid);

    
private:    
    static bool isBoundaryVertex(double angle);
    static bool isInnerVertex(double angle);
    static bool isCornerVertex(double angle);
    static bool isCurveCornerVertex(double angle);
};

#endif /* Classifier_hpp */
