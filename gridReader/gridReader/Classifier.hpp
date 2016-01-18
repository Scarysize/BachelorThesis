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

class Classifier {
public:
    static void classifyVertices(std::vector<Vertex*> *vertices, std::vector<Cell*> *cells);
    static bool isBoundaryEdge(Vertex *A, Vertex *B);
    static bool isInnerEdge(Vertex *A, Vertex *B);
    
private:
    static double calcSolidAngleSum(int seed, std::vector<Vertex*> *vertices ,std::vector<Cell*> *cells);
    
    static bool isBoundaryVertex(double angle);
    static bool isInnerVertex(double angle);
    static bool isCornerVertex(double angle);
    static bool isCurveCornerVertex(double angle);
};

#endif /* Classifier_hpp */
