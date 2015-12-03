//
//  Helper.hpp
//  gridReader
//
//  Created by Franz Neubert on 11/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#ifndef Helper_hpp
#define Helper_hpp

#include <stdio.h>
#include <vtkCell.h>

#endif /* Helper_hpp */

class Helper {
public:
    static bool edgesAreEqual(vtkCell *edgeA, vtkCell *edgeB);
};