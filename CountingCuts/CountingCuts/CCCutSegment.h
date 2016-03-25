//
//  CCCutSegment.h
//  CountingCuts
//
//  Created by Gaurav Joshi on 25/03/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#ifndef CCCutSegment_h
#define CCCutSegment_h

#include "CutSegment.h"

class CCCutSegment: public CutSegment {
  
public:
  
  CCCutSegment(int width, int height);
  virtual double gradient(double color1, double color2);
};

#endif /* CCCutSegment_h*/
