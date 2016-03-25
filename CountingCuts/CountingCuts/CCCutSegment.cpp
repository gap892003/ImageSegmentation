//
//  CCCutSegment.cpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 25/03/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#include "CCCutSegment.h"
#include "WeightFunction.h"

CCCutSegment::CCCutSegment(int width, int height):CutSegment(width,height){


}


double CCCutSegment::gradient(double color1, double color2){

  double diff = std::abs(color1-color2);
  return weightFunction( diff , this->getWidth(), this->getHeight());
}