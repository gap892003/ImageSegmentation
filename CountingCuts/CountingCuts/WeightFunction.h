//
//  WeightFunction.h
//  CountingCuts
//
//  Created by Gaurav Joshi on 25/03/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#ifndef WeightFunction_h
#define WeightFunction_h

#include <cmath>
#include "Constants.h"
#include "CutPlanarDefs.h"
/**
 * Weight function
 */
static WEIGHT_TYPE weightFunction ( double &luminance , int x, int y ){
  
  //  WEIGHT_TYPE newValue = 256 - luminance;
  //  WEIGHT_TYPE newValue = 100000000 - luminance*1000;
  //  should be inverse
  //  WEIGHT_TYPE newValue = pow (2, luminance+1);
  //  WEIGHT_TYPE newValue = exp (luminance+1);
  //  double temp  = ( 255 / log(luminance+2) );// x+2 as Dont want to deal with zeroes
  //  WEIGHT_TYPE newValue = temp + 1;
  
//    double temp  = 10000000000 - pow(luminance,4);
//    WEIGHT_TYPE newValue = temp + 1;
//    return newValue;
  
  /********Working 1 ******/
//    double temp  = pow ( 255-luminance, 8);// x+2 as Dont want to deal with zeroes
//    WEIGHT_TYPE newValue = temp + 1;
  /********Working 1 ******/
  
//  double temp =  ((double)1/(luminance+1))*1000;
//  WEIGHT_TYPE newValue = pow ( temp, 4) * (x*y);
  
//  double newValue = 1/(pow(2,luminance));
//  return newValue;
  
  double diff = luminance;
  diff = pow(diff, 2);
  double weight;
  weight = exp(-sqrt(diff));// + 0.0035;
  if (weight < EPSILON)
    weight = EPSILON;
  return weight;

}


#endif /* WeightFunction_h */
