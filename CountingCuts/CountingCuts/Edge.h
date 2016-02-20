//
//  Edge.hpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 12/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#ifndef Edge_h
#define Edge_h

#include <stdio.h>
#include "Constants.h"

class Edge{
  
public:
  int vertex1ID;
  int vertex2ID;
  
private:
  
  WEIGHT_TYPE weight; // Forward weight
  WEIGHT_TYPE residualWeight; // Backward weight
  
public:
  
  void setWeight( WEIGHT_TYPE weight );
  
  // getters and setters
  inline WEIGHT_TYPE getWeight( ){
    
    return this->weight;
  }
  
  inline WEIGHT_TYPE getResidualWeight( ){
    
    return this->residualWeight;
  }
  
  // method to subtract weight
  inline void setResidualWeight( WEIGHT_TYPE weight ){
    
    this->residualWeight = weight;
  }
  
  // method to subtract weight
  void subtractWeight( WEIGHT_TYPE weightToBeSubtracted );
  void subtractFromResidualWeight( WEIGHT_TYPE weightToBeSubtracted );
  
  // overloaded operator to subtract edges
  bool operator<(Edge& other){
    
    if ( this->weight < other.weight ){
      
      return true;
    }
    
    return false;
  }
};


#endif /* Edge_h */
