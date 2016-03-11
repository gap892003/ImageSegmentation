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
#include "LinkedList.h"

//typedef enum {
//
//  BOTTOM = 0,
//  RIGHT = 1,
//  TOP = 2,
//  LEFT = 3,
//}EdgeDirection;


class Edge{
  
public:
  int vertex1ID;
  int vertex2ID;
  
//  EdgeDirection direction; // dont need this 
  Node<Edge*> *nodeInVertex1AdjList;
  Node<Edge*> *nodeInVertex2AdjList;
  
private:
  
  WEIGHT_TYPE weight; // Forward weight
  WEIGHT_TYPE residualWeight; // Backward weight
  Edge* nonDualEdge; // used in case of dual to backtrack
  
public:
  
  Edge();
  void setWeight( WEIGHT_TYPE weight );

//  inline EdgeDirection getEdgeDirection( ){
//    
//    return direction;
//  }
//  
//  // method to subtract weight
//  inline void setEdgeDirection( EdgeDirection _direction ){
//    
//    this->direction = _direction;
//  }
  
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

  Edge* getNonDualEdge () {
    
    return nonDualEdge;
  }

  void setNonDualEdge ( Edge * _nonDualEdge) {
    
    nonDualEdge = _nonDualEdge;
  }

  // overloaded operator to subtract edges
  bool operator<(Edge& other){
    
    if ( this->weight < other.weight ){
      
      return true;
    }
    
    return false;
  }
};


#endif /* Edge_h */
