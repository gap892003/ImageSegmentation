//
//  Edge.cpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 12/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#include "Edge.h"

Edge::Edge(){

  nodeInVertex1AdjList = NULL;
  nodeInVertex2AdjList = NULL;
  nodeInMainList = NULL;
  nonDualEdge = NULL;
}

void Edge::setWeight( WEIGHT_TYPE weight ){
  
  // redundant but still here
  if ( weight < 0 ){
    
    this->weight = 0;
    this->residualWeight = -1*weight;
    return;
  }
  
  this->weight = weight;
  this->residualWeight = 0;
}

// method to subtract weight
void Edge::subtractWeight( WEIGHT_TYPE weightToBeSubtracted ){
  
  this->weight = this->weight - weightToBeSubtracted;
  this->residualWeight += weightToBeSubtracted;
}

void Edge::subtractFromResidualWeight( WEIGHT_TYPE weightToBeSubtracted ){
  
  this->weight = this->weight + weightToBeSubtracted;
  this->residualWeight -= weightToBeSubtracted;
}

Edge::~Edge(){

  if (nodeInVertex1AdjList != NULL){
  
    delete nodeInVertex1AdjList;
  }
  
  if (nodeInVertex2AdjList != NULL){
    delete nodeInVertex2AdjList;
  }
}