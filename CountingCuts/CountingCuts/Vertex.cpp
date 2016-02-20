//
//  Vertex.cpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 12/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#include "Vertex.h"
#include "Edge.h"

Vertex::Vertex(int totalVertices,int _id){
    
    adjacencyList = new Edge* [totalVertices];// since it can connect to all n-1 edges at max
    numberOfEdges = 0;
    id = _id;
    set = new int[totalVertices];
    size = 0;
    set[size++] = id;
    boss = id;
  }
  
  void Vertex::addEdge(Edge *e){
    
    adjacencyList[numberOfEdges++] = e;
  }
  
 Vertex::~Vertex(){
    
    //delete adjacencyList;
  }

