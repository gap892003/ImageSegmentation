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

// This function will be removed after introducting linkedLists
int Vertex::indexOfEdgeInList (int vertexID){
  
  for ( int i = 0 ; i < numberOfEdges ; ++i ){
    
    Edge* edge = adjacencyList[i];
    
    if (( ( (edge->vertex2ID == vertexID) && ( edge->vertex1ID == this->id ) ) )
      || (( (edge->vertex1ID == vertexID) && ( edge->vertex2ID == this->id ) ) )){
      
        return i;
    }
  }
  
  return -1;
}

void Vertex::insertEdgesInList(int index, Edge** edgesToAdd){


}

Vertex::~Vertex(){
  
  //delete adjacencyList;
}

