//
//  Vertex.hpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 12/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#ifndef Vertex_h
#define Vertex_h

#include <stdio.h>

class Edge;

class Vertex{
  
public:
  
  int id;
  Edge **adjacencyList;
  int numberOfEdges;
  
  //******* capstone specific
  int luminance;
  
  // union stucture
  int boss;
  int size;
  int *set;
  
  public:
  Vertex(int totalVertices,int _id);
  void addEdge(Edge *e);
  ~Vertex();
};

#endif /* Vertex_h */
