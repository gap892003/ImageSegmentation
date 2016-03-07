//
//  Vertex.hpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 12/02/16.
//  Copyright © 2016 Gj. All rights reserved.

#ifndef Vertex_h
#define Vertex_h

#include <stdio.h>
#include "LinkedList.h"
class Edge;

class Vertex{
  
public:
  
  int id;
  
  int numberOfEdges;
  LinkedList<Edge*> *adjacencyList;
  //******* capstone specific
  int luminance;
  
  // union stucture
  int boss;
  int size;
  int *set;
  
  public:
  Vertex(int totalVertices,int _id);
  Node<Edge*>* addEdge(Edge *e);
//  int indexOfEdgeInList (int vertexID); // not needed
  void insertEdgesInList(int index, Edge** edgesToAdd);
  void insertEdgesInList(Node<Edge*> *nodeToBeContractedPos1, Node<Edge*> *nodeToBeContractedPos2);
  void deleteEdge (Node<Edge*> *node) ;
  void printAdjacencyList();
  ~Vertex();
};

#endif /* Vertex_h */
