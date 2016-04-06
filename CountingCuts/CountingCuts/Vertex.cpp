//
//  Vertex.cpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 12/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#include "Vertex.h"
#include "Edge.h"
#include <iostream>
#include "Constants.h"

Vertex::Vertex(int totalVertices,int _id){
  
//  adjacencyList = new Edge* [totalVertices];// since it can connect to all n-1 edges at max
  adjacencyList = new LinkedList<Edge*>();
  numberOfEdges = 0;
  id = _id;
//  set = new int[totalVertices];
  size = 0;
//  set[size++] = id;
  boss = id;
}

Node<Edge*>* Vertex::addEdge(Edge *e){
  
//  adjacencyList[numberOfEdges++] = e;
  Node<Edge*> *nodeInLinkedList = adjacencyList->addValue(e);
  return nodeInLinkedList;
}

// This function will be removed after introducting linkedLists
//int Vertex::indexOfEdgeInList (int vertexID){
//  
//  for ( int i = 0 ; i < numberOfEdges ; ++i ){
//    
//    Edge* edge = adjacencyList[i];
//    
//    if (( ( (edge->vertex2ID == vertexID) && ( edge->vertex1ID == this->id ) ) )
//      || (( (edge->vertex1ID == vertexID) && ( edge->vertex2ID == this->id ) ) )){
//      
//        return i;
//    }
//  }
//  
//  return -1;
//}

void Vertex::insertEdgesInList(int index, Edge** edgesToAdd){

}

void Vertex::insertEdgesInList(Node<Edge*> *nodeToBeContractedPos1, Node<Edge*> *nodeToBeContractedPos2){

  adjacencyList->mergeLists(nodeToBeContractedPos1, nodeToBeContractedPos2);
}

void Vertex::deleteEdge (Node<Edge*> *node){

  adjacencyList->deleteNode(node);
}

void Vertex::printAdjacencyList(){
  
#ifdef DEBUG_ON
  std::cout << "********* Adj List" << id << "*********" << std::endl;
  for ( Edge *e = adjacencyList->beginIteration_debug(); e != NULL ;
       e = adjacencyList->getNextElement_debug()){
    
    std::cout << e->vertex1ID << "->" << e->vertex2ID << " W: " << e->getWeight() << " RW: " << e->getResidualWeight() << std::endl;
  }
  
  std::cout << "********* Adj List" << id << "*********" << std::endl;
#endif
}


Vertex::~Vertex(){
  
  // TODO: Deal with this
  // SHould be deleted when deleting edges 
//  delete adjacencyList;
}

