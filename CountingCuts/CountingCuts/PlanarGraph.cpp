//
//  PlanarGraph.cpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 24/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#include "PlanarGraph.h"
#include "Edge.h"
#include "Vertex.h"

// constructor
PlanarGraph :: PlanarGraph (int numberOfVertices, int numberOfEdges)
:Graph(numberOfVertices, numberOfEdges ){

}

/**
 * Recursive function to find out faces of graph
 */
void PlanarGraph::findFaces(){

  // pick a vertex
  // go through its adjacency list
  // for moving ahead you need to move only
  // in anticlockwise direction
  // an edge is done only when you are done exploring its 
  
  /**********TODO: Change this to iterator************************/
  // setting false for whether edge is done
  bool *seen = new bool [currentNumberOfVertices];
  
  for (int i = 0; i < currentNumberOfEdges; ++i ){
    
    edgesArray[i]->done = false;
  }
  
  for (int i = 0; i < currentNumberOfVertices; ++i ){
  
    seen[i] =false;
  }
  
  for (int i = 0; i < currentNumberOfVertices; ++i ){
    
    // get adjancency list
    // see edge by edge
    // if edge is done dont do anything
    // if it is not done then visit anticlockwise
    // i.e. call recursive function on it
    // IMP: Edges should only move in allowed direction
    // also face can not take an already done edge
    // when adjacency list is done mark vertex as done(seen)
    
    // catch here is seen[i] wont work, use seen [verticesArray [i]->id]
    if (!seen [verticesArray[i]->id]){
    
      Edge** edgesArray = verticesArray[i]->adjacencyList;
      
      for (int edgeNumber = 0 ; edgeNumber < (sizeof((edgesArray))/sizeof((edgesArray[0]))); ++edgeNumber){
      
        Edge *edge = edgesArray[edgeNumber];
        
        if ( !edge->done ){
          
          // recursively find all the faces that the edge is assosciated with
          std::list<Edge*> *path =  new std::list<Edge*>();
          path->push_back(edge);
          findFacesRec( path,  verticesArray[i] );
          delete path;
        }
      }
      
      seen[verticesArray[i]->id] = true;
    }
  }
}

void PlanarGraph::findFacesRec( std::list<Edge*> *path, Vertex *start ){

  // add all edges in same direction (or opposite if vertices are reversed)
  // in the face
  
  // may have to remove an edge from this path array after done with creating a face
  // as moving up that edge wont be in the path
  
  // base case
  // fetch last edge in path see if goes to startVertex
  
  //
}

Graph* PlanarGraph::calculateDual(){

  Graph *dualGraph = new Graph(numberOfFaces,this->currentNumberOfEdges);
  return dualGraph;
}