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
#include <iostream>
#include "Faces.h"


PlanarGraph::PlanarGraph ( int numberOfVertices ):Graph(numberOfVertices){

}

// constructor
PlanarGraph :: PlanarGraph (int numberOfVertices, int numberOfEdges)
:Graph(numberOfVertices, numberOfEdges ){

}

/**
 * Recursive function to find out faces of graph
 */
void PlanarGraph::findFaces(){

  // initialise list
  faces = new std::list<Faces*> ();
  
  // pick a vertex
  // go through its adjacency list
  // for moving ahead you need to move only
  // in anticlockwise direction
  // an edge is done only when you are done exploring its 
  
  /**********TODO: Change this to iterator************************/
  // setting false for whether edge is done
  bool *seen = new bool [currentNumberOfVertices];
  
  for (int i = 0; i < currentNumberOfEdges; ++i ){
    
    edgesArray[i]->doneFwd = false;
    edgesArray[i]->doneBakWd = false;
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
      int totalEdgesInAdj = verticesArray[i]->numberOfEdges;
      
      for (int edgeNumber = 0 ; edgeNumber < totalEdgesInAdj; ++edgeNumber){
      
        Edge *edge = edgesArray[edgeNumber];
        
        if ( !edge->doneBakWd || !edge->doneFwd )// I dont think I need this
        {
          // recursively find all the faces that the edge is assosciated with
          std::list<Edge*> *path =  new std::list<Edge*>();
          path->push_back(edge);
          if (findFacesRec( path,  verticesArray[i],  verticesArray[i] )){
            
            edge->doneFwd = true;
          }
//          delete path; // path will be deleted
        }
      }
      
      seen[verticesArray[i]->id] = true;
    }
  }
  
  delete seen;
}

bool PlanarGraph::findFacesRec( std::list<Edge*> *path, Vertex *start, Vertex *lastVertex ){

  // add all edges in same direction (or opposite if vertices are reversed)
  // in the face
  
  // may have to remove an edge from this path array after done with creating a face
  // as moving up that edge wont be in the path
  
  // base case
  // fetch last edge in path see if goes to startVertex
  
  // next edge should not between same vertices
  // while adding add all edges that are between same vertices
  
  Edge *lastEdgeAdded = path->back();
  Edge *nextEdge = NULL;
  bool lastEdgeWasFwd = lastEdgeAdded->vertex1ID == lastVertex->id;
  if ( (lastEdgeWasFwd && lastEdgeAdded->vertex2ID == start->id) ||
      ( !lastEdgeWasFwd && lastEdgeAdded->vertex1ID == start->id)){
    
    // create face here and move back
    Faces *newFace = new Faces(path);
    this->faces->push_back(newFace);
    return true;
  }else{
    
    int currentVertexID = lastEdgeWasFwd?lastEdgeAdded->vertex2ID:lastEdgeAdded->vertex1ID;
    Vertex *currentVertex = verticesArray[currentVertexID];
    
    // find this edge in adjancency list of currentVertex
    // TODO: Improve this
    Edge** adj = currentVertex->adjacencyList;
    int indexOfEdge = -1;
    
    // this will find edge as well it will find next edge to add
    // COOOOL !
    int totalEdgesInList = currentVertex->numberOfEdges;//sizeof(adj)/sizeof(adj[0]);
    for (int i = 0 ; i < totalEdgesInList; ++i) {
      
      Edge *edgeUnderQ = adj[i];
      
      // thsi will find all edges between 1 and 2 and set index to last
      // matching edge
//      if (((lastEdgeWasFwd && (edgeUnderQ->vertex1ID == lastVertex->id && edgeUnderQ->vertex2ID == currentVertexID))
//          ||
//          
//          (!lastEdgeWasFwd && (edgeUnderQ->vertex2ID == lastVertex->id && edgeUnderQ->vertex1ID == currentVertexID))) && (i != totalEdgesInList-1 )){
      if ((lastEdgeWasFwd && (edgeUnderQ->vertex1ID == lastVertex->id && edgeUnderQ->vertex2ID == currentVertexID))
           ||
           
           (!lastEdgeWasFwd && (edgeUnderQ->vertex2ID == lastVertex->id && edgeUnderQ->vertex1ID == currentVertexID))){
      
        indexOfEdge = i;
        nextEdge = adj[ (i+1) % totalEdgesInList ]; // id dont know if I need cyclic,
      }
      
      // this means we have skipped all similar edges and next edge is
      // our next edge in face
//      if (indexOfEdge != -1  && indexOfEdge != totalEdgesInList-1){
//        
////        nextEdge = adj[ (i+1) % totalEdgesInList ]; // cyclic
//          nextEdge = adj[ (i+1) ]; // id dont think I need cyclic
//      }
    }

    // check if nextEdge is fwd or backward
    bool nextEdgeFwd = currentVertex->id==nextEdge->vertex1ID;
    
    // if there is no path
    if ( (indexOfEdge == -1 || nextEdge == NULL )
        // next edge is fwd and it is done Fwd
        || ( nextEdgeFwd && nextEdge->doneFwd )
        // next edge is Bwd and it is done Bwd
        || ( !nextEdgeFwd && nextEdge->doneBakWd )
        ) {
      
      delete path; // to avoid leak
      return false;
    }
    
    path->push_back(nextEdge);
    if (findFacesRec(path, start ,currentVertex )){
    
      // mark that edge as done iff face was found
      nextEdge->doneFwd = lastEdgeWasFwd;
      nextEdge->doneBakWd = !lastEdgeWasFwd;
      return true;
    }else{
      
      return false;
    }
        
//    // move ahead
//    if (lastEdgeWasFwd){
//    
//      
//    }else{
//    
//      
//    }
  }
  
  return false;
}

void PlanarGraph::printFaces(){

  for ( std::list<Faces*>::iterator it = faces->begin(); it != faces->end(); ++it ){
    
      ((Faces*)*it)->printEdges();
  }
}

Graph* PlanarGraph::calculateDual(){

  Graph *dualGraph = new Graph( this->faces->size(),this->currentNumberOfEdges);
  return dualGraph;
}

PlanarGraph::~PlanarGraph(){

  faces->clear();
  
  // TODO: check if it deletes faces
  delete faces;
}


