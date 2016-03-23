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
#include "MyPlanarEdge.h"

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
  faces = new std::vector<Faces*> ();
  
  // pick a vertex
  // go through its adjacency list
  // for moving ahead you need to move only
  // in anticlockwise direction
  // an edge is done only when you are done exploring its 
  
  /**********TODO: Change this to iterator************************/
  // setting false for whether edge is done
  bool *seen = new bool [currentNumberOfVertices];
  
//  for (int i = 0; i < currentNumberOfEdges; ++i ){
  for ( Edge *edge = edgesArray->beginIteration(); edge != NULL ;
       edge = edgesArray->getNextElement()){
  
    ((MyPlanarEdge*)edge)->doneFwd = false;
    ((MyPlanarEdge*)edge)->doneBakWd = false;
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
    
//      Edge** edgesArray = verticesArray[i]->adjacencyList;
//      int totalEdgesInAdj = verticesArray[i]->numberOfEdges;
//      for (int edgeNumber = 0 ; edgeNumber < totalEdgesInAdj; ++edgeNumber){
      for ( Edge *e = verticesArray[i]->adjacencyList->beginIteration(); e != NULL ;
           e = verticesArray[i]->adjacencyList->getNextElement()){
      
        // check if next edge is same as previous
        Node<Edge*> *currentNode = verticesArray[i]->adjacencyList->getCurrentNode();
        Edge* nextEdge = currentNode->nextNode->val;
        
        // if it is same edge
        if (( (verticesArray[nextEdge->vertex1ID]->id == verticesArray[e->vertex1ID]->id )
            &&
            (verticesArray[nextEdge->vertex2ID]->id == verticesArray[e->vertex2ID]->id))
//            ||
//            ((verticesArray[nextEdge->vertex2ID]->id == verticesArray[e->vertex1ID]->id )
//             &&
//             (verticesArray[nextEdge->vertex1ID]->id == verticesArray[e->vertex2ID]->id))
            ) {
          
          continue;
        }
        
        
        MyPlanarEdge *edge = (MyPlanarEdge*)e;
        
        if ( !edge->doneBakWd || !edge->doneFwd )// I dont think I need this
        {
          // recursively find all the faces that the edge is assosciated with
          std::vector<Edge*> *path =  new std::vector<Edge*>();
          path->push_back(edge);
          int faceIDNew = findFacesRec( path,  verticesArray[i],  verticesArray[i] );
          bool fwdEdge = verticesArray[i]->id == edge->vertex1ID;
          
          if (faceIDNew !=-1){
            
            // add a check here
            if (fwdEdge){
            
              edge->doneFwd = true;
            }else{
            
              edge->doneBakWd = true;
            }
            
            edge->addFace(faceIDNew, fwdEdge);
          }
//          delete path; // path will be deleted
        }
      }
      
      seen[verticesArray[i]->id] = true;
    }
  }
  
  delete[] seen;
}


/**
 *  Faces are found out by traversing edges in anticlockwise direction
 *  (We get Outer plane anticlockwise, inner faces clockwise)
 */
int PlanarGraph::findFacesRec( std::vector<Edge*> *path, Vertex *start, Vertex *lastVertex ){

  // add all edges in same direction (or opposite if vertices are reversed)
  // in the face
  
  // may have to remove an edge from this path array after done with creating a face
  // as moving up that edge wont be in the path
  
  // base case
  // fetch last edge in path see if goes to startVertex
  
  // next edge should not between same vertices
  // while adding add all edges that are between same vertices
  
  Edge *lastEdgeAdded = path->back();
  MyPlanarEdge *nextEdge = NULL;
  bool lastEdgeWasFwd = verticesArray[lastEdgeAdded->vertex1ID]->id == lastVertex->id;
  if ( (lastEdgeWasFwd && verticesArray[lastEdgeAdded->vertex2ID]->id == start->id) ||
      ( !lastEdgeWasFwd && verticesArray[lastEdgeAdded->vertex1ID]->id == start->id)){
    
    int newFaceID = (int)this->faces->size();
    // create face here and move back
    Faces *newFace = new Faces(path, newFaceID);// this can cause precision loss, but wont happen in our case I guess as number of faces is unlikely to be greater than 2^32 (for bit machine)
    this->faces->push_back(newFace);
    return newFaceID;
  }else{
    
    int currentVertexID = lastEdgeWasFwd?verticesArray[lastEdgeAdded->vertex2ID]->id:verticesArray[lastEdgeAdded->vertex1ID]->id;
    Vertex *currentVertex = verticesArray[currentVertexID];
    
//    Edge** adj = currentVertex->adjacencyList;
//    int indexOfEdge = -1;
    // find this edge in adjancency list of currentVertex
    // TODO: Improve this
    // this will find edge as well it will find next edge to add
    // COOOOL !
//    int totalEdgesInList = currentVertex->numberOfEdges;//sizeof(adj)/sizeof(adj[0]);
//    for (int i = 0 ; i < totalEdgesInList; ++i) {
    
    for ( Edge *edgeUnderQ = currentVertex->adjacencyList->beginIteration(); edgeUnderQ != NULL ;
           edgeUnderQ = currentVertex->adjacencyList->getNextElement()){
              
      // this will find all edges between 1 and 2 and set index to last
      // matching edge
//      if (((lastEdgeWasFwd && (edgeUnderQ->vertex1ID == lastVertex->id && edgeUnderQ->vertex2ID == currentVertexID))
//          ||
//          
//          (!lastEdgeWasFwd && (edgeUnderQ->vertex2ID == lastVertex->id && edgeUnderQ->vertex1ID == currentVertexID))) && (i != totalEdgesInList-1 )){
      if ((lastEdgeWasFwd && ( verticesArray[edgeUnderQ->vertex1ID]->id == lastVertex->id && verticesArray[edgeUnderQ->vertex2ID]->id == currentVertexID))
           ||
           
           (!lastEdgeWasFwd && ( verticesArray[edgeUnderQ->vertex2ID]->id == lastVertex->id && verticesArray[edgeUnderQ->vertex1ID]->id == currentVertexID))){
      
//        indexOfEdge = i;
//        nextEdge = (MyPlanarEdge*)adj[ (i+1) % totalEdgesInList ]; // id dont know if I need cyclic,

        Node <Edge*>* currentNode = currentVertex->adjacencyList->getCurrentNode();
        nextEdge = (MyPlanarEdge*)currentNode->nextNode->val;
      }
      
      // this means we have skipped all similar edges and next edge is
      // our next edge in face
//      if (indexOfEdge != -1  && indexOfEdge != totalEdgesInList-1){
//        
////        nextEdge = adj[ (i+1) % totalEdgesInList ]; // cyclic
//          nextEdge = adj[ (i+1) ]; // id dont think I need cyclic
//      }
    }

    // this will happen if head = tail
    if ( nextEdge == NULL ) {
      
      delete path;
      return -1;
    }
    
    // check if nextEdge is fwd or backward
    bool nextEdgeFwd = currentVertex->id==nextEdge->vertex1ID;
    
    // if there is no path
    if ( (nextEdge == NULL )
        // next edge is fwd and it is done Fwd
        || ( nextEdgeFwd && nextEdge->doneFwd )
        // next edge is Bwd and it is done Bwd
        || ( !nextEdgeFwd && nextEdge->doneBakWd )
        ) {
      
      delete path; // to avoid leak
      return -1;
    }
    
    path->push_back(nextEdge);
    
    int faceID = findFacesRec(path, start ,currentVertex );
    
    if ( faceID != -1){
    
      // mark that edge as done iff face was found
      if (nextEdgeFwd){
        
        nextEdge->doneFwd = nextEdgeFwd;
      }else{
        
        nextEdge->doneBakWd = !nextEdgeFwd;
      }
//      nextEdge->doneFwd = lastEdgeWasFwd;
//      nextEdge->doneBakWd = !lastEdgeWasFwd;
      nextEdge->addFace(faceID, nextEdgeFwd);
      return faceID;
    }else{
      
      return -1;
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
  
  return -1;
}

Edge* PlanarGraph::insertEdgeInGraph(int idOfVertex1, int idOfvertex2, WEIGHT_TYPE weight, bool oneWay){

  MyPlanarEdge *edge = new MyPlanarEdge();
  return insertEdge (edge, idOfVertex1, idOfvertex2, weight, oneWay);
}

void PlanarGraph::printFaces(){

  for ( std::vector<Faces*>::iterator it = faces->begin(); it != faces->end(); ++it ){
    
      ((Faces*)*it)->printEdges();
  }
}

/**
 * This calculates dual in linear time 
 * But it does have extra edges which needs to be removed 
 * The edges which are to be excluded as a part of path finding 
 * from t to s
 */
Graph* PlanarGraph::calculateDual(){

  // go through all faces one by one
  // for each face go through all edges
  // vertex will be represented by same faceID
  // If this face id is same face1ID connect with Face2
  // *2 is accomodating for s-t path vertices
  Graph *dualGraph = new Graph( (int)this->faces->size() * 2 ,this->currentNumberOfEdges);
  int stPathEdgeNumber = 0; // this is to keep track of stPathedgenumber so that
                      // we can add new vertex with ID numberofFaces+stPathEdgeNumber
//  for (int i = 0 ; i < currentNumberOfEdges; ++i){
  for ( Edge *edge = edgesArray->beginIteration(); edge != NULL ;
       edge = edgesArray->getNextElement()){
  
    MyPlanarEdge* currentEdge = (MyPlanarEdge*)edge;
    
    // add edge between face and face 2
    // only if it not the pathe edge
    if (!currentEdge->belongsToStPath) {
      
      // weights wont matter
      if ( currentEdge->faceID2 != -1 && currentEdge->faceID1 != -1){
        
          Edge* newAddedEdge = dualGraph->insertEdgeInGraph( currentEdge->faceID2, currentEdge->faceID1, 1 );
        newAddedEdge->setNonDualEdge(currentEdge);
      }
    }else{

      // add new vertices corresponding to edge in st path and add edge to face 2
      // at the same time store it in pairs in map to get between which vertices
      // we want to find out paths
      // IMP : for  edges with same face1 do not add edge again
      // for this simply keep track of last face1
      // TODO : verify this we require face 1 or 2?
      int newFaceID = (int)faces->size() + stPathEdgeNumber;
      
      // this will happen only in case of only 1 cut
      if ( currentEdge->faceID1 != -1 ) {
        
        Edge* newAddedEdge = dualGraph->insertEdgeInGraph( newFaceID , currentEdge->faceID1, 1 );
        newAddedEdge->setNonDualEdge(currentEdge);
        
        if (currentEdge->faceID2 == -1) {
            // shouldnt come here
          std::cout << "something went wrong " << std::endl;
        }
        
        dualGraph->addVertexPair( newFaceID, currentEdge->faceID2 );
        ++stPathEdgeNumber;
      }
    }
  }
  
  return dualGraph;
}


void PlanarGraph::findAndMarkSTPath(){

  Edge** path = NULL;
  int edgesInPath = 0;
  path = BFS(edgesInPath, localSink, localsource);
  
  for ( int i = 0 ; i < edgesInPath; ++i ) {
    
    MyPlanarEdge* edge = (MyPlanarEdge*) path[i];
    edge->belongsToStPath = true;
  }
}

PlanarGraph::~PlanarGraph(){

  faces->clear();
  
  // TODO: check if it deletes faces
  delete faces;
}


