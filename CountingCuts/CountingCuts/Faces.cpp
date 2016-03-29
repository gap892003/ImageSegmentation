//
//  Faces.cpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 24/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#include "Faces.h"
#include "Edge.h"
#include <iostream>
#include "MyPlanarEdge.h"
#include "Constants.h"
#include "PlanarGraph.h"
#include "Vertex.h"
// constructor
Faces::Faces( std::vector<Edge*> *listOfEdges, int _faceID ): edgesInFace(listOfEdges),faceID(_faceID){

  // dont need this , adding Face id from outside
//  // add face identifier to all edges
//  for ( std::vector<Edge*>::iterator it = edgesInFace->begin(); it != edgesInFace->end(); ++it ){
//    
//      ((MyPlanarEdge*) *it)->addFace(_faceID);
//  }
  
  graph = NULL;
}

/**
 * Convinience function for debugging
 */
void Faces::printEdges(){

  if (graph == NULL) {
    
    std::cout << "Can not print edges " << std::endl;
  }
//#ifdef DEBUG_ON
  std::cout << "***** Face Starts here " << faceID << " ****" << std::endl;
  for ( std::vector<Edge*>::iterator it = edgesInFace->begin(); it != edgesInFace->end(); ++it ){
    
    MyPlanarEdge* edge = ((MyPlanarEdge*) *it);
    bool flag = edge->faceID1 == faceID;
//    std::cout << edge->vertex1ID << "->" << edge->vertex2ID << " Face1: "<< edge->faceID1 << " Face2: "<< edge->faceID2 << std::endl;
    std::cout << graph->verticesArray[edge->vertex1ID]->id << "->" << graph->verticesArray[edge->vertex2ID]->id << " Fwd: "<< flag << std::endl;

  }
    std::cout << "***** Face ends here ****" << std::endl;
//#endif
}

// destructor
Faces::~Faces(){
  
  edgesInFace->clear();
  delete edgesInFace;
}