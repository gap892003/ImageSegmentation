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

Faces::Faces( std::list<Edge*> *listOfEdges ): edgesInFace(listOfEdges){

}


void Faces::printEdges(){

  std::cout << "***** Face Starts here ****" << std::endl;
  for ( std::list<Edge*>::iterator it = edgesInFace->begin(); it != edgesInFace->end(); ++it ){
  
    std::cout << ((Edge*) *it)->vertex1ID << "->" << ((Edge*) *it)->vertex2ID << std::endl;
  }
    std::cout << "***** Face ends here ****" << std::endl;
}

Faces::~Faces(){
  
  edgesInFace->clear();
  delete edgesInFace;
}