//
//  Faces.h
//  CountingCuts
//
//  Created by Gaurav Joshi on 24/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#ifndef Faces_hpp
#define Faces_hpp
#include <stdio.h>
#include <list>
class Edge;

class Faces {

public:
  
  Faces( std::list<Edge*> *listOfEdges );
  ~Faces();
//  Edge **edgesInFace;
  std::list<Edge*> *edgesInFace;
  void printEdges();
};

#endif /* Faces_h */
