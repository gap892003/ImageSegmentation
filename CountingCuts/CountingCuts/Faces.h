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
#include <vector>

class Edge;

class Faces {

public:
  int faceID;
  Faces( std::vector<Edge*> *listOfEdges, int _faceID );
  ~Faces();
//  Edge **edgesInFace;
  std::vector<Edge*> *edgesInFace;
  void printEdges();
};

#endif /* Faces_h */
