//
//  PlanarEdge.h
//  CountingCuts
//
//  Created by Gaurav Joshi on 29/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#ifndef MyPlanarEdge_h
#define MyPlanarEdge_h

#include <stdio.h>
#include "Edge.h"

class Vertex;
class Faces;
class MyPlanarEdge : public Edge{

public:
  
  // specifically for planar graphs
  // to check whether the edge is done
  bool doneFwd = false ;
  bool doneBakWd = false ;
  
  // this two if set properly by outside code
  // will help find dual of quickely
  int faceID1 = -1; // face contains fwd
  int faceID2 = -1; // face contains backward
  
  // algorithm specific
  // this will help me identify that I need one
  // vertex per this types of edge and only connect to
  // one face
  bool belongsToStPath = false;
  
  inline void addFace (int faceID, bool fwd){
    
    if (fwd) {
      
      faceID1 = faceID;
    }else{
    
      faceID2 = faceID;
    }
  }
};
#endif /* MyPlanarEdge_hpp */
