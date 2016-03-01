//
//  PlanarEdge.h
//  CountingCuts
//
//  Created by Gaurav Joshi on 29/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#ifndef PlanarEdge_h
#define PlanarEdge_h

#include <stdio.h>
#include "Edge.h"

class Vertex;
class Faces;
class PlanarEdge : public Edge{

public:
  
  // specifically for planar graphs
  // to check whether the edge is done
  bool doneFwd = false ;
  bool doneBakWd = false ;
  int faceID1 = -1; // face contains fwd
  int faceID2 = -1; // face contains backward
  
  // algorithm specific
  bool belongsToStPath = false;
  
  inline void addFace (int faceID, bool fwd){
    
    if (fwd) {
      
      faceID1 = faceID;
    }else{
    
      faceID2 = faceID;
    }
  }
};
#endif /* PlanarEdge_hpp */
