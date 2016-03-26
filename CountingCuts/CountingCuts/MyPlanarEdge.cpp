//
//  PlanarEdge.cpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 29/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#include "MyPlanarEdge.h"

MyPlanarEdge::MyPlanarEdge():Edge(){

  doneFwd = false;
  doneBakWd = false;
  faceID1 = -1;
  faceID2 = -1;
  belongsToStPath = false;
}