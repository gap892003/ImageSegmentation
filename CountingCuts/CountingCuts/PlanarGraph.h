//
//  PlanarGraph.h
//  CountingCuts
//
//  Created by Gaurav Joshi on 24/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#ifndef PlanarGraph_h
#define PlanarGraph_h

#include <stdio.h>
#include "Graph.h"
#include <list>

class Vertex;
class PlanarGraph : public Graph{

private:
  int numberOfFaces;
  
public:
  PlanarGraph (int numberOfVertices, int numberOfEdges);
  void findFaces();
  void findFacesRec( std::list<Edge*> *path, Vertex* startVertex );
  Graph *calculateDual();
};


#endif /* PlanarGraph_h */
