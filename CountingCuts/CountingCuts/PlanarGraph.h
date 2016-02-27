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
class Faces;
class PlanarGraph : public Graph{

private:
  std::list <Faces*> *faces;
  
public:
  
  PlanarGraph ( int numberOfVertices );
  PlanarGraph (int numberOfVertices, int numberOfEdges);
  void findFaces();
  bool findFacesRec( std::list<Edge*> *path, Vertex* startVertex,Vertex *lastVertex);
  Graph *calculateDual();
  void printFaces();
  virtual  ~PlanarGraph();
};


#endif /* PlanarGraph_h */
