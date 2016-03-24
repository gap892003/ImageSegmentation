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
#include <vector>

class Vertex;
class Faces;
class PlanarGraph : public Graph{

private:
  std::vector <Faces*> *faces;
  int localsource, localSink;
public:
  
  PlanarGraph ( int numberOfVertices );
  PlanarGraph (int numberOfVertices, int numberOfEdges);
  void findFaces();
  int findFacesRec( std::vector<Edge*> *path, Vertex* startVertex,Vertex *lastVertex);
  Graph *calculateDual();
  void printFaces();
  virtual Edge* insertEdgeInGraph(int idOfVertex1, int idOfvertex2, WEIGHT_TYPE weight, bool oneWay = false);
  void findAndMarkSTPath();
  virtual  ~PlanarGraph();
  inline void setSource ( int source ){ localsource = source; }
  inline void setSink ( int sink ){ localSink = sink; }
  int findFacesRecNew( std::vector<Edge*> *path, Vertex *start, Vertex *lastVertex );
};


#endif /* PlanarGraph_h */
