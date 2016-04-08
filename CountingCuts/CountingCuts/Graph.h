//
//  Graph.hpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 12/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#ifndef Graph_h
#define Graph_h

#include <stdio.h>
#include "Constants.h"
#include <vector>
#include "LinkedList.h"

class Edge;
class Vertex;

class Graph{

public:
  Vertex **verticesArray;
  
protected:
  int totalVertices;
  int totalEdges;
//  Edge **edgesArray;
  int currentNumberOfVertices;
  int currentNumberOfEdges;
  int lastVertex;
//  int *visitedEdges;
  std::vector<int> *vertexPairArray; // specifically for counting cuts algorihm
  std::vector<long*> *vertexPairPathCounts; // will store array per vertexPair
  LinkedList<Edge*> *edgesArray;
public:
  
  // Constructor
  Graph (int numberOfVertices, int numberOfEdges);
  
  // constructor
  Graph ( int numberOfVertices );
  
  /**
   * Used internally to insert vertex in graph
   */
  void insertVertexInGraph( int idOfVertex );
  
  
  /**
   *  Function that takes vertex1 and 2 of the edge, creates those
   *  vertices and inserts it into the graph.
   *  returns created edge
   */
  virtual Edge* insertEdgeInGraph(int idOfVertex1, int idOfvertex2, WEIGHT_TYPE weight, bool oneWay = false);
  
  /**
   * Function to be used after max flow is found
   * It returns the matching vertex from other sets.
   */
  int getMatch(int i);

  /**
   *  Function used for debugging to print out edges in graph
   */
  void printEdges();
  
  void printEdgesBoss();
  /**
   *  Function used for debugging to print out edges passed
   */
  void printEdgesArray(Edge** array, int count);
  
  
  /**
   *  Ford Fulkerson algorithm implementation
   * @return Max flow value
   */
  WEIGHT_TYPE findMaxFlow(int s, int t);

  /**
   *  Ford Fulkerson algorithm implementation
   * @return Max flow value
   */
  bool* getMinCut(int s, int t);
  
  /**
   *  Function that finds and contracts a strongly connected 
   *  component
   */
  Graph* findAndContractSCC( int& source, int& sink );
  
  static void reverseEdges( Graph* graph );

  void deleteEdge(int id1, int id2);
  
  /**
   *
   */
  void addVertexPair (int vertex1, int vertex2);
 
  /**
   *  THE function
   **/
  long countMinCuts ();

  /**
   * Used for debugging
   */
  void printVertexPairArray();
  
  /**
   *   Sample a min cut
   */
  bool* sampleAMinCut();
  
protected:
  
  Edge** DFS( int& edgesInPath, int startIndex, int endIndex );
  /**
   *  DFS to find path
   */
  bool DFSRec( Vertex *ver, bool* seen , Edge** path, int &edgesInPath,  int endIndex );
  
  /**
   *  Function that returns path from sink to source
   */
  Edge** getOnePath(int sink, int source);
  
  
  /**
   *  BFS implementation for Edmonds Karp
   * @return Max flow value
   */
  Edge** BFS(int &edgesInPath, int s, int t );
  
  /**
   *  Finished Time version of DFS
   */
  void DFSTime(Vertex *ver, bool* seen ,int* discoveryTimeArray, int &time, std::vector<int> *verticesOrder);
  
  long countPaths (int source, int destination , bool *seen, std::vector<Edge*> &path
                  , long *countArray);
protected:
  
  Edge* insertEdge(Edge * edge, int idOfVertex1, int idOfvertex2, WEIGHT_TYPE weight, bool oneWay);
public:
  // destructor
  ~Graph();
};

#endif /* Graph_hpp */
