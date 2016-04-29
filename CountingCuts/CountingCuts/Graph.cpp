//
//  Graph.cpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 12/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//


#include "Edge.h"
#include "Vertex.h"
#include "Graph.h"
#include "Constants.h"
#include "PlanarGraph.h"
#include "MyPlanarEdge.h"
#include "CutPlanarDefs.h"
#include <iostream>
#include <list>
#include <vector>
#include <cmath>
#include <time.h>
#include <cstdlib>
using namespace std;

// Constructor
Graph::Graph (int numberOfVertices, int numberOfEdges):
totalVertices(numberOfVertices),totalEdges(numberOfEdges){
  
  verticesArray = new Vertex* [numberOfVertices];
  //  edgesArray = new Edge* [numberOfEdges];
  edgesArray = new LinkedList<Edge*> ();
  currentNumberOfEdges = 0;
  currentNumberOfVertices = 0;
  
  for ( int i = 0 ; i < numberOfVertices ; ++i){
    
    verticesArray[i] = NULL; // setting all to null
  }
  
  lastVertex = numberOfVertices-1;
//  visitedEdges = new int[numberOfEdges];
//  
//  for ( int i = 0 ; i < numberOfEdges ; ++i){
//    
//    visitedEdges[i] = 0;
//  }
  
  vertexPairArray = NULL;
  vertexPairPathCounts = NULL;
  collectedVerticesList = NULL;
}

// constructor
Graph::Graph ( int numberOfVertices ): totalVertices(numberOfVertices){
  
  verticesArray = new Vertex* [numberOfVertices+1]; // +1 as id start from 1
  edgesArray = NULL;
  currentNumberOfEdges = 0;
  currentNumberOfVertices = 0;
  for ( int i = 0 ; i < numberOfVertices ; ++i){
    
    verticesArray[i] = NULL; // setting all to null
  }
  
  vertexPairArray = NULL;
  vertexPairPathCounts = NULL;
  collectedVerticesList = NULL;
}

/**
 * Used internally to insert vertex in graph
 */
void Graph::insertVertexInGraph( int idOfVertex ){
  
  Vertex *vertex = new Vertex(totalVertices,idOfVertex);
  verticesArray[idOfVertex] = vertex;
  currentNumberOfVertices++;
}

/**
 *  Function that takes vertex1 and 2 of the edge, creates those
 *  vertices and inserts it into the graph.
 *  Default is false for one way operation
 */
Edge* Graph::insertEdgeInGraph(int idOfVertex1, int idOfvertex2, WEIGHT_TYPE weight, bool oneWay){
  
  Edge *edge = new Edge();
  return insertEdge( edge,idOfVertex1, idOfvertex2, weight, oneWay );
}


Edge* Graph::insertEdge(Edge * edge, int idOfVertex1, int idOfvertex2, WEIGHT_TYPE weight, bool oneWay){
  
  edge->vertex1ID = idOfVertex1;
  edge->vertex2ID = idOfvertex2;
  edge->setWeight( weight );
  
  if ( verticesArray[idOfvertex2] == NULL ) {
    
    insertVertexInGraph(idOfvertex2);
  }
  
  if ( verticesArray[idOfVertex1] == NULL ) {
    insertVertexInGraph(idOfVertex1);
  }
  
  Node<Edge*> *nodeIn1 =  verticesArray[idOfVertex1]->addEdge(edge);
  edge->nodeInVertex1AdjList = nodeIn1;
  
  /********** Analyse this *******************************/
  if (!oneWay){
    
    Node<Edge*> *nodeIn2 = verticesArray[idOfvertex2]->addEdge(edge);
    edge->nodeInVertex2AdjList = nodeIn2;
  }
  
  //  edgesArray[currentNumberOfEdges++] = edge;
  Node<Edge*> *mainNode = edgesArray->addValue(edge);
  edge->nodeInMainList = mainNode;
  return edge;
}

/**
 *  Function used for debugging to print out edges in graph
 */
void Graph::printEdges(){

//#ifdef DEBUG_ON
  for ( Edge *edge = edgesArray->beginIteration_debug(); edge != NULL ;
       edge = edgesArray->getNextElement_debug()){
    
    cout << verticesArray[edge->vertex1ID]->id << " " << verticesArray[edge->vertex2ID]->id << " " << edge->getWeight() << " " << edge->getResidualWeight() << endl  ;
  }
  
  cout <<endl;
//#endif
}

/**
 *  Function used for debugging to print out edges in graph
 */
void Graph::printEdgesBoss(){

//#ifdef DEBUG_ON
  for ( Edge *edge = edgesArray->beginIteration_debug(); edge != NULL ;
       edge = edgesArray->getNextElement_debug()){
    
    cout << verticesArray[edge->vertex1ID]->boss << " " << verticesArray[edge->vertex2ID]->boss << "  " << " W: "<< edge->getWeight() << "  RW: " << edge->getResidualWeight() << endl  ;
  }
  
  cout <<endl;
//#endif
  
}

/**
 *  Function used for debugging to print out edges passed
 */
void Graph::printEdgesArray(Edge** array, int count){

//#ifdef DEBUG_ON
  for (int i = 0 ; i < count ; ++i) {
    
    cout << array[i]->vertex1ID << " " << array[i]->vertex2ID << " " <<array[i]->getWeight() << endl  ;
  }
  
//#endif
}

Edge** Graph::DFS( int& edgesInPath, int startIndex,  int endIndex){
  
  Edge** path = new Edge*[totalVertices];
  bool *seen = new bool [totalVertices];
  
  for (int i = 0 ; i < totalVertices ; ++i){
    
    seen[i] = false;
  }
  
  DFSRec(verticesArray[startIndex], seen, path, edgesInPath, endIndex);
  delete[] seen;
  return path;
}


/**
 *  DFS to find path
 */
bool Graph::DFSRec( Vertex *ver, bool* seen , Edge** path, int &edgesInPath
                   ,int endIndex){
  
  // Here the basic logic is vertex1 is first end point and vertex2 is second endpoint, (Vertex 2 lies on the arrowhead)
  // So when we are traversing if we use vertex2 to vertex1 path we are using residual edge, if residual weight is greater than zero then only we can use that
  // when we are traversing if we use vertex1 to vertex2 path we are using normal edge and weight, if weight is greater than zero then only we can use that
  
  seen[ver->id] = true;
  
  if ( ver->id == endIndex ){
    
    return true;
  }
  
  // finding bottleneck
  //  Edge** array = ver->adjacencyList;
  //  int totalNeighbours = ver->numberOfEdges;
  //  for ( int i = 0; i < totalNeighbours; ++i) {
  //    Edge *edgeUnderQ = array[i] ;
  
  for ( Edge *edgeUnderQ = ver->adjacencyList->beginIteration(); edgeUnderQ != NULL ; edgeUnderQ = ver->adjacencyList->getNextElement()){
#ifdef DEBUG_ON
    
    cout << edgeUnderQ->vertex1ID << " " << edgeUnderQ->vertex2ID << " "
    << edgeUnderQ->getWeight() << " " << edgeUnderQ->getResidualWeight() << endl;
#endif
    
    // if we reach snk stop
    // Dont need this now
    //    if (
    //        ((edgeUnderQ->vertex1ID == lastVertex) && edgeUnderQ->getResidualWeight() != 0)
    //
    //        ||
    //
    //        ((edgeUnderQ->vertex2ID == lastVertex) && ( edgeUnderQ->getWeight() != 0))
    //        ){
    //
    //      if (path != NULL){
    //
    //        path[edgesInPath++] = edgeUnderQ;
    //      }
    //      return true;
    //
    //    }else
    
    if ( (!(seen[edgeUnderQ->vertex1ID] )) &&  isgreaterequal (edgeUnderQ->getResidualWeight(),EPSILON)){
      
      // using residual edge
      if ( DFSRec(verticesArray[edgeUnderQ->vertex1ID], seen, path, edgesInPath,endIndex) ){
        
        if (path != NULL){
          
          path[edgesInPath++] = edgeUnderQ;
        }
        return true;
      }
      
    }else if( (!(seen[edgeUnderQ->vertex2ID] )) && isgreaterequal(edgeUnderQ->getWeight(),EPSILON )){
      
      // using forward edge
      if (DFSRec(verticesArray[edgeUnderQ->vertex2ID], seen, path, edgesInPath,endIndex)){
        
        if (path != NULL){
          
          path[edgesInPath++] = edgeUnderQ;
        }
        return true;
      }
    }
  }
  
  return false;
}


/**
 *  Ford Fulkerson algorithm implementation
 *  s - Index of start point
 *  t - Index of endPoint point
 * @return Max flow value
 */
WEIGHT_TYPE Graph::findMaxFlow(int s, int t){
  
  int edgesInPath = 0;
  Edge **path = NULL;
#ifdef USE_EDMONDS_KARP
  
  path = BFS(edgesInPath, s,t);
#else
  
  path = DFS(edgesInPath, s,t);
#endif
  
  
  // edges that we get are in reverse order
  WEIGHT_TYPE maxFlow = 0 ;
  
  while ( edgesInPath != 0 ) {
    
#ifdef DEBUG_ON
    /***********REMOVE THIS************************************/
    cout << "Path is: " << endl;
    for (int i = 0 ; i < edgesInPath ; ++i){
      
      cout << path[i]->vertex1ID << "->" << path[i]->vertex2ID << " "<<  " W: "
      << path[i]->getWeight() << " RW: " << path[i]->getResidualWeight() << endl;
    }
    /***********REMOVE THIS************************************/
#endif
    
    int localLastVertex = t;
    WEIGHT_TYPE bottleNeckWeight = -1;
    
    // find out minimum weight in the path
    for ( int i = 0 ; i < edgesInPath ; ++i ){
      
      WEIGHT_TYPE weightToCompare = 0;
      
      // this means we have considered normal weight (forward edge)
      if (path[i]->vertex2ID == localLastVertex ){
        
        localLastVertex = path[i]->vertex1ID;
        weightToCompare = path[i]->getWeight();
        
      }else if (path[i]->vertex1ID == localLastVertex ){
        
        // this means we have considered residual weight (residual edge)
        localLastVertex = path[i]->vertex2ID;
        weightToCompare = path[i]->getResidualWeight();
      }
      
      if ( isless(weightToCompare,bottleNeckWeight) || bottleNeckWeight == -1) {
        
        bottleNeckWeight = weightToCompare;
      }
      
    }
    
    localLastVertex = t;
    // subtract that much
    for ( int i = 0 ; i < edgesInPath ; ++i ){
      
      // this means we have considered normal weight
      if (path[i]->vertex2ID == localLastVertex ){
        
        localLastVertex = path[i]->vertex1ID;
        path[i]->subtractWeight(bottleNeckWeight);
        
      }else if (path[i]->vertex1ID == localLastVertex ){
        
        // this means we have considered residual weight
        localLastVertex = path[i]->vertex2ID;
        path[i]->subtractFromResidualWeight(bottleNeckWeight);
      }
    }
    
    maxFlow += bottleNeckWeight;
    delete[] path;
    path = NULL;
    edgesInPath = 0;
#ifdef USE_EDMONDS_KARP
    
    path = BFS(edgesInPath, s,t);
#else
    
    path = DFS(edgesInPath, s,t);
#endif

  }
  
  if (path != NULL) {
    
    delete[] path;
    path = NULL;
  }
  
  return maxFlow;
}

/**
 *  Ford Fulkerson algorithm implementation
 *  calls findMaxFlow internally
 *
 *  @return Vertices that are in min cut
 */
bool* Graph::getMinCut(int s, int t){
  
  WEIGHT_TYPE maxFlow =  findMaxFlow(s, t);
  std::cout << "Max flow is : "<< maxFlow << std::endl;
  bool *seen = new bool [totalVertices];
  //Vertex** minCut = new Vertex* [totalVertices];
  
  for (int i = 0 ; i < totalVertices ; ++i){
    
    seen[i] = false;
  }
  
  int count = 0; // dummy wont be used
  DFSRec(verticesArray[s], seen, NULL, count, t);
  
  // go through all seen vertices and add them to
  // min cut Vertices array;
/*  count = 0;
  std::cout << "Min cut is { " ;
  for (int i = 0 ; i < totalVertices ; ++i){
    
    if (seen[i]){
      minCut[count++] = verticesArray[i];
      std::cout << i << ", " ;
    }
  }
  
  std::cout << " }"<< std::endl;
*/ 
//  delete[] seen;
//  seen = NULL;
  return seen;
//  return minCut;
}



/**
 *  Function that returns path from sink to source
 */
Edge** Graph::getOnePath(int sink, int source){
  
  Edge** path = new Edge*[totalVertices];
  bool *seen = new bool [totalVertices];
  int edgesCount = 0;
  
  for (int i = 0 ; i < totalVertices ; ++i){
    
    seen[i] = false;
  }
  
  // TODO: Change this
  DFSRec(verticesArray[sink], seen, path, edgesCount,source);
  
  //    printEdgesArray(path, edgesInPath);
  //    cout << endl;
  delete[] seen;
  return path;
}

/**
 * Function to be used after max flow is found
 * It returns the matching vertex from other sets.
 */
int Graph::getMatch(int i){
  
  // go through all neighbors to find whcih edge was finalized
  // to be part of maximum flow
  // We identify it by forward weight zero
  
  Vertex *ver = verticesArray[i] ;
  //  Edge** array = ver->adjacencyList;
  //  int totalNeighbours = ver->numberOfEdges;
  //  for ( int j = 0; j < totalNeighbours; ++j) {
  
  for ( Edge *edgeUnderQ = ver->adjacencyList->beginIteration(); edgeUnderQ != NULL ; edgeUnderQ = ver->adjacencyList->getNextElement()){
    
    if ((edgeUnderQ->vertex1ID == i) && edgeUnderQ->getWeight() == 0){
      
      return edgeUnderQ->vertex2ID;
    }
  }
  
  return -1;
}


/**
 * Implementing BFS for edmonds-Karp
 */
Edge** Graph::BFS(int &edgesInPath, int s, int t ){
  
  list<Vertex*> queue;
  Edge** path = new Edge*[totalVertices];
  Edge** traversedEdges = new Edge*[totalVertices];
  bool *seen = new bool [totalVertices];
  
  for (int i = 0 ; i < totalVertices ; ++i){
    
    seen[i] = false;
    traversedEdges[i] = NULL;
  }
  
  // put sink in queue
  queue.push_back(verticesArray[s]);
  
  while ( !queue.empty() ) {
    
    Vertex* vertexToVisit = queue.front();
    seen[vertexToVisit->id] = true;
    queue.pop_front();
    
    if (vertexToVisit->id == t){
      
      // done with BFS
      queue.clear();
      break;
    }
    
    //    Edge** array = vertexToVisit->adjacencyList;
    //    int totalNeighbours = vertexToVisit->numberOfEdges;
    //    for ( int i = 0; i < totalNeighbours; ++i) {
    //      Edge *edgeUnderQ = array[i] ;
    for ( Edge *edgeUnderQ = vertexToVisit->adjacencyList->beginIteration(); edgeUnderQ != NULL ; edgeUnderQ = vertexToVisit->adjacencyList->getNextElement()){
      
      // add vertices to queue if they are reachable
      // via residual or
      // verticesArray[edgeUnderQ->vertex1ID]->id is requred for planar graphs
      if ( (!(seen[verticesArray[edgeUnderQ->vertex1ID]->id] )) &&  isgreaterequal ( edgeUnderQ->getResidualWeight() , EPSILON)){
        
        if ( traversedEdges[verticesArray[edgeUnderQ->vertex1ID]->id] != NULL ){
          
          continue;
        }
        
        // using backward edge
        queue.push_back( verticesArray[edgeUnderQ->vertex1ID] );
        traversedEdges[verticesArray[edgeUnderQ->vertex1ID]->id] = edgeUnderQ;
        
      }else if( (!(seen[verticesArray[edgeUnderQ->vertex2ID]->id] )) && isgreaterequal( edgeUnderQ->getWeight() , EPSILON)){
        
        if ( traversedEdges[verticesArray[edgeUnderQ->vertex2ID]->id] != NULL ){
          
          continue;
        }
        
        // using forward edge
        queue.push_back( verticesArray[edgeUnderQ->vertex2ID] );
        traversedEdges[verticesArray[edgeUnderQ->vertex2ID]->id] = edgeUnderQ;
      }
    }
  }
  
  bool reached = false;
  int lastVertexId = verticesArray[t]->id;
  
  // could not reach sink
  if ( traversedEdges[t] == NULL) {
    
    delete [] traversedEdges;
    delete[] seen;
    delete[] path;
    return NULL;
  }
  
  
  while (!reached) {
    
    Edge *edgeInPath = traversedEdges[lastVertexId];
    path[edgesInPath++] = edgeInPath;
    
    // this means it was a forward edge(not residual edge) last vertex will be
    // vertex 1, vice versa
    if ( verticesArray[edgeInPath->vertex2ID]->id == lastVertexId ){
      
      lastVertexId = verticesArray[edgeInPath->vertex1ID]->id;
      // we will examine how we reached vertex 1
      // i.e. using which edge
    }else{
      
      // this means residual edge was used to reach vertex
      lastVertexId = verticesArray[edgeInPath->vertex2ID]->id;
    }
    
    if (lastVertexId == s){
      
      reached = true;
    }
  }
  
  delete[] seen;
  seen = NULL;
  delete[] traversedEdges;
  traversedEdges = NULL;
  return path;
}

/**
 *  Function finds Strongly connected Components, stores them
 *  then process ( contracts them )
 */
Graph* Graph::findAndContractSCC ( int &source, int& sink ){
  
  // Run DFS On residual graph
  bool *seen = new bool [totalVertices];
  int *finishTime =  new int[totalVertices];
  std::vector<int> *verticesOrder = new std::vector<int>();// this array is to store sorted
  // order of vertices during DFS
  int numberOfSCCFound = 0;
  collectedVerticesList = new vector<std::vector<int> *>();
  
  int time =  0 ;
  for (int i = 0 ; i < totalVertices ; ++i){
    
    finishTime[i] = -1;
    seen[i] = false;
  }
  
  for (int i = 0 ; i < totalVertices ; ++i){
    
    if (! seen[i]){
      
      DFSTime(verticesArray[i], seen, finishTime, time , verticesOrder);
    }
  }
  
  // reverse edges
  Graph::reverseEdges(this);
  
  // go through vertices order array
  // run DFS again
  for (int i = 0 ; i < totalVertices ; ++i){
    
    seen[i] = false;
  }
  
  // here it is reverse ordering of discovery time
  for (int i = totalVertices-1 ; i >= 0 ; --i){
    
    if (! seen[verticesOrder->at(i)]){
      
      int dummy = 0;
      
      // using int array here
      // dont want to use vector as , it can hamper
      // time complexity
      std::vector<int> *verticesCollected = new std::vector<int>();
      
      DFSTime( verticesArray[verticesOrder->at(i)], seen, NULL, dummy , verticesCollected);
      
      ++numberOfSCCFound;
      // do something with collcted vertices
      // right now just printing
#ifdef DEBUG_ON
      cout << "Strongly connected Component: Size:" << verticesCollected->size() << " Component: " ;
#endif
      // settting boss for each component
      // this will be useful while creating new graph
      int bossID =  verticesArray[verticesCollected->size()-1]->id;
      
      for (int k = 0 ; k < dummy ; ++k){
        
#ifdef DEBUG_ON
        cout << verticesCollected->at(k) << " ";
#endif
        verticesArray[verticesCollected->at(k)]->boss = bossID;
      }
      
#ifdef DEBUG_ON
      cout << endl;
#endif
      collectedVerticesList->push_back(verticesCollected);
      //      delete verticesCollected;
    }
  }
  
  // reverse again to make it normal
  Graph::reverseEdges(this);
  
  
  /********************************APPROACH 1*********************************/
  /*
   // create new graph here
   // while creating a new graph make sure that vertex ID of
   // all vertices in strongly connected component is same as
   // boss
   
   // using currentNumberOfvertices as vertices count as
   // there is a flaw that I am creating edges between "OLD"
   // ids. SO Graph should allocate exta space for that
   // 1____5____12 can be a case
   Graph *newGraph = new Graph(currentNumberOfVertices, currentNumberOfEdges);
   
   for ( int i = 0 ; i < currentNumberOfEdges ; ++i ){
   
   Edge* edge = edgesArray[i];
   cout << "Bosses:  " << verticesArray[edge->vertex1ID]->boss << "  "<< verticesArray[edge->vertex2ID]->boss<<endl;
   // see if their boss is same if yes then skip that edge
   if ( verticesArray[edge->vertex1ID]->boss == verticesArray[edge->vertex2ID]->boss ){
   
   continue;
   }else{
   
   // we can edge between two vertices. with vertex id as boss
   Edge* newEdge = newGraph->insertEdgeInGraph( verticesArray[edge->vertex1ID]->boss, verticesArray[edge->vertex2ID]->boss, edge->getWeight());
   newEdge->setResidualWeight(edge->getResidualWeight());
   }
   }
   */
  
  /********************************END OF APPROACH 1*************/
  
  /********************************APPROACH 2, 3 common part  *********************************/

  // create new replica of the graph with exact same bosses
  // vertices, then use edge contraction
  // IMP: while creating new graph, dont worry about residual graphs
  PlanarGraph *newGraph = new PlanarGraph (currentNumberOfVertices, currentNumberOfEdges*2);
  
  //currentNumberOfEdges*2 safer side to account for residual edges
  
  //  for ( int i = 0 ; i < currentNumberOfEdges ; ++i ){
  
  for ( Edge *edge = edgesArray->beginIteration(); edge != NULL ;
       edge = edgesArray->getNextElement()){
    
    //    Edge* edge = edgesArray[i];
    // We need double way adjacency list to find out faces of the graph    
    if ( isgreaterequal (edge->getWeight(),EPSILON) ){
      Edge* newEdge = newGraph->insertEdgeInGraph( edge->vertex1ID, edge->vertex2ID, 1 );
      newGraph->verticesArray[newEdge->vertex1ID]->boss = verticesArray[edge->vertex1ID]->boss;
      newGraph->verticesArray[newEdge->vertex2ID]->boss = verticesArray[edge->vertex2ID]->boss;
    }
    
    // if there is residual edge add actual edge
    if ( isgreaterequal (edge->getResidualWeight(),EPSILON) ){
      
      Edge* newEdgeRes = newGraph->insertEdgeInGraph( edge->vertex2ID, edge->vertex1ID, 1);
      newGraph->verticesArray[newEdgeRes->vertex1ID]->boss = verticesArray[edge->vertex2ID]->boss;
      newGraph->verticesArray[newEdgeRes->vertex2ID]->boss = verticesArray[edge->vertex1ID]->boss;
    }
  }

  
  // IMP:
  // CAN NOT DO FOLLOWING
  // EDGES CANT HAVE RESIDUAL WEIGHTS AFTER THIS POINT
  //PlanarGraph *newGraph = (PlanarGraph*)this;
  
  /**************************APPROACH 2, 3 common part ends  ***********************/
  
  
  /**************************APPROACH 2 ONLY ***********************/
  // Do edge contraction here,
  // go through graph one vertex at a time and then contract an
  // edge
  // this is poor version of this I will imporve on it later
  // IMP: Dont care about residual edges here
  
  // iterate through all strongly connected components collected above
  
  //  for ( int i = 0 ; i < newGraph->currentNumberOfEdges ; ++i ){
  
  /*  for ( int i = 0 ; i < collectedVerticesList.size() ; ++i ){
   
   int *verticesWithSameBoss = collectedVerticesList[i];
   
   int verticesCount = (sizeof((verticesWithSameBoss))/sizeof((verticesWithSameBoss[0])));
   
   // take vertex 1 and contract every vertex
   Vertex* ver1 = newGraph->verticesArray[verticesWithSameBoss[0]];
   
   for ( int j = 1 ; j < verticesCount ; ++j ){
   
   Vertex* ver2 = newGraph->verticesArray[verticesWithSameBoss[0]];
   
   // find edge with these two vertices and contract that
   // while contracting just change the pointer at
   // the location of vertexID, that would change pointers
   // of all other edges
   // so verticesArray[VertexID] will give a different vertex now
   // And we don't need to change pointers in incoming edges to that
   // vertex
   // IMP: Identify vertices by their boss names now
   // Awesome ! Right ?
   newGraph->verticesArray[ver2->id] = ver1;
   
   // now add all edges from ver2 to ver1's list
   // make sure the sequence of edges is same
   // TODO: writing it in a bad way, improve this
   int index = ver1->indexOfEdgeInList(ver2->id);
   
   // get edges array in order and then add it to the vertex
   ver1->insertEdgesInList( index , NULL);
   }
   }
   */
  /*******************APPROACH 2 ONLY ENDS********************************/
  
  /**********************APPROACH 3 ONLY Start ****************************/
  
  // go through all edges, see which needs to be contracted and contract it
  //  for ( Edge *edge = newGraph->edgesArray->beginIteration(); edge != NULL ;
  //       edge = newGraph->edgesArray->getNextElement()){
  
  for ( int i = 0 ; i < collectedVerticesList->size() ; ++i ){
    
    vector<int> *verticesWithSameBoss = collectedVerticesList->at(i);
    
    // take vertex 1 and contract every vertex
    Vertex* ver1 = newGraph->verticesArray[verticesWithSameBoss->at(verticesWithSameBoss->size()-1)];
    
    for ( int j = ( (int)verticesWithSameBoss->size()-2 ) ; j >= 0 ; --j ){
      
      Vertex* ver2 = newGraph->verticesArray[verticesWithSameBoss->at(j)];
      Edge *edge = NULL;
      
      // check if these are alreay merged
      // when we are moving according to verticesCollected list
      // ver2 == NULL should not ever happen
      // if it does then something is wrong
      if (ver1 == ver2 || ver2 == NULL) {
        
        cout << "Ignoring ver2" << endl;
        continue;
      }

      // find edge to contract
      // following is constant time as at the max there can be 8 edges
      for ( Edge *tempEdge = ver2->adjacencyList->beginIteration(); tempEdge != NULL ;
           tempEdge = ver2->adjacencyList->getNextElement()){
        
        // find matching edge
        if ( newGraph->verticesArray[tempEdge->vertex1ID]->id == ver1->id ||
            newGraph->verticesArray[tempEdge->vertex2ID]->id == ver1->id) {
          
          // this means there are duplicate edges and we should delete it
          if (edge != NULL){
            
            // delete duplicate edges
            if (newGraph->verticesArray[tempEdge->vertex1ID]->id == ver1->id &&  newGraph->verticesArray[tempEdge->vertex2ID]->id == ver2->id
                ) {
              
              ver1->deleteEdge(tempEdge->nodeInVertex1AdjList);
              ver2->deleteEdge(tempEdge->nodeInVertex2AdjList);
              
              // Backward edge
            }else if ( newGraph->verticesArray[tempEdge->vertex2ID]->id == ver1->id &&  newGraph->verticesArray[tempEdge->vertex1ID]->id == ver2->id){
              
              ver1->deleteEdge(tempEdge->nodeInVertex2AdjList);
              ver2->deleteEdge(tempEdge->nodeInVertex1AdjList);
              
            }
            
            tempEdge->nodeInVertex2AdjList = NULL;
            tempEdge->nodeInVertex1AdjList = NULL;
            newGraph->edgesArray->deleteNode(tempEdge->nodeInMainList);
            delete tempEdge;
          }else{
            
            edge = tempEdge;
          }
        }
      }
      
      if (edge == NULL) {
        
        std::cerr << "Something went wrong";
        continue;
      }
//    for ( Edge *edge = ver1->adjacencyList->beginIteration(); edge != NULL ;
//         edge = ver1->adjacencyList->getNextElement()){
//      
//      Vertex* ver2 = newGraph->verticesArray[edge->vertex1ID]->id == ver1->id ?  newGraph->verticesArray[edge->vertex2ID]: newGraph->verticesArray[edge->vertex1ID];
      
//      if ( ver1->boss == ver2->boss ){ // redundant in new approach
        
        // check if source or sink is being merged and change source and sink
        if ( ver1->id == source || ver2->id == source ){
          
          source = ver1->id;
        }
        
        if ( ver1->id == sink || ver2->id == sink){
          
          sink = ver1->id;
        }
        
#ifdef DEBUG_ON
        std::cout << " Before merging "<<  ver1->id << " " << ver2->id << std::endl;
        ver1->printAdjacencyList();
        ver2->printAdjacencyList();
#endif
        // node in ver1 list
        Node<Edge*> *mainNode  = newGraph->verticesArray[edge->vertex1ID]->id == ver1->id ? edge->nodeInVertex1AdjList:edge->nodeInVertex2AdjList;
        
        // node in vertex 2 list
        Node<Edge*> *secondNode  = newGraph->verticesArray[edge->vertex1ID]->id == ver1->id ? edge->nodeInVertex2AdjList:edge->nodeInVertex1AdjList;
/*
        // before calling insert we need to find all edges between these two
        // vertices for maintaining consistency
        // we will select last node from list2 and first node from list 1
        // will have to do this for forward as well reverse direction
        Node<Edge*> *tempNode = NULL;
        
        if ( mainNode != NULL ) {
          
          tempNode = mainNode->nextNode;
        }
        
        //deleting from list 1
        while (true && tempNode != NULL) {
          
          Edge* tempEdge = tempNode->val;
          
          // if edge is between same vertices then remove it from both lists
          // we will remove it next iteration
          // FWD edge
          if (newGraph->verticesArray[tempEdge->vertex1ID]->id == ver1->id &&  newGraph->verticesArray[tempEdge->vertex2ID]->id == ver2->id
              ) {
            
            ver1->deleteEdge(tempNode);
            ver2->deleteEdge(tempEdge->nodeInVertex2AdjList);
           
            // Backward edge
          }else if ( newGraph->verticesArray[tempEdge->vertex2ID]->id == ver1->id &&  newGraph->verticesArray[tempEdge->vertex1ID]->id == ver2->id){
            
            ver1->deleteEdge(tempNode);
            ver2->deleteEdge(tempEdge->nodeInVertex1AdjList);
            
          }else{
            
            break;
          }
          
          tempEdge->nodeInVertex2AdjList = NULL;
          tempEdge->nodeInVertex1AdjList = NULL;
          newGraph->edgesArray->deleteNode(tempEdge->nodeInMainList);
          
          if ( mainNode != NULL && mainNode->val != NULL ){
            
            tempNode = mainNode->nextNode;
          }else{
            
            tempNode = NULL;
          }
        }
        
        tempNode = NULL;
        
        if (  mainNode != NULL ) {
          
          tempNode = mainNode->prevNode;
        }
        
        //deleting from list 1
        while (true && tempNode!= NULL) {
          
          Edge* tempEdge = tempNode->val;
          // if edge is between same vertices then remove it from both lists
          // we will remove it next iteration
          if (newGraph->verticesArray[tempEdge->vertex1ID]->id == ver1->id &&  newGraph->verticesArray[tempEdge->vertex2ID]->id == ver2->id
              ) {
            
            ver1->deleteEdge(tempNode);
            ver2->deleteEdge(tempEdge->nodeInVertex2AdjList);
            
          }else if (newGraph->verticesArray[tempEdge->vertex2ID]->id == ver1->id &&  newGraph->verticesArray[tempEdge->vertex1ID]->id == ver2->id){
            
            ver1->deleteEdge(tempNode);
            ver2->deleteEdge(tempEdge->nodeInVertex1AdjList);
            
          }else{
            
            break;
          }
          
          tempEdge->nodeInVertex2AdjList = NULL;
          tempEdge->nodeInVertex1AdjList = NULL;
          newGraph->edgesArray->deleteNode(tempEdge->nodeInMainList);
          if ( mainNode != NULL && mainNode->val != NULL ){
            
            tempNode = mainNode->prevNode;
          }else{
            
            tempNode = NULL;
          }
        }
*/
        // TODO: change following to mainNode and Second Node
        if (edge->nodeInVertex1AdjList != NULL && edge->nodeInVertex2AdjList != NULL){
          
          // merge lists
          // second node will be deleted inside
          ver1->insertEdgesInList( mainNode, secondNode);
          
          // very important not to delete this edge
          // it will be deleted afterwards, if you delete here
          // in case this is first edge , iterator will move back
          //  and the way iterator code is written
          // in next call for getNextElement iteration will stop.
          
          ver1->deleteEdge(mainNode);
          edge->nodeInVertex1AdjList = NULL;
          edge->nodeInVertex2AdjList = NULL;
          newGraph->edgesArray->deleteNode( edge->nodeInMainList );
          delete edge;
//          if (newGraph->verticesArray[edge->vertex1ID]->id == ver1->id){
//            
//            edge->nodeInVertex2AdjList = NULL;
//          }else{
//          
//            edge->nodeInVertex1AdjList = NULL;
//          }
        }
        
        // contract his edge
        // while contracting just change the pointer at
        // the location of vertexID, that would change pointers
        // of all other edges
        // so verticesArray[VertexID] will give a different vertex now
        // And we don't need to change pointers in incoming edges to that
        // vertex
        // IMP: Identify vertices by id
        // Awesome ! Right ?
        newGraph->verticesArray[ver2->id] = ver1;

#ifdef DEBUG_ON
        std::cout << " after merging "<<  ver1->id << " " << ver2->id << std::endl;
        ver1->printAdjacencyList();
        LinkedList<Edge*> *adjacencyList = ver1->adjacencyList;
        for ( Edge *e = adjacencyList->beginIteration_debug(); e != NULL ;
             e = adjacencyList->getNextElement_debug()){
          
          std::cout << newGraph->verticesArray[e->vertex1ID]->id << "->" << newGraph->verticesArray[e->vertex2ID]->id << " W: " << e->getWeight() << " RW: " << e->getResidualWeight() << std::endl;
        }
#endif
        delete ver2;
        ver2 = NULL;
//      }
    }
    
  }
  
#ifdef DEBUG_ON
  cout << "******************************************" << endl;
  cout << " New graph edges" << endl;
  newGraph->printEdgesBoss();
  cout << "******************************************" << endl;
#endif
  // go through list again and remove all useless edges
  Edge* lastEdge = NULL;
  for ( Edge *edge = newGraph->edgesArray->beginIteration(); edge != NULL ;){
    
    if ((newGraph->verticesArray[edge->vertex1ID]->id == newGraph->verticesArray[edge->vertex2ID]->id)
        ||
        (edge->getWeight() == 0 && edge->getResidualWeight()==0 )
        ||
        // also checking if edge is similar to last undeleted edge
        // if it is we dont need that 
        (lastEdge != NULL && newGraph->verticesArray[lastEdge->vertex1ID]->id == newGraph->verticesArray[edge->vertex1ID]->id && newGraph->verticesArray[lastEdge->vertex2ID]->id == newGraph->verticesArray[edge->vertex2ID]->id)
        ) {
      
      // get current node
      Node<Edge*> *nodeToBeDeleted = newGraph->edgesArray->getCurrentNode();
      
      Node<Edge*> *node = NULL;
      
      // see which node is part of which list
      // wont go insde it for edges which are deleted from adjacency list
      // but only for edges which has weight 0;
      if ( edge->nodeInVertex1AdjList == NULL
          && edge->nodeInVertex2AdjList != NULL) {
        
        node = edge->nodeInVertex2AdjList;
        edge->nodeInVertex2AdjList = NULL;
      }else if ( edge->nodeInVertex2AdjList == NULL
                && edge->nodeInVertex1AdjList != NULL) {
        
        node = edge->nodeInVertex1AdjList;
        edge->nodeInVertex1AdjList = NULL;
      }else if (edge->nodeInVertex2AdjList != NULL
                && edge->nodeInVertex1AdjList != NULL){
        
        newGraph->verticesArray[edge->vertex1ID]->deleteEdge(edge->nodeInVertex1AdjList);
        newGraph->verticesArray[edge->vertex2ID]->deleteEdge(edge->nodeInVertex2AdjList);
        edge->nodeInVertex1AdjList = NULL;
        edge->nodeInVertex2AdjList = NULL;
      }
      
      newGraph->verticesArray[edge->vertex1ID]->printAdjacencyList();
      
      // ideally this shouldnt be required
      // as we are deleting unnecessary edges from adjacency list
      // in the above code
      // if it comes here , there is some error
      // this part of code should strictly be to remove edges
      // from main list as we cant delete it above
      if (node != NULL){
        
        newGraph->verticesArray[edge->vertex1ID]->deleteEdge(node);
      }
      
      newGraph->verticesArray[edge->vertex1ID]->printAdjacencyList();
      Edge* tempEdge = edge;
      // need to move ahead before deleting node from edgesArray
      edge = newGraph->edgesArray->getNextElement();
      
      // delete node from main list and vertex's adj list
      newGraph->edgesArray->deleteNode( nodeToBeDeleted );
      delete tempEdge;
    }else {
      
      lastEdge = edge; // last undeleted edge
      edge = newGraph->edgesArray->getNextElement();
    }
  }
  
  newGraph->setSource( source );
  newGraph->setSink ( sink );
  
#ifdef DEBUG_ON
  cout << "******************************************" << endl;
  cout << "After Deleting graph edges" << endl;
//  newGraph->printEdgesBoss();

//  for ( int  i = 0 ; i < newGraph->currentNumberOfVertices; ++i){
//    
//    cout << newGraph->verticesArray[i]->id << " " ;
//  }

    for ( int  i = 0 ; i < collectedVerticesList->size(); ++i){
      
      std::vector <int> *verWithSameBoss = collectedVerticesList->at(i);
      cout << verWithSameBoss->at(verWithSameBoss->size()-1) << " " ;
    }
  cout << endl <<"NUMBER OF VERTICES: " << collectedVerticesList->size() << endl;
//  LinkedList<Edge*> *adjacencyList = newGraph->verticesArray[0]->adjacencyList;
//  std::cout << std::endl << "********* Adj List 0  " << "*********" << std::endl;
//  for ( Edge *e = adjacencyList->beginIteration_debug(); e != NULL ;
//       e = adjacencyList->getNextElement_debug()){
//    
//    std::cout << newGraph->verticesArray[e->vertex1ID]->id << "->" << newGraph->verticesArray[e->vertex2ID]->id << " W: " << e->getWeight() << " RW: " << e->getResidualWeight() << std::endl;
//  }
//
  cout << endl;
  cout << "******************************************" << endl;
#endif
  
  /***************************APPROACH 3 ONLY ENDS*****************************/
  // clearing memory
  delete[] seen;
  seen = NULL;
  delete[] finishTime;
  delete verticesOrder;
  return newGraph;
}

/**
 *  Finished Time version of DFS
 */
void Graph::DFSTime(Vertex *ver, bool* seen ,int* discoveryTimeArray, int &time, std::vector<int> *verticesOrder){
  
  seen[ver->id] = true;
  //  Edge** array = ver->adjacencyList;
  //  int totalNeighbours = ver->numberOfEdges;
  //  for ( int i = 0; i < totalNeighbours; ++i) {
  //    Edge *edgeUnderQ = array[i] ;
  for ( Edge *edgeUnderQ = ver->adjacencyList->beginIteration(); edgeUnderQ != NULL ; edgeUnderQ = ver->adjacencyList->getNextElement()){
    
    if ( (!(seen[edgeUnderQ->vertex1ID] )) &&  !(edgeUnderQ->getResidualWeight() < EPSILON)){
      
      // using residual edge
      DFSTime(verticesArray[edgeUnderQ->vertex1ID], seen, discoveryTimeArray, time , verticesOrder);
      
    }else if( (!(seen[edgeUnderQ->vertex2ID] )) && !(edgeUnderQ->getWeight() < EPSILON) ){
      
      // using forward edge
      DFSTime(verticesArray[edgeUnderQ->vertex2ID], seen, discoveryTimeArray, time , verticesOrder);
    }
  }
  
  // we add it when it is DONE
  if ( discoveryTimeArray != NULL ){
    
    discoveryTimeArray[ver->id] = time; //will remove if not needed
    //in future
  }
  
  //  verticesOrder[time] = ver->id;
  verticesOrder->push_back( ver->id );
  ++time;
}


/**
 * Function to delete edge in the graph
 * Use only for planar graphs. (For images max number of edges is 8).
 */
void Graph::deleteEdge(int id1, int id2){
  
  
}


/**
 * Function to reverse edges
 */
void Graph::reverseEdges( Graph* graph ){
  
  //  for ( int i = 0 ; i < graph->currentNumberOfEdges ; ++i){
  for ( Edge *edge = graph->edgesArray->beginIteration(); edge != NULL ;
       edge = graph->edgesArray->getNextElement()){
    
    //    Edge* edge = graph->edgesArray[i];
    int tempId = edge->vertex2ID;
    edge->vertex2ID = edge->vertex1ID;
    edge->vertex1ID = tempId;
    
    // other way to reverse edges
//    WEIGHT_TYPE tempWeight = edge->getWeight();
//    edge->setWeight(edge->getResidualWeight());
//    edge->setResidualWeight(tempWeight);
  }
}

/**
 *  Vertex pair array if not initialized then initialize it
 *  and add pair to the array
 */
void Graph::addVertexPair ( int vertex1, int vertex2 ){
  
  if ( vertexPairArray == NULL ){
    
    vertexPairArray = new std::vector<int>();
  }
  
  vertexPairArray->push_back(vertex1);
  vertexPairArray->push_back(vertex2);
}

/**
 * THE function of the algorithm which calculates
 */
double Graph::countMinCuts (){
  
  double minCutsCount = 0;
  vertexPairPathCounts = new std::vector<double*>();
  
  // if following is the case that means no faces existed
  // which means only one min cut was there
  if (vertexPairArray == NULL) {
    
    return 1;
  }
  
  // go through each vertex pair
  // find number of paths in between those vertices
  for ( int i = 0 ; i < vertexPairArray->size(); i += 2 ){
    
    // create a visited edges array
    // mark each edge after visiting it
    // if that edge can be taken then only take that path
    // this is only for Testing
    // In DAG we should not have such case
    // TODO: Remove this after successful test
    bool *seen = new bool [currentNumberOfVertices];
    
    for (int j = 0 ; j < currentNumberOfVertices; ++j ){
      
      seen[j] = false;
    }
    
    vector<Edge*> path;
    double *countArray =  new double[currentNumberOfVertices];
    
    // intialize to -1
    for (int i = 0 ; i < currentNumberOfVertices ; ++i){
      
      countArray[i] = -1;
    }
    
    // TODO: Remove path after testing
    double temp = countPaths( vertexPairArray->at(i), vertexPairArray->at(i+1), seen, path, countArray);
    minCutsCount += temp;
    countArray[vertexPairArray->at(i)] = temp;
    
//    for (int i = 0 ; i < currentNumberOfVertices ; ++i){
//      
//      std::cout <<  i<<": "<< countArray[i] << std::endl;
//    }

    vertexPairPathCounts->push_back(countArray);
    delete [] seen;
//    delete [] countArray;// will b
  }
  
  return minCutsCount;
}

double Graph::countPaths (int source, int destination , bool *seen, std::vector<Edge*> &path, double *countArray){
  
  if (source == destination){

#ifdef DEBUG_ON
    std::cout << "Path START: "<< std::endl;
    
    for (int i = 0 ; i < path.size(); ++i){
      
      Edge *edge = path.at(i);
      std::cout << edge->vertex1ID << " -> " << edge->vertex2ID << std::endl;
      
      // see if we can print dual
      Edge *nonDualEdge = edge->getNonDualEdge();
      if (nonDualEdge != NULL) {
       
        std::cout << "Non dual edge: " <<  nonDualEdge->vertex1ID << " -> " << nonDualEdge->vertex2ID << std::endl;
      }
    }
    
    std::cout << "Path END: "<< std::endl;
#endif
    return 1;
  }
  
  // get adjacency list of source
  seen[source] = true;
  
  // check only forward edges here
  //  Edge** adjList = verticesArray[source]->adjacencyList;
  //  int numberOfEdgesInList = verticesArray[source]->numberOfEdges;
  //  for ( int i = 0 ; i < numberOfEdgesInList ; ++i){
  
  if (countArray[source] != -1){
    
    return countArray[source];
  }
  
  double pathCount = 0;
  for ( Edge *edge = verticesArray[source]->adjacencyList->beginIteration(); edge != NULL ; edge = verticesArray[source]->adjacencyList->getNextElement()){
    
    //    Edge *edge = adjList[i];
    
    // move ahead only if its a fwd edge and next vertex is not seen
    // seen check is necessary as we dont want to move backwards
    if ( edge->vertex1ID == source && !seen[edge->vertex2ID] ){
      
      //path.push_back(edge);
        double pathsFromV2 = countPaths(edge->vertex2ID, destination, seen, path, countArray);
        countArray[edge->vertex2ID] = pathsFromV2;
        pathCount = pathCount + pathsFromV2;
      //path.pop_back();
      seen[edge->vertex2ID] = false;
    }
  }
  
  return pathCount;
}


void Graph::printVertexPairArray(){

  // this means there were only two strongly connected components
  // or this is not a dual graph
  if (vertexPairArray == NULL){
    
    std::cerr << "Nothing to print here" << endl;
    return;
  }
  
  std::cout << "************* Vertex pairs **********" << std::endl;
  for (int i = 0; i < vertexPairArray->size()-1; i+=2 ){
  
    std::cout << " Pair " << ((i/2) + (1)) << ": " << vertexPairArray->at(i) << " " <<vertexPairArray->at(i+1) << std::endl;
  }
  std::cout << "************* Vertex pairs **********" << std::endl;
}


/**
 *   Sample a min cut
 */
/*std::set<int>* Graph::sampleAMinCut(int seed ){

//  bool *sampledMinCut = new bool[currentNumberOfVertices];
//  
//  for (int i = 0 ; i < currentNumberOfVertices; ++i ){
//    
//    sampledMinCut[i] = false;
//  }
//  
  // pick a vertex pair at random
  // consider path count array for that
  // start at end vertex pick a path from
  // destination to source,by going through adjacency list of selected
  // vertex, while picking vertex make sure its path count is not zero
  // beacuse in that case it will not be included in the path
  // aftter picking that vertex traverse back
  // store selected edges in an array.
  // go through that list in reverse order, see non dual edge, collect vertices
  // in a set
  if (seed == 0) {
    
    srand (time(NULL));
  }else{
    
    srand(seed);
  }

  int pairNumber  =  rand()%(vertexPairArray->size()/2);
  
  int sourceSelected = vertexPairArray->at( pairNumber );
  int destinationSelected = vertexPairArray->at( pairNumber+1 );
  double *pathCountArray = vertexPairPathCounts->at(pairNumber);
  std::vector<Edge*> pathChosen;
  
  for (int i = 0 ; i < currentNumberOfVertices; ++i ){
  
    if (destinationSelected == sourceSelected) {
      
      break;
    }
    
    Vertex *localDestination = verticesArray[destinationSelected];
    std::vector<Edge*> edgeLottery; // old_implementation
#ifdef WEIGHTED_PROBABILITY
    std::vector<double> edgeProbability;
    double lastProbab = 0.0;
#endif
    // go through its edge list, where it is the second endpoint
    // and put other points in an array as many times as they have count of
    // path, so that vertices with higher number of paths gets picked
    // with higher probability
    for ( Edge *edge = localDestination->adjacencyList->beginIteration(); edge != NULL ; edge = localDestination->adjacencyList->getNextElement()){

      // edge is incoming to this vertex
      if ( edge->vertex2ID == localDestination->id ){
      
//        if (edge->vertex1ID ==  sourceSelected){
//        
//          // no lottery here , pick this edge
//          pathChosen.push_back(edge);
//          break;
//        }
        
        // check if we can actually pick this vertex
        if (pathCountArray[edge->vertex1ID] != 0){
          
          // put it in lottery as many times as its count
          // so that it has chance proportional to number of paths to it
#ifdef WEIGHTED_PROBABILITY
            
            edgeLottery.push_back(edge);
            double probab = lastProbab + pathCountArray[edge->vertex1ID] / pathCountArray[edge->vertex2ID];
          
            if (probab > 1.0){
          
              cout << "Something went wrong: Probability can not be greater than 1" << endl  ;
              return NULL;
            }
          
            edgeProbability.push_back(probab);
            lastProbab = probab;
#else
          
           edgeLottery.push_back(edge);
#endif
        }
      }
    }
    
    // pick a vertex from vertex lottery
    int indexOfNextChosenEdge = 0;
    
#ifdef WEIGHTED_PROBABILITY
    
    double randomGen = rand()%10000; // allows 4 digit precision after decimal pt
    randomGen = randomGen/10000;
    
    for (int i = 0 ; i < edgeProbability.size() ; ++i) {
      
      if (isless(randomGen, edgeProbability.at(i))){
        
        indexOfNextChosenEdge = i;
        break;
      }
    }
    
#else
    indexOfNextChosenEdge = rand()%(edgeLottery.size());
#endif
    
    
    Edge *selectedEdge = edgeLottery.at(indexOfNextChosenEdge);
    pathChosen.push_back(selectedEdge);
    destinationSelected = selectedEdge->vertex1ID;// because we are going reverse
  }
  
  // now we have a chosen path. Put all the vertices of non dual edges
  // of these path edges in a set, vertices above the source is our
  // sampled min cut, i.e. we have a forward cut in the graph in which
  // SCCs are contracted. So all verticesArray[edge->vertex1ID]->id are vertices
  // in our min cut. How to get vertices in SCCs ?
  set<int> *minCut = new set<int>;
  
  for ( int j = 0 ; j < pathChosen.size() ; ++j ){
  
    Edge *edge = pathChosen.at(j);
    Edge* nonDualEdge = edge->getNonDualEdge();
    minCut->insert(nonDualEdge->vertex1ID);
  }
  
#ifdef DEBUG_ON
  std::cout << "********************************************" << std::endl;
  std::cout << "Sampled Cut:" << std::endl;
  std::set<int>::iterator it;
  for ( it = minCut->begin(); it != minCut->end(); ++it ){
    
    int bossId = *it;
    std::cout << bossId << " ";
  }
  std::cout << std::endl;
  std::cout << "********************************************" << std::endl;
#endif
  
  return minCut;
}*/


/**
 *   Sample a min cut
 */
std::set<int>* Graph::sampleAMinCut(int seed ){
  
  //  bool *sampledMinCut = new bool[currentNumberOfVertices];
  //
  //  for (int i = 0 ; i < currentNumberOfVertices; ++i ){
  //
  //    sampledMinCut[i] = false;
  //  }
  //
  // pick a vertex pair at random
  // consider path count array for that
  // start at end vertex pick a path from
  // destination to source,by going through adjacency list of selected
  // vertex, while picking vertex make sure its path count is not zero
  // beacuse in that case it will not be included in the path
  // aftter picking that vertex traverse back
  // store selected edges in an array.
  // go through that list in reverse order, see non dual edge, collect vertices
  // in a set
  if (seed == 0) {
    
    srand (time(NULL));
  }else{
    
    srand(seed);
  }
  
  int pairNumber  =  rand()%(vertexPairArray->size()/2);
  
  int sourceSelected = vertexPairArray->at( pairNumber );
  int destinationSelected = vertexPairArray->at( pairNumber+1 );
  double *pathCountArray = vertexPairPathCounts->at(pairNumber);
  std::vector<Edge*> pathChosen;
  
  for (int i = 0 ; i < currentNumberOfVertices; ++i ){
    
    if (destinationSelected == sourceSelected) {
      
      break;
    }
    
    Vertex *localSource = verticesArray[sourceSelected];
    std::vector<Edge*> edgeLottery; // old_implementation
#ifdef WEIGHTED_PROBABILITY
    std::vector<double> edgeProbability;
    double lastProbab = 0.0;
#endif
    // go through its edge list, where it is the second endpoint
    // and put other points in an array as many times as they have count of
    // path, so that vertices with higher number of paths gets picked
    // with higher probability
    for ( Edge *edge = localSource->adjacencyList->beginIteration(); edge != NULL ; edge = localSource->adjacencyList->getNextElement()){
      
      // edge is incoming to this vertex
      if ( edge->vertex1ID == localSource->id ){
        
        //        if (edge->vertex1ID ==  sourceSelected){
        //
        //          // no lottery here , pick this edge
        //          pathChosen.push_back(edge);
        //          break;
        //        }
        
        // check if we can actually pick this vertex
        if (pathCountArray[edge->vertex2ID] != 0){
          
          // put it in lottery as many times as its count
          // so that it has chance proportional to number of paths to it
#ifdef WEIGHTED_PROBABILITY
          
          edgeLottery.push_back(edge);
          double probab = lastProbab + (pathCountArray[edge->vertex2ID] / pathCountArray[edge->vertex1ID]);
          
//          if ( isgreater(probab , 1.0) == 1){
//            
//            cout << "Something went wrong: Probability can not be greater than 1: " << probab <<  endl  ;
//            return NULL;
//          }
          
          edgeProbability.push_back(probab);
          lastProbab = probab;
#else
          
          edgeLottery.push_back(edge);
#endif
        }
      }
    }
    
    // pick a vertex from vertex lottery
    int indexOfNextChosenEdge = 0;
    
#ifdef WEIGHTED_PROBABILITY
    
    double randomGen = rand()%10000; // allows 4 digit precision after decimal pt
    randomGen = randomGen/10000;
    
    for (int i = 0 ; i < edgeProbability.size() ; ++i) {
      
      if (isless(randomGen, edgeProbability.at(i))){
        
        indexOfNextChosenEdge = i;
        break;
      }
    }
    
#else
    indexOfNextChosenEdge = rand()%(edgeLottery.size());
#endif
    
    
    Edge *selectedEdge = edgeLottery.at(indexOfNextChosenEdge);
    pathChosen.push_back(selectedEdge);
    sourceSelected = selectedEdge->vertex2ID;// because we are going reverse
  }
  
  // now we have a chosen path. Put all the vertices of non dual edges
  // of these path edges in a set, vertices above the source is our
  // sampled min cut, i.e. we have a forward cut in the graph in which
  // SCCs are contracted. So all verticesArray[edge->vertex1ID]->id are vertices
  // in our min cut. How to get vertices in SCCs ?
  set<int> *minCut = new set<int>;
  
  for ( int j = 0 ; j < pathChosen.size() ; ++j ){
    
    Edge *edge = pathChosen.at(j);
    Edge* nonDualEdge = edge->getNonDualEdge();
    minCut->insert(nonDualEdge->vertex1ID);
  }
  
#ifdef DEBUG_ON
  std::cout << "********************************************" << std::endl;
  std::cout << "Sampled Cut:" << std::endl;
  std::set<int>::iterator it;
  for ( it = minCut->begin(); it != minCut->end(); ++it ){
    
    int bossId = *it;
    std::cout << bossId << " ";
  }
  std::cout << std::endl;
  std::cout << "********************************************" << std::endl;
#endif
  
  return minCut;
}


bool* Graph::getMaskingForSet(std::set<int> *minCut){

  if (collectedVerticesList == NULL) {
    
    cout << "INVALID OPERATION " << endl;
    return NULL;
  }
  
  bool *sampledMinCut = new bool[currentNumberOfVertices];
  
  for (int i = 0 ; i < currentNumberOfVertices; ++i ){
  
    sampledMinCut[i] = false;
  }
  
  std::set<int>::iterator it;
  
  for ( it = minCut->begin(); it != minCut->end(); ++it ){
    
    int bossId = verticesArray[*it]->boss;
    
    // search this boss in collectedVerticesList
    // assign all of the values to true;
    
    for (int i = 0 ; i < collectedVerticesList->size(); ++i ) {
      
      vector<int> *list = collectedVerticesList->at(i);
      int listBossID = list->at(list->size()-1);
      
      if (listBossID == bossId){
      
        for (int j = 0 ; j < list->size(); ++j) {
        
          sampledMinCut[list->at(j)] = true;
        }
        
        break;
      }
    }
  }
  
  return sampledMinCut;
}


// destructor
Graph::~Graph(){
  
  std::set<Vertex*> *uniqueVer = new std::set<Vertex*>();
  
  for (int i = currentNumberOfVertices-1 ; i >=0 ; --i) {
    
    Vertex *temp = verticesArray[i];
//    
//    for (int j = i-1 ; j >=0 ; --j) {
//        
//        if ( temp == verticesArray[j]) {
//          
//          verticesArray[j] = NULL;
//        }
//      }
//      
//      delete temp ;
//      temp = NULL;
    
    uniqueVer->insert(temp);
  }
  
  std::set<Vertex*>::iterator it;
  for (it = uniqueVer->begin(); it != uniqueVer->end(); ++it)
  {
    Vertex* ver = *it;
    delete ver;
  }
  
  delete uniqueVer;
  //  for (int i = 0 ; i < currentNumberOfEdges; ++i) {
  //    
  //    if (edgesArray[i] != NULL) {
  //     
  //      delete edgesArray[i];
  //    }
  //  }
  
  delete edgesArray;
  edgesArray = NULL;
  
  if (vertexPairArray != NULL){
    
    vertexPairArray->clear();
    delete vertexPairArray;
    vertexPairArray = NULL;
  }
  
  
  if (vertexPairPathCounts != NULL){
  
      for (int i = 0 ; i < vertexPairPathCounts->size(); ++i) {
        
        double *arr = vertexPairPathCounts->at(i);
        delete arr;
      }
  }
  
  delete vertexPairPathCounts;
  
  if (collectedVerticesList != NULL ){
  
    collectedVerticesList->clear();
    delete collectedVerticesList;
  }
}
