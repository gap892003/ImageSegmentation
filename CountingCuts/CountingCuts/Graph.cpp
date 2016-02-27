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
#include <iostream>
#include <list>
#include <vector>
using namespace std;

// Constructor
Graph::Graph (int numberOfVertices, int numberOfEdges):
totalVertices(numberOfVertices),totalEdges(numberOfEdges){
  
  verticesArray = new Vertex* [numberOfVertices];
  edgesArray = new Edge* [numberOfEdges];
  currentNumberOfEdges = 0;
  currentNumberOfVertices = 0;
  
  for ( int i = 0 ; i < numberOfVertices ; ++i){
    
    verticesArray[i] = NULL; // setting all to null
  }
  
  lastVertex = numberOfVertices-1;
  visitedEdges = new int[numberOfEdges];
  
  for ( int i = 0 ; i < numberOfEdges ; ++i){
    
    visitedEdges[i] = 0;
  }
  
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
  edge->vertex1ID = idOfVertex1;
  edge->vertex2ID = idOfvertex2;
  edge->setWeight( weight );
  
  if ( verticesArray[idOfvertex2] == NULL ) {
    
    insertVertexInGraph(idOfvertex2);
  }
  
  if ( verticesArray[idOfVertex1] == NULL ) {
    insertVertexInGraph(idOfVertex1);
  }
  
  verticesArray[idOfVertex1]->addEdge(edge);
  
  /********** Analyse this *******************************/
  if (!oneWay){
  
    verticesArray[idOfvertex2]->addEdge(edge);
  
  }
  
  edgesArray[currentNumberOfEdges++] = edge;
  return edge;
}

/**
 *  Function used for debugging to print out edges in graph
 */
void Graph::printEdges(){
  
  for (int i = 0 ; i < currentNumberOfEdges; ++i) {
    
    cout << edgesArray[i]->vertex1ID << " " << edgesArray[i]->vertex2ID << "  " << " W: "<< edgesArray[i]->getWeight() << "  RW: " << edgesArray[i]->getResidualWeight() << endl  ;
  }
  
  cout <<endl;
}

/**
 *  Function used for debugging to print out edges passed
 */
void Graph::printEdgesArray(Edge** array, int count){
  
  for (int i = 0 ; i < count ; ++i) {
    
    cout << array[i]->vertex1ID << " " << array[i]->vertex2ID << " " <<array[i]->getWeight() << endl  ;
  }
}

Edge** Graph::DFS( int& edgesInPath, int startIndex,  int endIndex){
  
  Edge** path = new Edge*[totalVertices];
  bool *seen = new bool [totalVertices];
  
  for (int i = 0 ; i < totalVertices ; ++i){
    
    seen[i] = false;
  }
  
  DFSRec(verticesArray[startIndex], seen, path, edgesInPath, endIndex);
  delete seen;
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
  
  Edge** array = ver->adjacencyList;
  int totalNeighbours = ver->numberOfEdges;
  
  // finding bottleneck
  for ( int i = 0; i < totalNeighbours; ++i) {
    
    Edge *edgeUnderQ = array[i] ;
    
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
    
    if ( (!(seen[edgeUnderQ->vertex1ID] )) &&  edgeUnderQ->getResidualWeight() !=0  ){
      
      // using residual edge
      if ( DFSRec(verticesArray[array[i]->vertex1ID], seen, path, edgesInPath,endIndex) ){
        
        if (path != NULL){
          
          path[edgesInPath++] = edgeUnderQ;
        }
        return true;
      }
      
    }else if( (!(seen[array[i]->vertex2ID] )) && edgeUnderQ->getWeight() != 0 ){
      
      // using forward edge
      if (DFSRec(verticesArray[array[i]->vertex2ID], seen, path, edgesInPath,endIndex)){
        
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
    
    
    /***********REMOVE THIS************************************/
    cout << "Path is: " << endl;
    for (int i = 0 ; i < edgesInPath ; ++i){
      
      cout << path[i]->vertex1ID << "->" << path[i]->vertex2ID << " "<<  " W: "
      << path[i]->getWeight() << " RW: " << path[i]->getResidualWeight() << endl;
    }
    /***********REMOVE THIS************************************/
    
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
      
      if ( weightToCompare < bottleNeckWeight || bottleNeckWeight == -1) {
        
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
    delete path;
    path = NULL;
    edgesInPath = 0;
    path = DFS(edgesInPath,s,t);
  }
  
  if (path != NULL) {
    
    delete path;
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
Vertex** Graph::getMinCut(int s, int t){
  
  WEIGHT_TYPE maxFlow =  findMaxFlow(s, t);
  std::cout << "Max flow is : "<< maxFlow << std::endl;
  bool *seen = new bool [totalVertices];
  Vertex** minCut = new Vertex* [totalVertices];
  
  for (int i = 0 ; i < totalVertices ; ++i){
    
    seen[i] = false;
  }
  
  int count = 0; // dummy wont be used
  DFSRec(verticesArray[s], seen, NULL, count, t);
  
  // go through all seen vertices and add them to
  // min cut Vertices array;
  count = 0;
  std::cout << "Min cut is { " ;
  for (int i = 0 ; i < totalVertices ; ++i){
    
    if (seen[i]){
      minCut[count++] = verticesArray[i];
      std::cout << i << ", " ;
    }
  }
  
  std::cout << " }"<< std::endl;
  delete seen;
  seen = NULL;
  return minCut;
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
  delete seen;
  return path;
}

/**
 * Function to be used after max flow is found
 * It returns the matching vertex from other sets.
 */
int Graph::getMatch(int i){
  
  Vertex *ver = verticesArray[i] ;
  Edge** array = ver->adjacencyList;
  int totalNeighbours = ver->numberOfEdges;
  
  // go through all neighbors to find whcih edge was finalized
  // to be part of maximum flow
  // We identify it by forward weight zero
  for ( int j = 0; j < totalNeighbours; ++j) {
    
    if ((array[j]->vertex1ID == i) && array[j]->getWeight() == 0){
      
      return array[j]->vertex2ID;
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
    
    Edge** array = vertexToVisit->adjacencyList;
    int totalNeighbours = vertexToVisit->numberOfEdges;
    
    for ( int i = 0; i < totalNeighbours; ++i) {
      
      // add vertices to queue if they are reachable
      // via residual or
      Edge *edgeUnderQ = array[i] ;
      
      if ( (!(seen[edgeUnderQ->vertex1ID] )) &&  edgeUnderQ->getResidualWeight() !=0  ){
      
        if ( traversedEdges[verticesArray[edgeUnderQ->vertex1ID]->id] != NULL ){
          
          cout << " Vertex already there not adding to queue" << endl;
          continue;
        }

        // using backward edge
        queue.push_back( verticesArray[edgeUnderQ->vertex1ID] );
        traversedEdges[verticesArray[edgeUnderQ->vertex1ID]->id] = edgeUnderQ;
        
      }else if( (!(seen[array[i]->vertex2ID] )) && edgeUnderQ->getWeight() != 0 ){

        if ( traversedEdges[verticesArray[edgeUnderQ->vertex2ID]->id] != NULL ){
          
          cout << " Vertex already there not adding to queue" << endl;
          continue;
        }

        // using forward edge
        queue.push_back( verticesArray[edgeUnderQ->vertex2ID] );
        traversedEdges[verticesArray[edgeUnderQ->vertex2ID]->id] = edgeUnderQ;
      }
    }
  }
  
  bool reached = false;
  int lastVertexId = t;
  
  while (!reached) {
    
    Edge *edgeInPath = traversedEdges[lastVertexId];
    path[edgesInPath++] = edgeInPath;
    
    // this means it was a forward edge(not residual edge) last vertex will be
    // vertex 1, vice versa
    if ( edgeInPath->vertex2ID == lastVertexId ){
      
      lastVertexId = edgeInPath->vertex1ID;
      // we will examine how we reached vertex 1
      // i.e. using which edge
    }else{
      
      // this means residual edge was used to reach vertex
      lastVertexId = edgeInPath->vertex2ID;
    }
    
    if (lastVertexId == s){
      
      reached = true;
    }
  }
  
  delete seen;
  seen = NULL;
  delete traversedEdges;
  traversedEdges = NULL;
  return path;
}

/**
 *  Function finds Strongly connected Components, stores them 
 *  then process ( contracts them )
 */
Graph* Graph::findAndContractSCC ( int source, int sink ){
  
  // Run DFS On residual graph
  bool *seen = new bool [totalVertices];
  int *finishTime =  new int[totalVertices];
  int *verticesOrder = new int[totalVertices];// this array is to store sorted
                                              // order of vertices during DFS
  int numberOfSCCFound = 0;
  vector<int*> collectedVerticesList;
  
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
    
    if (! seen[verticesOrder[i]]){

      int dummy = 0;
      
      // using int array here
      // dont want to use vector as , it can hamper
      // time complexity
      int *verticesCollected = new int[totalVertices];

      DFSTime( verticesArray[verticesOrder[i]], seen, NULL, dummy , verticesCollected);
      
      ++numberOfSCCFound;
      // do something with collcted vertices
      // right now just printing
      cout << "Strongly connected Component: " ;

      // settting boss for each component
      // this will be useful while creating new graph
      for (int k = 0 ; k < dummy ; ++k){
        
        cout << verticesCollected[k] << " ";
        verticesArray[verticesCollected[k]]->boss = verticesCollected[0];
      }
      
      cout << endl;
      collectedVerticesList.push_back(verticesCollected);
//      delete verticesCollected;
    }
  }

  // reverse again to make it normal
  Graph::reverseEdges(this);

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
  
  
  // create new replica of the graph with exact same bosses
  // vertices, then use edge contraction
  // IMP: while creating new graph, dont worry about residual graphs
  PlanarGraph *newGraph = new PlanarGraph (currentNumberOfVertices, currentNumberOfEdges*2);

  //currentNumberOfEdges*2 safer side to account for residual edges
  
  for ( int i = 0 ; i < currentNumberOfEdges ; ++i ){
    
    Edge* edge = edgesArray[i];
    
    //*****************Analyse this
    // We need double way adjacency list to find out faces of the graph
    Edge* newEdge = newGraph->insertEdgeInGraph( edge->vertex1ID, edge->vertex2ID, edge->getWeight() );
    newGraph->verticesArray[newEdge->vertex1ID]->boss = verticesArray[edge->vertex1ID]->boss;
    newGraph->verticesArray[newEdge->vertex2ID]->boss = verticesArray[edge->vertex2ID]->boss;
    
    // if there is residual edge add actual edge
    if (edge->getResidualWeight() > 0){
    
      Edge* newEdgeRes = newGraph->insertEdgeInGraph( edge->vertex2ID, edge->vertex1ID, edge->getResidualWeight());
      newGraph->verticesArray[newEdge->vertex1ID]->boss = verticesArray[edge->vertex2ID]->boss;
      newGraph->verticesArray[newEdge->vertex2ID]->boss = verticesArray[edge->vertex1ID]->boss;

    }
  }
  
  cout << "******************************************" << endl;
  cout << " New graph edges" << endl;
  newGraph->printEdges();
  cout << "******************************************" << endl;
  
  // Do edge contraction here,
  // go through graph one vertex at a time and then contract an
  // edge
  // this is poor version of this I will imporve on it later
  // IMP: Dont care about residual edges here

  // iterate through each all strongly connected components collected above

//  for ( int i = 0 ; i < newGraph->currentNumberOfEdges ; ++i ){
  
  for ( int i = 0 ; i < collectedVerticesList.size() ; ++i ){
    
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
  
  // clearing memory
  delete seen;
  seen = NULL;
  delete finishTime;
  delete verticesOrder;
  return newGraph;
}

/**
 *  Finished Time version of DFS
 */
void Graph::DFSTime(Vertex *ver, bool* seen ,int* discoveryTimeArray, int &time, int *verticesOrder){

  seen[ver->id] = true;
  Edge** array = ver->adjacencyList;
  int totalNeighbours = ver->numberOfEdges;
  
  for ( int i = 0; i < totalNeighbours; ++i) {
    
    Edge *edgeUnderQ = array[i] ;

    if ( (!(seen[edgeUnderQ->vertex1ID] )) &&  edgeUnderQ->getResidualWeight() !=0  ){
      
      // using residual edge
      DFSTime(verticesArray[array[i]->vertex1ID], seen, discoveryTimeArray, time , verticesOrder);
      
    }else if( (!(seen[array[i]->vertex2ID] )) && edgeUnderQ->getWeight() != 0 ){
      
      // using forward edge
      DFSTime(verticesArray[array[i]->vertex2ID], seen, discoveryTimeArray, time , verticesOrder);
    }
  }
  
  // we add it when it is DONE
  if ( discoveryTimeArray != NULL ){

    discoveryTimeArray[ver->id] = time; //will remove if not needed
    //in future
  }
  
  verticesOrder[time] = ver->id;
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

  for ( int i = 0 ; i < graph->currentNumberOfEdges ; ++i){
    
    Edge* edge = graph->edgesArray[i];
    int tempId = edge->vertex2ID;
    edge->vertex2ID = edge->vertex1ID;
    edge->vertex1ID = tempId;
  }
}

// destructor
Graph::~Graph(){
  
  for (int i = 0 ; i < currentNumberOfVertices; ++i) {
    
    delete verticesArray[i];
  }
  
  for (int i = 0 ; i < currentNumberOfEdges; ++i) {
    
    delete edgesArray[i];
  }
}
