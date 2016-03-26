//
//  main.cpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 09/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#include "Constants.h"
#include <iostream>
#include <fstream>
#ifdef OPEN_IMAGE_IO_AVAILABLE
  #include <imageio.h>
#endif

#include <math.h>
#include <string>
#include <sstream>
#include <cstdlib>

#include "WeightFunction.h"
#include "Graph.h"
#include "Vertex.h"
#include "Edge.h"
#include "PlanarGraph.h"
#include "LinkedList.h"
#include "CutPlanar.h"
#include "CutPlanarDefs.h"
#include "CCCutSegment.h"

#ifdef OPEN_IMAGE_IO_AVAILABLE
  using namespace OpenImageIO;
#endif

using namespace std;

/**
 *   Function reference: OpenImageIO 1.7 Programmer Documentation
 */

#ifdef OPEN_IMAGE_IO_AVAILABLE
void readImageAndCreateGraph ( Graph *graph ){
  
  //  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/colorCircle.jpg");
    ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/blackCircleSmall.jpg");
//    ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/blackCircleSmall2.jpg");
  //  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/blackCircle.jpg");
  //  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/sample1.jpg");
  //  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/square.jpg");
//  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/lena_color_small.png");
//  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/lena_bw_small2.jpg");
  
  if (!imageFile){
    
    return;
  }
  
  const ImageSpec &imageSpecs = imageFile->spec();
  int xResolution = imageSpecs.width;
  int yResolution = imageSpecs.height;
  int channels = imageSpecs.nchannels;
  int pixelsArraySize = ( xResolution * yResolution * channels ) ;
  //  char *pixels = new char [ pixelsArraySize ];
  std::vector<unsigned char> pixels (pixelsArraySize);
  imageFile->read_image ( TypeDesc::UINT8, &pixels[0] );
  
  std::cout << "Pixel data: " << pixelsArraySize << std::endl;
  std::cout << "Pixel channels: " << imageSpecs.nchannels << std::endl;
  std::cout << "Resolution : " << xResolution << " * " << yResolution << std::endl;
  
  // create Graph here, (edges *4 , safer side they will be less than that)
  
  graph = new Graph (  xResolution * yResolution, xResolution * yResolution*4 );
  int row = 0 ;
  int column = 0 ;
  int *luminanceArray = new int[ xResolution * yResolution ];
  
  // TODO: add a check for number of channels
  for (int i = 0; i < pixelsArraySize ; i+=imageSpecs.nchannels){
    
    int R = (int)pixels[i];
    int G = (int)pixels[i+1];
    int B = (int)pixels[i+2];
    
    // relative luminance
    float luminance = (0.2126*R + 0.7152*G + 0.0722*B);
    
#ifdef DEBUG_ON
    std::cout <<  "R:" << R << " G:" << G << " B: " << B << std::endl;
    std::cout <<  row*yResolution + column << " " << luminance << std::endl;
#endif
    // put luminance values and make edges
    luminanceArray[ row*yResolution + column ] = luminance;
    ++column;
    
    // reset column to 0
    if ( column == yResolution ) {
      
      column = 0;
      ++row;
    }
    
    // unnecessary but added for testing
    if (row == xResolution){
      
      break;
    }
  }
  
  //  for (int i = 0; i < spec.channelnames.size() ; ++i){
  //
  //    std::cout << spec.channelnames[i] << std::endl;
  //  }
  
  
  // connect all vertices now
  // edge weights are subtraction of luminance values
  
  for (int  i = 0 ; i < xResolution ; ++i ){
    
    for (int  j = 0 ; j < yResolution ; ++j ){
      
      //      int currentPixelIndex = i * xResolution + j * xResolution;
      //      int rightPixelIndex = i * xResolution + (j+1) * xResolution;
      //      int bottomPixelIndex = (i+1) * xResolution + j * xResolution;
      
      int currentPixelIndex = i * yResolution + j;
      int rightPixelIndex = i * yResolution + (j+1);
      int bottomPixelIndex = (i+1) * yResolution + j;
      
      // insert edges
      // NOTE: DO NOT check of left and up, as will be taken
      // care by following two cases.
      // if right side exists
      
      // storing edges in a sequence where it will be useful for
      // finding faces
      // add a check for bottom pixel
      if (i != xResolution - 1 ){
        
        // if there is a pixel at bottom
        double bottomWeight  = abs(luminanceArray[currentPixelIndex] - luminanceArray[bottomPixelIndex]);
        WEIGHT_TYPE newWeight = weightFunction ( bottomWeight, xResolution, yResolution );
        
#ifdef USE_BIDIRECTIONAL_EDGES
        
        Edge* newEdge = graph->insertEdgeInGraph(currentPixelIndex, bottomPixelIndex,newWeight);
        newEdge = graph->insertEdgeInGraph(bottomPixelIndex, currentPixelIndex, newWeight);
#else
        
        Edge* newEdge = graph->insertEdgeInGraph(currentPixelIndex, bottomPixelIndex,newWeight,true);
        newEdge = graph->insertEdgeInGraph(bottomPixelIndex, currentPixelIndex, newWeight,true);
        
#endif
      }
      
#ifdef TRY_DIAGONAL_EDGES
      if (i != xResolution - 1 && j != yResolution-1 ){
        
        int diagonalPixelIndex = (i+1) * yResolution + j + 1;
        // if there is a pixel at bottom
        int weight  = abs(luminanceArray[currentPixelIndex] - luminanceArray[diagonalPixelIndex]);
        WEIGHT_TYPE newWeight = weightFunction ( weight, xResolution, yResolution );
        
        Edge* newEdge = graph->insertEdgeInGraph(currentPixelIndex, diagonalPixelIndex,newWeight);
        newEdge = graph->insertEdgeInGraph(diagonalPixelIndex, currentPixelIndex, newWeight);
      }
      
      
#endif
      // add a check for right side
      if ( j != yResolution-1 ){
        
        double rightWeight  = abs(luminanceArray[currentPixelIndex] - luminanceArray[rightPixelIndex]);
        WEIGHT_TYPE newWeight = weightFunction ( rightWeight, xResolution, yResolution );

#ifdef USE_BIDIRECTIONAL_EDGES
        Edge* newEdge = graph->insertEdgeInGraph(currentPixelIndex, rightPixelIndex,newWeight);
        
        newEdge = graph->insertEdgeInGraph(rightPixelIndex, currentPixelIndex, newWeight);
#else 
        
        Edge* newEdge = graph->insertEdgeInGraph(currentPixelIndex, rightPixelIndex,newWeight,true);
        newEdge = graph->insertEdgeInGraph(rightPixelIndex, currentPixelIndex, newWeight,true);
#endif
      }
    }
  }
  
#ifdef PRINT_GRAPH
  graph->printEdges();
#endif
  
  // close and destroy the file
  imageFile->close ();
  ImageInput::destroy (imageFile);
  
  int source = 2;
  //  int sink = 1000;
  int sink =floor(xResolution * yResolution/2);
  //  int sink = (xResolution * yResolution) - 1;
  
  // find min cut value
  graph->getMinCut( source, sink);
  
  // contract strongly connected components here
  Graph *graphDash = graph->findAndContractSCC( source, sink );
  ((PlanarGraph*)graphDash)->findFaces();
  ((PlanarGraph*)graphDash)->findAndMarkSTPath();
  Graph *dualGraph = ((PlanarGraph*)graphDash)->calculateDual();
  
#ifdef DEBUG_ON
  std::cout << "************* DUAL GRAPH **********" << std::endl;
  dualGraph->printEdges();
  std::cout << "************* DUAL GRAPH **********" << std::endl;
#endif
  
  std::cout << "Number of min cuts: " <<   dualGraph->countMinCuts() << std::endl;
  delete graphDash;
  delete dualGraph;
}

void testCountingCuts (){
  
  // reading image and creating a planar graph
  Graph *planarGraph;
  readImageAndCreateGraph( planarGraph );
  
  //   delete graph after done
  delete planarGraph;
  planarGraph = NULL;
}
#endif

void testPlanarGraphs(){
  
  std::ifstream arq(getenv("GRAPH"));
  std::cin.rdbuf(arq.rdbuf());
  
  int numberOfVertices = 0 ;
  std::cin >> numberOfVertices;
  PlanarGraph *planarGraph = new PlanarGraph(numberOfVertices,numberOfVertices);
  int id1;
  std::cin >> id1;
  
  while ( id1 != -1) {
    
    int id2;
    std::cin >> id2;
    planarGraph->insertEdgeInGraph(id1, id2, 1);
    std::cin >> id1;
  }
  
  planarGraph->findFaces();
  planarGraph->printFaces();
  planarGraph->setSource( 4 );
  planarGraph->setSink( 0 );
  planarGraph->findAndMarkSTPath();
  Graph *dual = planarGraph->calculateDual();
  dual->printEdges();
  delete dual;
  delete planarGraph;
}

void testCountingPaths (){
  
  std::ifstream arq(getenv("GRAPH"));
  std::cin.rdbuf(arq.rdbuf());
  
  int numberOfVertices = 0 ;
  std::cin >> numberOfVertices;
  Graph *planarGraph = new Graph(numberOfVertices,numberOfVertices);
  int id1;
  std::cin >> id1;
  
  while ( id1 != -1) {
    
    int id2;
    std::cin >> id2;
    planarGraph->insertEdgeInGraph(id1, id2, 1);
    std::cin >> id1;
  }
  
  planarGraph->addVertexPair(0, 4);
  std::cout << "Number of paths " << planarGraph->countMinCuts() << std::endl;
}

void testLinkedList(){
  
  LinkedList<int> *list = new LinkedList<int>();
  
  for (int i = 0 ; i < 10 ; ++i ){
    
    Node<int> *node1 = new Node<int>();
    node1->val = i;
    list->addNodeAtEnd( node1 );
  }
  
  LinkedList<int> *list2 = new LinkedList<int>();
  
  for (int i = 0 ; i < 4 ; ++i ){
    
    Node<int> *node1 = new Node<int>();
    node1->val = i*10;
    list2->addNodeAtEnd( node1 );
  }
  
  std::cout <<  " Begins "<<  std::endl;
  int val = list->beginIteration();
  
  for (int i = val; i < list->getListSize()-1 ; i = list->getNextElement()) {
    
    std::cout <<  i << "  " ;
  }
  
  std::cout <<  std::endl;
  std::cout <<  " ends "<<  std::endl;
  list->mergeLists( list->getHeadNode(), list2->getHeadNode() );
  list->deleteNode(list->getHeadNode());
  list->printList();
  list->printListRev();
  delete list;
  delete list2;
}

void testCountingOnSchmidtGraph(){
  
  std::ifstream arq(getenv("SchmidtGraph"));
  std::cin.rdbuf(arq.rdbuf());
  
  int numberOfVertices = 0 ;
  std::cin >> numberOfVertices;
  
  int numberOfEdges = 0 ;
  std::cin >> numberOfEdges;
  Graph *planarGraph = new Graph(numberOfVertices,numberOfEdges);
  int id1;
  std::cin >> id1;
  double weight;
  double revWeight;
  while ( id1 != -1) {
    
    int id2;
    std::cin >> id2;
    std::cin >> weight;
    std::cin >> revWeight;

    // I need to compare here with EPSILON as I dont want
    // very small values to be included
    double actualWeight  = 1 - isless( weight, EPSILON);
    double actualRevWeight  = 1 - isless( revWeight, EPSILON);
    
    Edge *edge =  planarGraph->insertEdgeInGraph(id1, id2, actualWeight);
    edge->setResidualWeight(actualRevWeight);
    std::cin >> id1;
  }
  
  // find min cut value
  int source, sink;
  
  std::cin >> source;
  std::cin >> sink;
  
  // min cut will be calculated from Schimdt segmentation
  //  planarGraph->getMinCut( source , sink );
  //  planarGraph->printEdges();
  Graph *graphDash = planarGraph->findAndContractSCC( source, sink );
  ((PlanarGraph*)graphDash)->findFaces();
  ((PlanarGraph*)graphDash)->printFaces();
  ((PlanarGraph*)graphDash)->findAndMarkSTPath();
  Graph *dualGraph = ((PlanarGraph*)graphDash)->calculateDual();
  
  std::cout << "************* DUAL GRAPH **********" << std::endl;
  dualGraph->printEdges();
  std::cout << "************* DUAL GRAPH **********" << std::endl;
  
  std::cout << "Number of min cuts: " << dualGraph->countMinCuts() << std::endl;
  delete graphDash;
  delete dualGraph;
  delete planarGraph;
}


void testCountingOnGraph(){
  
  std::ifstream arq(getenv("GRAPH5"));
  std::cin.rdbuf(arq.rdbuf());
  
  int numberOfVertices = 0 ;
  std::cin >> numberOfVertices;
  Graph *planarGraph = new Graph(numberOfVertices,numberOfVertices);
  int id1;
  std::cin >> id1;
  int weight;
  
  while ( id1 != -1) {
    
    int id2;
    std::cin >> id2;
    std::cin >> weight;
    planarGraph->insertEdgeInGraph(id1, id2, weight);
    std::cin >> id1;
  }
  
  // find min cut value
  int source = 0 , sink = 7 ;//sink = numberOfVertices-1;
  planarGraph->getMinCut( source , sink );
  planarGraph->printEdges();
  Graph *graphDash = planarGraph->findAndContractSCC( source, sink );
  ((PlanarGraph*)graphDash)->findFaces();
  ((PlanarGraph*)graphDash)->printFaces();
  ((PlanarGraph*)graphDash)->findAndMarkSTPath();
  Graph *dualGraph = ((PlanarGraph*)graphDash)->calculateDual();
  
  std::cout << "************* DUAL GRAPH **********" << std::endl;
  dualGraph->printEdges();
  std::cout << "************* DUAL GRAPH **********" << std::endl;
  
  std::cout << "Number of min cuts: " << dualGraph->countMinCuts() << std::endl;
  delete graphDash;
  delete dualGraph;
  delete planarGraph;
}

/******************************************************************
 INTEGRATION WITH SCHIMDTSEGMENTATION
 
 Citation :     [1] Efficient Planar Graph Cuts with Applications in Computer Vision.
 F. R. Schmidt, E. Töppe, D. Cremers,
 IEEE CVPR, Miami, Florida, June 2009
 
 ******************************************************************/

unsigned char *loadSimplePPM(int &w, int &h, const string &filename) {
  
  char line[1000];
  int depth = 0;
  unsigned char *rgb = 0;
  long lastpos;
  
  /*  streampos lastpos;
   ifstream ifs(filename.c_str(), ios_base::binary);*/
  
  FILE *fh = fopen(filename.c_str(), "rb");
  
  w = 0, h = 0;
  
  if (!fgets(line, 1000, fh))
    return NULL;
  
  if (strcmp(line, "P6\n")) {
    //cerr << filename << " is no PPM-Datei\n";
    return NULL;
  }
  
  while (!feof(fh)) {
    
    lastpos = ftell(fh);
    
    if (!fgets(line, 1000, fh))
      return NULL;
    
    if (line[0] == '#') {
      //      cout << "Comment: " << line;
    } else if (!w) {
      if (sscanf(line, "%d %d", &w, &h) < 2) {
        cerr << "error while reading the file " << filename;
        cerr << " expected width and height of image\n";
        return NULL;
      }
    } else if (!depth) {
      if (sscanf(line, "%d", &depth) < 1) {
        cerr << "error while reading the file " << filename;
        cerr << " expected color depth\n";
        return NULL;
      }
    } else {
      rgb = new unsigned char[w*h*3];
      fseek(fh, lastpos, SEEK_SET);
      if (fread(rgb, 1, w*h*3, fh) != size_t(w*h*3)) {
        fclose(fh);
        return NULL;
      }
      break;
    }
    
  }
  
  fclose(fh);
  
  return rgb;
  
}

bool saveSimplePPM(unsigned char *rgb, int w, int h, const string &filename) {
  
  ofstream fos(filename.c_str(), ios_base::binary);
  ostringstream ost;
  string s;
  
  if (!fos)
    return false;
  
  fos << "P6" << endl;
  
  ost << w << " " << h << endl;
  
  fos << ost.str();
  fos << "255" << endl;
  
  fos.write((const char*)rgb, w*h*3);
  
  fos.close();
  
  return true;
  
}


unsigned char *RGBDataToGrey(unsigned char *rgb, int w, int h) {
  
  unsigned char *pic = new unsigned char[w*h];
  int i;
  
  for (i=0; i<w*h; i++)
    pic[i] = rgb[i*3];
  
  return pic;
  
}


unsigned char *GreyDataToRGB(unsigned char *pic, int w, int h) {
  
  unsigned char *rgb = new unsigned char[w*h*3];
  int i;
  
  for (i=0; i<w*h; i++)
    rgb[i*3] = rgb[i*3+1] = rgb[i*3+2] = pic[i];
  
  return rgb;
  
}


unsigned char *SegMaskAndGreyDataToRGB(CutPlanar::ELabel *mask,
                                       unsigned char *pic,
                                       int w, int h) {
  
  unsigned char *rgb = new unsigned char[w*h*3];
  int i;
  
  for (i=0; i<w*h; i++) {
    
    rgb[i*3] = (mask[i]==CutPlanar::LABEL_SINK) ? (unsigned char)(pic[i*3] / 255.f * 200.f) : 255;
    rgb[i*3+1] = (unsigned char)(pic[i*3+1] / 255.f * 200.f);
    rgb[i*3+2] = (mask[i]==CutPlanar::LABEL_SINK) ? 255 : (unsigned char)(pic[i*3+2] / 255.f * 200.f);
    
  }
  
  return rgb;
  
}

void countingCutsThroughSchmidt ( std::string picName, bool useCustomWeightFunction, int sinkRow = 0, int sinkColumn = 0, int sourceRow = 0 , int sourceColumn = 0 ){
  
  unsigned char* rgbData = NULL;
  int xResolution, yResolution;
  //  string picName = "/Users/Gaurav/Documents/STudies/Capstone/blackCircleSmall.ppm";
  //  string picName = "/Users/Gaurav/Documents/STudies/Capstone/lena_bw_small2.ppm";
#ifdef OPEN_IMAGE_IO_AVAILABLE
  ImageInput *imageFile = ImageInput::open(picName);
  
  if (!imageFile){
    
    return;
  }
  
  const ImageSpec &imageSpecs = imageFile->spec();
  xResolution = imageSpecs.width;
  yResolution = imageSpecs.height;
  int channels = imageSpecs.nchannels;
  std::string format(imageFile->format_name());
  
  if (format.find("pnm") != -1) {
    
    // this is a ppm file follow different flow.
    rgbData = loadSimplePPM(xResolution, yResolution, picName);
  }else{
    
    int pixelsArraySize = ( xResolution * yResolution * channels ) ;
    //  char *pixels = new char [ pixelsArraySize ];
    std::vector<unsigned char> pixels (pixelsArraySize);
    imageFile->read_image ( TypeDesc::UINT8, &pixels[0] );
    
    std::cout << "Pixel data: " << pixelsArraySize << std::endl;
    std::cout << "Pixel channels: " << imageSpecs.nchannels << std::endl;
    std::cout << "Resolution : " << xResolution << " * " << yResolution << std::endl;
    // close and destroy the file
    imageFile->close ();
    ImageInput::destroy (imageFile);
    rgbData =  pixels.data();
  }
  
  // This means code will work only for ppm file
#else
  
  // this is a ppm file follow different flow.
  rgbData = loadSimplePPM(xResolution, yResolution, picName);
#endif
  
  
  if (rgbData == NULL) {
    
    std::cout << " ERROR could not read image!" << std::endl;
    exit(1);
  }
  
  // CutSegment initialize
  //perform segmentation task
  unsigned char *grey = RGBDataToGrey(rgbData, xResolution, yResolution);
  CutSegment *sc;
  
  if (useCustomWeightFunction){
  
    sc = new CCCutSegment ( xResolution, yResolution);
  }else{
  
    sc = new CutSegment( xResolution, yResolution);
  }
  
  sc->setImageData(grey);
  
  if (sinkRow == 0 && sinkColumn == 0) {
    
    sinkRow = yResolution/2;
    sinkColumn = xResolution/2;
  }
  
  int source[2] = {sourceRow,sourceColumn};
  int sink[2] = {sinkRow,sinkColumn};
  
  sc->setSourceSink(NULL, source, sink);
  cout << "Cut: " << sc->segment() << "\n";
  CutPlanar::ELabel *mask = new CutPlanar::ELabel[xResolution*yResolution];
  sc->getLabels(mask);
  
  //read out segmentation result and save to disk
  unsigned char *rgbNew = SegMaskAndGreyDataToRGB(mask, rgbData, xResolution,yResolution);
  saveSimplePPM(rgbNew, xResolution, yResolution, string("result.ppm"));
  cout << "\nSegmentation result written to result.ppm'\n\n";
  
//  for (int  i = 0 ; i < xResolution*yResolution ; ++i ){
//    
//    std::cout << mask[i] << " " ;
//  }
  
  PlanarEdge* changed_Edges = sc->getEdges();
  int numberOfEdges = sc->getNumberOfEdges();
  int numberofHorizontalEdges =  (yResolution - 1)  *  xResolution;
  int numberOfVertices = xResolution*yResolution;
  
  int sourceToWrite = sourceRow*xResolution + sourceColumn;
  int sinkToWrite = sinkRow*xResolution + sinkColumn;

  Graph *planarGraph = new Graph(numberOfVertices,numberOfVertices*2);
  
  try{
    
    // Write graph to file
    ofstream myfile;
    myfile.open ("/Users/Gaurav/Documents/STudies/Capstone/SchimdtOutput.txt");
    myfile << xResolution*yResolution << " " <<  numberOfEdges << std::endl;
    
    for (int i = 0 ; i < numberOfEdges/2 ; ++i ){
      
      // bottom edge
      PlanarEdge *edge2 = &changed_Edges[i+numberofHorizontalEdges];

      // top edge
      PlanarEdge *edge = &changed_Edges[i];
      
      double actualWeight  = 1 - isless( edge2->getCapacity(), EPSILON);
      double actualRevWeight  = 1 - isless( edge2->getRevCapacity(), EPSILON);

      Edge *newEdge = planarGraph->insertEdgeInGraph((int)edge2->getTail()->vertexID, (int)edge2->getHead()->vertexID, actualWeight);
      newEdge->setResidualWeight(actualRevWeight);

      actualWeight  = 1 - isless( edge->getCapacity(), EPSILON);
      actualRevWeight  = 1 - isless( edge->getRevCapacity(), EPSILON);
      
      newEdge = planarGraph->insertEdgeInGraph((int)edge->getTail()->vertexID, (int)edge->getHead()->vertexID, actualWeight);
      newEdge->setResidualWeight(actualRevWeight);
      
#ifdef WRITE_INTERMEDIATE_TO_FILE
      myfile << edge2->getTail()->vertexID << " " << edge2->getHead()->vertexID << " " << edge2->getCapacity() << " " << edge2->getRevCapacity() << "\n";
      
      myfile << edge->getTail()->vertexID << " " << edge->getHead()->vertexID << " " << edge->getCapacity() << " " << edge->getRevCapacity() << "\n";
#endif
    }

#ifdef WRITE_INTERMEDIATE_TO_FILE
    myfile << "-1\n" ;
    myfile << sourceToWrite << " " << sinkToWrite << "\n" ;
    myfile.close();
#endif
    delete sc;
  }catch (exception e){
    
    cout << "Something went wrong while writing file" << endl;
    delete rgbNew;
    delete mask;
    exit(1);
  }
  
  delete[] rgbNew;
  delete mask;

  Graph *graphDash = planarGraph->findAndContractSCC( sourceToWrite, sinkToWrite );
  ((PlanarGraph*)graphDash)->findFaces();
  ((PlanarGraph*)graphDash)->printFaces();
  ((PlanarGraph*)graphDash)->findAndMarkSTPath();
  Graph *dualGraph = ((PlanarGraph*)graphDash)->calculateDual();
  
  std::cout << "************* DUAL GRAPH **********" << std::endl;
  dualGraph->printEdges();
  std::cout << "************* DUAL GRAPH **********" << std::endl;
  
  std::cout << "Number of min cuts: " << dualGraph->countMinCuts() << std::endl;
  delete graphDash;
  delete dualGraph;
  delete planarGraph;

}

void countCutsWithArguments(int argc, const char * argv[]){

  //  countingCutsThroughSchmidt("/Users/Gaurav/Documents/STudies/Capstone/enso1.ppm",40,40);
  //  countingCutsThroughSchmidt("/Users/Gaurav/Documents/STudies/Capstone/blackCircleSmall.ppm",true);

  if (argc == 0) {
    
    cerr << "Invalid Input !" << endl;
  }
  
  string picName = argv[1];
  
  if (picName.rfind("ppm") == string::npos){
  
    cerr << "Image should be .ppm" << endl;
    exit(1);
  }
  
  bool useCustomWeightFunction = false;
  int sinkRow = 0;
  int sinkColumn = 0;
  int sourceRow = 0;
  int sourceColumn = 0;
  
  if ( argc >= 3 ){
  
    useCustomWeightFunction = atoi(argv[2]);
  }

  if ( argc >= 4 ){
    
    sinkRow = atoi(argv[3]);
  }

  if ( argc >= 5 ){
    
    sinkColumn = atoi(argv[4]);
  }

  if ( argc >= 6 ){
    
    sourceRow = atoi(argv[5]);
  }

  if ( argc >= 7 ){
    
    sourceRow = atoi(argv[6]);
  }

  countingCutsThroughSchmidt( picName, useCustomWeightFunction,sinkRow,sinkColumn,sourceRow,sourceColumn );
}


int main(int argc, const char * argv[]) {
  
//    testCountingCuts();
//    testPlanarGraphs();
//    testCountingPaths();
  //  testLinkedList();
//    testCountingOnGraph();
  //  testCountingOnSchmidtGraph();
  //  countingCutsThroughSchmidt("/Users/Gaurav/Documents/STudies/Capstone/simmons.ppm",false);
  // countingCutsThroughSchmidt("/Users/Gaurav/Documents/STudies/Capstone/simmons2_small.ppm", false, 38, 162 , 35,70);
  //  countingCutsThroughSchmidt("/Users/Gaurav/Documents/STudies/Capstone/simmons2_small.ppm", 55,70 );
  
  //  countingCutsThroughSchmidt("/Users/Gaurav/Documents/STudies/Capstone/simmons2_small2.ppm",6,21);
  countCutsWithArguments(argc, argv);  
}