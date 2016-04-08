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


template <typename T>
std::string to_string(T value)
{
  std::ostringstream os ;
  os << value ;
  return os.str() ;
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
/******************************************************************
 INTEGRATION WITH SCHIMDTSEGMENTATION
 
 Citation :     [1] Efficient Planar Graph Cuts with Applications in Computer Vision.
 F. R. Schmidt, E. Töppe, D. Cremers,
 IEEE CVPR, Miami, Florida, June 2009
 
 ******************************************************************/


Graph *createGraph (int xResolution, int yResolution, int *luminanceArray){
  
  Graph* graph = new PlanarGraph (  xResolution * yResolution, xResolution * yResolution*4 );
  // connect all vertices now
  // edge weights are subtraction of luminance values
  int totalEdgesAdded = 0;
  for (int  i = 0 ; i < yResolution ; ++i ){
    
    for (int  j = 0 ; j < xResolution ; ++j ){
      
      //      int currentPixelIndex = i * xResolution + j * xResolution;
      //      int rightPixelIndex = i * xResolution + (j+1) * xResolution;
      //      int bottomPixelIndex = (i+1) * xResolution + j * xResolution;
      
      int currentPixelIndex = i * xResolution + j;
      int rightPixelIndex = i * xResolution + (j+1);
      int bottomPixelIndex = (i+1) * xResolution + j;
      
      // insert edges
      // NOTE: DO NOT check of left and up, as will be taken
      // care by following two cases.
      // if right side exists
      
      // storing edges in a sequence where it will be useful for
      // finding faces
      // add a check for bottom pixel
      if (i < yResolution-1 ){
        
        // if there is a pixel at bottom
        double bottomWeight  = abs(luminanceArray[currentPixelIndex] - luminanceArray[bottomPixelIndex]);
        WEIGHT_TYPE newWeight = weightFunction ( bottomWeight, xResolution, yResolution );
        
#ifdef USE_BIDIRECTIONAL_EDGES
        
        Edge* newEdge = graph->insertEdgeInGraph(currentPixelIndex, bottomPixelIndex,newWeight);
        newEdge = graph->insertEdgeInGraph(bottomPixelIndex, currentPixelIndex, newWeight);
        totalEdgesAdded = totalEdgesAdded+2;
#else
        
        Edge* newEdge = graph->insertEdgeInGraph(currentPixelIndex, bottomPixelIndex,newWeight);
        //        newEdge = graph->insertEdgeInGraph(bottomPixelIndex, currentPixelIndex, newWeight,true);
        newEdge->setResidualWeight(newWeight);
        ++totalEdgesAdded;
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
        
        totalEdgesAdded = totalEdgesAdded +2;
      }
      
      
#endif
      // add a check for right side
      if ( j < xResolution-1 ){
        
        double rightWeight  = abs(luminanceArray[currentPixelIndex] - luminanceArray[rightPixelIndex]);
        WEIGHT_TYPE newWeight = weightFunction ( rightWeight, xResolution, yResolution );
        
#ifdef USE_BIDIRECTIONAL_EDGES
        Edge* newEdge = graph->insertEdgeInGraph(currentPixelIndex, rightPixelIndex,newWeight);
        
        newEdge = graph->insertEdgeInGraph(rightPixelIndex, currentPixelIndex, newWeight);
        totalEdgesAdded = totalEdgesAdded + 2;
#else
        
        Edge* newEdge = graph->insertEdgeInGraph(currentPixelIndex, rightPixelIndex,newWeight);
        //newEdge = graph->insertEdgeInGraph(rightPixelIndex, currentPixelIndex, newWeight,true);
        
        newEdge->setResidualWeight(newWeight);
        ++totalEdgesAdded;
#endif
      }
    }
  }
  
  cout << "Totaledges = " << totalEdgesAdded << endl;
  return graph;
}

void calculateCuts( Graph *graph, int source, int sink, int xresol, int yresol, string name, unsigned char* rgbData = NULL ){
  
  // contract strongly connected components here
  Graph *graphDash = graph->findAndContractSCC( source, sink );
  ((PlanarGraph*)graphDash)->printEdges();
  ((PlanarGraph*)graphDash)->findFaces();
  ((PlanarGraph*)graphDash)->printFaces();
  ((PlanarGraph*)graphDash)->findAndMarkSTPath();
  Graph *dualGraph = ((PlanarGraph*)graphDash)->calculateDual();
  dualGraph->printVertexPairArray();
  
#ifdef DEBUG_ON_2
  std::cout << "************* DUAL GRAPH **********" << std::endl;
  dualGraph->printEdges();
  std::cout << "************* DUAL GRAPH **********" << std::endl;
#endif
  
  std::cout << "Number of min cuts: " <<   dualGraph->countMinCuts() << std::endl;
  
  for (int  i = 0 ; i < 5 ; ++i){
    
    std::set<int> *minCutSet = dualGraph->sampleAMinCut(i*10);
    bool *minCut = graph->getMaskingForSet(minCutSet);
    
    if (minCut == NULL) {
      
      cout << "Could not sample" << endl;
      return;
    }
    
    int numberOfVertices = graphDash->getNumberOfVertices();
    CutPlanar::ELabel *mask = new CutPlanar::ELabel[numberOfVertices];
    
    for (int i = 0; i < numberOfVertices ; ++i){
      
      if ( !minCut[i] ) {
        
        // belongs to source
        mask[i] = CutPlanar::LABEL_SOURCE;
        
      }else{
        
        // belongs to sink
        mask[i] = CutPlanar::LABEL_SINK;
      }
    }
    
    string str("sampledCut");
    str = str + to_string(i) + "_";
    string newName = str + name;
    unsigned char *rgbNew = SegMaskAndGreyDataToRGB( mask, rgbData, xresol, yresol);
    saveSimplePPM(rgbNew, xresol, yresol, newName);
    delete minCut;
    delete[] mask;
    delete[] rgbNew;
  }
  delete graphDash;
  delete dualGraph;
}


/**
 *   Function reference: OpenImageIO 1.7 Programmer Documentation
 */

#ifdef OPEN_IMAGE_IO_AVAILABLE
Graph* readImageAndCreateGraph (int &xres, int &yres){
  
  //  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/colorCircle.jpg");
//    ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/simmons2_small2.jpg");
//    ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/blackCircleSmall2.jpg");
  //  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/blackCircle.jpg");
  //  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/sample1.jpg");
  //  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/square.jpg");
  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/lena_color_small.png");
//  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/lena_bw_small2.jpg");
  
  if (!imageFile){
    
    return NULL;
  }
  
  const ImageSpec &imageSpecs = imageFile->spec();
  int xResolution = imageSpecs.width;
  int yResolution = imageSpecs.height;
  xres = xResolution;
  yres = yResolution;
  int channels = imageSpecs.nchannels;
  int pixelsArraySize = ( xResolution * yResolution * channels ) ;
  //  char *pixels = new char [ pixelsArraySize ];
  std::vector<unsigned char> pixels (pixelsArraySize);
  imageFile->read_image ( TypeDesc::UINT8, &pixels[0] );
  
  std::cout << "Pixel data: " << pixelsArraySize << std::endl;
  std::cout << "Pixel channels: " << imageSpecs.nchannels << std::endl;
  std::cout << "Resolution : " << xResolution << " * " << yResolution << std::endl;
  
  // create Graph here, (edges *4 , safer side they will be less than that)
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
    luminanceArray[ row*xResolution + column ] = luminance;
    ++column;
    
    // reset column to 0
    if ( column == xResolution ) {
      
      column = 0;
      ++row;
    }
    
    // unnecessary but added for testing
    if (row == yResolution){
      
      break;
    }
  }
  
  //  for (int i = 0; i < spec.channelnames.size() ; ++i){
  //
  //    std::cout << spec.channelnames[i] << std::endl;
  //  }
  
  Graph *graph =  createGraph(xResolution, yResolution, luminanceArray);
  // close and destroy the file
  imageFile->close ();
  ImageInput::destroy (imageFile);
  
  return graph;
}

void testCountingCuts (){
  
  // reading image and creating a planar graph
  int xRes;
  int yRes;
  Graph* graph = readImageAndCreateGraph(xRes,yRes);
  
  int source = 0;
  int sink = floor(xRes* yRes/2);
  //sink = 20*xRes + 19; // simmons2_small2
  sink = 20*xRes + 20;
  cout << "Sink: " << sink << endl;
  // find min cut value
  graph->getMinCut( source, sink);
  calculateCuts(graph, source, sink);
  delete graph;
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
  
  std::ifstream arq(getenv("GRAPH2"));
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
  int source = 0 ;
  //int sink = 7 ;
  int sink = numberOfVertices-1;
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

void countingCutsThroughSchmidt ( std::string picName, bool useCustomWeightFunction,bool useSchmidt = true,   int sinkRow = 0, int sinkColumn = 0, int sourceRow = 0 , int sourceColumn = 0 ){
  
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
  
  if (sinkRow == 0 && sinkColumn == 0) {
    
    sinkRow = (int)(yResolution/2);
    sinkColumn = (int)(xResolution/2);
  }
  
  int source[2] = {sourceRow,sourceColumn};
  int sink[2] = {sinkRow,sinkColumn};
  int sourceToWrite = sourceRow*xResolution + sourceColumn;
  int sinkToWrite = sinkRow*xResolution + sinkColumn;
  cout << "Source to write " << sourceToWrite << endl;
  cout << "Sink to write " << sinkToWrite << endl;
  
  // CutSegment initialize
  //perform segmentation task
  unsigned char *grey = RGBDataToGrey(rgbData, xResolution, yResolution);
  
  if ( !useSchmidt ){
  
    int *pic = new int[xResolution*yResolution];
    
    for ( int i=0; i < xResolution*yResolution; ++i){
      pic[i] = (int)grey[i];
    }
    
    Graph *g = createGraph(xResolution,yResolution,pic);
    bool* minCut = g->getMinCut (sourceToWrite,sinkToWrite);
    CutPlanar::ELabel *mask = new CutPlanar::ELabel[xResolution*yResolution];
    
    for (int i = 0; i < xResolution*yResolution; ++i){
      
      if ( minCut[i] ) {
        
        // belongs to source
        mask[i] = CutPlanar::LABEL_SOURCE;
        
      }else{

        // belongs to sink
        mask[i] = CutPlanar::LABEL_SINK;
      }
    }
    
    string name("result_");
    
    if (picName.rfind("/") == string::npos){
      
      name = name + picName.substr(0, picName.rfind(".ppm"));
    }else{
    
      name = name + picName.substr(picName.rfind("/")+1, picName.rfind(".ppm"));

    }
    
    name = name + to_string(sourceRow) + "," + to_string(sourceColumn) + "_" + to_string(sinkRow) + "," + to_string(sinkColumn);
    name.append(".ppm");

    unsigned char *rgbNew = SegMaskAndGreyDataToRGB( mask, rgbData, xResolution,yResolution );
    saveSimplePPM(rgbNew, xResolution, yResolution, name);
    cout << "\nSegmentation result written to " << name << "\n\n";
    delete[] minCut;
    delete[] mask;
    delete[] rgbNew;
    calculateCuts(g,sourceToWrite,sinkToWrite,xResolution,yResolution, name ,rgbData);
    delete[] pic;
    delete[] grey;
    delete[] rgbData;
    delete g;
    return;
  }
  
  CutSegment *sc;
  
  if (useCustomWeightFunction){
  
    sc = new CCCutSegment ( xResolution, yResolution);
  }else{
  
    sc = new CutSegment( xResolution, yResolution);
  }
  
  sc->setImageData(grey);
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
  int numberofHorizontalEdges =  (xResolution - 1)  *  yResolution;
  int numberOfVertices = xResolution*yResolution;
  
  Graph *planarGraph = new PlanarGraph(numberOfVertices, numberOfVertices*2);
  
  try{
    
    // Write graph to file
    ofstream myfile;
    myfile.open ("/Users/Gaurav/Documents/STudies/Capstone/SchimdtOutput.txt");
    myfile << xResolution*yResolution << " " <<  numberOfEdges << std::endl;
    int horizontalEdgeIndex = 0;
    int verticalEdgeIndex = 0;
    
//    for (int i = 0 ; i < numberOfEdges/2 ; ++i ){
    for (int i = 0 ; i < numberOfVertices ; ++i){
      
      // bottom edge
//      PlanarEdge *edge2 = &changed_Edges[i+numberofHorizontalEdges];
      PlanarEdge *edge2 = NULL;
      
      // dont need to add any vertical edge for last row
      if ( (verticalEdgeIndex+numberofHorizontalEdges) < numberOfEdges ) {
        
        edge2 = &changed_Edges[verticalEdgeIndex+numberofHorizontalEdges];
        ++verticalEdgeIndex;
        double actualWeight  = 1 - isless( edge2->getCapacity(), EPSILON);
        double actualRevWeight  = 1 - isless( edge2->getRevCapacity(), EPSILON);
        
        Edge *newEdge = planarGraph->insertEdgeInGraph((int)edge2->getTail()->vertexID, (int)edge2->getHead()->vertexID, actualWeight);
        newEdge->setResidualWeight(actualRevWeight);
      }
      // top edge
//      PlanarEdge *edge = &changed_Edges[i];
      PlanarEdge *edge = NULL;
      
      // dont need to add horizontal edge for last vertex of every row
      if ( (i+1)%xResolution != 0 ) {
        
        edge = &changed_Edges[horizontalEdgeIndex++];
        double actualWeight  = 1 - isless( edge->getCapacity(), EPSILON);
        double actualRevWeight  = 1 - isless( edge->getRevCapacity(), EPSILON);
        
        Edge * newEdge = planarGraph->insertEdgeInGraph((int)edge->getTail()->vertexID, (int)edge->getHead()->vertexID, actualWeight);
        newEdge->setResidualWeight(actualRevWeight);
      }
      
#ifdef WRITE_INTERMEDIATE_TO_FILE
      if (edge2 != NULL ){
        
        myfile << edge2->getTail()->vertexID << " " << edge2->getHead()->vertexID << " " << (1 - isless( edge2->getCapacity(), EPSILON)) << " " << (1 - isless( edge2->getRevCapacity(), EPSILON)) << "\n";
      }
      
      if (edge != NULL ){
        
        myfile << edge->getTail()->vertexID << " " << edge->getHead()->vertexID << " " << (1 - isless( edge->getCapacity(), EPSILON)) << " " << (1 - isless( edge->getRevCapacity(), EPSILON)) << "\n";
      }
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
    delete[] mask;
    exit(1);
  }
  
  delete[] mask;
  delete[] rgbNew;
  delete[] grey;
  delete[] rgbData;

  Graph *graphDash = planarGraph->findAndContractSCC( sourceToWrite, sinkToWrite );
  graphDash->printEdges();
  std::cout << "Source: " << ((PlanarGraph*)graphDash)->getSource() << " Sink: " << ((PlanarGraph*)graphDash)->getSink() << std::endl;
  
  ((PlanarGraph*)graphDash)->findFaces();
  ((PlanarGraph*)graphDash)->printFaces();
  ((PlanarGraph*)graphDash)->findAndMarkSTPath();
  Graph *dualGraph = ((PlanarGraph*)graphDash)->calculateDual();  
  dualGraph->printVertexPairArray();
  
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
  bool useSchmidt = true;

  int sinkRow = 0;
  int sinkColumn = 0;
  int sourceRow = 0;
  int sourceColumn = 0;
  
  if ( argc >= 3 ){
  
    useCustomWeightFunction = atoi(argv[2]);
  }
  
  if ( argc >= 4 ){
    
    useSchmidt = atoi(argv[3]);
  }


  if ( argc >= 5 ){
    
    sinkRow = atoi(argv[4]);
  }

  if ( argc >= 6 ){
    
    sinkColumn = atoi(argv[5]);
  }

  if ( argc >= 7 ){
    
    sourceRow = atoi(argv[6]);
  }

  if ( argc >= 8 ){
    
    sourceColumn = atoi(argv[7]);
  }

  cout << "Sourcerow: "<< sourceRow << " Sourcecolumn:" << sourceColumn << endl;
  cout << "SinkRow: "<< sinkRow << " SinkColumn:" << sinkColumn << endl;
  countingCutsThroughSchmidt( picName, useCustomWeightFunction, useSchmidt, sinkRow,sinkColumn,sourceRow,sourceColumn );
}


int main(int argc, const char * argv[]) {
  
    //testCountingCuts();
//    testPlanarGraphs();
//    testCountingPaths();
//    testLinkedList();
//    testCountingOnGraph();
  //  testCountingOnSchmidtGraph();
//    countingCutsThroughSchmidt("/Users/Gaurav/Documents/STudies/Capstone/lena_bw_Small2.ppm",false);
  // countingCutsThroughSchmidt("/Users/Gaurav/Documents/STudies/Capstone/simmons2_small.ppm", false, 38, 162 , 35,70);
  //  countingCutsThroughSchmidt("/Users/Gaurav/Documents/STudies/Capstone/simmons2_small.ppm", 55,70 );
  
  //  countingCutsThroughSchmidt("/Users/Gaurav/Documents/STudies/Capstone/simmons2_small2.ppm",6,21);
    countCutsWithArguments(argc, argv);
  //countingCutsThroughSchmidt("/Users/Gaurav/Documents/STudies/Capstone/enso1.ppm",false);
}