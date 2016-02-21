//
//  main.cpp
//  CountingCuts
//
//  Created by Gaurav Joshi on 09/02/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#include <iostream>
#include <imageio.h>
#include "Graph.h"
#include "Vertex.h"
#include "Edge.h"
#include <math.h>
#include "Constants.h"

using namespace OpenImageIO;

/**
 * Weight function
 */
WEIGHT_TYPE weightFunction ( int &luminance ){

//  WEIGHT_TYPE newValue = 256 - luminance;
//  WEIGHT_TYPE newValue = 100000000 - luminance*1000;
  
  //  should be inverse
//  WEIGHT_TYPE newValue = pow (2, luminance+1);
//  WEIGHT_TYPE newValue = exp (luminance+1);
  
  double temp  = ( 255 / log(luminance+2) );// x+2 as Dont want to deal with zeroes
  WEIGHT_TYPE newValue = temp + 1;
  return newValue;
}

/**
 *   Function reference: OpenImageIO 1.7 Programmer Documentation
 */
void readImageAndCreateGraph ( Graph *graph ){
  
//  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/colorCircle.jpg");
//  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/blackCircleSmall.jpg");
//  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/blackCircle.jpg");
//  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/sample1.jpg");
  ImageInput *imageFile = ImageInput::open("/Users/Gaurav/Documents/STudies/Capstone/square.jpg");

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
    
#ifdef PRINT_RGB
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
      
      // add a check for right side
      if ( j != yResolution-1 ){
        
        int rightWeight  = abs(luminanceArray[currentPixelIndex] - luminanceArray[rightPixelIndex]);
        WEIGHT_TYPE newWeight = weightFunction ( rightWeight );
        
        graph->insertEdgeInGraph(currentPixelIndex, rightPixelIndex,newWeight);
        graph->insertEdgeInGraph(rightPixelIndex, currentPixelIndex, newWeight);
      }
      
      
      // add a check for bottom pixel
      if (i != xResolution - 1 ){
      
        // if there is a pixel at bottom
        int bottomWeight  = abs(luminanceArray[currentPixelIndex] - luminanceArray[bottomPixelIndex]);
        WEIGHT_TYPE newWeight = weightFunction ( bottomWeight );
        
        graph->insertEdgeInGraph(currentPixelIndex, bottomPixelIndex,newWeight);
        graph->insertEdgeInGraph(bottomPixelIndex, currentPixelIndex, newWeight);

      }
    }
  }
 
#ifdef PRINT_GRAPH
  graph->printEdges();
#endif

    // close and destroy the file
  imageFile->close ();
  ImageInput::destroy (imageFile);
  
  // find min cut value
  graph->getMinCut( 2, floor(xResolution * yResolution/2));
  
  // contract strongly connected components here
  Graph *graphDash = graph->findAndContractSCC();
  delete graphDash;
}

int main(int argc, const char * argv[]) {
  
  // reading image and creating a planar graph
  Graph *planarGraph;
  readImageAndCreateGraph( planarGraph );

  // delete graph after done
  delete planarGraph;
  planarGraph = NULL;
  return 0;
}