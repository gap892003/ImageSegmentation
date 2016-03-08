//
//  LinkedList.h
//  CountingCuts
//
//  Created by Gaurav Joshi on 04/03/16.
//  Copyright © 2016 Gj. All rights reserved.
//

#ifndef LinkedList_h
#define LinkedList_h

#include <stdio.h>
#include <iostream>
#include "Constants.h"

template <class T>
struct Node{
  
  T val;
  Node* nextNode;
  Node* prevNode;
};

template <class T>
class LinkedList{

protected:
  Node<T> *headNode;
  Node<T> *tailNode;
  Node<T> *iteratorNode;
  Node<T> *debugIteratorNode;
  int sizeOfList; // not realiable , part of design , just for debugging
  
public:
  
  LinkedList(){
  
    headNode = NULL;
    tailNode = NULL;
    iteratorNode = NULL;
    debugIteratorNode = NULL;
  }
  
  Node<T>* addValue ( T val ){
  
    Node<T> *node = new Node<T>();
    node->val = val;
    addNodeAtEnd(node);
    return node;
  }
  
  void addNodeAtEnd( Node<T> * nodeToAdd ){
  
    if (headNode == NULL) {
      
      headNode = nodeToAdd;
      headNode->nextNode = NULL;
      headNode->prevNode = NULL;
      tailNode = nodeToAdd;
    }else{
      
      nodeToAdd->nextNode = headNode;
      headNode->prevNode = nodeToAdd;
      
      nodeToAdd->prevNode = tailNode;
      tailNode->nextNode = nodeToAdd;
      
      tailNode = nodeToAdd;
    }
    
    ++sizeOfList;
  }
  
  void addNodeAtStart( Node<T> *nodeToAdd ){
  
    if (headNode == NULL) {
      
      headNode = nodeToAdd;
      headNode->nextNode = NULL;
      headNode->prevNode = NULL;
      tailNode = nodeToAdd;
    }else{
      
      nodeToAdd->nextNode = headNode;
      headNode->prevNode = nodeToAdd;
      
      nodeToAdd->prevNode = tailNode;
      tailNode->nextNode = nodeToAdd;
      
      headNode = nodeToAdd;
    }
    
    ++sizeOfList;
  }
  
  void addListAfterNode( Node<T> *nodeInList1, LinkedList *listToAdd){
  
    if (headNode == NULL) {
      
      headNode = listToAdd->getHeadNode();
      tailNode = listToAdd->getTailNode();
    }else{
      
      Node<T> *temp = nodeInList1->nextNode;
      nodeInList1->nextNode = listToAdd->getHeadNode();
      Node<T> *lastNode = listToAdd->getTailNode();
      lastNode->nextNode = temp;
    }
    
    sizeOfList += listToAdd->getListSize();
  }

  // both node will have same value here but different next and prev
  void mergeLists(Node<T> *nodeInList1, Node<T> *nodeInList2){
  
    // dont think head and tail will need to be changed
    // we dont want nodeInList2 in our linkedList
    Node<T> *temp = nodeInList1->nextNode;
    
    if (nodeInList2->nextNode != NULL){
      
      nodeInList1->nextNode = nodeInList2->nextNode;
      nodeInList2->nextNode->prevNode = nodeInList1;
    }
    
    if (nodeInList2->prevNode != NULL){
      
      nodeInList2->prevNode->nextNode = temp;
      temp->prevNode = nodeInList2->prevNode;
    }
    
    nodeInList2->nextNode = NULL;
    nodeInList2->prevNode = NULL;
    delete nodeInList2;
    nodeInList2 = NULL;
  }
  
  inline Node<T> * getHeadNode(){
  
    return headNode;
  }
  
  inline Node<T>* getTailNode(){
    
    return tailNode;
  }
  
  inline int getListSize(){
    
    return sizeOfList;
  }
  
  void deleteNode( Node<T> *node){
    
    if ( node  == NULL) return;
    
    if (node->prevNode != NULL && node->nextNode != NULL){

      if ( headNode == node ){
        
        headNode = node->nextNode;
      }
      
      if ( tailNode == node ){
        
        tailNode = node->nextNode;
      }

      node->prevNode->nextNode = node->nextNode;
      node->nextNode->prevNode = node->prevNode;
      
      // check if iteration is going on , delete that node
      // and move pointer back , so that iterator is not
      if (iteratorNode != NULL && iteratorNode == node) {
        iteratorNode = node->prevNode;
      }
    }else {
    
      // this means its only node in list
      headNode = NULL;
      tailNode = NULL;
      iteratorNode = NULL;
    }
    
    --sizeOfList;
    node->nextNode = NULL;
    node->prevNode = NULL;
    delete node;
    node = NULL;
  }

  void printList(){

#ifdef DEBUG_ON
    Node<T> *iNode = headNode;
    
    do {
      
      std::cout<< iNode->val << " ";
      iNode = iNode->nextNode;
    } while (iNode != tailNode);
    
    std::cout<< tailNode->val << " ";
    std::cout<< std::endl;
#endif
    
  }
  
  
  void printListRev(){
    
#ifdef DEBUG_ON
    Node<T> *iNode = tailNode;
    
    do {
      
      std::cout<< iNode->val << " ";
      iNode = iNode->prevNode;
    } while (iNode != headNode);
    
    std::cout<< headNode->val << " ";
    std::cout<< std::endl;
#endif
    
  }

  void clearList(){
  
    Node<T> *iNode = headNode;
    headNode = NULL;
    
    while (iNode != tailNode) {
      
      Node<T> *temp = iNode;
      iNode = iNode->nextNode;
      if (iNode != NULL){
      
        delete temp;
        temp = NULL;
      }
    }
    
    if (tailNode != NULL){
      
      delete tailNode;
      tailNode = NULL;
    }
  }

  /**
   * Begin iteration
   */
  T beginIteration(){
    
    if (headNode != NULL) {
      
      iteratorNode = headNode;
      return iteratorNode->val;
    }
    
    return NULL;
  }
  
  T beginIteration_debug(){
    
    if (headNode != NULL) {
      
      debugIteratorNode = headNode;
      return debugIteratorNode->val;
    }
    
    return NULL;
  }

  
  /**
   *  Used for iteration
   */
  T getNextElement(){
    
    iteratorNode = iteratorNode->nextNode;
    
    // this indicates end iteration
    // iteratorNode == NULL , ensures begin is called
    // also if only one node was there in linkedList and it was
    // deleted during iteration it takes care of that.
    if ( iteratorNode == headNode || iteratorNode == NULL){
    
      iteratorNode = tailNode; // this is just to keep iteratorNode at last node
                               // in the ring
      return NULL;
    }
    
    return iteratorNode->val;
  }
  
  T getNextElement_debug(){
    
    debugIteratorNode = debugIteratorNode->nextNode;
    
    // this indicates end iteration
    // iteratorNode == NULL , ensures begin is called
    // also if only one node was there in linkedList and it was
    // deleted during iteration it takes care of that.
    if ( debugIteratorNode == headNode || debugIteratorNode == NULL){
      
      debugIteratorNode = tailNode; // this is just to keep iteratorNode at last node
      // in the ring
      return NULL;
    }
    
    return debugIteratorNode->val;
  }
  
  /**
   * TO be used during iteration
   */
  
  Node<T> *getCurrentNode (){
  
    return iteratorNode;
  }
  
  /**
   * Begin iteration
   */
  void endIteration(){
    
    
  }

  
  ~LinkedList(){
  
    if (headNode != NULL) {
      
         clearList();
    }
  }
};

#endif /* LinkedList_h */
