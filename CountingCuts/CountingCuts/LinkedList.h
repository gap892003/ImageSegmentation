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
  int sizeOfList; // not realiable , part of design , just for debugging
  
public:
  
  LinkedList(){
  
    headNode = NULL;
    tailNode = NULL;
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

    node->prevNode->nextNode = node->nextNode;
    node->nextNode->prevNode = node->prevNode;
    --sizeOfList;
    delete node;
  }

  void printList(){
  
    Node<T> *iNode = headNode;
    
    do {
      
      std::cout<< iNode->val << " ";
      iNode = iNode->nextNode;
    } while (iNode != tailNode);
    
    std::cout<< tailNode->val << " ";
    std::cout<< std::endl;
  }
  
  
  void printListRev(){
    
    Node<T> *iNode = tailNode;
    
    do {
      
      std::cout<< iNode->val << " ";
      iNode = iNode->prevNode;
    } while (iNode != headNode);
    
    std::cout<< headNode->val << " ";
    std::cout<< std::endl;
  }

  void clearList(){
  
    Node<T> *iNode = headNode;
    headNode = NULL;
    
    do {
      
      Node<T> *temp = iNode;
      iNode = iNode->nextNode;
      delete temp;
      temp = NULL;
    } while (iNode != tailNode);
    
    if (tailNode != NULL){
      
      delete tailNode;
      tailNode = NULL;
    }
  }

  /**
   * Begin iteration
   */
  T beginIteration(){
  
    iteratorNode = headNode;
    return iteratorNode->val;
  }
  
  /**
   *  Used for iteration
   */
  T getNextElement(){
    
    iteratorNode = iteratorNode->nextNode;
    
    // this indicates end iteration
    if ( iteratorNode == headNode ){
    
      return NULL;
    }
    
    return iteratorNode->val;
  }
  
  /**
   * Begin iteration
   */
  void endIteration(){
    
    
  }

  
  ~LinkedList(){
  
    clearList();
  }
};

#endif /* LinkedList_h */
