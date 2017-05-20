//
//  Trees.h
//  PlayArea
//
//  Created by zayan on 2/15/16.
//  Copyright (c) 2016 Eftiquar. All rights reserved.
//

#ifndef __PlayArea__Trees__
#define __PlayArea__Trees__

#include <stdio.h>
#include "Algorist.h"
#include <limits>

template <typename T >
class TRNode {
    public:

    T m_data;
    TRNode* m_Left;
    TRNode* m_Right;
    TRNode(T data);
};



template<typename T>
TRNode<T>::TRNode(T data):m_data(data)
{
    
}
template <typename T, typename TIter>
TRNode<T>* BuildTreebalanced(TIter begin,TIter end)
{
    
    if(begin == end)
        return nullptr;
    
    auto mid = begin + (end- begin) /2;
    TRNode<T>* root = new TRNode<T>(*mid);
    root->m_Left = BuildTreebalanced<T>(begin,mid);
    root->m_Right = BuildTreebalanced<T>(mid+1, end);
    
    return root;
}
template <typename T>
bool IsBST(TRNode<T> * root,T upperBound, T lowerBound)
{
    if(!root)
        return true;
    if(root->m_data > upperBound || root->m_data < lowerBound)
      return false;
    
    if(root->m_Left && !IsBST(root->m_Left,root->m_data,lowerBound))
    {
        return false;
    }
    if(root->m_Right)
    {
        return IsBST(root->m_Right,upperBound,root->m_data);
    }
    return true;

   
}
template <typename T>
void OnAddNodeSL(TRNode<T>* elem,TRNode<T>*& head)
{
    elem->m_Right = head;
    head = elem;
}
template <typename T>
void FlattenSL(TRNode<T>* root,TRNode<T>*& head)
{
   if(root == nullptr)
       return;

    FlattenSL<T>(root->m_Right,head);
    OnAddNodeSL(root,head);
    FlattenSL<T>(root->m_Left,head);
    root->m_Left = nullptr;
}
template <typename T>
void OnAddNodeDL(TRNode<T>* elem,TRNode<T>*& head)
{
    elem->m_Right = head;
    if(head)
    {
        head->m_Left = elem;
    }
    head = elem;
}
template <typename T>
void FlattenDL(TRNode<T>* root,TRNode<T>*& head)
{
    if(root == nullptr)
        return;
    
    FlattenDL<T>(root->m_Right,head);
    OnAddNodeDL(root,head);
    FlattenDL<T>(root->m_Left,head);

}

template <typename T>
void ByOrder(TRNode<T>* node,vector<T>& output)
{
    output.push_back(node->m_data);
}
template <typename T>
void Preorder(TRNode<T>* root,vector<T>& output)
{
  if(root == nullptr)
      return ;
    
    ByOrder(root,output);
    Preorder(root->m_Left, output);
    Preorder(root->m_Right,output);
}

template <typename T>
void Inorder(TRNode<T>* root,vector<T>& output)
{
    if(root == nullptr)
        return ;
    Inorder(root->m_Left, output);
    ByOrder(root,output);
    Inorder(root->m_Right,output);
}
template <typename T>
TRNode<T>* TreeFromPreAndInOrder(const vector<T>& inorder,const vector<T>& preorder )
{
    
    auto preorderFirst = preorder.cbegin();
    auto preorderLast = preorder.cend();
    
    if(preorderFirst == preorderLast)
        return nullptr;
    
    auto inorderFirst = inorder.cbegin();
    auto inorderLast = inorder.cend();
    
    if(inorderFirst == inorderLast)
        return nullptr;
    
    TRNode<T>* root = new TRNode<T>(*preorderFirst);
    
    using NodeParentPair = std::pair<TRNode<T>*, TRNode<T>*>;
    
    queue<NodeParentPair> predecessors;
    predecessors.push( {root,nullptr});
    
    
    while(!predecessors.empty())
    {
        {
            auto pred = predecessors.front();
            predecessors.pop();
            auto preorderPos = std::find(preorderFirst,preorderLast,pred.first->m_data);
            auto inorderParentPos = std::find(inorderFirst,inorderLast,pred.first->m_data);

            if(preorderPos != preorderLast)
            {
                //advance to potential left child
                ++preorderPos;
                if(preorderPos == preorderLast)
                    break;
    
                {
                    auto inorderLeftPos = std::find(inorderFirst,inorderLast,*preorderPos);
                    if(inorderLeftPos < inorderParentPos)
                    {
                        //confirmed left child
                        pred.first->m_Left = new TRNode<T>(*inorderLeftPos);
                        predecessors.push({pred.first->m_Left,pred.first});
                        ++preorderPos;
                    }
                    
                }
                //advance to right child
                
                {
                    
                    for(;preorderPos != preorderLast;++preorderPos)
                    {
                        if(*preorderPos > pred.first->m_data)
                            break;
                    }
                    if(preorderPos == preorderLast)
                        continue;
                    
                    auto inorderRightPos = std::find(inorderFirst,inorderLast,*preorderPos);
                    if(inorderRightPos > inorderParentPos)
                    {
                       if(pred.second)
                       {
                           auto upperBoundPos = std::find(inorderFirst,inorderLast,pred.second->m_data);
                           if(inorderRightPos < upperBoundPos)
                           {
                               pred.first->m_Right = new TRNode<T>(*preorderPos);
                               predecessors.push( {pred.first->m_Right,pred.first});
                               
                           }
                       }
                        else
                        {
                            pred.first->m_Right = new TRNode<T>(*preorderPos);
                            predecessors.push( {pred.first->m_Right,pred.first});
                            
                        }
                    }
                }
                
                
            }
        }
    }
    
    return root;
}

template <typename T>
TRNode<T>* TreeFromPreOrder(const vector<T>& preorder )
{
    
    auto preorderFirst = preorder.cbegin();
    auto preorderLast = preorder.cend();
    
    if(preorderFirst == preorderLast)
        return nullptr;
    
    
    TRNode<T>* root = new TRNode<T>(*preorderFirst);
    
    using NodeParentPair = std::pair<TRNode<T>*, TRNode<T>*>;
    
    queue<NodeParentPair> predecessors;
    predecessors.push( {root,nullptr});
    
    
    while(!predecessors.empty())
    {
        {
            auto pred = predecessors.front();
            predecessors.pop();
            auto thisRootPos = std::find(preorderFirst,preorderLast, pred.first->m_data);
            if(/*thisRootPos != preorderLast &&*/ thisRootPos + 1 != preorderLast)
            {
                //check for left
                auto tryLeftPos = thisRootPos + 1;
                if(*tryLeftPos < pred.first->m_data)
                {
                    pred.first->m_Left = new TRNode<T>(*tryLeftPos);
                    predecessors.push({pred.first->m_Left,pred.first});
                }

                auto tryRightPos = std::find_if(thisRootPos,preorderLast, [thisRootPos](T key)
                {
                    return *thisRootPos < key;
                });
                if(tryRightPos != preorderLast)
                {
                    //if(pred.second && *tryRightPos < pred.second->m_data)
                    if(pred.second)
                    {
                        if(*tryRightPos < pred.second->m_data)
                        {
                        pred.first->m_Right = new TRNode<T>(*tryRightPos);
                        predecessors.push({pred.first->m_Right,pred.first});
                        }
                        
                    }
                    else
                    {
                        pred.first->m_Right = new TRNode<T>(*tryRightPos);
                        predecessors.push({pred.first->m_Right,pred.first});
                    }
                }
            }
            
        }
    }
    
    return root;
}

template <typename T>
void OnLevel(vector<TRNode<T>*> & nodes)
{
    nodes.clear();
}
template <typename T>
void LevelWalk(TRNode<T> * root)
{
    if(!root)
        return;
    
    queue<TRNode<T>*> parents;
    queue<TRNode<T>*> children;
    parents.push(root);
    vector<TRNode<T>*> result;
    while(!parents.empty())
    {
        auto thisParent = parents.front();
        result.push_back(thisParent);
        parents.pop();
        if(thisParent->m_Left)
            children.push(thisParent->m_Left);
        if(thisParent->m_Right)
            children.push(thisParent->m_Right);
        

        
        if(parents.empty())
        {
            OnLevel(result);
            swap(parents,children);
        }
        
    }
}

template <typename RandIter>
void Rotate(RandIter first,RandIter last)
{
    if(first == last)
        return;
    
    for( auto next = first+1; next != last;next ++)
    {
        std::swap(*next,*first);
    }
}

template <typename RandIter>
void Rotate(RandIter first,RandIter nFirst,RandIter last)
{
    
    for( auto next = nFirst; first != next;)
    {
        std::swap(*first++,*next++);
        if(next == last)
            next = nFirst;
        else if(first == nFirst)
            nFirst = next;
    }
}
template <typename T>
void FindKHelper(TRNode<T>* root,std::vector<TRNode<T>*>& result,size_t target)
{
    if(root == nullptr)
        return ;
    FindKHelper(root->m_Right,result,target);
    result.push_back(root);
    if((target +1) == result.size())
    {
        return;
    }
    FindKHelper(root->m_Left,result,target);

    
}

template <typename T>
TRNode<T>* FindKLargest(TRNode<T>* root, size_t k)
{
    if(root == nullptr)
        return nullptr;

    std::vector<TRNode<T>*> result;

    FindKHelper(root,result,k);
    
}
template <class T>
TRNode<T>* ReverseIt(TRNode<T>* list)
{
    
    TRNode<T> * reverse = nullptr;
    while (list)
    {
        auto next = list->m_Right;
        list->m_Right = reverse;
        reverse = list;
        list = next;
    }
    return reverse;
}

#endif /* defined(__PlayArea__Trees__) */
