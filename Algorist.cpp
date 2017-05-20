//
//  Algorist.cpp
//  PlayArea
//
//  Created by zayan on 8/29/15.
//  Copyright (c) 2015 Eftiquar. All rights reserved.
//

#include "Algorist.h"
#include <vector>
#include <map>
using std::map;
using std::vector;
//knapsack

class KnapSackItem
{
public:
    KnapSackItem(size_t weight = 0, size_t value= 0):m_Weight(weight), m_Value(value){}
    
    size_t Weight()const {return  m_Weight;}
    
    size_t Value()const {return m_Value;}
    bool operator == (const KnapSackItem& item) const
    {
        return m_Weight == item.m_Weight && m_Value == item.m_Value;
    }
private:
    
    size_t m_Weight;
    size_t m_Value;
};

using KnapSackItems = vector<KnapSackItem>;
void DoKnapSack(const KnapSackItems& items, size_t capacity,vector<KnapSackItems>& finalSolutionVector)
{
    
    vector<KnapSackItems> interimSolutionVector(capacity);
    
    for(const auto& itemCurser: items)
        {
            size_t capacityCursor = 1;
            
            finalSolutionVector = interimSolutionVector;
            interimSolutionVector.assign(capacity, KnapSackItems());

            for( auto& valueCursor: finalSolutionVector)
            {

                if ( itemCurser.Weight() < capacityCursor)
                {
                    interimSolutionVector[capacityCursor-1] = valueCursor;
                }

                else // check if this item can be part of solution
                {
                        size_t prevBestVal = 0;
                        for(auto value : valueCursor)
                        {
                            prevBestVal += value.Value();
                        }
                    
                        size_t indexMinusThisItem = capacityCursor == itemCurser.Weight() ? itemCurser.Weight() - 1 : capacityCursor- itemCurser.Weight() -1;
                    
                        auto solMinusThisItem = finalSolutionVector[indexMinusThisItem];
                    
                        size_t valMinusThis = 0;
                        size_t weightMinusThis = 0;
                 
                        for(auto prevSol:solMinusThisItem)
                        {
                            valMinusThis += prevSol.Value();
                            weightMinusThis = prevSol.Weight();
                        }

                        if( (capacityCursor >= (weightMinusThis + itemCurser.Weight()) ) &&
                           ((valMinusThis + itemCurser.Value()) > prevBestVal))
                        {
                            solMinusThisItem.push_back(itemCurser);
                            interimSolutionVector[capacityCursor-1] = solMinusThisItem;
                        }
                        else if(itemCurser.Value() > prevBestVal)
                        {
                            solMinusThisItem.clear();
                            solMinusThisItem.push_back(itemCurser);
                            interimSolutionVector[capacityCursor-1] = solMinusThisItem;
                        }
                        else
                        {
                            interimSolutionVector[capacityCursor-1] = valueCursor;

                        }
                    
                    }
                ++capacityCursor;
            }
    }
    
    finalSolutionVector = interimSolutionVector;
}
size_t ATOI(const wchar_t* numberString)
{
    size_t res =0;
    while(*numberString)
    {
        res = res*10 + (*numberString++ - L'0');
    }
    return  res;
}

int main_back(int argc, const char * argv[])
{
    
    size_t res = ATOI(L"20100");
    KnapSackItems items
    {
        {4,3},
        {3,2},
        {2,4},
        {3,4}
    };

    
    
    map<size_t,size_t> myTree;
    myTree[0] = 1;
    myTree[1] = 2;
    myTree[2] = 3;
    
    
    
    auto mapIter = myTree.begin();
    mapIter--;
    while(mapIter != myTree.end())
    {
        mapIter++;
        
    }
    mapIter++;
    vector<KnapSackItems> finalSolutionVector;
    DoKnapSack(items,6,finalSolutionVector);
    return 0;
}
/*
 // Returns:  pointer to the left-most node under __x.
 // Precondition:  __x != nullptr.
 template <class _NodePtr>
 inline _LIBCPP_INLINE_VISIBILITY
 _NodePtr
 __tree_min(_NodePtr __x) _NOEXCEPT
 {
 while (__x->__left_ != nullptr)
 __x = __x->__left_;
 return __x;
 }

 // Returns:  pointer to the right-most node under __x.
 // Precondition:  __x != nullptr.
 template <class _NodePtr>
 inline _LIBCPP_INLINE_VISIBILITY
 _NodePtr
 __tree_max(_NodePtr __x) _NOEXCEPT
 {
 while (__x->__right_ != nullptr)
 __x = __x->__right_;
 return __x;
 }
 
 // Returns:  pointer to the next in-order node after __x.
 // Precondition:  __x != nullptr.
 template <class _NodePtr>
 _NodePtr
 __tree_next(_NodePtr __x) _NOEXCEPT
 {
 if (__x->__right_ != nullptr)
 return __tree_min(__x->__right_);
 while (!__tree_is_left_child(__x))
 __x = __x->__parent_;
 return __x->__parent_;
 }
 
 // Returns:  pointer to the previous in-order node before __x.
 // Precondition:  __x != nullptr.
 template <class _NodePtr>
 _NodePtr
 __tree_prev(_NodePtr __x) _NOEXCEPT
 {
 if (__x->__left_ != nullptr)
 return __tree_max(__x->__left_);
 while (__tree_is_left_child(__x))
 __x = __x->__parent_;
 return __x->__parent_;
 }
 
 // Returns:  pointer to a node which has no children
 // Precondition:  __x != nullptr.
 template <class _NodePtr>
 _NodePtr
 __tree_leaf(_NodePtr __x) _NOEXCEPT
 {
 while (true)
 {
 if (__x->__left_ != nullptr)
 {
 __x = __x->__left_;
 continue;
 }
 if (__x->__right_ != nullptr)
 {
 __x = __x->__right_;
 continue;
 }
 break;
 }
 return __x;
 }

 */
#include <list>
#include <array>
using std::vector;
using std::array;

template <typename InputIterator, typename OutputItrator, typename Transformer>
OutputItrator Ordered_transformer(InputIterator first,InputIterator last,OutputItrator out, Transformer trfr)
{
    OutputItrator result = out;
    while(first != last)
    {
        *out++ = trfr(*first++);
    }
    
    return result;
}

size_t ComputePartitions( const vector<size_t>& input,const vector<vector<size_t>>& offsets, size_t partitionNo, size_t elements, vector<vector<size_t>>& output)

{
    if(partitionNo == 1)
    {
        output.push_back(vector<size_t>(input.begin() ,input.begin() + elements));
        return 0;
    }
   
  
    ComputePartitions(input, offsets, partitionNo -1,offsets[elements-1][partitionNo-1] +1, output);
    output.push_back(vector<size_t>(input.begin() + offsets[elements-1][partitionNo-1] + 1 ,input.begin() + elements ));
    
    return 0;
}
vector<vector<size_t>> LinearPartition(const vector<size_t>& input, size_t partitions)
{
    vector<vector<size_t>> result ;
    vector<size_t> cumulativesums;
    vector<size_t> offsets(input.size());
    vector<vector<size_t>> distMap;
    
    size_t sum = 0;
    Ordered_transformer(input.begin(), input.end(), std::back_inserter(cumulativesums), [&sum](size_t elem)
                        {
                            return sum += elem;
                        });
    
    for(auto row =0; row < input.size();row++)
    {
        result.emplace_back(partitions,0);
        result[row][0] = cumulativesums[row];
    }
    for(auto col = 1; col < partitions;col++)
    {
        result[0][col] = cumulativesums[0];
    }
    for(auto row =0; row < input.size();row++)
    {
        distMap.emplace_back(partitions,0);
    }

   
    for(auto elemRow = 1;elemRow < input.size(); elemRow++)
    {
        for(auto partitionCol = 1; partitionCol < partitions;partitionCol++)
        {
            //x is where partitionCol - 1 partitions on left side end, from x + 1 to elemRow we insert last partition
            //we exhaustively calculate cost of each partition placement from zeroth position to one position less than elem row
            //meaning ,
            //distance map is not computed for zeroth partition as there is no partition left of it
            for(auto x = 0;x < elemRow; x++)
            {
            auto maxCost = std::max(result[x][partitionCol-1] , cumulativesums[elemRow] - cumulativesums[x] );
            if(result[elemRow][partitionCol] == 0)
            {
                result[elemRow][partitionCol] = maxCost;
                //offsets[elemRow] = x;
                distMap[elemRow][partitionCol] = x;
            }
            else if(result[elemRow][partitionCol] > maxCost)
            {
                result[elemRow][partitionCol] = maxCost;
                distMap[elemRow][partitionCol] = x;
            }
            }
        }
    }
    vector<vector<size_t>> output;

    
    ComputePartitions( input, distMap, partitions, input.size(), output);
    return result;
    
}