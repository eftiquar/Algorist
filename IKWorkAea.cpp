//
//  IKWorkAea.cpp
//  PlayArea
//
//  Created by zayan on 5/6/17.
//  Copyright © 2017 Eftiquar. All rights reserved.
//

#include "IKWorkAea.hpp"
std::vector<int>::iterator partition(std::vector<int>::iterator begin,std::vector<int>::iterator end, int bitmask)
{
    auto pivot = find_if(begin,end,
                         [bitmask](const int elem){
                             return !(elem & bitmask);
                         }) ;
    if(pivot == end )
        return begin;
    
    iter_swap(pivot, begin);
    pivot =  begin;
    ++pivot;
    for(auto next = pivot ; next != end;++next)
    {
        if( !((*next) & bitmask))
        {
           
            iter_swap(pivot++,next);
           
        }
        
    }
    return pivot;
}
int FindMissing(std::vector<int> & input, size_t bitcount)
{
    int result = 0 ;
    int bitmask = 1;
    auto begin = input.begin();
    auto end =  input.end();
    for(auto i =0; i < bitcount; ++i)
    {
        if(begin == end)
            break;
      
        auto group1 =  partition(begin,end, bitmask);
        auto sizeGroup0 = distance(begin,group1);
        auto sizeGroup1 = distance(group1,end);
        if(sizeGroup0 < sizeGroup1)
        {
            //missing number has bitmask off
            end = group1;
        }
        else
        {
            result|= bitmask;
            begin =group1;
        }
        bitmask <<= 1;;
    }
    //let us count int's ending with '0'
    //and those ending with '1'
    //whichever is lesser will have missing interger
    
    return result;
}
