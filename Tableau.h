//
//  Tableau.h
//  PlayArea
//
//  Created by zayan on 8/15/15.
//  Copyright (c) 2015 Eftiquar. All rights reserved.
//

#ifndef __PlayArea__Tableau__
#define __PlayArea__Tableau__

#include <stdio.h>
#include <utility>



namespace Tableau

{

        // IMPORTANT -- please read the behaviour of BSearch
        //@brief : Returns lower bound for key in the sorted array of T ranged as [begin,end).
        //if no key is found, it returns where the key should be placed to maintain sorted-ness
        //returns end if the key is greater than the largest element.
        //returns begin if key is smaller than the smallest element, or key is equal to the smallest element
        //end is the position immediately past the last element and should never be dereferenced.
    
        //begin + (end - begin)/2 is preferred over (begin + end) /2 to avoid range overflow errors, mathematically the  expressions evalaute to
        //same result except under overflow conditions.
        //more importantly, mid is bounded as [begin,end), thus mid is never an invalid index, hence no checks are needed to validate it
        //why ?
        //begin < end ==> begin + end < end + end
        //therefore (begin + end) < 2*end
        //(begin + end)/2 < end
        //hence mid < end
        // begin + (end - begin)/2;
        //if end - begin == 1, then mid = begin, so mid can equal begin, but can never equal end
        //our loop sets either begin = mid +1 (potentially making begin = end) or end = mid so the inherent asymmetry between begin being valid position and end one off the last
        //element helps us converge begin and end

        //also we are using just one comparison - strictly greater than. T just needs to provide > opererator.
    
        template <typename T>
        T* BSearch(T* begin,T* end,const T& key)
        {

            while (begin != end)
            {
                T* mid = begin + (end - begin)/2;
        
                if(key > *mid) // key is strictly greater than mid, begin from next element onwards
                    begin = mid + 1;
                else // key is either less than mid or equal to mid, exclude this element from search, let us see if there is duplicate key else we will converge here if no duplicate exists. If key is less then end needs to be excluded to look for bounds lower than end
                    end = mid;
            }
    
            return begin;

        }
    
    

        // @brief swaps elements lhs and rhs

 
        template <typename T>
        inline void Swap(T& lhs, T& rhs)
        {
            T tmp = std::move(lhs);
            lhs = std::move(rhs);
            rhs = std::move(tmp);
        }
    
        //@brief reverse elements in [begin,end). In-place.
        template <typename T>
        void InplaceReverse(T * begin, T* end)
        {
            while(begin < end)
            {
                Swap(*begin++,*--end); // economy of expression, thy name is C 
            }
        }
    
};

#endif /* defined(__PlayArea__Tableau__) */