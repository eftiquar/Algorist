//
//  IKWorkAea.hpp
//  PlayArea
//
//  Created by zayan on 5/6/17.
//  Copyright © 2017 Eftiquar. All rights reserved.
//

#ifndef IKWorkAea_hpp
#define IKWorkAea_hpp

#include <stdio.h>
#include <vector>
#include <iterator>
#include <algorithm>
#include <map>
using namespace std;
int FindMissing(std::vector<int> & input, size_t bitcount);
namespace IK
{
    
    template <class RandIter>
    RandIter Partition(RandIter first, RandIter last)
    {
        
        if(first == last)
            return last;
        
        RandIter mid = std::next(first, distance(first,last)/2);
        
        iter_swap(first, mid);
        RandIter pivotIter = find_if(first +1 , last, [first](const typename iterator_traits<RandIter>::value_type & elem)
                                         {
                                             return (*first < elem  );
                                         });
        
        if(pivotIter == last)
        {
            iter_swap(first, pivotIter-1);
            return last;
        }
        for( auto next = pivotIter ; next != last;++next)
        {
            if(!( *next > *first))
            {
                iter_swap(pivotIter++, next);
            }
        }
        iter_swap(first,pivotIter -1);
        return pivotIter;
    }
    
    template <class RandIter>
    RandIter QuickSelect(RandIter first, RandIter last, size_t rank)
    {
      
        if(first == last)
            return last;
        
        RandIter result = last;
        RandIter begin = first;
        while(first != last)
        {
           
           result = Partition(first, last);
            auto dist = distance(begin,result);
            if(dist < rank)
                first  = result ;
            
            else if(rank < dist)
                last = result -1;
            
            else
                break;
        }
        return result;
    }
//  space o(n) - n-1 + n-1
    // time - degree pow height
    //2 to pow 2(n-1)
    // tree with height = 2(n-1) which is space complexity
    //degree = 2
    // leaves = 2 pow 2(n-1) which is time complexity
#if 0
    int MaxPath(std::vector<vector<int>> grid, int x, int y, int endx, int endy)
    {
        
        if(x == endx && y == endy)
               return grid[x][y];
        
        if(x+1 == endx)
            return MaxPath(grid, x , y +1 , endx, endy) + grid[x][y];
        
        if(y+1 == endy)
            return MaxPath(grid, x +1, y , endx, endy) + grid[x][y];
        
        int right = MaxPath(grid, x +1, y, endx, endy);
        int down = MaxPath(grid, x , y +1, endx, endy);
        
        return std::max(right,down) + grid[x][y];
        
    }
   
    int MaxPath0(std::vector<vector<int>> grid)
    {
        int result = 0;
        
        return result;
    }
 #endif
    
    template<class T>
    class pqueue
    {
    public:
        pqueue(const vector<T>& in):m_queue(in.size() +1)
        {
            m_begin = m_queue.begin()+1;
            copy(in.cbegin(),in.cend(),m_begin);
        
            m_end  = m_queue.end();
            m_root  = m_queue.begin();
            
            heapify();
        }
        
        
        const T& top()const
        {
            if(size() == 0)
                throw "Underflow";
            return m_begin[0];
        }
        T pop()
        {
            auto res = top();
            iter_swap(m_begin, --m_end);
            sink(m_begin);
            
        }
        void push(const T & elem)
        {
            m_queue.push_back(elem);
            m_end = m_queue.end();
            swim(m_end-1);
        }
        
        bool empty()const
        {
            return  m_end == m_begin;
        }
    private:
        using qiterator = typename vector<T>::iterator ;
        void heapify()
        {
            auto back = firstnonleaf();
            while(back != m_root)
            {
                sink(back);
                back = parent(back);
            }
        }
        void sink(qiterator iter)
        {
            while(!isleaf(iter))
            {
                auto leftIter = leftchild(iter);
                auto rightIter = rightchild(iter);
                
                auto min = leftIter;
                
                if(rightIter != m_end && *rightIter < *min)
                    min = rightIter;
                if(*min < *iter)
                {
                    iter_swap(iter,min);
                    iter = min;
                }
                else
                {
                    break;
                }
                
            }
        }
        void swim(qiterator iter)
        {
            while(iter != m_root)
            {
                auto parentIter = parent(iter);
                if(*iter < *parentIter)
                {
                    swap_iter(iter,parentIter);
                    iter = parentIter;
                }
                else
                {
                    break;
                }
            }
            
        }
        qiterator leftchild(qiterator iter)
        {
            auto dist = std::distance(m_root, iter);
            
            if(size() < dist)
                return m_end;
            
            auto leftchild = dist*2;
            
            if(size() < leftchild)
                return m_end;
            
            return m_root + leftchild;
            
        }
        qiterator rightchild(qiterator iter)
        {
            auto leftchildIter = leftchild(iter);
            
            return leftchildIter != m_end ? leftchildIter + 1 : m_end;
            
        }
        qiterator parent(qiterator iter)
        {
            auto dist = std::distance(m_root, iter);
            dist /=2;
            return m_root + dist;
        }
        
        bool isleaf(qiterator iter)
        {
            return leftchild(iter) == m_end;
        }
        qiterator firstnonleaf()
        {
            return m_root + size()/2;
            
        }
        size_t size()const
        {
            return distance(m_begin, m_end);
        }
        qiterator m_begin;
        qiterator m_end;
        qiterator m_root;
       
        vector<T> m_queue;
        
        
    };
    
    template <class RandIter>
    void revvup(RandIter first, RandIter last)
    {
        while(first != last && first != --last)
        {
            //if is vowel
            iter_swap(first++,last);
            //else ++ first;
                
        }
    }

    /*
     _,_,_,_,
     1,2,3_ _
     _,_,1,2,3
     
     */
    template <typename RandIter>
    class CircularCacheLru
    {
        
    using Val_Type = typename std::iterator_traits<RandIter>::value_type;
    public:CircularCacheLru(RandIter first, RandIter last):
        m_first(first), m_last(last),m_back(last),m_front(last)
        {
        }
      void push_back(RandIter elem)
        {
            if(m_back == m_front)
            {
                m_front=  m_back = m_first;
                *m_back++= *elem;
                return;
            }
            if(m_back > m_front)//all elements lie between (front,back], no wrap around has happened yet
            {
                while(--m_back != m_front && *m_back<*elem)
                {
                    
                }
                
                if(m_back == m_front  )
                {
                    if(*m_back<*elem)
                    {
                        *m_back++ = *elem;
                    }
                    else if(*elem < *m_back)
                    {
                        if(++m_back != m_last)
                           *++m_back = *elem;
                        else
                           {//reset circular queue
                               *m_first = *m_front;
                               m_front = m_first;
                               m_back = m_front +1;
                               *m_back++ = *elem;
                           }
                    }
                    else
                    {
                        ++m_back;
                    }
                }
                else //back != front but NOT *m_back<*elem, also back > front so no wrap around has happened
                {
                    ++m_back;
                    if(m_back == m_last)//wrap around now
                    {
                        m_back = m_first;
                        *m_back++ = *elem;
                        if(m_front == m_first)
                            ++m_front;//evict first
                        
                    }
                }
                
            }
            else //wrapped around
            {
                
                while(--m_back != m_first && *m_back<*elem)
                {
                    if(m_back == m_first)
                    {
                        if(*m_back<*elem)//overwrite
                        {
                            *m_back++ = *elem;
                        }
                        else //need to insert
                        {
                           if(m_back +1 == m_front)
                           {
                               
                           }
                        }
                    }
                }
            }
           }
        void pop_front()
        {
        }
        Val_Type front()const
        {
            if(m_front != m_back)
            return *m_front;
            
            return Val_Type(0);//throw error
        }
        
    private:
        RandIter m_back;
        RandIter m_front;
        
        RandIter m_first;
        RandIter m_last;
        
    };
    template<class RandIter>
    vector<typename std::iterator_traits<RandIter>::value_type> MAxInSlidingWindow(RandIter first, RandIter last, size_t windowSize)
    {
        
    }
    //pqueue<int> pq(10);
    
 
    //1  printpermutations - space complexity - O(n), time complexity - n * n!
    //2 given an array , print subset
    //3 recursion - do stuff at leaf level instead of returning results back
    // 4 edit distance
    
    //1 3 4
    /*000
      001
      010
      011
     for subset 
     
     space - O(n) 
     time - 2 pow n
     
     Printsubsets (input, r+1, out, w)
     out[w] = inp[r]
     Printsubsets (input, r+1, out, w+1)
     
     N queens problem
     */
    
    template <class T>
    class TrieNode
    {
    public:
        TrieNode(T ch);
        T m_ch;
        std::map<T,TrieNode*> m_children;
    };
    
    
}

#endif /* IKWorkAea_hpp */
