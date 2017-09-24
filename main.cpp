//
//  main.cpp
//  PlayArea
//
//  Created by zayan on 6/20/15.
//  Copyright (c) 2015 Eftiquar. All rights reserved.
//
/*
 6
 /   | \
 1    2  5
 /| \  |
 0 1  2 1
 |
 7
 
 6 
 
 
 long myFindSum(Node *node long parentsum) {
 
 if(node-> children -> empty)
 cout >>  parentsum *10 + node -> data;
 
 if(node)
 {
 
 for each child i
 {
 //return parentum*10 + myfindsum(childhood , node-data * 10)
   myfindsum(i , parentsum *10 + node-data )
 }
 }
 }
 */
#include <iostream>
#include <memory>
#include "assert.h"
#include <functional>
#include <algorithm>
#include <map>
#include <stack>
#include "play.h"
#include "Trees.h"

#include "IKWorkAea.hpp"


constexpr wchar_t NULL_CHAR = L'\0';
constexpr wchar_t SUFFIX_CHAR = L'$';
using std::string;

string sortCharacters(string inString);
template < typename T>
size_t RankIT(T in[], size_t end, const T& val );

bool RegExpCore(const wchar_t* inputString,const wchar_t* regExp);

bool WildCardMatch(const wchar_t* inputString,const wchar_t* regExp, wchar_t prefix)
{

    do
    {
         if(RegExpCore(inputString,regExp))
             return true;
    }    while( *inputString && (*inputString++ == prefix || prefix == L'.'));
    
    return false;
}
bool RegExpCore(const wchar_t* inputString,const wchar_t* regExp)
{
    if(!*inputString  && !*regExp )
        return true;
    if(*regExp == L'$' && *(regExp +1) == NULL_CHAR )
        return *inputString == NULL_CHAR;
    
    if((*(regExp +1)) == L'*')
    {
        return WildCardMatch(inputString, regExp +2, *regExp);
    }
    
    if(*inputString && (*inputString == *regExp || *regExp == L'.'))
        return RegExpCore(++inputString, ++regExp);
    
    return false;
}

bool RegExp(const wchar_t* inputString,const wchar_t* regExp)
{
    if(*regExp == SUFFIX_CHAR)
        return RegExpCore(inputString, ++regExp);
    
    do
    {
        if(RegExpCore(inputString, regExp))
            return true;
    }while(*inputString++);
    return false;
}
template <typename T>
size_t FindTurningPoint(T in[], size_t end)
{
    size_t begin = 0;
    while(begin != end)
    {
        size_t mid = begin + ( end - begin)/2;
        //note that the aim is to find original begin
        //if mid is above rotated begin (in) , it cannot be turning point
        if(in[mid] > *in)
            begin = mid +1;
        else
            end = mid;
    }
    return begin;
}

//1,2,3,4,5,5,6,6

// 56612345
//123455555

//55123455

/*
copy the set of points 
 1. PX = sorted by x co-ordinates - XSort
 2. PY = sorted by y co-ordinates   YSort
 
 Q = left half of PX, PY
 >>QX, QY
 R = right half of PX, PY
 >> RX, RY
 
 D & C ( divide and conquer) recursively find P1Q1 and P2Q2
 P1Q1 = ClosestPair(QX,QY)
 P2Q2 = ClosestPair(RX,RY)
 
 delta = min(P1Q1,P2Q2)
 P3Q3 = ClosestSplitPair(PX,PY,delta)
 xbar = |PX|/2 element, means largest x co-ordinate in left half of PX
 SY = [xbar-delta, xbar + delta] sorted by y co-ordinates
 the shortest pair is within 7 positions in the SY strip
 
 return min(p1Q1,P2Q2,P3,Q3)
 
 */
//fib n = fib(n-1) + fib(n-2)

void Justify(size_t width, const vector<std::string>& tokens)
{
    if(!width || !tokens.size())
        return;
    
    size_t lengthSofar = 0;
    auto first = tokens.begin();
    for(auto token = first;token != tokens.end(); ++token )
    {
        if(lengthSofar + token->length() < width)
        {

            lengthSofar += (token->length() +1);
            if(lengthSofar < width)
            {
                continue;
            }

        }
        
        if(lengthSofar == token->length() == width)
        {
         
            std::copy(first,token + 1,std::ostream_iterator<std::string>(std::cout," "));
            std::cout << "\r\n";
            first = token +1;

            lengthSofar = 0;
        }
        else
        {
            std::copy(first,token ,std::ostream_iterator<std::string>(std::cout," "));
            std::cout << "\r\n";
            
            //token -= 1;
            first = token;
            
            lengthSofar = 0;
        }
    }
    
    if(first != tokens.end())
    {
            std::copy(first,tokens.end() ,std::ostream_iterator<std::string>(std::cout," "));
    }
}


template <typename T>
class Node
{
// friend
public:
    using NodePtr = Node<T>* ;
    Node<T>(T* key):m_key(key),m_pNext(nullptr)
    {}
    
private:
    T* m_key;
    Node<T>* m_pNext;
    
};
size_t Factorial_RT(size_t partial_result,size_t factorialcursor)
{
    return factorialcursor < 2 ? partial_result : Factorial_RT(partial_result*factorialcursor,factorialcursor-1);
}
size_t Factorial(size_t input)
{
    return Factorial_RT(1, input);
}
size_t power_rt(size_t partial_result,size_t base, size_t pow)
{
    if(pow ==1)
        return partial_result;
    return power_rt(partial_result*base,base,pow -1);
    
}
size_t Power(size_t base,size_t pow)
{
    assert(base != 0 || pow !=0);
    
    if (pow ==0)
    {
        return 1;
    }
    if(base == 0)
        return 0;
    
    return power_rt(1, base, pow);
}

size_t PowerBasic(size_t base,size_t pow)
{
    assert(pow != 0 || base !=0);
    
    size_t result = 1;
    
    while(pow)
    {
        if(pow %2 ==0)
        {
            base *= base;
            pow/=2;
        }
        
        result *= base;
        pow--;
        
    }
    return  result;
}

template <typename T>
Node<T>* ReverseTR(Node<T> * head,Node<T>*& reverse)
{
    if(head == nullptr)
        return reverse;
    
    Node<T>* next = head->m_pNext;
    head->m_pNext = reverse;
    reverse = head;
    return ReverseTR(next, reverse);
}

template <typename T>
Node<T>* Reverse(Node<T> * head)
{
    Node<T>* reverseList = nullptr;
  
    while (head)
    {
        Node<T>* next = head->m_pNext;
        head->m_pNext = reverseList;
        reverseList = head;
        head = next;
    }
    return  reverseList;
}

//n is used as counter
size_t fib_aux(size_t fib0,size_t fib1, size_t n )
{
    if(n ==0)
        return fib0;

    if(n==1)
        return fib1;
    
    return fib_aux(fib1, fib0 + fib1, n-1);
    
}
using namespace std;
// 1,2,3,4,5,<end>
size_t FibRt(size_t n )
{
    return fib_aux(0, 1, n);
}
// fib(n-1) = fib(n-2) + fib(n-3)
size_t Fib(size_t input)
{
    size_t prev0 = 0;
    size_t prev1 = 1;
    size_t thisFib = 0;
    if(input == 0 || input ==1)
        return 1;
    
    for (size_t cursor = 2;cursor!= input;cursor++)
    {
        thisFib = prev0 + prev1;
        prev0 = prev1;
        prev1 = thisFib;
    }
    
    return thisFib;
}
//2223
//2189

constexpr unsigned MAX_SIZE = 255;
size_t fibTable[MAX_SIZE]= {0};
size_t FibDynamic(size_t in)
{

    if(in <2)
    {
        return fibTable[in] = in;
    }
    if(fibTable[in])
        return fibTable[in];
    
    return fibTable[in] = FibDynamic(in-1) + FibDynamic(in-2);

}

template <typename T >
void RShift(T array[],size_t fromPosition, size_t nPositions, size_t endPos)
{
    size_t source = endPos;
    size_t destination = endPos + nPositions;
    size_t end = fromPosition;
    while(source != end)
    {
        array[--destination] = array[--source];
    }
}

size_t CountOnes(size_t inputVal)
{
    size_t nOnes = 0;
    while(inputVal)
    {
        /*
        if(inputVal & 0x1)
            ++nOnes;
        
        inputVal /= 2;
         */
        if(inputVal & ~(inputVal-1))
        {
            ++nOnes;
        }
        inputVal &= (inputVal-1);
    }
    
    return nOnes;
}



#include <random>
size_t Random(size_t min, size_t max)
{
    
    std::random_device rd;
    std::uniform_int_distribution<size_t> dis(min, max);
    return dis(rd);
    
}
//when we are seeking upper bound, for duplicate keys, end will get locked at the element after first
//non-duplicate key or at the end if no such key exists
//we will approach the upper bound by setting begin = mid. as mid is always < end, it will never become
// ==end so begin will never become ==end, instead if(mid == begin && key == in[mid]) will be exit condition
//if the key is non existent, we either set end = mid or begin = mid +1 which eventually converge and
// begin != end becomes false.why ? because  end = mid ensures that end approaches begin as mid is in range
//[begin, end) and mid approaches begin or begin = mid + 1 ensures begin approahes end
template <typename T >
size_t UpperBound(const T in[],size_t end, T key)
{
    size_t begin = 0;
    while(begin != end)
    {
        size_t mid = begin + (end - begin)/2;
        
        
        if(mid == begin)
            break;
        
        if( key < in[mid])
            end = mid;
        else
            begin = mid;
        
        /*else if(key == in[mid]) // this condition makes converging hard when we have duplicates as mid can never equal end, hence begin can never equal end, but mid will equal begin when convergence happens
         begin = mid ;*/

        
    }
    
    return begin ;
}

//this is clone of  bin search, as bin search returns index of first occurence of matching key
//or
template <typename T >
size_t LowerBound(const T in[],size_t end, T key)
{
    size_t begin = 0;
    while(begin != end)
    {
        size_t mid = begin + (end - begin)/2;
        
        if( in[mid] < key)
            begin = mid +1;
        else
            end = mid;
        
    }
    
    return begin ;
}

template <typename T >
size_t BinSearch(const T in[],size_t end, T key)
{
    size_t begin = 0;
    while(begin != end)
    {
        size_t mid = begin + (end - begin)/2;
        
        if( key > in[mid])
            begin = mid +1;
        else
            end = mid;
        
    }
    
    return begin;
}
template <typename T>
void Swap(T in[], size_t lPos, size_t rPos)
{
    T temp = in[lPos];
    in[lPos] = in[rPos];
    in[rPos] = temp;
}
template <typename T >
void Shuffle( T in[],size_t end)
{
    for(size_t begin =0; begin != end; begin++)
    {
        size_t random = Random(begin,end);
        Swap<T>(in, begin,random);
        
    }
}

template <typename T>
void Dijkstra(T in[],size_t begin, size_t end)
{
  
    size_t whiteBegin = begin;
    size_t blueBegin = end;
    size_t xBegin = begin;
    
    T white = in[begin + (end - begin)/2];
    
    while(xBegin != blueBegin)
    {
        if(in[xBegin] < white) // it's red
        {
            Swap<T>(in,whiteBegin++,xBegin++);
            
        }
        else if(in[xBegin] == white)
        {
            xBegin++;
        }
        else // it's blue
        {
            Swap(in,xBegin,--blueBegin);
        }
    }
}
/*
 mid = begin + (end - begin)/2 = (2*begin + end - begin) /2 = (begin + end)/2
 begin < end ==> begin + end < end + end ==> begin + end < 2*end
 hence (begin + end) /2 < end , hence mid < end
 */
template <typename T>
void QSortDjikstra( T in[],size_t begin, size_t end)
{
    if(begin == end)
        return;
    
    size_t whiteBegin = begin;
    size_t blueBegin = end;
    size_t xBegin = begin;
    
    T white = in[begin + (end - begin)/2];
    //white is the pivot element, this algorithm accounts for duplicate
    // pivot entries and places them in right place
    //if
    while(xBegin != blueBegin)
    {
        if(in[xBegin] < white) // it's red
        {
            Swap<T>(in,whiteBegin++,xBegin++);
            
        }
        else if(in[xBegin] == white)
        {
            xBegin++;
        }
        else // it's blue
        {
            Swap(in,xBegin,--blueBegin);
        }
    }
    QSortDjikstra<T>(in,begin,whiteBegin);
    QSortDjikstra<T>(in,blueBegin,end);

}

template <typename T>
void QSort( T in[],size_t begin, size_t end)
{
    if(begin == end)
        return;

    Swap<T>(in, begin, begin + (end - begin)/2);
    size_t partitionIndex = begin ;
    
    for(size_t cursor = partitionIndex + 1;cursor < end;cursor ++ )
    {
        if(in[cursor] < in[begin])
        {
            Swap<T>(in,cursor,++partitionIndex);
        }
    }
    Swap<T>(in,begin,partitionIndex);
    
    QSort(in, begin,partitionIndex );
    
    QSort(in, partitionIndex +1, end);
    
}
// 0314
//1
template <typename T>
size_t Partition(T in[] , size_t begin, size_t end)
{
    size_t mid = begin + (end - begin)/2;
    size_t pivotIndex = begin;
    Swap<T>(in, begin, mid );
    for(size_t cursor = begin +1;cursor != end; cursor++)
    {
        if(in[cursor]< in[begin])
        {
            Swap(in, cursor, ++pivotIndex);
        }
    }
    Swap<T>(in,begin,pivotIndex);
    return pivotIndex;
}
template <typename T>

size_t QSelect(T in[],size_t end, T key)
{
    size_t begin = 0;
    if(begin == end)
        return begin;
    size_t pivot = 0;
    while(begin != end)
    {
        pivot = Partition(in,begin, end);
        
        if(in[pivot] == key)
        {
            return pivot;
        }
        if(key > in[pivot] )
        {
            begin = pivot +1;

        }
        else
        {
            end = pivot ;
        }
        
    }
    return begin;
}

//0,1,2,3
// mid = 0 + 4-0 / 2 >> 2
//0,1,2,3,4 >> 0 + 0 + 5-0/2 = 2 >> 0,1 && 2,3,4

template <typename T>
void InsertionSort(T* inputArray, size_t nLength)
{
    for(size_t begin =1; begin != nLength; begin++)
        for(size_t cursor = begin; cursor != 0 && inputArray[cursor]< inputArray[cursor-1]; cursor--)
            Swap(inputArray,cursor, cursor-1);
}

template <typename T >
class MergerSorter {
public:
    MergerSorter(T * inputArray, size_t nLength):m_array(inputArray),m_Length(nLength)
    {
        m_auxiallary = new T[nLength];
    }
    ~MergerSorter()
    {
        delete [] m_auxiallary;
    }
    
    void Sort()
    {
        MergeSort(0, m_Length );
    }
    
    void BottomUp()
    {
        for(size_t subArrayLength =1; subArrayLength < m_Length; subArrayLength*=2)
        {
            for(size_t subArrayBegin = 0;subArrayBegin < m_Length - subArrayLength; subArrayBegin += 2*subArrayLength)
                MergePartitions(subArrayBegin, subArrayBegin + subArrayLength,subArrayBegin + subArrayLength,
                                subArrayBegin + subArrayLength + subArrayLength > m_Length ? m_Length :   subArrayBegin + 2*subArrayLength);
        }
    }
void    MergePartitions(size_t lhsBegin,size_t lhsEnd,size_t rhsBegin, size_t rhsEnd)
    {
        size_t auxBegin = lhsBegin;
        while(auxBegin < rhsEnd)
        {
            m_auxiallary[auxBegin] = m_array[auxBegin];
            auxBegin++;
            
        }

        auxBegin = lhsBegin;
        while(auxBegin < rhsEnd)
        {
            if(lhsBegin == lhsEnd)
                m_array[auxBegin++] = m_auxiallary[rhsBegin++];
            
            else if(rhsBegin == rhsEnd)
                m_array[auxBegin++] = m_auxiallary[lhsBegin++];
            
            else if(m_auxiallary[rhsBegin] < m_auxiallary[lhsBegin])
                m_array[auxBegin++] = m_auxiallary[rhsBegin++];
            
            else
                m_array[auxBegin++] = m_auxiallary[lhsBegin++];
            
        }
        
    }
    
private:
void    MergeSort(size_t begin, size_t end)
    {
        if(end - begin < 2)
            return ;
        
        size_t mid = begin + (end-begin)/2;
        MergeSort(begin, mid );
        MergeSort(mid, end);
        MergePartitions(begin,mid , mid , end);
    }
    T* m_array;
    size_t m_Length;
    T* m_auxiallary;
};

#include "Maze.h"

void MazeMain(const char* mazeFile)
{
    CUndirectedGraph maze;
    ifstream in(mazeFile);
    in >> maze;
    
}
template <typename T>
void ReverseInplace(T in[], size_t end)
{
    assert(end >0);
    size_t begin = 0;
    while (begin < end)
    {
        Swap(in, begin++, --end);
      
    }
    
}
inline size_t Parent(size_t thisIndex)
{
    return (thisIndex + 1)/2 -1;
}
inline size_t LeftIndex(size_t thisIndex)
{
    return (thisIndex+1)*2 -1;
}

inline size_t RightIndex(size_t thisIndex)
{
    return (thisIndex +1)*2;
}

template<typename T>
void Heapify(T v[],size_t index, size_t length)
{
    size_t max = index;
    size_t left = LeftIndex(index);
    size_t right = RightIndex(index);
    
    if(left < length && v[left] > v[max])
        max = left;
    if(right < length && v[right] > v[max] )
        max  = right;
    
    if(max != index)
    {
        Swap(v,max,index);
        Heapify(v,max,length);
    }
}

template<typename T>
void BuildHeap(T v[],size_t length)
{
    for(size_t begin = (length/2) -1;begin !=0; begin--)
        Heapify(v,begin, length);
}
template<typename T>
void HeapSort(T v[],size_t length)
{
    BuildHeap(v,length);
    for(size_t begin = length -1;begin !=0;begin--)
    {
        Swap(v,0,--length);
        Heapify(v,0,length);
        
    }
}
template<typename T>
struct Less
{
    bool operator()( T l, T r)
    {
        return l < r;
    }
};

template <typename T, size_t N =10,typename C  = Less<T>>
class PriorityQueue
{
public:
    PriorityQueue(T* v ):m_queue(v),m_nIndex(N)
    {
 
    }
    void Insert(const T & key)
    {
        C less;
        m_queue[++m_nIndex] = key;        //index 0 is sentinel , items are one based, please note
        size_t index = m_nIndex;
        while(index >1 && less(m_queue[index],m_queue[index/2]))
        {
            Swap(m_queue,index, index/2);
            index = index/2;
        }
    }
    void Swim(size_t index)
    {
        C less;
        //index ==1 implies we are at root
        while(index > 1 && less (m_queue[index],m_queue[index/2]))
              {
                  Swap(m_queue, index,index/2);
                  index = index/2;
              }
    }
    void Heapify()
    {
        for(size_t begin = m_nIndex/2;begin >=1;begin--)
        {
            Sink(begin);
        }
    }
    void Sort()
    {
        Heapify();
        while(m_nIndex >1)
        {
            Swap(m_queue,1,m_nIndex--);
            Sink(1);
        }
    }
    
/*
 
 Children 
 
 (n +1)*2 -1 
 (n +1)*2 
 
 (n +1)/2 -1 
 
     0
  1     2
 3  4  5  6
 
 */
    void Sink(size_t index)
    {
        C less;

        while(index*2 <= m_nIndex)
        {
            size_t child = index*2;
            if(child < m_nIndex && less(m_queue[child],m_queue[child +1]))
                child ++;
            if(!less(m_queue[index],m_queue[child]))
                break;
            Swap(m_queue,index, child);
            index = child;
        }
        
    }
    T Delete()
    {
        size_t delCursor = 1;
        T ret = m_queue[delCursor];
        Swap(m_queue,1,m_nIndex--);
        
        C less;
        while(delCursor*2 <= m_nIndex)
        {
            size_t j = delCursor*2;
            if( j < m_nIndex && less(m_queue[j],m_queue[j+1]))
                j++;
            if(!less(delCursor,j))
                break;
            Swap(m_queue,j,delCursor);
            delCursor = j;
        
        }
    }
private:
    size_t m_nIndex = 0;
    T * m_queue;
    
};

template <typename K = size_t, typename V = size_t >
class TNode
{
    using NodePtr = TNode<K,V>*;
public:
    TNode(const K & key = K{}, const V & val = V{}):m_key(key), m_val(val)
    {}
    
//private:
    NodePtr m_left, m_right,m_parent;
    K m_key;
    V m_val;
};
template <typename K , typename V >
using NodePtr = TNode<K,V>*;

//thisNode is assumed to be non-null
template <typename K,typename V>
NodePtr<K,V> TreeMin(NodePtr<K,V> thisNode)
{
    while(thisNode->m_left)
    {
        thisNode = thisNode->m_left;
    }
    return thisNode;
}
template <typename K,typename V>
NodePtr<K,V> TreeMax(NodePtr<K,V> thisNode)
{
    while (thisNode->m_right)
    {
        thisNode = thisNode->m_right;
    }
    
    return thisNode;
}
template <typename K,typename V>
bool TreeIsLeftChild(NodePtr<K,V> thisNode)
{
    return thisNode == thisNode->m_parent->m_left;
}


template <typename K,typename V>
NodePtr<K,V> NextNode(NodePtr<K,V> thisNode)
{
    if(thisNode->m_right)
        return TreeMin(thisNode->m_right);
    
    while(!TreeIsLeftChild(thisNode))
        thisNode = thisNode->m_parent;
    
    return thisNode->m_parent;
}


template <typename K = size_t, typename V = size_t >
class Tree
{
public:

    Tree()
    {
        if(m_endMarker.m_nodePtr)
        {
        m_endMarker.m_nodePtr->m_left = m_root;
        }
    }
    void Insert(const K& key, const V & val)
    {

        if(!m_root)
        {
            m_root= new TNode<K,V>(key,val);
            return;
        }
        
        NodePtr<K,V> cursor = m_root;
        while(cursor)
        {
            if(key < cursor->m_key)
            {
                if(cursor->m_left)
                {
                    cursor = cursor->m_left;
                    continue;
                }
                else
                {
                    cursor->m_left = new TNode<K,V>(key,val);
                    cursor->m_left->m_parent = cursor->m_right;
                    break;
                }
            }
            
            else if(cursor->m_Val > key )
                {
                    if(cursor->m_right)
                    {
                        cursor = cursor ->m_right;
                        continue;
                    }
                    else
                    {
                        cursor->m_right = new TNode<K,V>(key,val);
                        cursor->m_right->m_right->m_parent = cursor;
                        break;
                    }
                }
        }
    }
    class Iterator
    {
        friend class Tree<K,V>;
    public:

        Iterator& operator ++()
        {
            if(m_nodePtr)
            {
                m_nodePtr = NextNode<K,V>(m_nodePtr);
            }
            return *this;
        }
        Iterator operator++(int)
        {
            Iterator current = *this;
            ++(*this);
            return current;
        }
        TNode<K,V>& operator *()
        {
            if(!m_nodePtr)
                throw L"Null pointer derefeenced";
            
            return *m_nodePtr;
        }
        NodePtr<K,V> operator->()
        {
            if(!m_nodePtr)
                throw L"Null pointer derefeenced";

            return m_nodePtr;
        }
    private:
        Iterator(NodePtr<K,V> node = nullptr):m_nodePtr(node)
        {}
        NodePtr<K,V> m_nodePtr;
    };
    
    Iterator end()
    {
        if(m_root)
        {
            m_endMarker.m_nodePtr->m_left = m_root;
        }
        return Iterator(m_endMarker);
    }
    Iterator begin()
    {
        return m_root == nullptr ? Iterator(m_endMarker) :  Iterator(TreeMin(m_root));
    }
    
private:
    Iterator m_endMarker;
    NodePtr<K,V> m_root;
};
pair<int, int> FindMaximumSubarray(const vector<int>& A)
{
    // A[range.first : range.second - 1] will be the maximum subarray.
    
    pair<int, int> range(0, 0);
    int min_idx = -1, min_sum = 0, sum = 0, max_sum = 0;
    for (int i = 0; i < A.size(); ++i)
    {
        sum += A[i];
        if (sum < min_sum)
        {
            min_sum = sum, min_idx = i;
        }
        if (sum - min_sum > max_sum)
        {
            max_sum = sum - min_sum, range = {min_idx + 1, i + 1};
            
        }
    }
    return range;
}
constexpr char SumMatrix[10][10] =
{
    {0,1,2,3,4,5,6,7,8,9},
    {1,2,3,4,5,6,7,8,9},
    {2,3,4,5,6,7,8,9},
    {3,4,5,6,7,8,9},
    {4,5,6,7,8,9},
    {5,6,7,8,9},
    {6,7,8,9},
    {7,8,9},
    {8,9},
    {9}
};

class LargeDigit
{
public:
    LargeDigit(const char* digit):m_digit(digit)
    {}
    const std::string Get()const
    {
        return m_digit;
    }
private:
    const std::string m_digit;
};
std::string operator + (const LargeDigit& lhs, const LargeDigit& rhs)
{
    std::string result;
    auto sumCursorLhs = lhs.Get().rbegin();
    auto sumCursoRhs = rhs.Get().rbegin();
    for(auto lhsDig = sumCursorLhs;lhsDig!= lhs.Get().rend();lhsDig++)
    {
        
    }
    return result;
}
#include <vector>
#include <string>
#include <iterator>
using std::wstring;
using std::vector;
using STRLIST = vector<wstring>;
/*
 A 
 AB 
 BA
 
 C
 CAB, 
 ACB,
 ABC
 CBA,
 BCA,
 BAC
 
 D
 DCAB,CDAB,CADB,CABD
 DACB,ADCB,ACDB,ACBD,
 cabernet
 */

STRLIST Permutate(wstring ch, wstring strIn)
{
    STRLIST strPerms;
    auto pos = strIn.begin();

    do
    {
        wstring strThisPerm = strIn;
        strThisPerm.insert(std::distance(strIn.begin(), pos),ch);
        strPerms.push_back(strThisPerm);
        
    }while(pos++ != strIn.end());
    
    return strPerms;
}
STRLIST Permutate(wchar_t ch, wstring strIn)
{
    STRLIST strPerms;
    auto pos = strIn.begin();
    
    do
    {
        wstring strThisPerm = strIn;
        strThisPerm.insert(pos,1,ch);
        strPerms.push_back(strThisPerm);
        
    }while(pos++ != strIn.end());
    
    return strPerms;
}

STRLIST Permutate(wstring strIn)
{
    STRLIST strPerms;

    for(auto elem :strIn)
    {
        
        if(strPerms.empty())
        {
            strPerms.push_back(wstring(1,elem));
            continue;
        }
        STRLIST interim;
        for(auto perm : strPerms)
        {
            STRLIST strPerm = Permutate(elem, perm);
            std::copy(strPerm.begin(), strPerm.end(),back_insert_iterator<STRLIST>(interim));
        }
        strPerms = interim;
    }
    
#if 0

    if(strIn.length() == 1)
    {
        strPerms.push_back(strIn);
        return strPerms;
    }
    wstring strFirst = strIn.substr(0,1);
    strIn.erase(0,1);
    STRLIST permsSecond = Permutate(strIn);
    for(auto strPermElemsSecond : permsSecond)
    {
        STRLIST permSub = Permutate(strFirst, strPermElemsSecond);
        std::copy(permSub.begin(),permSub.end(),back_insert_iterator<STRLIST>(strPerms));
    }
    

#endif
        return strPerms;
}
/*
 if node X is right child of Y, then all elements in (Y +1,X) in the inorder sequence should be to the right of Y in preorder sequence
 if node X is left child of Y, then all elements in (X,Y) in inorder sequence should be to the left of X in preorder sequence
 
 if left child, then it should precede parent in in-order traversal, if right child, should succeed in in-order, but be lesser than parent
 
 std::queue rootQ;
 rootQ.push(preorder[0]);
 
 while(!rootQ.empty)
 {
 auto root = rootQ.front();
 rootQ.pop;
 auto left = root->getLeft();
 if(left)
 {
 root->left = left;
rootQ.push(left);
 }
 
 auto right = root->right;
 if(right)
 {
 root->right = right;
 rootQ.push(right);
 }
 }
 
 */
/*
        A
   B         C
 D   F     G    E
 
in order D,B,F,A,G,C,E
pre order A,B,D,F,C,G,E
 
      A
   B
     D
 
B,D,A
A,B,D
 
   A
      B
    C
   D
  E
A,E,D,C,B
A,B,C,D,E

    A
   B
  D
D,B,A
A,B,D
 
 root = NUll;
 
 for(cursor : preorder)
 {
   if(root == nullptr)
    root = cursor;
   if(hasleft(cursor)
      cursor->left = *(cursor +1);
   if(hasright(cursor)
     push(cursor)
 
 
 }
 
 if this is right child and its predecessor != parent
 or left child and successor != parent 
 //we know node is right or left child, because when we traverse from root, we push the root 
 //and subsequent item's position relative to root tells us if it is left or right
 bool hasLeftChild( node)
 {
   if(isroot)
    return pos(node,inorder) >0;
 
   if(isLeft(node, parent.top)
   {
     return pos(node, in order) !=0
   }
   return pos (predecesor (node) > pos(parent) ;
 
 
 }
 bool hasrightchild(node)
 {
   if(isroot(node)
  inorder.end() - return pos(node,inorder) > 1;
 
    if(isLeftNode)
    return successor < parent 
   return pos(node,inorder) != inorder.end -1;
 }
 
 
 
 bool is parent(this)
 {
 (leftchild  and succesor = parent  and no pred)
   false; it is root
 
 (right child  and pred == parent and no succesor)
   return false
 
  return true
 }
 
 stack parents;
 auto root = nullptr;
 for (rootcursor : pre-Order)
 {

 
    if(root == nullptr)
  {
        root = rootCursor;
        parents.push( rootcursor);
  }
    else
     {
        rootparent = parents.top ;
        bool bLeft = isLeft(rootparent); // pos child < pos parent in in-order vector
        if(bLeft)
        {
            rootparent->left = rootcursor);
            if rootparent has no right child
            {
                inorder.erase
            }
        }
        else
        {
             rootParent->right = rootcursor;
            inorder.erase rootParent;
            parents.pop;
        }
        if(is parent(rootcursor)) //pos rootcursor in inorder != 0
        parents.push( rootcursor);
        else
        {
          inorder.erase rootCursor
        }
    }
 
 }
               A
        AB            BA
  ABC    ACB CAB   BAC BCA CBA
 ABCD ABDC
 
 ABCD
 ABDC
 ADBC
 DABC
 
 ACBD
 ACDB
 ADCB
 DACB
 
 BAC
 
 BCA
 
 CAB
 
 CBA
 
 123
 132
 213
 231
 312
 321
 
 */
bool isLeftChild(vector<int>& inorder, int elem, int parent)
{
    
    auto elempos = find(inorder.begin(),inorder.end(),elem);
    auto parentpos = find(inorder.begin(),inorder.end(),parent);
  
    return elempos < parentpos;
}

using std::queue;
void BuildTree(vector<int> inorder,vector<int> preorder)
{
    int root =-1;
    std::stack<int> parents;
    for(auto rootCursor : preorder)
    {
        if(root == -1)
        {
            root = rootCursor;
            parents.push(rootCursor);
        }
        else
        {
            int parent = parents.top();
            bool isLeft = isLeftChild(inorder, rootCursor, parent);
            if(isLeft)
            {
                //parent ->left = rootCursor;
            }
            else
            {
                //parent-> right = rootCursor;
                
            }
        }
    }
}

template <typename FwdIter>
pair<FwdIter,FwdIter> Kadane(FwdIter first, FwdIter last)
{
    pair<FwdIter,FwdIter> result = {last,last};
    pair<FwdIter,FwdIter> interimResult = result;
    
    auto sum = 0;
    auto interimSum = 0;
    
    for(auto next = first;next != last;++next)
    {
        
        if(interimSum == 0 && *next >0)
        {
            interimResult.first =   next;
            
        }
        if(interimSum + *next >0)
        {
            interimSum = interimSum + *next;
            interimResult.second = next +1;
            if(interimSum > sum)
            {
                sum = interimSum;
                result = interimResult;
            }
        }
        else /*if(localSum <= 0)*/
        {
            
            interimSum = 0;
            interimResult = {last,last};
        }
    }
    
    return result;
}
template <typename FwdIter>
pair<FwdIter,FwdIter> MaxProfit(FwdIter first, FwdIter last)
{
    pair<FwdIter,FwdIter> result = {first,first};
    pair<FwdIter,FwdIter> interimResult = result;
    
    for(auto next = first;++next != last;++first)
    {
       if(*next > *interimResult.second)
       {
           interimResult.second = next;
       }
       else if( *next < *interimResult.first) // end of subsequence corresponding to local minimum
       {
          if(* interimResult.second - *interimResult.first > *result.second-*result.first)
              result = interimResult;
           
           interimResult.first = interimResult.second = next;
       }
    }

    if(* interimResult.second - *interimResult.first > *result.second-*result.first)
        result = interimResult;

    return result;
}


template <typename FwdIter>
pair<FwdIter,FwdIter> MonoTonic(FwdIter first, FwdIter last)
{
    pair<FwdIter,FwdIter> result = {first,first};
    
    pair<FwdIter,FwdIter> interimResult = {first,first};
    
    for(auto next = first;++next != last;++first)
    {
        if(*next < *first) // found a dip
        {
            interimResult.second = next;
            
            if(std::distance(interimResult.first, interimResult.second) > std::distance(result.first, result.second))
                result = interimResult;
            
            interimResult.first = interimResult.second = next;
            
        }
        else
        {
            interimResult.second = next + 1;
        }
    }
        
    if(std::distance(interimResult.first, interimResult.second) > std::distance(result.first, result.second))
        result = interimResult;

    
    return result;
}
/*
 
 int KadaneSeq[] = {-2, 1, -3, 4, -1, 2, 1, -5, 4,-3,100};
 size_t KIndex[] = {10,9,8,7,6,5,4,3,2,1,0};
 
 

 auto list = {3,2,0,1};
for( auto idx : list)
 {
 autp tmp = list[idx];
 list[idx] = list[list[idx];
 list[list[tmp]] = list[list[list[tmp]]];
 //arr 0 = arr[3], arr[3] = arr[1]
 }

 */
template <typename RandIterSeq, typename RandIterIndex>
void AlignIndex(RandIterSeq first,RandIterSeq last,RandIterIndex firstIndex, RandIterIndex lastIndex)
{
    auto elemCount = std::distance(first ,last);

    for(auto idx = 0;idx != elemCount;idx++)
    {

        while(idx != firstIndex[idx])
        {
            std::swap(first[idx],first[firstIndex[idx]]);
            std::swap(firstIndex[idx],firstIndex[firstIndex[idx]]);

        }

    }
    
}
/*template <typename T>
std::pair<T*,T*> Monotonic(T* begin, T* end)
{
    T* gBegin = begin;
    T* gEnd = begin;
    
    T* lBegin = begin;
    T* lEnd = begin;
    auto iter = begin;
    while (iter != end && iter +1 != end)
    {
        if(*(iter +1) < *iter ) // found break in sequence
        {
            lEnd = iter + 1;
            if(lEnd - lBegin > gEnd - gBegin)
            {
                gBegin = lBegin;
                gEnd = lEnd;
            }
            
            lBegin = lEnd;
            
        }
        else
        {
            lEnd = iter +2;
        }
        
        iter ++;
    }
    if(lEnd - lBegin > gEnd - gBegin)
    {
        gBegin = lBegin;
        gEnd = lEnd;
    }
    
    return std::make_pair(gBegin,gEnd);
}*/

using mystring =std::wstring;
//    wstring result = AddString(L"11001", L"00111");
//    11001
//   000111

template <size_t base>
mystring AddString(mystring lhs,mystring rhs)
{
   
    auto op1rbegin = lhs.rbegin();
    auto op1rEnd   = lhs.rend();
    auto op2rBegin = rhs.rbegin();
    auto op2rEnd = rhs.rend();
    
    if(lhs.length() < rhs.length())
    {
        swap(op1rbegin,op2rBegin);
        swap(op1rEnd,op2rEnd);
        
    }
    auto carryFLag = 0;
    
    mystring result;
    transform(op1rbegin, op1rEnd, back_inserter(result),
              [&carryFLag,&op2rBegin,&op2rEnd](wchar_t bit)
               {
  
                  auto lhsbit = bit - L'0';
                  auto rhsbit = op2rBegin != op2rEnd ? *op2rBegin++ - L'0': 0;
                  auto result = lhsbit + rhsbit + carryFLag;
                  carryFLag = result/base;
                  result %=base;

                  return result + L'0';
              });

    if(carryFLag)
        result += L'1';
    std::reverse(result.begin(),result.end());


  
    return result;
}
template <typename FwdIter>
std::pair<FwdIter,FwdIter> MaxMin(FwdIter first, FwdIter last)
{
    if(first == last)
        return {last,last};
    std::pair<FwdIter,FwdIter> result{first,first};
    
    for(auto next = first+1; next != last;++next)
    {
        if(  *next < *result.first)
        {
            result.first = next;
        }
        else if(*result.second < *next)
        {
            result.second = next;
        }
    }
    
    return result;
}

template <typename T>
pair<T*,T*> MaxMin(T in[],size_t sz)
{
    pair<T*,T*> result = {in,in};
    
    if(sz==0)
        return result;
    
    auto begin = in;
    auto end = in + sz;
    auto next = begin +1;
    
    for(;;)
   {
        
        if(next == end)
        {
           if(*begin < *(result.first))
               result.first = begin;
            else if(*begin > *(result.second))
                result.second = begin;
            return result;
        }
        if(*begin < *next)
        {
            if(*begin < *(result.first))
                result.first = begin;
            if(*next > *(result.second))
                result.second = next;
                
        }
        else
        {
            if(*begin > *(result.second))
                result.second = begin;
            if(*next < *(result.first))
                result.first = next;
        }


        if(++next == end)
            break;
        
        begin = next;
        ++next;
    }
    
    return result;
}
template<class InputIt, class UnaryPredicate>
InputIt find_if_not_x(InputIt first, InputIt last, UnaryPredicate q)
{
    while(first != last)
    {
        if(q(*first))
            break;
        first++;
    }
    return first;
}
auto lamb = [](size_t t) -> bool
{
    return t == 4;
};
template<typename T, typename U>
void PrintRecursive(T begin, T end, T cursor, U output )
{
    if(cursor == end )
        return;
    if(*cursor == L' ')
    {
        *output = *cursor;
    }
}



const wchar_t* text =

       L"RULE 230 of the International Association of Athletics Federations (IAAF) rulebook is not common knowledge to most fans" L"watching the Olympics, but it is likely to be among the most frequently invoked regulation of this Olympiad. Paragraph two of the laws of" L"race walking define the sport as “a progression of steps so taken that the walker makes contact with the ground, so that no visible (to the" L"human eye) loss of contact occurs”. But try walking at any pace and you quickly realise that keeping one foot in touch with the ground is" L"easier said than done. Competitors who break this rule are first cautioned, then eventually disqualified if they repeatedly infringe it.";
std::pair<size_t,size_t> Stocked(const size_t* begin, const size_t* end);

class IUnknown
{
public:
    void Release()
    {
        delete this;
    }
};

template <typename T>
auto RelOnDel = [](T * ptr)
{
    ptr->Release();
    
};

template <typename T>
using RelDeleter = decltype(RelOnDel<T>);
#include <functional>

STRINGS DialPad(const wstring& number );
vector < int > MergeSort(vector < int > intArr);
int findMinimum(vector < int > arr);
#include <array>

void RLE(string input);
int getMaxArea(int hist[], int n);
string moveLetters(string input);
void PartitonInplace(std::vector<int>& input);
vector<size_t> MaxAreaHistogram(vector<size_t>& bars);
void StablePart(std::vector<int>& input);

size_t WaterinBar(vector<size_t> bars);
size_t NapSak(vector<pair<size_t,size_t>> items, size_t Cap);
template<typename T>
class Link
{
public:
    Link(const T & data):m_data(data),m_pNext(nullptr)
    {}
//private:
    T m_data;
    Link<T>* m_pNext;
};
//a->b->c

template<typename T>
Link<T>* Reverse(Link<T>* link,Link<T>*& head )
{

    if(link == nullptr)
        return link;
    
    {
        auto next = link->m_pNext;
        link->m_pNext = nullptr;
        while(next)
        {
            auto cursor = next->m_pNext;
            next->m_pNext = link;
            link = next;
            next = cursor;
        }
    }
    
    auto prefix = Reverse(link->m_pNext,head);
    if(prefix)
    {
        link->m_pNext = prefix->m_pNext;
        prefix->m_pNext = link;
    }
    else
    {   head = link;
    }
    return link;
    
}

template <typename RandIter>
RandIter LowerBound(RandIter begin, RandIter end, const typename iterator_traits<RandIter>::value_type& val)
{
    while (begin!= end)
    {
        auto mid = begin + (end - begin)/2;
        if( *mid < val)
        {
            begin = mid +1;
        }
        else
            end = mid;
    }
    
    return begin;
}

template <typename RandIter, typename Pred>
RandIter Bound(RandIter begin, RandIter end, Pred predicate)
{
    while (begin!= end)
    {
        auto mid = begin + (end - begin)/2;
        if(! predicate(*mid))
        {
            begin = mid +1;
        }
        else
            end = mid;
    }
    
    return begin;
}

pair<Operation,int> ParseOperation(string& op);
void superStack(vector < string > operations);
void GoSpiral(const vector<vector<size_t>>& grid );
void DetectIntersection();
bool IsRotatedPalindrome(string input);
void Permuda(string& input);
void SSet(string&);
void ParaSynthesize(size_t open, size_t close);
void InsertS(IntIter begin, IntIter end);
void NChooseK(const wstring& in, size_t k);
void Find(wstring& strText, wstring& strPat);
void Test();
size_t LCSLite(const wchar_t* lbegin, const wchar_t* lend, const wchar_t* rbegin, const wchar_t* rend);
int main(int argc, const char * argv[])

{
    JustifyText(text,57);
    wstring strText = L"  Hello World Hello";
    wstring pattern = L"Hello";
    auto lcs = LCSLite(strText.c_str(), strText.c_str() + strText.size(), pattern.c_str(),pattern.c_str() + pattern.size());
    
    Test();
    
    Find(strText,pattern);
    vector<size_t> binsort {0,1,2,3,4};
    NChooseK(L"123", 3);
    auto lwrBound = [](size_t elem) ->bool
    {
        return !(elem < 2);
    };
    auto lbu = Bound(binsort.begin(), binsort.end(),lwrBound);
    
    vector<int> srt= {6,7};
    InsertS(srt.begin(),srt.end());
    //vector
    size_t array1[]= {0,3,2,9,4,7,8,1,6};
    //InsertS(array1,array1 + sizeof(array1)/sizeof(array1[0]));
    
    std::sort(array1, array1 + sizeof(array1)/sizeof(array1[0]));
    auto myTree = BuildTreebalanced<size_t,size_t*>(array1,array1 + sizeof(array1)/sizeof(array1[0]));
    
    {
        //ParaSynthesize(4,4);
        //string p("1221");
        //Permuda(p);
        //SSet(p);
        //auto pali = IsRotatedPalindrome("hjijhgi");
        DetectIntersection();
        auto firstList = nullptr;
        for(int i = 0; i!= 5; i++)
        {
        }
        vector<vector<size_t>> spiral =
        {
            {0,1,2},
            {3,4,5},
            {6,7,8},
            {9,10,11}
        };
        GoSpiral(spiral);
        
        vector<string> inputops{
            {"push -36"},
            {"push 20"},
            {"pop "},
            {"push -9"},
            {"pop "},
            {"push -53"},
            {"pop "},
            {"inc 1 -17"}
            };
        superStack(inputops);
        
        
        Link<int>* head = new Link<int>(5);
        head->m_pNext = new Link<int>(4);
        head->m_pNext->m_pNext =new Link<int>(3);
        Link<int>* rev;
        Reverse(head,rev);
        
        vector<pair<size_t,size_t>> items =
        {{20,65},{8,35},{60,245},{55,195},{40,65},{70,150},{85,275},{25,155},
            {30,120},{65,320},{75,75},{10,40},{95,200},{50,100},{40,220},{10,99}
        };
        NapSak(items, 130);
        

        
        vector<int> fb = {2,3,-4,-9,-1,-7,1,-5,-6};
        StablePart(fb);
        vector<size_t> bars ={6,2,5,4,6,5,6};
        auto water = WaterinBar(bars);
        MaxAreaHistogram(bars);
        string letter = "111xyz";
        auto moved =moveLetters(letter);
        
        auto maxsub =  MaxSubVector(fb);
        vector<int> sub(maxsub.first,maxsub.second);
        

        RLE("aaabcccd");
        vector<int> inorder {2, 2, 2, 2, 2, 2, 2, 2, 0,2,2,2};
        auto res =  findMinimum(inorder);
        vector<int> preorder {50,30,10,40,70,60,90};
        //20,65;8,35;60,245;55,195;40,65;70,150;85,275;25,155
        //30,120;65,320;75,75;10,40;95,200;50,100;40,220;10,99
        
        
    }
    
    {
        //vector<string> parans =  GenerateNestedBrackets(4);
        GenerateParanthesis(3);
        
    }
    {
        vector<int> ip{4,2,3,4,5,6,0};
       // FindMissing(ip,3);
        IK::pqueue<int> pq(ip);
    }
    std::vector<size_t> points{9,10,100,0,5,22,1,4};
    auto rangerank = IK::QuickSelect(points.begin(), points.end(), 7);
    
    auto resort= MergeSort(vector < int >{12,4,9,-1});
    
    auto resstr = sortCharacters("bca000011111hello");
    Perms(6);
    std::vector<size_t> vec;
    myrandom_shuffle(vec.begin(),vec.end());
    
    auto powres = PowerX(3, 5);
    auto r = Multfast(13,3);
    {
        unique_ptr<IUnknown,RelDeleter<IUnknown>> ptr(new IUnknown{},RelOnDel<IUnknown>);
        unique_ptr<IUnknown,std::function<void(IUnknown*)>> iunk(new IUnknown{}, [](IUnknown* ptr){if(ptr) ptr->Release();});
    }
    size_t missing[] = {5,4,3,2,4,5,5,5,8};
    Insert_Sort(missing,missing + sizeof(missing)/sizeof(missing[0]));
    
    BiPartition(missing, missing + sizeof(missing)/sizeof(missing[0]));
    
    auto missed =  FindDuplicateOrMissing(missing,missing + sizeof(missing)/sizeof(size_t));

    size_t ages[] = {5,2,8,9,0};
    auto rk = RankIT<size_t>(ages, 5,10);
    
    std::vector<size_t> rangei{320,1200};
    
    std::vector<wstring> ranges{L"320",L"1200",L"1523", L"1999"};
    
    std::sort(rangei.begin(),rangei.end());
    
    std::sort(ranges.begin(),ranges.end());
    
    //std::binary_search(<#_ForwardIterator __first#>, <#_ForwardIterator __last#>, <#const _Tp &__value_#>, <#_Compare __comp#>)
    
    std::vector<size_t> skis{10,5,30,32};
    std::vector<size_t>skier{10,35};
    auto res = Skifall(skis,skier);
    
    
    const wchar_t* lhs = L"zayan";
    const wchar_t* rhs = L"ayan";
    auto maxlcs =  LCS (lhs,lhs+wcslen(lhs),rhs,rhs+wcslen(rhs));
    
    std::vector<std::pair<size_t,size_t>> intervals{{1,2},{2,4}, {0,1},{5,10}};
    MergeIngervals(intervals);
    
    auto den=  GetDenoes(5, std::vector<size_t>{1,3});

    const size_t stockst[] = {10, 7, 5, 8, 11, 9};
    std::pair<size_t,size_t> max = Stocked(stockst, stockst + sizeof(stockst)/sizeof(stockst[0]));
                                      
    
    
    
    size_t rods[]= {0,3,5,10,12,14};
    auto mm = MaxMin(rods,sizeof(rods)/sizeof(rods[0]));
    
    MaxRod(rods,rods + sizeof(rods)/sizeof(rods[0]));
    vector<std::pair<size_t,size_t> >values {{10,5},{40,4},{30,6},{50,3}};
    KnapTheSack(values, 10);;
    KnapsackFull(values, 10);

    
    DialPad(L"0321");
   
    using std::vector;
    vector<size_t> solutions;
    STRINGS input = {L"abc", L"DEF",L"hij",L"klm"};
    vector<string> ju = {"abc", "DEF","hij","klm"};
    Justify(6,ju);
    //getchar();
    //BackTrack(solutions, 0, input[0].size(),  input);
    //AlignIndex
    PowerSet(L"hello");
    GetPossiblePositions(2,4);
    int rotate[] = {0,1,2,3,4};
    Rotate(rotate, rotate + 3,rotate+sizeof(rotate)/sizeof(int));
    
    //auto pt = FindTurningPoint(rotate,sizeof(rotate)/sizeof(int));
    
    int KadaneSeq[] = {-2, 1, -3, 4, -1, 2, 1, -5, 4, -3, 100};
    auto lis = LIS(std::vector<int>({1,2,3,0,0,100,101,102,103,110}));
    auto zig = MaxZigZAg(std::vector<int>({1,2,3,0,0,100,101,102,103,110}));
    
   // auto just = Justify(L"hello world");
    
    
    size_t KIndex[] = {10, 9,  7, 8,  6, 5, 4,  3, 2,  0, 1};
    std::vector<size_t> LinearPart = {1,2,3,4,5,6,7,8,9};
    LinearPartition(LinearPart, 4);
    AlignIndex(KadaneSeq, KadaneSeq + sizeof(KadaneSeq)/sizeof(int),KIndex, KIndex + sizeof(KIndex)/sizeof(KIndex[0]));
    
    auto kresult = MaxProfit(KadaneSeq, KadaneSeq + sizeof(KadaneSeq)/sizeof(int));
    
    
    //Kadane
 
    size_t stocks[] = {200, 200,50,100,200,10,500};
    auto profit = MaxProfit(stocks, stocks + sizeof(stocks)/sizeof(stocks[0]));
    int array[] = {1,1,5,5,5,5,5,5,6,6,6,6,8,9,10,12};
    size_t pivot = QSelect(array, 1, -1);

    
    //LevelWalk(myTree);
    {
        vector<size_t> inorder;
        Inorder(myTree,inorder);
        vector<size_t> preorder;
        Preorder(myTree,preorder);
        TreeFromPreOrder(preorder);
        TreeFromPreAndInOrder(inorder,preorder);
    }
    TRNode<size_t>* head = nullptr;
    //FlattenDL(myTree,head);
    //FlattenSL(myTree,head);
    
    
    bool isBst = IsBST<size_t>(myTree,std::numeric_limits<size_t>::max(),std::numeric_limits<size_t>::min());

    
   // auto its = find_if_not_x< >(array1, array1 + sizeof(array1)/sizeof(array1[0] , lamb   ));
    auto its = find_if_not_x(array1, array1 + sizeof(array1)/sizeof(array1[0]), lamb);
    
    wstring perm(L"tta");
    STRINGS permutes =  PermutateF(perm);
    
    map<wchar_t,wstring> neighBours = { {L'a',L"abcd"},{L'e',L"efgh"},{L'i',L"ijkl"},{L'm',L"mnop"} };
    INPUTLIST il =  DialThis(L"aem",neighBours);
    
    FindMin({0,0}, {3,4});
    
    size_t sizeArray [] = {1,2,3,2,5,6,7,8,9};
    wstring result = AddString<10>(L"9", L"10009");
    
    auto mono = MonoTonic(sizeArray, sizeArray +9);
    
    STRLIST perms = Permutate(L"ABCDEFG");
    vector<int> maxarray = {-1, 3 ,4,-10,50,1000};
    pair<int, int> range = FindMaximumSubarray(maxarray);
    
    
    //std::map<size_t,size_t> myMap;
    //myMap[0] =1;
    //myMap[1] =2;
    using std::sort;
    using size_pair = std::pair<size_t,size_t>;
    vector<size_pair> myMap;
    sort(myMap.begin(),myMap.end(), [](const size_pair& p1, const size_pair& p2)
         {
             return p1.first < p2.first;
         });
    //Tree<> t;
    //Tree<size_t,size_t>::Iterator iter = t.begin();
   // const std::function<bool( size_t,size_t )>& comparator = LessComparator<size_t>;
    size_t array0[]= {0,3,2,9,4,7};
    
    
    
    PriorityQueue<size_t,6> pqueue(array0);
    pqueue.Sort();

    
   // assert(argc == 2);
    
    //MazeMain(argv[1]);
    
    size_t pow = Power(4,5);
    size_t fac = Factorial(5);
    size_t fib = Fib(6);
    
    
    //int array[] = {0,0,0,0,0,0,0};

    

    constexpr size_t arrayLength = sizeof(array)/sizeof(int);
    ReverseInplace(array,arrayLength);
    FindTurningPoint(array,arrayLength);
    //Dijkstra<int>(array,0,arrayLength);
    QSortDjikstra<int>(array,0,arrayLength);
    size_t lb = LowerBound(array,arrayLength,1);
    size_t ub = UpperBound(array,arrayLength,0);
    
    
   // for(;;)
    {
    
    Shuffle<int>(array, arrayLength);
    
    //QSort<int>(array,0,arrayLength);
        
    //InsertionSort(array,arrayLength);
    MergerSorter<int> merge(array,arrayLength);
//    merge.BottomUp();
    merge.Sort();
    //size_t x = BinSearch<int>(array, 7, 7000);
    
    //if(x)
    {
    }
    }
    RShift<int>(array, 2, 3, 5);
    wchar_t charArray[30] = {L"Hello world"};
    
    wchar_t* begin = charArray;
    while(*begin)
    {
        if((*begin) == L' ')
        {
            RShift(charArray, begin - charArray, 2, wcslen(charArray)+1);
            *begin++ = L'%';
            *begin++ = L'2';
            *begin++ = L'0';
        }
        begin ++;
    }
    
    
    
    return 0;
}
template <typename  T >
void Merge(T lhs[],size_t lhsrbegin, T rhs[], size_t rhsrbegin, size_t end)
{
    size_t lhsrend = 0;
    size_t rhsrend = 0;
    
    size_t rbegin = end;
    while(rbegin--)
    {

        if(rhsrbegin && lhsrbegin && (rhs[rhsrbegin-1] > lhs[lhsrbegin-1]))
                rhs[rbegin] = rhs[--rhsrbegin];
        
         else  if(lhsrbegin && rhsrbegin && (lhs[lhsrbegin-1] > rhs[rhsrbegin-1] ))
                rhs[rbegin] = lhs[--lhsrbegin];
        
         else if(lhsrbegin)
             rhs[rbegin] = lhs[--lhsrbegin];
        
        else if(lhsrbegin)
           rhs[rbegin] = lhs[--lhsrbegin];
        
        
    }
}

// need to revise this, use breadth first search algorithm, sort of dijkstra min path

template <typename K = size_t, typename V = size_t>
class XNode
{
    vector<XNode*> FindMin( V& val)
    {
        vector<XNode*> minPath;
        V minVal = m_Val + val;
        
        for (auto aChildNode : m_Children)
        
        {
            V thisval  = minVal;
            vector<XNode*> thisminPath =    aChildNode.FindMin(thisval);

            if(minPath.empty())
            {
                minPath = thisminPath;
                minVal = thisval;
            }
            else if(thisval < minVal)
            {
                minPath = thisminPath;
                minVal = thisval;
                
            }
        }
        val = minVal;
        minPath.push_back(this);
        return minPath;
        
        
    }
private:
    vector<XNode*> m_Children;
    K m_Key;
    V m_Val;
    };
using std::shared_ptr;;

template <typename T>
class TreeNode
{
public:
    TreeNode(T data):m_data(data){}
    
    shared_ptr<T> m_pRight,m_pLeft;
    
private:
    T m_data;
};


/*
   0,1,2,3,4,5,6,7,8,9
 0 0,1,2,3,4,5,6,7,8,9
 1 1,2,3,4,5,6,7,8,9
 2 2,3,4,5,6,7,8,9
 3 3,4,5,6,7,8,9
 4 4,5,6,7,8,9
 5 5,6,7,8,9
 6 6,7,8,9
 7 7,8,9
 8 8,9
 9 9
 */
/*
 increasing sub sequence
 
 */
//
template <typename T>
size_t PartIt(T* in, size_t end)
{
   //auto mid = in +
    auto first = in;
    auto last = first + end;
    auto mid = first + (last - first)/2;
    auto pivotIndex = first;
    
    std:iter_swap(pivotIndex,mid);
    
    for(auto next = first +1;next != last;++next )
    {
        if(*next < *first)
            std::iter_swap(++pivotIndex,next);
    }
    
    std::iter_swap(first,pivotIndex);
    
    return pivotIndex - first;
}

template < typename T>
size_t RankIT(T in[], size_t end, const T& val )
{
    auto pivot = PartIt(in,end);
    
    while(pivot != end && in[pivot] != val)
    {
      if(in[pivot] < val)
      {
          pivot = PartIt<size_t>(in +pivot +1, end -(pivot+1));
      }
      else
      {
          pivot = PartIt(in,pivot);
      }
    }
    return pivot;
}
template <typename T>
struct TXNode
{
    TXNode *m_pLeft, *m_pRight;
    
    T m_data;
};
template <typename T>
TXNode<T> * Successor(TXNode<T> * root, TXNode<T>* node)
{
    TXNode<T>* succPrev = nullptr;
    TXNode<T>* succCurr = root;
    while(succCurr != node)
    {
        if(node->m_data <  succCurr->m_data)
        {
            succPrev = succCurr;
            succCurr = succCurr->m_pLeft;
        }
        else
        {
            succCurr = succCurr->m_pRight;
        }
    }
    return succPrev;
}

template <typename T>
TXNode<T> * Predecessor(TXNode<T> * root, TXNode<T>* node)
{
    TXNode<T>* pred = nullptr;
    TXNode<T>* ancestor = root;
    
    while(ancestor != node)
    {
        if(node->m_data > ancestor->m_data)
        {
            pred = ancestor;
            ancestor = ancestor->m_pRight;
        }
        else
        {
            ancestor = ancestor->m_pLeft;
        }
    }
}
template <typename T>
TXNode<T> * FindMin(TXNode<T> * root)
{
    if(root == nullptr)
        return root;
    
    while(root->m_pLeft)
    {
        root = root->m_pLeft;
    }
    return root;
}

template <typename T>
TXNode<T> * DeleteNode(TXNode<T> * root, T data)
{
    if(root == nullptr)
        return root;
    
    if(data < root->m_data)
        root->m_pLeft =  DeleteNode(root->m_pLeft, data);
    
    else if(data > root->m_data)
        root->m_pRight =  DeleteNode(root->m_pRight, data);
    
    else if(root->m_pLeft && root->m_pRight)
    {
        auto successor = FindMin(root->m_pRight);
        root->m_data = successor->m_data;
        root->m_pRight = DeleteNode(root->m_pRight, successor->m_data);
        
    }
    else if(root->m_pLeft || root->m_pRight)
    {
        auto newroot = root->m_pLeft ? root->m_pLeft : root->m_pRight;
        delete root;
        root = newroot;
    }
    else
    {
        delete root;
        root = nullptr;
    }
    return root;
}

/*
 #include <bits/stdc++.h>
 
 using namespace std;
 class LinkedListNode{
 public:
 int val;
 LinkedListNode *next;
 
 LinkedListNode(int node_value) {
 val = node_value;
 next = NULL;
 }
 };
 
 LinkedListNode* _insert_node_into_singlylinkedlist(LinkedListNode *head, LinkedListNode *tail, int val){
 if(head == NULL) {
 head = new LinkedListNode(val);
 tail = head;
 }
 else {
 LinkedListNode *node = new LinkedListNode(val);
 tail->next = node;
 tail = tail->next;
 }
 return tail;
 }
 */
#include <unordered_map>
class TrieNode
{
    friend class Trie;
public:
    TrieNode():m_bLeaf(false)
    {}
    
private:
    unordered_map<wchar_t,TrieNode*> m_Next;
    bool m_bLeaf;
};
using StrIter = const wstring::const_iterator;
class Trie
{
public:
    Trie():m_root(nullptr)
    {
    }
    void Insert(TrieNode* root,StrIter begin, StrIter end )
    {
        if(begin == end)
            return;
        
        if(!root)
        {
            root = new TrieNode{};
        }
        
        if(!root->m_Next[*begin])
        {
            root->m_Next[*begin] = new TrieNode{};
            
            if(next(begin) == end)
            {
                root->m_Next[*begin]->m_bLeaf = true;
                return;
            }
        }
        Insert(root->m_Next[*begin], next(begin), end);
        
    }
    void LongestPrefix(StrIter begin, StrIter end, wstring& result )
    {
        if(begin == end)
            return;
        
        if(m_root->m_Next[*begin] && m_root->m_Next.size() ==1)
        {
            result.push_back(*begin);
            LongestPrefix(next(begin), end, result);
        }
    }
    wstring LongestPrefix(const wstring& input)
    {
        wstring result;
        LongestPrefix(input.cbegin(),input.cend(), result);
        return result;
    }
    
    void Insert(const wstring& input)
    {
        if(!m_root)
            m_root = new TrieNode{};
        Insert(m_root,input.cbegin(),input. cend());
    }
private:
    TrieNode* m_root;
};
