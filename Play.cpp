//
//  Play.cpp
//  PlayArea
//
//  Created by zayan on 2/13/16.
//  Copyright (c) 2016 Eftiquar. All rights reserved.
//

#include "Play.h"
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <ctime>
//find min path in a grid
// input source pair - x , y , dest pair x, y
//return list of pair



PATH FindMinHelper(POINT start, POINT destination)
{

    PATH result ;
    if(start == destination)
        return {start};
    
    if (destination.second == start.second && destination.first == 1)
    {
        return {start};
    }
    if(destination.first == start.first && destination.second ==1)
    {
        return {start};
    }
    auto minXPrev = FindMinHelper(start,{destination.first > 1 ? destination.first -1 :0, destination.first ? destination.second : destination.second -1});
    auto minYPrev = FindMinHelper(start,{destination.second ? destination.first : destination.first-1,destination.second > 1 ?destination.second-1 :0});
    

    
    if(minXPrev.size() < minYPrev.size())
    {
        if(destination.first == 0 && destination.second ==0)
        {
            minXPrev.push_back(destination);
        }
        else
        {
        minXPrev.push_back({destination.first > 1 ? destination.first -1 :0, destination.first ? destination.second : destination.second -1});
        }
        return minXPrev;
    }
    else
    {
        if(destination.first == 0 && destination.second ==0)
        {
            minXPrev.push_back(destination);
        }
        else
        {
        minYPrev.push_back({destination.second ? destination.first : destination.first-1,destination.second > 1 ?destination.second-1 :0});
        }
        return minYPrev;
    }
    

}


void FindMin(POINT start, POINT destination)
{
    auto path =    FindMinHelper(start,destination);
    path.push_back(destination);
    
    for(auto segment :path)
    {
        std::cout << segment.first << segment.second << L"\n";
    }
}

INPUTLIST DialThis(wstring::iterator cursor,wstring::iterator end,  map<wchar_t,wstring> & neighbours)
{
    if(cursor == end)
    {
        return {L""};
    }
    if(cursor + 1 == end)
    {
        INPUTLIST thisResult;
        //thisResult.push_back({wstring(1,*cursor)});
        std::transform(neighbours[*cursor].begin(),neighbours[*cursor].end(),std::back_inserter(thisResult), [](wchar_t ch){return wstring(1,ch);});
        return thisResult;
        
    }
    INPUTLIST result = DialThis(cursor+1, end, neighbours);
    INPUTLIST thisResult;
    
    for(auto subResultCursor:result)
    {
        //thisResult.push_back(*cursor + subResultCursor);
        std::transform(neighbours[*cursor].begin(),neighbours[*cursor].end(),std::back_inserter(thisResult),
                       [&subResultCursor](wchar_t ch)
        {
            
            return ch + subResultCursor;
        });
    }
    return thisResult;
}
INPUTLIST DialLoop(wstring::iterator cursor,wstring::iterator end,  map<wchar_t,wstring> & neighbours);
INPUTLIST DialThis(wstring input,map<wchar_t,wstring> & neighbours)
{
 //return DialThis(input.begin(),input.end(),neighbours);
    return DialLoop(input.begin(), input.end(), neighbours);

}

/*
 a    d     g
 b    e     h
 c    f     i
 
 000
 001
 002
 
 010
 011
 012
 
 020
 021
 022
 
 100
 101
 102
 
 110
 111
 112
 
 120
 121
 122
 
 200
 201
 202
 
 210
 211
 212
 
 220
 221
 222
 
 
 b b
   c
 ab
 ac
 
 bc
 
 a b c
 abc
 bc
 
 a b
 b c
 c
 
 b c
 c
 
 ab
 ac
 bc
 
 helto
 h  e l
 e  l t
 l
 t
 0 1 2 
 0 1 3
 0 1 4
 
 0 2 3
 0 2 4
 
 0 3 4
 
 1 2  3
   3
   4
 
 1 3 4
 
 
*/


INPUTLIST DialLoop(wstring::iterator cursor,wstring::iterator end,  map<wchar_t,wstring> & neighbours)
{
    INPUTLIST result;
    wstring base(cursor,end);

    auto distanceEnd = neighbours[*cursor].length();
    
    vector< wstring::iterator::difference_type> cursors;
    std::transform(cursor, end, std::back_inserter(cursors),[](wchar_t ch)
                   {
                       return 0;
                   }
                   );
    for(;;)
    {
        wstring thisResult;
        auto distance = 0;
        std::transform(cursor, end, std::back_inserter(thisResult),
                       [&cursors, &distance,&neighbours](wchar_t ch)
                       {
                           
                           return neighbours[ch][cursors[distance++]];
                       });
        
        result.push_back(thisResult);
        
        bool carryFlag = 0;
        for(auto rcursor = cursors.rbegin();rcursor!=cursors.rend();rcursor++)
        {
                if(*rcursor +1 == distanceEnd)
                {
                    *rcursor = 0;
                    carryFlag = true;
                }
                else
                {
                    //carry absorbed, carry on
                    *rcursor +=1;
                    carryFlag = false;
                }
              if(!carryFlag)
                  break;
        }
        //overflow occured, quit
        if(carryFlag)
            break;
    }
    return result;
}


/*
 7 / -3 = -2*-3 + 1
 -2/-3 = 0 , -2 
 
 -2 = -3*0 - 2
    = -3*(0 + 1 - 1) - 2
    =
 a/b = q, 
 a = bq - r
 a = bq +b -b -r
 bq + b -(b+r)
 q(b+1) - (b+r)
 
 r >=0 && r < b 
 if(r <0) 
 a = b(q +1 -1) + r
 a = b(q+1) -b + r
 a = b(q-1) +b +r 
 
 
 */
/*
 123 
 132
 231 >> 213 
 231
 321 >> 312
 321
 
 */
bool PermutateF(wstring::iterator first,wstring::iterator end)
{
    if(first == end)
        return false;
    
    auto lastCursor = end;
    
    if(first == --lastCursor)
        return false;

    
    for(;;)
    {

        auto thisLast = lastCursor;
        if(*--lastCursor < *thisLast)
        {
            auto dipPoint = lastCursor;
            auto endCursor = end;
            //while dip point is not less than endcursor
            while(!(*dipPoint < *--endCursor))
                ;
            std::swap(*dipPoint, *endCursor);
            std::reverse(thisLast,end);
            return true;
        }
        else if(lastCursor == first)
        {
            return false;
        }
    }

}
STRINGS PermutateF(wstring in)
{
    STRINGS res;
    std::sort(in.begin(),in.end());
    
    auto first = in.begin();
    auto last = in.end();

    //base case
    res.push_back(wstring(first,last));
    
   while(PermutateF(first,last))
   {
       res.push_back(wstring(first,last));
   }
    
    return res;
}
std::map<size_t,size_t> GenerateLogVector()
{
    vector<size_t>  logs(32);
    size_t n = 0;
    std::map<size_t,size_t> logMap;
    std::generate(logs.begin(),logs.end(), [&n,&logMap]()
                  {
                      size_t powThis = 1 << n;
                      logMap[powThis] = n;
                      return n++;
                  });
    return logMap;
}

STRINGS PowerSet(wstring in)
{
    STRINGS powerset;
    std::map<size_t,size_t>  logsMap = GenerateLogVector();
    
    size_t end = 1 << in.length();
    size_t strlen = in.length();
    for(size_t begin = 0; begin !=end;++begin)
    {
       wstring thisSubset;
       auto masks = begin;
       while(masks)
       {
           size_t mask = masks & (~ (masks -1));
           thisSubset += in[strlen - 1 - logsMap[mask]];
           masks &= (masks-1);
       }
        std::reverse(thisSubset.begin(),thisSubset.end());
        powerset.emplace_back(thisSubset);
    }
    return powerset;
}
size_t StoI(const wchar_t* in)
{
    size_t result = 0;
    while(*in)
    {
        result*= 10;
        result += *in++ - L'0';
        
    }
    return result;
}
auto constexpr hashMod = 997;
auto constexpr base = 26;
bool RabinCarp(wstring target, wstring in)
{
    size_t targetHash = 0;
    size_t inHash = 0;
    size_t baseMaxPower = 0;
    std::transform(target.begin(),target.end(),in.begin(),in.begin(), [&targetHash,&inHash,&baseMaxPower] (wchar_t t, wchar_t s)
              {
                  baseMaxPower = baseMaxPower ? (baseMaxPower* base) % hashMod : 1;
                  targetHash = (targetHash * baseMaxPower + t - L'0')/hashMod;
                  inHash = (inHash * baseMaxPower + s - L'0')/hashMod;
                  
                  return s;
              }
              );
    for(auto cursor = target.length(); in.length() - cursor < target.length();cursor ++)
    {
        if(targetHash == inHash && !in.compare(cursor - target.length(),cursor,target))
            return cursor;
            
        
    }
      return false;
}

vector<size_t> GetPossiblePositions(size_t npos, int nDays)
{
    vector<size_t> positions;
    
    auto start = 1 << npos;
    positions.push_back(start);
    
    auto next = 0;
    
    for(;;)
    {
        auto rightMost = start & ~ (start -1);
        if(rightMost >1 /*&& rightMost < (1 << 31)*/)
        {
            next |= rightMost >> 1;
            next |= rightMost << 1;
        }
        if(rightMost ==1 )
        {
            next |= rightMost <<1;
        }
        else if(rightMost == (1 <<31))
        {
            next |= rightMost >> 1;
        }
        start &= (start -1);
        
        if(!start)
        {
            positions.push_back(next);
            if(positions.size() == nDays)
                break;
            start = next;
            next = 0;
        }
    }
    
    return positions;
}
using std::vector;

void BuildCandidates(vector<size_t>& solutions,const STRINGS& input, size_t level,  vector<size_t>& candidates)
{
    vector<size_t> vecCandidates;
    vecCandidates.resize(input[0].size());
    for(auto levelC=0;levelC< level-1;levelC++)
    {
        vecCandidates[solutions[levelC]] = 1;
    }
    
    for(auto candidate= 0;candidate < input[0].size();++candidate)
    {
        if(!vecCandidates[candidate])
        {
            candidates.push_back(candidate);
        }
    }
    
    
#if 0
    auto index = 0;
    for(auto character :input[level-1])
    {
        candidates.push_back(index++);
    }
#endif
}
bool IsSolution(vector<size_t>& solutions, size_t level, size_t target, const STRINGS& input)
{
    return level == target;
}
void PrintSolution(vector<size_t> solutions, const STRINGS& input)
{
    wstring strResult;
    auto index = 0;
    std::transform(solutions.begin(),solutions.end(),std::back_inserter(strResult),[&index,&input](size_t ch)
                   {
                       return input[0][ch];
                       
                   });
}
void BackTrack(vector<size_t> solutions, size_t level, size_t target, const STRINGS& input)
{
    vector<size_t> candidates;
    if(IsSolution(solutions,level,target,input))
    {
        PrintSolution(solutions,input);
        solutions.clear();
    }
    
    else
    {
        level ++;
        BuildCandidates(solutions,input,level,candidates);
        solutions.resize(level);

        for(auto candidate:candidates)
        {
            solutions[level-1] = candidate;
            BackTrack(solutions, level, target, input);
        }
        
    }
}
//123
//
bool PermutateQ(std::wstring::iterator first, std::wstring::iterator last)
{
    using FwdIter = std::wstring::iterator;
    if(first == last)
        return false;
    
    FwdIter endCursor = last;
    
    if(first == --endCursor)
        return false;
    
    for(;;)
    {
        auto lastCursor = endCursor;
        if(*--endCursor < *lastCursor)
        {
            auto dipPoint = endCursor;
            auto nextIter = last;

            while(!(*dipPoint < *--nextIter))
                ;
            std::swap(*dipPoint,*nextIter);
            std::reverse(lastCursor,last);
            return true;
        }
        else if(endCursor == first)
            return false;
    }
}

/*
 LIS[j] = max LIS[i] , 0 <=i < j for each sequence[i] < sequence[j]  
          if(seq[j] > seq[i] && lis[j] <= lis[i]
           lis[j] = lis[i] + 1.
 */
vector<int> LIS( const vector<int> & sequence)
{
    if(sequence.size() <2)
        return sequence;
    
    vector<size_t> listIndices(sequence.size(),1);
    vector<size_t> predecessors(sequence.size());
    
    size_t max = 1;
    size_t max_cursor = 0;
    for(auto outerLoop = sequence.cbegin() + 1;outerLoop != sequence.cend();++outerLoop)
        for(auto innerLoop = sequence.cbegin();innerLoop < outerLoop;++innerLoop)
        {
          if(*innerLoop < *outerLoop)
          {
              if(listIndices[innerLoop - sequence.cbegin()] +1 > listIndices[outerLoop - sequence.cbegin()])
              {
                  listIndices[outerLoop - sequence.cbegin()] = listIndices[innerLoop - sequence.cbegin()] +1;
                  predecessors[outerLoop - sequence.cbegin()] = innerLoop - sequence.cbegin();
                  if(listIndices[outerLoop - sequence.cbegin()] > max)
                  {
                      max = listIndices[outerLoop - sequence.cbegin()];
                      max_cursor = outerLoop - sequence.cbegin();
                  }
              }
                  
          }
        }
    
    std::vector<int> result;

     result.push_back(sequence[max_cursor]);
     while (--max !=0)
     {
         result.push_back(sequence[predecessors[max_cursor]]);
         max_cursor = predecessors[max_cursor];
     }
    
    return result;
}
/*
 offsets[zig][i] = max(offsets[zag][j]) + 1, such that sequence[j] < sequence[i] for each j< i
 offsets[zag][i]
 */

vector<int> MaxZigZAg(const vector<int> & sequence)
{
    vector<int> result;
    enum  ZigZag{zig,zag};
    vector<int> offsets[2] =
    {
        vector<int>(sequence.size(),1),
        vector<int>(sequence.size(),1)
        
    };
    
    vector<size_t> predecessors(sequence.size(),0);
    size_t maxLength = 1;
    size_t maxCursor = 0;
    
    if(sequence.size() <2)
        return sequence;
    
    for(auto outerLoop = sequence.cbegin() + 1;outerLoop != sequence.end();outerLoop ++ )
        for(auto innerLoop = sequence.cbegin();innerLoop < outerLoop;innerLoop ++)
        {
            if(*outerLoop < *innerLoop)
            {
                if(offsets[zig][innerLoop-sequence.cbegin()] +1 > offsets[zag][outerLoop - sequence.cbegin()])
                {
                    offsets[zag][outerLoop-sequence.cbegin()] = offsets[zig][innerLoop-sequence.cbegin()] +1;
                    predecessors[outerLoop- sequence.cbegin()] = innerLoop-sequence.cbegin();
                }
            }
            else if(  *innerLoop < *outerLoop)
            {
                if(offsets[zag][innerLoop - sequence.cbegin()] + 1 > offsets[zig][outerLoop - sequence.cbegin()])
                {
                    offsets[zig][outerLoop - sequence.cbegin()] =offsets[zag][innerLoop - sequence.cbegin()] + 1 ;
                    predecessors[outerLoop- sequence.cbegin()] = innerLoop - sequence.cbegin();
                }
            }
            auto maxLocal = offsets[zig][outerLoop - sequence.cbegin()];
            if(maxLocal < offsets[zag][outerLoop- sequence.cbegin()])
                maxLocal = offsets[zag][outerLoop- sequence.cbegin()];
            
            if(maxLength < maxLocal)
            {
                maxLength = maxLocal;;
                maxCursor = outerLoop - sequence.cbegin();
            }
        }

    result.push_back(sequence[maxCursor]);
    while (--maxLength !=0)
    {
        result.push_back(sequence[predecessors[maxCursor]]);
        maxCursor = predecessors[maxCursor];
    }
    
    return result;

}
#include <sstream>


#include <array>

/*
 0000
 
 */

STRINGS DialPad(const wstring& number )
{
    const std::array<const wchar_t*,9> keyPad = {L"0",L"1",L"abc",L"def",L"hij",L"klm",L"nop",L"qrs",L"tuv"};
    std::vector<size_t> indexedLength;

    std::transform(keyPad.cbegin(),keyPad.cend(),std::back_inserter(indexedLength),[](const wchar_t* key)
                   {
                       return wcslen(key);
                   });
    
    std::vector<size_t> cursors(number.length(),0);
    std::vector<size_t> keys;
    std::transform(number.cbegin(), number.cend(),std::back_inserter(keys),
                   [](const wchar_t key)
                   {
                      return  key - L'0';
                   });
    
    STRINGS result;
    
    wstring entry;
    std::vector<size_t> nextCursors;
    
    for(;;)
    {

            for(const auto keyIndex : keys)
            {
              auto  ch =  keyPad[keyIndex][cursors[keyIndex]];
                entry += ch;
            }
        
        result.push_back(entry);
        entry.clear();
        
           size_t carry = 1;
           size_t index = keys.size();
           std::transform(cursors.rbegin(),cursors.rend(),std::back_inserter(nextCursors),
                       [&indexedLength,&carry,&index,&keys](const size_t cursor) ->size_t
                       {
                           --index;
                           if(carry)
                           {
                               if(cursor + 1 == indexedLength[keys[index]] || cursor == indexedLength[index])
                               {
                                   carry =1;
                                   return 0;
                               }

                               else
                               {
                                   carry = 0;
                                   if(cursor == indexedLength[index])
                                       carry =1;
                                   return cursor == indexedLength[index] ? 0 : cursor +1;
                               }
                               //return cursor +1;
                           }
                           return cursor;
                       });

        if(carry) // overflowed
            break;
        
        std::reverse(nextCursors.begin(),nextCursors.end());
        cursors = nextCursors;
        nextCursors.clear();


    }
    
    return result;
}
/*
 length	0	1	2	3	4	5
 pi	    0	3	5	10	12	14
 ri	0	/	/	/	/	/
 si	/	/	/	/	/	/
 */

//Rev[i]

size_t MaxRod(const size_t* begin, const size_t* end)
{
    size_t result = 0;
    auto elemcount = end - begin;
    if (elemcount == 0 )
        return result;

    std::vector<size_t> partitions(elemcount);
    
    std::vector<size_t> revenue(elemcount);
    revenue[0] = begin[0];

    
    for(auto outer = 1;outer < elemcount;outer++)
    {
        result = begin[outer];
        partitions[outer] = outer;
        
        for(auto inner = 1; inner < outer ;inner++)
        {
            if(result < begin[inner] + revenue[outer - inner])
            {
                result = begin[inner] + revenue[outer - inner];
      
                partitions[outer] = inner;
            }
        }
        revenue[outer] = result;
        
    }
    std::vector<size_t> solution;
    size_t length = elemcount-1;
    while(length)
    {
        solution.push_back(partitions[length]);
        length -=partitions[length];
    }
    return  revenue[elemcount-1];
}
size_t KnapTheSack(vector<std::pair<size_t,size_t> >values, size_t capacity)
{
    
    vector<size_t> knapSackMatrix(capacity + 1);
    
    //outer loop is
    for(size_t item = 0;item != values.size();item ++  )
    {
        //need to walk backwards else we will count same item twice
        //for(size_t capacityCursor = 0;capacityCursor < capacity + 1;++capacityCursor)
        for(size_t capacityCursor = capacity;capacityCursor >=values[item].second;--capacityCursor )
        {
            knapSackMatrix[capacityCursor] = std::max(knapSackMatrix[capacityCursor-values[item].second] + values[item].first,knapSackMatrix[capacityCursor]);
            
        }

    }
    return  knapSackMatrix[capacity];
}
/*
size_t KnapTheSack(vector<std::pair<size_t,size_t> >values, size_t capacity)
{

    vector<vector<size_t>> knapSackMatrix(values.size(),vector<size_t>(capacity +1,0) );
    
    for(size_t item = 0;item < values.size();item ++  )
        for(size_t capacityCursor = 1;capacityCursor <= capacity;capacityCursor++ )
    {
        size_t localMax = 0;
        if(values[item].second <= capacityCursor)
        {
            //take this item
            size_t maxWOthisItem = item ? knapSackMatrix[ item-1][capacityCursor] : 0;
            size_t maxWthisItem =  item ? knapSackMatrix[item-1][capacityCursor- values[item].second] + values[item].first :values[item].first ;
            localMax = maxWOthisItem > maxWthisItem ? maxWOthisItem : maxWthisItem;
        }
        else
            localMax = knapSackMatrix[item >0 ? item-1:0][capacityCursor];
        
        knapSackMatrix[item][capacityCursor] = localMax;

    }
    return  knapSackMatrix[values.size()-1][capacity];
}
*/
/*
 c[j] = c[i-1] + lc[i,j]; if lc[i,j] != inf
        = 0 if j == n 
 
 */

std::vector<wstring> JustifyText(const wstring& text, int lineWidth)
{
    std::vector<wstring> result;
    std::wistringstream input(text);
    std::vector<wstring> tokens;
    
    wstring token;
    while(input >> token)
        tokens.push_back(token);

    size_t maxTokensPerLine = lineWidth %2 ? lineWidth/2 +1: lineWidth/2;
    std::vector<std::vector<int>> lineCost(tokens.size(),std::vector<int>(maxTokensPerLine,-1));
    
    for(int tokenIndex = 0;tokenIndex != tokens.size();++tokenIndex)
    {
        int spaceRemaining = lineWidth;
        for(size_t tokenCursor = tokenIndex; spaceRemaining > 0;tokenCursor++)
        {
            int thisCost = static_cast<int>(spaceRemaining - tokens[tokenCursor].size() -1);
            
            if(thisCost == -1) // gobble space, we found exact fit
                thisCost = 0;
            if(tokenCursor == tokens.size()-1)// we are on last line
                thisCost = 0;
            if(thisCost >=0)
            {
                lineCost[tokenIndex][tokenCursor-tokenIndex] = thisCost*thisCost;
                spaceRemaining = thisCost;
            }
            else break;
        }
    }
    std::vector<int> aggregateCost(tokens.size(),-1);
   
    std::vector<int> preds(tokens.size(),-1);
    //walk from each token, backwards
    for(int tokenIndex = 0;tokenIndex != tokens.size();++tokenIndex)
    {
        int minCost = -1;
        for(int tokenCursor = tokenIndex;tokenCursor>=0;--tokenCursor)
        {
            int delta = tokenIndex - tokenCursor;
            // cost upto tokencursor -1
            int prevCost = tokenCursor ? aggregateCost[tokenCursor -1]: 0;
            //line from tokencursor to tokenindex
            int thislineCost = lineCost[tokenCursor][delta];
            
            if(thislineCost < 0)
                 break;

            int thisCost = prevCost + thislineCost;
            if(minCost == -1 || thisCost < minCost)
            {
                minCost = thisCost;
                //predecessor is start of line at token cursor which ends at token index
                preds[tokenIndex] = tokenCursor;
            }
        }
        aggregateCost[tokenIndex] = minCost;
    }
    for(size_t lineIndex = preds.size()-1;;)
    {
        size_t lineStart = preds[lineIndex];
           
        wstring line;
        for(;lineStart < lineIndex +1; lineStart++)
        {
            line += tokens[lineStart];
            if(lineStart != lineIndex)
                line += L" ";
        }
        result.push_back(line);
        
        if(!preds[lineIndex])
            break;
        lineIndex = preds[lineIndex] - 1;
    }
    std::reverse(result.begin(), result.end());
    return result;
    
}
    
/*
 
 shuffle
 grufle
 
    g  r  u f l e
 s  0  0  0 0 0 0
 h  0  0  0 0 0 0
 u  0  0  1 1 1 1
 f  0  0  1 2 2 2
 f  0  0  1 2 2 2
 l  0  0  1 2 3 3
 e  0  0  1 2 3 4
 
 */
/*
 rod of length L 
 rod 1 = 5
 rod  2 = 7
 rod  3 = 8
 
 5 
 1,4
 2,3
 3,2
 4,1
 
 */
template <size_t M, size_t N>
size_t MaxGrid(size_t grid[M][N])
{
    size_t solutionGrid[M][N] = grid;
    
    for(auto ycursor = N;ycursor;)
    {
        --ycursor;
        for(auto xcursor = M;xcursor;)
        {
            --xcursor;
            size_t right = xcursor + 1 == M ? 0 : grid[xcursor + 1][ycursor];
            size_t bottom = ycursor + 1 == N ? 0 : grid[xcursor][ycursor +1];
            solutionGrid[xcursor][ycursor] = std::max(right,bottom) + grid[xcursor][ycursor];
        }
    }
    return solutionGrid[0][0];
}
size_t LCS (const wchar_t* lhsbegin,const wchar_t* lhsend,const wchar_t* rhsbegin, const wchar_t* rhsend )
{
    
    auto lhsElems = lhsend - lhsbegin;
    auto rhsElems = rhsend -rhsbegin;
    
    if(lhsElems == 0 || rhsElems ==0)
        return 0;
    
    std::vector<size_t> result;
    
    std::vector<std::vector<size_t>> lcsMatrix(lhsElems,std::vector<size_t>(rhsElems));
    
    for(size_t lhsCursor = 0; lhsCursor != lhsElems; ++lhsCursor)
        for(size_t rhsCursor =0; rhsCursor != rhsElems; ++rhsCursor)
        {
            if(lhsbegin[lhsCursor] == rhsbegin[rhsCursor])
            {
                size_t prevLength = (rhsCursor && lhsCursor) ?  lcsMatrix[lhsCursor -1][rhsCursor -1] : 0;
                lcsMatrix[lhsCursor][rhsCursor] = prevLength +1;
            }
            else
            {
                size_t left = rhsCursor ? lcsMatrix[lhsCursor][rhsCursor -1]:0;
                size_t top  = lhsCursor ? lcsMatrix[lhsCursor -1][rhsCursor]:0;
                size_t max = std::max(left,top);
                lcsMatrix[lhsCursor][rhsCursor] = max;
                
            }
        }
    
    std::wstring lcs;
    size_t rhsCursor = rhsElems;
    size_t lhsCursor = lhsElems;
    
    while(rhsCursor && lhsCursor && lcsMatrix[lhsCursor-1][rhsCursor-1])
    {
        size_t rhsElem = rhsCursor -1;
        size_t lhsElem = lhsCursor -1;
        
        size_t thisElem = lcsMatrix[lhsElem][rhsElem];
        size_t topElem = lhsElem ? lcsMatrix[lhsElem -1][rhsElem]:0;
        size_t leftElem = rhsElem ? lcsMatrix[lhsElem][rhsElem-1]:0;
        
        if(thisElem == topElem)
            lhsCursor--;
        else if(thisElem == leftElem)
            rhsCursor--;
        else
        {
            lcs.push_back(lhsbegin[lhsElem]);
            lhsCursor--;
            rhsCursor--;
        }
    }
    
    std::reverse(lcs.begin(),lcs.end());
                 
    return lcsMatrix[lhsElems-1][rhsElems-1];
    
#if 0
    auto lhsElems = lhsend - lhsbegin;
    auto rhsElems = rhsend -rhsbegin;
    
    if(lhsElems == 0 || rhsElems ==0)
        return 0;
    
    if( rhsElems > lhsElems)
    {
        std::swap(lhsbegin, rhsbegin);
        std::swap(lhsend, rhsend);
        std::swap(lhsElems, rhsElems);
    }
    //previous == lhsindex
    vector<size_t> previousIndices(lhsElems);

    for(auto lhsIndex = 0;lhsIndex != lhsElems;++lhsIndex)
    {
        //current = rhs
        vector<size_t> currentIndices(rhsElems);
        for(auto rhsIndex = 0;rhsIndex != rhsElems;++rhsIndex)
        {
            if(lhsbegin[lhsIndex] == rhsbegin[rhsIndex])
            {
                currentIndices[rhsIndex] = rhsIndex ? previousIndices[rhsIndex-1] + 1 : 1;
            }
            else
            {
                auto zig = previousIndices[rhsIndex];
                auto zag = rhsIndex ? currentIndices[rhsIndex-1]:0;
                currentIndices[rhsIndex] = std::max(zig,zag);
            }
        }
        std::swap(previousIndices,currentIndices);
    }
        
    return  previousIndices[rhsElems-1];
#endif
}

using std::pair;
pair<size_t,size_t> Stocked(const size_t* begin, const size_t* end)
{

    if(begin == end)
        return pair<size_t,size_t>{0,0};

    pair<size_t,size_t> result{*begin,*begin};
    pair<size_t,size_t> interimresult{result};
    
    for(auto next = begin +1;next != end; ++next)
    {
        if(*next > interimresult.second)
            interimresult.second = *next;
        else if(*next < interimresult.first)
        {
            if(interimresult.second - interimresult.first >
               result.second - result.first)
                result = interimresult;
            interimresult = pair<size_t,size_t>(*next,*next);
        }
        
    }
    if(interimresult.second - interimresult.first >
       result.second - result.first)
        result = interimresult;
    
    return result;
}
/*
 
 1 -  {1}
 2  - {2},  {1,1}
 3 -  {3},  {1,2}, {1,1,1}
 4 -  {1,3},{2,2}, {1,1,2}, {1,1,1,1}

 A(0) = 1
 A(n) - ways to change ncents using 1c coin
 B(n) - ways to change n cents using 1c and 2c coins -> An + B(n-2)
 */

size_t GetDenoes(const size_t targetAmount, const std::vector<size_t>& denominations)
{
    size_t result =1;
    if(targetAmount ==0)
        return result;
        
    std::vector<size_t> subSolutions(targetAmount +1);
    subSolutions[0] = 1;
    for(const size_t denom: denominations)
    {
        for(size_t amtCursor = denom; amtCursor < targetAmount +1;++amtCursor )
        subSolutions[amtCursor] += subSolutions[amtCursor - denom ];
    }
    return subSolutions[targetAmount];
}

std::vector<size_t> KnapsackFull(std::vector<std::pair<size_t,size_t>> items,size_t capacity )
{
    std::vector<size_t> result;
    struct itemcost
    {
        itemcost(size_t value = 0, bool keep= false):
        m_value(value),  m_keep(keep)
        {
        }
        size_t m_value;

        bool m_keep;
    };
    
    std::vector<std::vector<itemcost>> costMatrix(capacity +1,std::vector<itemcost>(items.size()));
    for(size_t nItem = 0; nItem < items.size();++nItem)
    {
        
        for(size_t capacityCursor = capacity; items[nItem].second <= capacityCursor; --capacityCursor)
        {
            size_t valueWItem = costMatrix[capacityCursor- items[nItem].second][nItem ? nItem-1:nItem].m_value + items[nItem].first;
            size_t valueWOItem = nItem ? costMatrix[capacityCursor][nItem-1].m_value :0;
            if(valueWItem > valueWOItem)
            {
                costMatrix[capacityCursor][nItem] = itemcost{valueWItem,true};
            }
            else
            {
                costMatrix[capacityCursor][nItem] = itemcost{valueWOItem,false};
            }
            
        }
    }
    

    size_t ncapactiy = capacity;
    
    for(size_t nItem = items.size(); nItem;--nItem)
    {
        //if(costMatrix[ncapactiy][nItem -1].m_keep)
        
        size_t thisItem = nItem-1;
        size_t valWOthisItem = thisItem ? costMatrix[ncapactiy][thisItem-1].m_value : 0;
        if(costMatrix[ncapactiy][thisItem].m_value != valWOthisItem)
        
        {
            result.push_back(items[thisItem].second);
            ncapactiy-= items[nItem-1].second;
        }
    }
    
    return  result;
}

std::vector<std::pair<size_t,size_t>> MergeIngervals(std::vector<std::pair<size_t,size_t>> intervals)
{
    std::vector<std::pair<size_t,size_t>> result;
    if(intervals.size() == 0)
        throw L"invalid arg";
    
    std::sort(intervals.begin(),intervals.end(),[](std::pair<size_t,size_t> left, std::pair<size_t,size_t> right)
              {return left.first< right.first;});
    

    
    std::pair<size_t,size_t> currentInterval = intervals[0];
    
    using intervalIter = std::vector<std::pair<size_t,size_t>>::const_iterator;
    
    for (intervalIter next = intervals.cbegin(); ++next != intervals.cend(); )
    {
         if(next->first <= currentInterval.second)
         {
             if(next ->second > currentInterval.second)
                 currentInterval.second = next->second;
         }
        else
        {
            result.push_back(currentInterval);
            currentInterval = *next;
        }
        
    }
    
    result.push_back(currentInterval);
    return result;
}

/*
 
    s0 s1 h0 h1
 
        skier
      0 1 2 3
 ski
  -1  i i i i
   0  x
   1  0 x
   2  0 0 x
   3  0 0 0 x
   4
 
 */
// skiers < skis else solution is trivial, just sort the two and allot each pair in order - ski 0 -> skier 0, ski 1 - skier 1 and so on
// skicost [i,j] = min (skicost[i-1,j], skicost[i-1,j-1] + |i-j|)

std::vector<std::pair<size_t, size_t>> Skifall( std::vector<size_t>& skis, std::vector<size_t> skiers)
{
    
    std::vector<std::pair<size_t, size_t>> result;
    std::sort(skis.begin(), skis.end());
    std::sort(skiers.begin(),skiers.end());
    
    std::vector<size_t> baseResult;
    size_t accumulate = 0;
    std::transform(skiers.cbegin(), skiers.cend(),skis.cbegin(), std::back_inserter(baseResult),[&accumulate] (size_t skier, size_t ski)
                   {
                       return accumulate += (ski < skier ? skier -ski : ski - skier);
                   });
                   
                   
    
    std::vector<std::vector<size_t>> costMatrix(skis.size(), std::vector<size_t>(skiers.size()));
    
    costMatrix[0] = baseResult;
    
    for(size_t ski = 1;ski != skis.size(); ++ski)
    {
        for(size_t skier =0; skier != skiers.size();++skier)
        {
            auto costWOthisSki = costMatrix[ ski-1][skier];
            auto thisCost = skiers[skier] < skis[ski] ? skis[ski] - skiers[skier] : skiers[skier] - skis[ski];
            size_t costtWthisSki = 0;
            if( skier > 0)
            {
                costtWthisSki =  costMatrix[ski-1][skier-1] + thisCost;
            }
            else
            {
                costtWthisSki =  std::min(thisCost,costMatrix[ski-1][skier/*0*/]);
            }
            costMatrix[ski][skier] = std::min(costWOthisSki,costtWthisSki);
        }
    }
    
    size_t ski = skis.size();
    for(size_t skierCursor = skiers.size(); skierCursor!=0;)
    {
        auto thisSkier = skierCursor-1;
        auto thisSki = ski -1;
        auto topVal = thisSki ? costMatrix[thisSki-1][thisSkier]:costMatrix[thisSki][thisSkier];
        //auto topLeft = costMatrix[thisSki ? thisSki-1:thisSki][thisSkier ? thisSkier-1:0];
        auto thisVal = costMatrix[thisSki][thisSkier];
        if(thisVal == topVal)
        {
            --ski;
        }
        else
        {
            result.push_back(std::pair<size_t, size_t>(thisSkier,thisSki));
            --ski;
            --skierCursor;

        }
        
    }
    return result;
}


/*
 LCS fast - only compute the length
 
 */

size_t LCSFast(const wchar_t* lbegin, const wchar_t* lend, const wchar_t* rbegin, const wchar_t* rend)
{

    auto lengthLeft = lend - lbegin;
    auto lengthRight = rend- rbegin;

    if(lengthRight ==0 || lengthLeft ==0)
        return 0;
    
    auto minLength = std::min(lengthLeft,lengthRight);
    const wchar_t* xbegin = nullptr;
    const wchar_t* xend = nullptr;
    const wchar_t* ybegin = nullptr;
    const wchar_t* yend = nullptr;
    
    size_t xlength = 0;
    size_t ylength = 0;
    if(minLength == lengthRight)
    {
        xbegin = rbegin;
        xend = rend;
        ybegin = lbegin;
        yend = lend;
        
        xlength = lengthRight;
        ylength = lengthLeft;
        
    }
    else
    {
        xbegin = lbegin;
        xend = lend;
        ybegin = rbegin;
        yend = rend;
        
        xlength = lengthLeft;
        ylength = lengthRight;
    }
    
    std::vector<size_t> lcsyCost(xlength);
    for(auto ycursor = 0;ycursor != ylength;++ycursor)
    {
        std::vector<size_t> lcsxCost(xlength);
        for(auto xcursor = 0; xcursor != xlength;++xcursor)
        {
            auto prevXcost = xcursor ? lcsxCost[xcursor-1] :0;
            auto prevYcost = ycursor ?lcsyCost[xcursor] :0;
            if(xbegin[xcursor] == ybegin[ycursor])
            {
                lcsxCost[xcursor] = std::max(prevXcost,prevYcost) + 1;
            }
            else
            {
                lcsxCost[xcursor] = std::max(prevXcost,prevYcost);
            }
        }
        lcsyCost = std::move(lcsxCost);
    }
    
    return lcsyCost[xlength-1];
}
template<typename randomiter>
size_t PartitionRandom(randomiter base, size_t left, size_t right)
{
    auto pivot = left + (right - left)/2;
    std::swap(base[left],base[pivot]);
    pivot = left;
    for(auto next = left +1;next != right; ++ next)
    {
        if(base[next] < base[left] )
            std::swap(base[next],base[++pivot]);
    }
    std::swap(base[left],base[pivot]);
    
    return  pivot;
}
/*
 0,1,2,3,4,5
 0,1,1,2,3,4
 
 0,1,1 >> dup 
 
 3 = 2 + missing - dup
 
 miising - dup = -1 
 mising = dup -1 
 dup = missing + 1
 
 0,1,3 >> missing
 
 */
template<typename randomiter >
size_t FindDuplicateOrMissing(randomiter begin, randomiter end)
{
    
    auto elemcount = end - begin;
    size_t first = 0;
    auto last = elemcount;
    while(first != last)
    {
        auto pivot = PartitionRandom(begin, first,last);
        if(begin[pivot] != pivot)
        {
            last = pivot;
        }
        else
        {
            first = pivot + 1;
        }
    }
    return first;
}
size_t FindDuplicateOrMissing(size_t* begin, size_t* end)
{
    return FindDuplicateOrMissing<size_t*>(begin, end);
}
template <class RandIter>
void BiPartition(RandIter first, RandIter last)
{
    if(first == last)
        return;
    auto pivot = first + (last - first)/2;
    std::iter_swap(first,  pivot);
    pivot = first;
    auto next = first;
    
    while (next != last)
    {
        while(++next != last)
        {
            if(*next > *pivot)
                break;
        
        }
        while(last != next)
        {
            if(*--last < *pivot)
                break;
        }
        if(last != next)
        {
            std::iter_swap(last, next);
        }
    }
   
    std::iter_swap(--next, pivot);
   
}

void BiPartition(size_t* begin, size_t* end)
{
    BiPartition<size_t*>(begin, end);
}

void MakeRandom(size_t* begin, size_t* end)
{
    
    std::srand(std::time(0)); // use current time as seed for random generator
    auto  random_variable = std::rand();
    auto elemCount = std::distance(begin,end);
    
    for(;elemCount > 1;--elemCount )
    {
        
    }
}

size_t Multfast(size_t multiplier, size_t multiplicand)
{
    size_t result = 0;
    
    while(true)
    {
        if(multiplier %2 ==1)
        {
            
            result += multiplicand;
        }
        multiplier /= 2;
        
        if(multiplier ==0)
            return result;
        
        multiplicand += multiplicand;
   
    }
    
    return result;
}

size_t Power_Accumulate(size_t base, size_t exp, size_t result)
{
    
    while(true)
    {
        if(exp %2 ==1)
        {
            result *= base;
        }
        exp /=2;
        if(exp == 0) return result;
        base *= base;
    }
}
size_t PowerX(size_t base, size_t exp)
{
    if(exp ==0) return 1;
    while(exp %2 == 0)
        {
            base *= base;
            exp /=2;
        }
    if(exp ==1) return base;
    
    exp /=2;
    size_t result = base;
    
    base *=base;
    
    return  Power_Accumulate(base,exp,result);
}



template< class RandomIt >
void random_shuffle( RandomIt first, RandomIt last )
{
    typename std::iterator_traits<RandomIt>::difference_type i, n;
    n = last - first;
    for (i = n-1; i > 0; --i) {
        using std::swap;
        swap(first[i], first[std::rand() % (i+1)]);
        // rand() % (i+1) isn't actually correct, because the generated number
        // is not uniformly distributed for most values of i. A correct implementation
        // will need to essentially reimplement C++11 std::uniform_int_distribution,
        // which is beyond the scope of this example.
    }
}

int FindDuplicate(std::__1::vector<int> input)
{
    typedef typename std::__1::vector<int>::difference_type diff_t ;
    diff_t first, last;
    first = 0;
    last = input.size();
    
    while(first != last)
    {
        diff_t mid = first + (last - first)/2;
        size_t count = 0;
        for(auto elem :input)
        {
            if(elem >= first && elem < mid)
                ++count;
        }
        if(count == mid - first)
        {
            first = mid;
        }
        else
        {
            last = mid;
        }
            
    }
    return first;
}
template <typename T>
T exchange(T& l, T&&r)
{
    T temp = l;
    l = std::forward<T>(r);
    return temp;
}

size_t Fibonacci(size_t in)
{
    
    if(in < 2)
        return in;
    
    size_t fib0 = 0, fib1 = 1;
    
    for(size_t first = 0;first != (in -2 +1); ++ first)
    {
        fib0 = exchange(fib1, fib0 + fib1);
    }
    return fib1;
}

template < typename RandIter>
void Insert_SortT(RandIter first, RandIter last)
{
    typename std::iterator_traits<RandIter>::difference_type elems, indexOuter,indexInner ;
    elems = last - first;
    
    for( indexOuter = 1;indexOuter < elems; ++indexOuter)
    {
        auto newElem = first[indexOuter];
        for(indexInner = indexOuter; indexInner >0 && newElem < first[indexInner -1];--indexInner)
            first[indexInner] = first[indexInner -1];
        first[indexInner] = newElem;
    }
}

void Insert_Sort(size_t* begin, size_t* end)
{
    Insert_SortT(begin, end);
}

template <typename FwdIter>
FwdIter Unique(FwdIter first, FwdIter last)
{
    if(first == last)
        return first;
    
    for(FwdIter next = first +1; next != last; ++next)
    {
        if(*first != *next && ++first != next)
        {
            *first = std::move(*next);
        }
    }
    return ++first;
}
template <typename FwdIter, typename V>
FwdIter Remove(FwdIter first, FwdIter last,const V & v )
{
    FwdIter result = std::find(first,last,v);
    
    if(result != last)
        for(FwdIter next = result +1; next != last; ++next)
        {
            if(!(*next == v))
            {
                *result++ = std::move(*next);
            }
        }
    return result;
}
std::vector<std::vector<size_t>> Perms(size_t len)
{
    std::vector<std::vector<size_t>> permresult;
    for(size_t count = 0; count != len; ++count)
    {
        if(count > 0)
        {
            std::vector<std::vector<size_t>> res;
            for( auto seedcopy:permresult)
            {
                for(size_t thisperm = 0; thisperm != seedcopy.size()+1;++thisperm)
                {
                    std::vector<size_t> interim_result;
                    using std::transform;
                    transform(seedcopy.begin(), seedcopy.end(), std::back_inserter(interim_result),
                              [&thisperm](size_t elem)
                              {
                                  return elem < (thisperm + 1) ? elem : elem +1;
                              });
                    interim_result.push_back(thisperm+1);
                    res.push_back(interim_result);
                }
            }
            permresult = res;
        }
        else
        {
            permresult.push_back({count+1});
        }
    }
    return permresult;
}

// 1,2,3,4
//  2,5

template <typename FwdIter, typename OutputIter>
void Unionize (FwdIter first1, FwdIter last1,FwdIter first2, FwdIter last2, OutputIter out)
{
    while(first1 != last1)
    {
        if(first2==last2)
        {
            std::copy(first1,last1, out);
            return;
        }
        if(*first2 < *first1)
        {
            *out++ = *first2++;
        }
        else
        {
            *out++ = *first1++;
            if(!(*first1 < first2))
               {
                   ++first2;
               }
        }
        
    }
    std::copy(first2,last2, out);
}

template <typename FwdIter, typename OutputIter>
void InterSect (FwdIter first1, FwdIter last1,FwdIter first2, FwdIter last2, OutputIter out)
{
    while(first1 != last1)
    {
        if(first2==last2)
        {
            return;
        }
        if(*first2 < *first1)
        {
           ++first2;
        }
        else
        {
            if(!(*first1 < first2))
            {
                *out++ = *first1++;
                ++first2;
            }
            else
            {
                ++first1;
            }
            
        }
        
    }
 
}

template <typename FwdIter, typename OutputIter>
void DIfference (FwdIter first1, FwdIter last1,FwdIter first2, FwdIter last2, OutputIter out)
{
    while(first1 != last1)
    {
        if(first2==last2)
        {
            std::copy(first1, last1, out);
            return;
        }
        
        if(*first1 < *first2)
        {
            *out++ = *first1++;
        }
        else
        {
            if(!(*first2 < *first1))
            {
                ++first1;
            }
            ++first2;
        }
    }
    
}

template <typename NutIter, typename BoltIter>
void NutsNBolts(NutIter firstNut, NutIter lastNut,BoltIter firstBolt, BoltIter lastBolt  )
{
    if(firstNut == lastNut)
        return;
    
    if(std::distance(firstNut, lastNut) != std::distance(firstBolt, lastBolt))
        return;
    
    
    
    auto midNut1 = std::next(firstNut, std::distance(firstNut, lastNut)/2);
    auto midBolt1 = std::partition(firstBolt, lastBolt,
                                   [&midNut1](const typename std::iterator_traits<BoltIter>::value_type& bolt)
                                   {
                                       return bolt < *midNut1;
                                   });
    auto midBolt2 = std::partition(midBolt1, lastBolt,
                                   [&midNut1](const typename std::iterator_traits<BoltIter>::value_type& bolt)
                                    {
                                        return !(*midNut1 < bolt );
                                    });
    
    midNut1 = std::partition(firstNut, lastNut,
                             [&midBolt1](const typename std::iterator_traits<NutIter>::value_type& nut)
                             {
                                 return nut < *midBolt1;
                             });
    
    auto midNut2  = std::partition(midNut1, lastNut,
                                   [&midBolt1](const typename std::iterator_traits<NutIter>::value_type& nut)
                                   {
                                       return !(*midBolt1 < nut);
                                   });
    NutsNBolts(firstNut, midNut1, firstBolt, midBolt1);
    NutsNBolts(midNut2,lastNut,midBolt2,lastBolt);
    
}


template <typename FwdIter>
size_t Rank(FwdIter begin, FwdIter end, size_t rank)
{
    
    auto pivot = begin + (end-begin)/2;
    std::iter_swap(begin, begin + pivot);
    pivot = begin;
    for(auto next = begin+1; next != end;next ++)
    {
        if(*next < *begin)
        {
            std::iter_swap(++pivot,next);;
        }
    }
    std::iter_swap(pivot,begin);
    return pivot - begin;
}

template <typename FwdIter>
FwdIter Rank(FwdIter begin, FwdIter end, size_t rank)
{
    if(begin == end)
        return begin;
    
}
template< typename T, typename outIter>
void MergeSort(T begin, T end, outIter out)
{
    auto mid = begin + (end - begin)/2;
    
    if(begin == end)
        return;
    
}

template <typename FwdIter>
void Heapify(FwdIter begin, FwdIter end)
{
    auto distance = end - begin;
    auto thisRoot = distance/2;
    
}

//heap sort stable
//256 GB with 64 GB
// 2 sets - is one a s=ubset of other
// home work - given a list of documents and a list of words with their positions  in each of the documents, rank these docs based on the distance between //words w1 and w2 in the doc
//
constexpr unsigned char ASCIIMAX = 255;
using std::string;

string sortCharacters(string inString) {
    
    vector<int> hash(256);
    for(auto thisChar : inString)
    {
        unsigned char size = thisChar;
        if(size <= ASCIIMAX)
            ++hash[thisChar];
    }
    string result;
    for(unsigned char code = 0;;++code )
    {
        
        while(hash[code])
        {
            result.push_back(char(code));
            --hash[code];
        }
        if(code == ASCIIMAX)
            break;
        
    }
    return result;
}

using VecIter = std::vector < int >::iterator ;
using RevVecIter = vector < int >::reverse_iterator;


template <class RevIter>
void MergeReverse(RevIter rbegin0,RevIter rend0,RevIter rbegin1,RevIter rend1, RevIter out )
{
    while(rbegin1 != rend1)
    {
        if(rbegin0== rend0)
        {
            while(rbegin1 != rend1)
            {
                *--out = *--rend1;
            }
            return;
        }
        
        if(*--rend1 < *--rend0)
        {
            *--out = *rend0;
            ++rend1;
        }
        else
        {
            *--out = *rend1;
            ++rend0;
        }
    }
    
}

void MergeSortHelper(VecIter begin,VecIter end, VecIter aux )
{
    auto dist = distance(begin,end);
    
    if(dist < 2)
    return;
    
    dist /=2;
    
    auto mid = begin + dist;
    
    MergeSortHelper(begin,mid,aux);
    
    MergeSortHelper(mid,end,aux);
    std::copy(mid,end,aux);
    
    MergeReverse(begin,mid,aux, aux + distance(mid,end),end);
}
vector < int > MergeSort(vector < int > intArr) {
    if(intArr.size() <2 )
    {
        return intArr;
    }
    auto dist = std::distance(intArr.begin(),intArr.end());
    vector < int > aux(dist/2+1);
    
    MergeSortHelper(intArr.begin(),intArr.end(),aux.begin());
    
    return intArr;
}

  /* 0
   1   2 
  3 4 5 6
   
   */

template <typename T>
class TNode
{
public:
    TNode(T data):m_data(data)
    {
        
    }
    void Traverse()
    {
        if(m_pLeft)
            m_pLeft->Traverse();
        //cout << m_data << endl;
        if(m_pRight)
            m_pRight->Traverse();
    }
    T m_data;
    TNode* m_pLeft, * m_pRight;
};

template <typename T>
TNode<T>* BuildTree(vector < int >::const_iterator begin, vector < int >::const_iterator end)
{
    if(begin == end)
        return nullptr;
    auto mid = begin + distance(begin,end)/2;
    auto root = new TNode<T>(*(mid));
    root->m_pLeft =  BuildTree<int>(begin,mid);
    root->m_pRight =  BuildTree<int>(mid+1, end);
    return root;
}




