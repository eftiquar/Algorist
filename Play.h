//
//  Play.h
//  PlayArea
//
//  Created by zayan on 2/13/16.
//  Copyright (c) 2016 Eftiquar. All rights reserved.
//

#ifndef __PlayArea__Play__
#define __PlayArea__Play__

#include "algorist.h"

INPUTLIST DialThis(wstring input,map<wchar_t,wstring> & neighbours);
STRINGS CombineThem(const wstring& in,size_t n );
STRINGS PermutateF(wstring in);
STRINGS PowerSet(wstring in);
vector<size_t> GetPossiblePositions(size_t npos, int nDays);
void BackTrack(vector<size_t> solutions, size_t level, size_t target, const STRINGS& input);

vector<int> LIS( const vector<int> & sequence);
vector<int> MaxZigZAg(const vector<int> & sequence);
vector<wstring> Justify(const wstring & textInput, const size_t length);
size_t MaxRod(const size_t* begin,const size_t* end);
size_t KnapTheSack(vector<std::pair<size_t,size_t> >values, size_t capacity);
std::vector<wstring> JustifyText(const wstring& text, int lineWidth);
size_t LCS (const wchar_t* lhsbegin,const wchar_t* lhsend,const wchar_t* rhsbegin,const wchar_t* rhsend );
size_t GetDenoes(const size_t targetAmount, const std::vector<size_t>& denominations);
std::vector<size_t> KnapsackFull(std::vector<std::pair<size_t,size_t>> items,size_t capacity );
std::vector<std::pair<size_t,size_t>> MergeIngervals(std::vector<std::pair<size_t,size_t>> intervals);
std::vector<std::pair<size_t, size_t>> Skifall( std::vector<size_t>& skis,  std::vector<size_t> skiers);


size_t FindDuplicateOrMissing(size_t* begin, size_t* end);
void BiPartition(size_t* begin, size_t* end);
size_t Multfast(size_t multiplier, size_t multiplicand);
size_t PowerX(size_t base, size_t exp);

#include <random>
template< class RandomIt >
void myrandom_shuffle( RandomIt first, RandomIt last )
{
    auto max = PTRDIFF_MAX;
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
void Insert_Sort(size_t* begin, size_t* end);
std::vector<std::vector<size_t>> Perms(size_t len);

#endif /* defined(__PlayArea__Play__) */
