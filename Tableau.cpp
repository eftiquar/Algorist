//
//  Tableau.cpp
//  PlayArea
//
//  Created by zayan on 8/15/15.
//  Copyright (c) 2015 Eftiquar. All rights reserved.
//

#include "Tableau.h"
#include "assert.h"
#include <cwchar>
using std::wcscmp;
void TestBinSearch()
{
    // floating point comparison is fraught business, the rounding errors need to be accounted for
    // "Being approximately right is better than being exactly wrong". This adage fits like a t or should it be 'f' for floating points ?
    // For customizing float1 > float2 operations, a class can be implemented with operator > overloaded
    // one way could be to promote both floating operands to bouble and then compare them for extra precision
    //or a defined delta threshold which defines when a rounding should be down or up
    
    // our BSearch function is a template and relies on operator > to perform search
    
    float floaters[] = {0.0, 0.4,0.5,0.5,0.5,0.78,10.11};
    constexpr size_t arrayLength = sizeof(floaters)/sizeof(float);
    
    float* begin = floaters;
    float* end = floaters + arrayLength;
    
    
    
    //search for an arbitrary key at position 5
    float* key = floaters + 5;
    float* searchResult = Tableau::BSearch<float>(begin, end, *key);
    assert(searchResult == key);
    
    
    //search an arbitrary key that does not exist, but its value is  between begin and begin + 1, so expected result is begin + 1
    //the key does not exist, but its position would be begin +1 if we were to insert it and maintain sortedness
    searchResult = Tableau::BSearch<float>(begin, end, 0.0000001);
    assert(searchResult == begin +1);

    //search for duplicate elements must return position of first occurence of duplicated element
    searchResult = Tableau::BSearch<float>(begin, end, 0.5);
    assert(searchResult == begin +2);

    //boundary condition - lower bound, search for zeroth element
    searchResult = Tableau::BSearch<float>(begin,end, *begin);
    assert(searchResult == begin);
    //non existent element which is smaller than smallest
    searchResult = Tableau::BSearch<float>(begin,end, -0.0);
    assert(searchResult == begin);
    

    //boundary condition - upper bound

    searchResult = Tableau::BSearch<float>(begin,end, *(end -1));
    assert(searchResult == end -1);
    
    //search for an element greater than greatest element , so its position should be immediately past the last element, aka end
    searchResult = Tableau::BSearch<float>(begin,end, 11);
    assert(searchResult == end);
    
    //omitting tests for even length and odd length arrays. But they must be considered for completeness
}
void TestReverse()
{
    //basic test on arbitrary string
    //non palindrome string
    wchar_t testString[] =       L"Hello World";
    wchar_t shadowTestString[] = L"Hello World";
    
    // it is paramount that we do not consider the null terminator for reversing, hence deducting 1
    size_t stringLength = sizeof(testString)/sizeof(wchar_t) - 1;
    
    wchar_t* begin = testString;
    wchar_t* end = begin + stringLength;
    
    Tableau::InplaceReverse(begin,end );
    Tableau::InplaceReverse(begin, end);
    assert(wcscmp(testString, shadowTestString) == 0);
    
    //string with even length
    wchar_t evenString[] = L"Even";
    wchar_t evenShadow[] = L"Even";
    
    stringLength = sizeof(evenString)/sizeof(wchar_t) -1;
    begin = evenString;
    end = begin + stringLength;
    
    Tableau::InplaceReverse(begin,end );
    Tableau::InplaceReverse(begin, end);
    assert(wcscmp(evenString, evenShadow) == 0);
    
    //string with odd length
    wchar_t oddString[] = L"Odd";
    wchar_t oddShadowString[] = L"Odd";
    
    stringLength = sizeof(oddString)/sizeof(wchar_t) -1;
    begin = oddString;
    end = begin + stringLength;
    
    Tableau::InplaceReverse(begin,end );
    Tableau::InplaceReverse(begin, end);
    assert(wcscmp(oddString, oddShadowString) == 0);
    
    
    //test for palindrome, should be unchanged after reversal
    wchar_t palindromeString[] = L"level";
    wchar_t shadowPalindrome[] = L"level";
    
    stringLength = sizeof(palindromeString)/sizeof(wchar_t) - 1;
    begin = oddString;
    end = begin + stringLength;
    
    Tableau::InplaceReverse(begin,end );

    assert(wcscmp(palindromeString, shadowPalindrome) == 0);
    
}
int tabmain()
{
    TestBinSearch();
    TestReverse();
    return 0;
}