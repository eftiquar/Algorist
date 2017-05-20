//
//  Algorist.h
//  PlayArea
//
//  Created by zayan on 8/29/15.
//  Copyright (c) 2015 Eftiquar. All rights reserved.
//

#ifndef __PlayArea__Algorist__
#define __PlayArea__Algorist__

#include <stdio.h>
#include <stdio.h>
#include <tuple>
#include <map>
#include <vector>
#include <algorithm>
#include <queue>

using POINT = std::pair<size_t,size_t>;
using PATH = std::vector<POINT>;
using DISTANCE = size_t;
void FindMin(POINT start, POINT destination);
using std::vector;
std::vector<vector<size_t>> LinearPartition(const vector<size_t>& input, size_t partitions);
using std::vector;
using std::wstring;
using std::map;
using std::queue;

using INPUTLIST = vector<std::wstring>;
using STRINGS = vector<wstring>;


#endif /* defined(__PlayArea__Algorist__) */
