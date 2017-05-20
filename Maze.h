//
//  Maze.h
//  PlayArea
//
//  Created by zayan on 8/1/15.
//  Copyright (c) 2015 Eftiquar. All rights reserved.
//

#ifndef __PlayArea__Maze__
#define __PlayArea__Maze__

#include <stdio.h>

#include <fstream>

#include <vector>
#include <map>
#include <set>
#include <queue>
#include <algorithm>
#include <string>
#include <sstream>

using std::string;
using std::ifstream ;
using std::queue;
using std::vector;
using std::map;
using std::set;
using std::istream_iterator;


using GraphElementIterator = istream_iterator<size_t>;
using OnPath = void (*) (const vector<size_t>& );


class CUndirectedGraph
{
    friend ifstream& operator >> (ifstream & in, CUndirectedGraph & graph );
public:
    CUndirectedGraph();
    const set<size_t>& AdjacentList(size_t nVertex)const;
    size_t VerticesCount()const;
private:
    void LoadGraph(GraphElementIterator begin, GraphElementIterator end);
    
    void  InsertEdge(size_t V, size_t W);
    size_t m_nVertices;
    size_t m_nEdges;
    
    map<size_t,set<size_t>> m_mapAdjacencyList;
    
};

class CBFSearch
{
public:
    
    void GetShortestPath(size_t destVertex,vector<size_t>& shortestPath);
    CBFSearch(const CUndirectedGraph& uGraph, int argc , const char* argv[]);
    void DoBFS(const OnPath& onExitFound,const OnPath& onNoExitFound);
    
private:
    const CUndirectedGraph& m_Graph;
    size_t m_nSourceVertex;
    mutable vector<bool> m_bVisited;
    map<size_t,size_t> m_mapPathElements;
    vector<size_t> m_exitVertices;
};
#endif /* defined(__PlayArea__Maze__) */
