//
//  Maze.cpp
//  PlayArea
//
//  Created by zayan on 8/1/15.
//  Copyright (c) 2015 Eftiquar. All rights reserved.
//
///Users/zayan/Documents/Projects
#include "Maze.h"
#include <iterator>
#include <assert.h>
#include <algorithm>
#include <iostream>
using std::copy;
using std::ostream_iterator;
void PrintUsage()
{
    std::cout << "Maze <input file> <entry> <exit0> <exit1> \r\n";
    std::cout << "input file describes maze in terms of vertices and edges."
                 "\r\n<entry> is vertex id of starting point and exitn are vertex ids of available exits";
    
}

void OnExitPathFound(const vector<size_t>& shortestPath)
{
    std::cout << "Found Exit \r\n";
    
    copy(shortestPath.rbegin(),shortestPath.rend(),ostream_iterator<size_t>(std::cout,"-->"));
    
    std::cout << "\r\n";
    std::cout.flush();
    
}
void OnNoExitPathFound(const vector<size_t>& unfoundExits)
{
    std::cout << "No exit path found for following exit vertices that were specified: \r\n";
    
    copy(unfoundExits.begin(),unfoundExits.end(),ostream_iterator<size_t>(std::cout,",\r\n"));
    
    std::cout.flush();
}

int Mazemain(int argc, const char * argv[])
{
    try
    {

        if(argc < 3)
        {
            PrintUsage();
            return -1;
        }
        
        CUndirectedGraph maze;
        ifstream in(argv[1]);
        if(!in)
            throw std::invalid_argument("Unable to open input file. Terminating");
        in >> maze;
   
        string strMazeData = argv[2];
        
        CBFSearch bfs(maze,argc,argv);
        
        //bfs.DoBFS([](const vector<size_t>& shortestPath){PrintPath(shortestPath);});
        bfs.DoBFS(OnExitPathFound, OnNoExitPathFound);

    }
    catch(std::invalid_argument ex)
    {
        std::cerr << ex.what();
    }
    return 0;
}
CUndirectedGraph::CUndirectedGraph():m_nEdges(0),m_nVertices(0)
{
}
size_t CUndirectedGraph::VerticesCount()const
{
    return m_nVertices;
}

ifstream& operator >> (ifstream & in, CUndirectedGraph & graph )
{
    GraphElementIterator begin = GraphElementIterator(in);
    GraphElementIterator end = std::istream_iterator<size_t>();
    
    //read number of vertices
    
    graph.m_nVertices = *begin++;
    
    assert(begin!= end);
    if(begin == end)
        throw std::invalid_argument("inadequate data: unable to read number of edges");
    
    //read number of edges
    graph.m_nEdges = *begin++;
    
    
    assert(begin!= end);
    if(begin == end)
        throw std::invalid_argument("inadequate data: no edges or vertices found");
 
    //read source
    
    graph.LoadGraph(begin, end);
    
    return in;
}

void CUndirectedGraph::LoadGraph(GraphElementIterator begin, GraphElementIterator end)
{
    size_t edges = 0;
  
    while(begin != end)
    {
        edges++;
        size_t edgeFrom  = *begin++;
        assert(begin!= end);
        
        if(begin == end)
            throw std::invalid_argument("inadequate data: unable to read graph");

        size_t edgeTo = *begin++;

        assert(edgeFrom < m_nVertices);
        assert(edgeTo < m_nVertices);
        
        if(edgeTo >= m_nVertices || edgeFrom >= m_nVertices)
            throw std::invalid_argument("Vertex id cannot exceed the number of vertices. Vertex id range- [0,total Vertices) " );
        
        m_mapAdjacencyList[edgeFrom].insert(edgeTo);
        m_mapAdjacencyList[edgeTo].insert(edgeFrom);
    }
    if(edges != m_nEdges)
            throw std::invalid_argument("Number of edges found does not match the specified number of edges");
}
 const set<size_t>& CUndirectedGraph::AdjacentList(size_t nVertex)const

{
    const auto& adjList = m_mapAdjacencyList.find(nVertex);
    
    if(adjList == m_mapAdjacencyList.end())
       throw std::invalid_argument("Invalid argument to AdjacentList ");
       
    return adjList->second;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
CBFSearch::CBFSearch(const CUndirectedGraph& uGraph, int argc , const char* argv[]):m_Graph(uGraph)
{
    m_bVisited.resize(uGraph.VerticesCount());
    int argIndex = 2;
    m_nSourceVertex = atoi(argv[argIndex++]);
    
    while(argIndex < argc)
    {
        m_exitVertices.push_back(atoi(argv[argIndex++]));
    }
    
}

void CBFSearch:: DoBFS(const OnPath& onExitFound,const OnPath& onNoExitFound)
{
    m_bVisited[m_nSourceVertex] = true;
    queue<size_t> bfsQueue;
    bfsQueue.push(m_nSourceVertex);
    
    while (!bfsQueue.empty())
    {
        size_t parentVertex = bfsQueue.front();
        bfsQueue.pop();
        for( auto adjacentVertex : m_Graph.AdjacentList(parentVertex))
        {
            if(!m_bVisited[adjacentVertex])
            {
                m_mapPathElements[adjacentVertex] = parentVertex;
                m_bVisited[adjacentVertex] = true;
                bfsQueue.push(adjacentVertex);
                
                auto removedElement = remove(m_exitVertices.begin(),m_exitVertices.end(),adjacentVertex);
                
                if(removedElement != m_exitVertices.end())
                {
                    vector<size_t> shortestPath;
                    GetShortestPath(adjacentVertex,shortestPath);
                    onExitFound(shortestPath);
                    
                    m_exitVertices.erase(removedElement);
                    if(m_exitVertices.empty())
                        return;
                }
            }
        }
    }
    if(!m_exitVertices.empty())
    {
        onNoExitFound(m_exitVertices);
    }
    
}

void CBFSearch::GetShortestPath(size_t destVertex,vector<size_t>& shortestPath)
{

    auto pathIter = m_mapPathElements.find(destVertex);
    
    if(pathIter == m_mapPathElements.end())
        return ;

    //backtracking from dest vertex, so push it first, walk up from here to m_nSourceVertex
    shortestPath.push_back(destVertex);
    
    for(;pathIter != m_mapPathElements.end() && pathIter->second != m_nSourceVertex; pathIter = m_mapPathElements.find(pathIter->second) )
    {
        shortestPath.push_back(pathIter->second);
    }
    shortestPath.push_back(m_nSourceVertex);
}

