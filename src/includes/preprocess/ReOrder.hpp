/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2020/11/20
 */

#ifndef IOE_SORW_REORDER_HPP
#define IOE_SORW_REORDER_HPP

#include "engine/Settings.hpp"
#include <fstream>
#include <algorithm>

class Node{
public:
    vid_t vertex;
    Node *next;
    explicit Node(vid_t _vertex){
        this->vertex = _vertex;
        next = nullptr;
    }
};

class ParTbl{
private:
    std::vector<Node *> subGraph;
    std::vector<Node *> tail;
public:
    explicit ParTbl(bid_t nBlocks);
    void addVertex(vid_t vertex, bid_t block);
    Node *getSubGraphFirstNode(bid_t block);
};

ParTbl::ParTbl(bid_t nBlocks): subGraph(nBlocks, nullptr), tail(nBlocks, nullptr) {}

void ParTbl::addVertex(vid_t vertex, bid_t block) {
    auto pNode = new Node(vertex);
    if (!subGraph[block]){
        subGraph[block] = pNode;
        tail[block] = pNode;
    }else{
        tail[block]->next = pNode;
        tail[block] = pNode;
    }
}

Node *ParTbl::getSubGraphFirstNode(bid_t block) {
    return subGraph[block];
}

class Edge{
public:
    vid_t from;
    vid_t to;
public:
    Edge(vid_t _from, vid_t _to);
};

Edge::Edge(vid_t _from, vid_t _to) {
    this->from = _from;
    this->to = _to;
}

class ReOrder{
private:
    ParTbl parTbl;
    std::vector<vid_t> newId;
    std::vector<Edge> edgeList;
    vid_t nVertices;
    static bool cmpEdge(Edge edge1, Edge edge2);
public:
    ReOrder(const std::string& graphFile, const std::string& parFile, bid_t nBlocks);
    void genReOrderFile(const std::string& fileName);
    vid_t getStartVertex(bid_t block);
    vid_t getNVertex();
};

bool ReOrder::cmpEdge(Edge edge1, Edge edge2) {
    return edge1.from < edge2.from;
}

ReOrder::ReOrder(const std::string& graphFileName, const std::string& parFileName, bid_t nBlocks): parTbl(nBlocks) {
    vid_t vertex = 0;
    bid_t block;
    std::ifstream parFile(parFileName);
    std::string line;
    while (parFile.good()){
        getline(parFile, line);
        if (line.empty())
            continue;
        std::istringstream iss(line);
        iss >> block;
        parTbl.addVertex(vertex, block);
        vertex++;
    }
    parFile.close();
    nVertices = vertex;

    newId.assign(nVertices, 0xffffffff);
    vertex = 0;
    for(bid_t b = 0; b < nBlocks; b++){
        auto pNode = parTbl.getSubGraphFirstNode(b);
        while (pNode){
            newId[pNode->vertex] = vertex;
            pNode = pNode->next;
            vertex++;
        }
    }

    std::ifstream graphFile(graphFileName);

    vid_t sourceVertex, destVertex;
    while (graphFile.good()){
        getline(graphFile, line);
        if (line.empty() || line.at(0) == '#' || line.at(0) == '%'){
            continue;
        }
        std::istringstream iss(line);
        iss >> sourceVertex;
        iss >> destVertex;
        if (sourceVertex == destVertex)
            continue;
        edgeList.emplace_back(newId[sourceVertex], newId[destVertex]);
    }
    graphFile.close();
    std::sort(edgeList.begin(), edgeList.end(), cmpEdge);
}

void ReOrder::genReOrderFile(const std::string& fileName) {
    std::ofstream reOrderFile(fileName);
    for (auto & edge: edgeList){
        reOrderFile << edge.from << '\t' << edge.to << std::endl;
    }
    reOrderFile.close();
}

vid_t ReOrder::getStartVertex(bid_t block) {
    return newId[parTbl.getSubGraphFirstNode(block)->vertex];
}

vid_t ReOrder::getNVertex() {
    return nVertices;
}

#endif //IOE_SORW_REORDER_HPP
