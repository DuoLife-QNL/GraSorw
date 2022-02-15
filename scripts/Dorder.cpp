/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/6/10 11:41
 * 给定子图数量，输出子图划分使得子图的总度数、总节点数都相近
 * 输出：一个按照新的顺序重新编号后的边列表文件
 * 独立执行没有意义，需要得到执行后的结果后再进行顺序子图划分
 * 输入的原始图文件需要是稠密的，即编号从0开始，没有孤立节点
 */

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sstream>
#include <vector>
#include <cassert>
#include <bits/stdc++.h>
#include "engine/Settings.hpp"

typedef struct Node {
    vid_t vertexId;
    eid_t outDegree;
}Node;

bool compareNodes(Node n1, Node n2){
    return n1.outDegree > n2.outDegree;
}

int main (int argc, char *argv[]){
    int opt = 0;
    std::string baseFileName = " ";
    int nBlocks = 0;
    while ((opt = getopt (argc, argv, "f:p:")) != -1){
        switch (opt) {
            case 'f':
                baseFileName = optarg;
                break;
            case 'p':
                nBlocks = atoi(optarg);
                break;
        }
    }
    if (baseFileName == " "){
        std::cout << "No file name!" << std::endl;
        abort();
    }else if(!nBlocks){
        std::cout << "Partition number should be given!" << std::endl;
        abort();
    }
    std::string outFileName = baseFileName + ".dorder.part." + std::to_string(nBlocks);

    /* 读原始图文件并直接统计度信息 */
    std::ifstream inFile(baseFileName);
    vid_t maxOriId = 0, sourceVertex, destVertex;
    std::string line;
    std::vector<Node> nodes;
    std::vector<vid_t> oriSrcV, oriDestV;
    vid_t currentVertexId = 0;
    vid_t countNeighbors = 0;
    while (inFile.good()){
        getline(inFile, line);
        if (line.empty()){
            continue;
        }
        if (line.at(0) == '#' || line.at(0) == '%'){
            continue;
        }
        std::istringstream iss(line);
        iss >> sourceVertex;
        iss >> destVertex;
        if (sourceVertex > maxOriId)
            maxOriId = sourceVertex;
        if (destVertex > maxOriId)
            maxOriId = destVertex;
        oriSrcV.push_back(sourceVertex);
        oriDestV.push_back(destVertex);
        if (sourceVertex == currentVertexId){
            countNeighbors++;
        }else{
            assert(sourceVertex == currentVertexId + 1);
            Node node = {.vertexId = currentVertexId, .outDegree = countNeighbors};
            nodes.push_back(node);
            currentVertexId++;
            countNeighbors = 1;
        }
    }
    inFile.close();
    Node node = {.vertexId = currentVertexId, .outDegree = countNeighbors};
    nodes.push_back(node);

    /* 将原始ID边列表整理成neighbor数组*/
    std::vector<std::vector<vid_t>> neighbors;
    neighbors.resize(maxOriId + 1);
    for (eid_t e = 0; e < oriSrcV.size(); e++){
        sourceVertex = oriSrcV[e];
        destVertex = oriDestV[e];
        neighbors.at(sourceVertex).push_back(destVertex);
    }

    /* 根据节点的度降序排序 */
    std::sort(nodes.begin(), nodes.end(), compareNodes);

    /* 将节点分配到子图中 */
    std::vector<std::vector<vid_t>> blockVs;
    blockVs.resize(nBlocks);
    vid_t i = 0;
    while (i < nodes.size()){
        for (bid_t b = 0; b < nBlocks; b++){
            if (i == nodes.size()){
                break;
            }
            blockVs.at(b).push_back(nodes.at(i).vertexId);
            i++;
        }
        for (bid_t b = nBlocks - 1; b >= 0 && b < nBlocks; b--){
            if (i == nodes.size()){
                break;
            }
            blockVs.at(b).push_back(nodes.at(i).vertexId);
            i++;
        }
    }

    /* 生成newId和oldId双向查找表 */
    std::vector<vid_t> newId, oldId;
    newId.resize(maxOriId + 1);
    oldId.resize(maxOriId + 1);
    vid_t currentId = 0;
    for (bid_t b = 0; b < nBlocks; b++){
        for (const auto vid: blockVs.at(b)){
            newId.at(vid) = currentId;
            oldId.at(currentId) = vid;
            currentId++;
        }
    }

    /* 输出 */
    std::ofstream reorderFile(outFileName);
    for (vid_t v = 0; v < newId.size(); v++){
        vid_t source_new = v;
        vid_t source_old = oldId.at(source_new);
        for (const auto & dest_old: neighbors.at(source_old)){
            vid_t dest_new = newId.at(dest_old);
            reorderFile << source_new << " " << dest_new << std::endl;
        }
    }
    reorderFile.close();
}