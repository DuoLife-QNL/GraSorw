/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/12/13 16:30
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cassert>

#include "engine/Settings.hpp"
#include "metrics/GetMemoryInfo.hpp"

const vid_t nVertices_CW = 3563602788;
//const vid_t nVertices_CW = 100000000;
typedef std::vector<std::vector<vid_t>> graph_t;
const std::string CRAWLWEB_GRAPH = "/mnt/data/hongzheng/data/crawlweb/crawlweb.txt";
const std::string TEST_GRAPH = "/home/hongzheng/data/testData/testGraph.txt";

class graph{
public:
    graph_t g;

    vid_t src_index = 0, dst_index = 0;

    void getFirstEdge(vid_t &src, vid_t &dst, bool &successful){
        if (g.at(0).size()){
            src = 0;
            dst = g.at(0).at(0);
        }else{
            getNextEdge(src, dst, successful);
        }
    }

    void getNextEdge(vid_t &src, vid_t &dst, bool &successful){
        dst_index ++;
        if (dst_index >= g.at(src_index).size()){
            src_index ++;
            dst_index = 0;
        }
        while (src_index < g.size() && g.at(src_index).size() == 0){
            src_index ++;
            dst_index = 0;
        }
        if (src_index >= g.size()){
            successful = false;
        }else{
            src = src_index;
            dst = g.at(src_index).at(dst_index);
            successful = true;
        }
    }
};

/* get the edges which need to be added into the new file */
graph getOneWayEdges(std::string fileName, vid_t nVertices){
    graph addingEdges;
    addingEdges.g.resize(nVertices);
    std::ifstream edgeListFile(fileName);
    vid_t src, dst;
    while (edgeListFile.good()){
        edgeListFile >> src >> dst;
        if (src % 1000 == 0){
            std::cout << "Scanning source: " << src << std::endl;
        }
        auto pos = std::find(addingEdges.g[src].begin(), addingEdges.g[src].end(), dst);
        if (pos != addingEdges.g[src].end()){
            /* is not one-way edge */
            addingEdges.g[src].erase(pos);
            addingEdges.g[src].shrink_to_fit();
        }else{
            addingEdges.g[dst].push_back(src);
        }
    }
    edgeListFile.close();
    return addingEdges;
}

/* return True if <src1, dst1> should appear before <src2, dst2> */
bool edgeFirstThan(const vid_t &src1, const vid_t &dst1, const vid_t &src2, const vid_t &dst2){
    if (src1 < src2){
        return true;
    }
    if (src1 == src2 && dst1 < dst2){
        return true;
    }
    return false;
}

void insertEdges(std::string fileName, graph &addingEdges){
    std::string fileName_unDir = fileName + ".unDir";
    std::ifstream originFile(fileName);
    std::ofstream newFile(fileName_unDir);
    vid_t source_file, dest_file;
    vid_t source_insert, dest_insert;
    bool moreEdgesToAdd = true;
    originFile >> source_file >> dest_file;
    addingEdges.getFirstEdge(source_insert, dest_insert, moreEdgesToAdd);
    bool lineLeft = false;
    while (moreEdgesToAdd && originFile.good()){
        if (edgeFirstThan(source_file, dest_file, source_insert, dest_insert)){
            newFile << source_file << "\t" << dest_file << std::endl;
            originFile >> source_file >> dest_file;
            lineLeft = true;
        }else{
            newFile << source_insert << "\t" << dest_insert << std::endl;
            lineLeft = false;
            addingEdges.getNextEdge(source_insert, dest_insert, moreEdgesToAdd);
        }
    }
    while (moreEdgesToAdd){
        newFile << source_insert << "\t" << dest_insert << std::endl;
        addingEdges.getNextEdge(source_insert, dest_insert, moreEdgesToAdd);
    }
    while (originFile.good()){
        newFile << source_file << "\t" << dest_file << std::endl;
        originFile >> source_file >> dest_file;
        lineLeft = true;
    }
    if (lineLeft){
        newFile << source_file << "\t" << dest_file << std::endl;
    }
    originFile.close();
    newFile.close();
}

void memoryTest(){
    graph testG;
    std::vector<vid_t> neighbors;
    neighbors.assign(10, 5828);
    std::cout << "size of neighbor list: " << sizeof(neighbors[0]) * neighbors.capacity() << std::endl;
//    std::cout << "size of test graph: " << (double)(sizeof(testG.g[0]) * testG.g.capacity() + sizeof(neighbors[0]) * neighbors.capacity() * nVertices_CW) / (1024 * 1024 * 1024) << "GB" << std::endl;
    testG.g.assign(nVertices_CW, neighbors);
    std::cout << "size of test graph: " << (sizeof(testG.g[0]) * testG.g.capacity() + sizeof(neighbors[0]) * neighbors.capacity() * nVertices_CW) << "B" << std::endl;
    double memoryUsage;
    memoryUsage = get_memory();
    std::cout << "Memory usage now: " << memoryUsage << std::endl;
    for (vid_t src = 0; src < testG.g.size(); src++){
        testG.g.at(src).erase(testG.g.at(src).begin() + 2, testG.g.at(src).begin() + 6);
        testG.g.at(src).shrink_to_fit();
    }
    std::cout << "delete the 8th element of all sources, and shrink to fit" << std::endl;
    std::cout << "size of a single neighbor list now: " << sizeof(testG.g.at(0).at(0)) * testG.g.at(0).capacity() << std::endl;
//    std::cout << "size of the graph now: " << (double)(sizeof(testG.g[0]) * testG.g.capacity() + sizeof(neighbors[0]) * neighbors.capacity() * nVertices_CW) / (1024 * 1024 * 1024) << "GB" << std::endl;
    std::cout << "size of the graph now: " << (sizeof(testG.g[0]) * testG.g.capacity() + sizeof(testG.g[0][0]) * testG.g[0].capacity() * nVertices_CW) << "B" << std::endl;
    memoryUsage = get_memory();
    std::cout << "Memory usage now: " << memoryUsage << std::endl;
    for (vid_t src = 0; src < testG.g.size(); src++){
        testG.g.at(src).push_back(1);
        testG.g.at(src).shrink_to_fit();
    }
    std::cout << "push back a new element to each neighbor list, and shrink to fit" << std::endl;
    std::cout << "size of a single neighbor list now: " << sizeof(testG.g.at(0).at(0)) * testG.g.at(0).capacity() << std::endl;
//    std::cout << "size of the graph now: " << (double)(sizeof(testG.g[0]) * testG.g.capacity() + sizeof(neighbors[0]) * neighbors.capacity() * nVertices_CW) / (1024 * 1024 * 1024) << "GB" << std::endl;
    std::cout << "size of the graph now: " << (sizeof(testG.g[0]) * testG.g.capacity() + sizeof(testG.g[0][0]) * testG.g[0].capacity() * nVertices_CW) << "B" << std::endl;
    memoryUsage = get_memory();
    std::cout << "Memory usage now: " << memoryUsage << std::endl;

    for (vid_t src = 0; src < testG.g.size(); src++){
        testG.g.at(src).clear();
        testG.g.at(src).shrink_to_fit();
    }
    std::cout << "clear all neighbor list, and shrink to fit" << std::endl;
    std::cout << "size of a single neighbor list now: " << sizeof(testG.g.at(0).at(0)) * testG.g.at(0).capacity() << std::endl;
//    std::cout << "size of the graph now: " << (double)(sizeof(testG.g[0]) * testG.g.capacity() + sizeof(neighbors[0]) * neighbors.capacity() * nVertices_CW) / (1024 * 1024 * 1024) << "GB" << std::endl;
    std::cout << "size of the graph now: " << (sizeof(testG.g[0]) * testG.g.capacity() + sizeof(testG.g[0][0]) * testG.g[0].capacity() * nVertices_CW) << "B" << std::endl;
    memoryUsage = get_memory();
    std::cout << "Memory usage now: " << memoryUsage << std::endl;

}

void reverseEdges(std::string originalFileName, std::string reverseFileName){
    std::ifstream originalGraph(originalFileName);
    std::ofstream reverseGraph(reverseFileName);
    vid_t source, dest;
    vid_t lastInputSource, lastInputDest;
    while (originalGraph.good()){
        originalGraph >> source >> dest;
        if (source == lastInputSource && dest == lastInputDest)
            continue;
        reverseGraph << dest << "\t" << source << std::endl;
        lastInputSource = source;
        lastInputDest = dest;
    }
    originalGraph.close();
    reverseGraph.close();
}

class Edge{
public:
    vid_t source;
    vid_t dest;

    bool operator < (const Edge &b) const
    {
        if      (source < b.source)  return true;
        else if (source > b.source)  return false;
        // we get here when source are the same. now sort on dest
        if      (dest < b.dest)  return true;
        else if (dest >= b.dest) return false;
    }

    bool operator != (const Edge &b) const
    {
        if (source == b.source && dest == b.dest){
            return false;
        }else{
            return true;
        }
    }

    bool operator == (const Edge &b) const
    {
        if (source == b.source && dest == b.dest){
            return true;
        }else{
            return false;
        }
    }

    // overload the << operator for writing an Edge
    friend std::ostream& operator<<(std::ostream &os, const Edge &b)
    {
        os  << b.source  << "\t"
            << b.dest;
        return os;
    }
    // overload the >> operator for reading into an Edge
    friend std::istream& operator>>(std::istream &is, Edge &b)
    {
        is  >> b.source
            >> b.dest;
        return is;
    }
};


void mergeTwoFiles(const std::string &originalFileName, const std::string &reverseFileName, const std::string &outputFileName){
    Edge edge1, edge2, lastOutputEdge;
    lastOutputEdge.source = INVALID_VID;
    lastOutputEdge.dest = INVALID_VID;
    std::ifstream originalFile(originalFileName);
    std::ifstream reverseFile(reverseFileName);
    std::ofstream outputFile(outputFileName);

    /* read the first line of two files */
    originalFile >> edge1;
    reverseFile >> edge2;

    while (originalFile.good() && reverseFile.good()){
        while (edge1 < edge2 && originalFile.good()){
            if (edge1 != lastOutputEdge){
                outputFile << edge1 << std::endl;
                lastOutputEdge = edge1;
            }
            originalFile >> edge1;
        }
        while (!(edge1 < edge2) && reverseFile.good()){
            if (edge2 != lastOutputEdge){
                outputFile << edge2 << std::endl;
                lastOutputEdge = edge2;
            }
            reverseFile >> edge2;
        }
    }
    if (edge1 < edge2){
        if (edge1 != lastOutputEdge){
            outputFile << edge1 << std::endl;
            lastOutputEdge = edge1;
        }
        if (edge2 != lastOutputEdge){
            outputFile << edge2 << std::endl;
            lastOutputEdge = edge2;
        }
        while (reverseFile >> edge2 && lastOutputEdge != edge2){
            outputFile << edge2 << std::endl;
            lastOutputEdge = edge2;
        }
    }else{
        if (edge2 != lastOutputEdge){
            outputFile << edge2 << std::endl;
            lastOutputEdge = edge2;
        }
        if (edge1 != lastOutputEdge){
            outputFile << edge1 << std::endl;
            lastOutputEdge = edge1;
        }
        while (originalFile >> edge1 && lastOutputEdge != edge1){
            outputFile << edge1 << std::endl;
            lastOutputEdge = edge1;
        }
    }
}

int main(int argc, char *argv[]){
//    graph CW = getOneWayEdges(CRAWLWEB_GRAPH, nVertices_CW);
//    insertEdges(CRAWLWEB_GRAPH, CW);
//    memoryTest();
    std::string inFileName = argv[1];
    std::string outFileName = argv[2];
    reverseEdges(inFileName, outFileName);
//    mergeTwoFiles(TEST_GRAPH, TEST_GRAPH + ".reverse", TEST_GRAPH + ".unDir");
//    std::string inFile1 = argv[1];
//    std::string inFile2 = argv[2];
//    std::string outFile = argv[3];
//    mergeTwoFiles(inFile1, inFile2, outFile);
    return 0;
}