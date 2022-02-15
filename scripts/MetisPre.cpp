/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2020/11/13
 */
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sstream>
#include <vector>
#include <map>
#include "engine/Settings.hpp"
using namespace std;
class VertexNode{
public:
    vid_t id;
    VertexNode *next;
    explicit VertexNode(vid_t _id) {
        id = _id;
        next = nullptr;
    }
};
class GraphTbl{
private:
    vid_t nVertex;
    vector<VertexNode *> tbl;
    vector<VertexNode *> tail;
    eid_t nEdges;
public:

    void addEdge(vid_t source, vid_t dest){
        if (tbl[source] == nullptr){
            auto *pNode = new VertexNode(dest);
            tbl[source] = pNode;
            tail[source] = pNode;
            nEdges++;
        }else{
            if (!checkVertex(source, dest)){
                auto *pNode = new VertexNode(dest);
                tail[source]->next = pNode;
                tail[source] = pNode;
                nEdges++;
            }
        }
    }

    explicit GraphTbl(vid_t _nVertex): tbl(_nVertex, nullptr), tail(_nVertex, nullptr){
        nVertex = _nVertex;
        nEdges = 0;
    }

    bool checkVertex(vid_t source, vid_t dest){
        VertexNode *pNode = tbl[source];
        while (pNode){
            if (pNode->id == dest){
                return true;
            }
            pNode = pNode->next;
        }
        return false;
    }

//    void addNeighbor(vid_t vertex1, vid_t vertex2){
//        addEdge(vertex1, vertex2);
//        addEdge(vertex2, vertex1);
//    }

    void genMetisFile(const string& fileName){
        ofstream outFile(fileName);
        outFile << nVertex << '\t' << nEdges / 2;
        VertexNode *pNode;
        for (vid_t v = 0; v < nVertex; v++){
            outFile << endl;
            pNode = tbl[v];
            outFile << pNode->id + 1;
            pNode = pNode->next;
            while (pNode){
                outFile << " " << pNode->id + 1;
                pNode = pNode->next;
            }
        }
    }
};

int main (int argc, char *argv[])
{
    int opt;
    string baseFileName;
    opt = getopt (argc, argv, "f:");
    if(opt){
        baseFileName = optarg;
    }else{
        cout << "Use -f to declare input file name!" << endl;
        return 0;
    }
    string outFileName = baseFileName;
    outFileName.append(".toMetis");
    ifstream inFile(baseFileName);
    vid_t maxOriId = 0, sourceVertex, destVertex;
    string line;
    vector<vid_t> oriSrcV, oriDestV;
    while (inFile.good()){
        getline(inFile, line);
        if (line.empty()){
            continue;
        }
        if (line.at(0) == '#' || line.at(0) == '%'){
            continue;
        }
        istringstream iss(line);
        iss >> sourceVertex;
        iss >> destVertex;
//        if (sourceVertex == destVertex){
//            continue;
//        }
        if (sourceVertex > maxOriId)
            maxOriId = sourceVertex;
        if (destVertex > maxOriId)
            maxOriId = destVertex;
        oriSrcV.push_back(sourceVertex);
        oriDestV.push_back(destVertex);
    }
//    vid_t nVertex = 0;
//    vector<vid_t> newId(maxOriId + 1, -1);
//    vid_t newId[maxOriId + 1];
//    for (vid_t v = 0; v <= maxOriId; v++){
//        newId[v] = -1;
//    }
//    GraphTbl graphTbl(maxOriId + 1);
    std::vector<std::vector<vid_t>> neighbors;
    neighbors.resize(maxOriId + 1);
//    for (eid_t e = 0; e < oriSrcV.size(); e++){
//        sourceVertex = oriSrcV[e];
//        destVertex = oriDestV[e];
////        if (newId[sourceVertex] == -1){
////            newId[sourceVertex] = nVertex;
////            nVertex++;
////        }
////        if (newId[destVertex] == -1){
////            newId[destVertex] = nVertex;
////            nVertex++;
////        }
//        graphTbl.addEdge(sourceVertex, destVertex);
//
//    }
    for (eid_t e = 0; e < oriSrcV.size(); e++){
        sourceVertex = oriSrcV[e];
        destVertex = oriDestV[e];
        neighbors.at(sourceVertex).push_back(destVertex);
    }
//    for (eid_t e = 0; e < oriSrcV.size(); e++){
//        sourceVertex = newId[oriSrcV[e]];
//        destVertex = newId[oriDestV[e]];
//        graphTbl.addNeighbor(sourceVertex, destVertex);
//    }
//    graphTbl.genMetisFile(outFileName);
    ofstream outFile(outFileName);
    outFile << maxOriId + 1 << '\t' << oriSrcV.size() / 2;
    for (vid_t v = 0; v <= maxOriId; v++){
        eid_t n = 0;
        outFile << std::endl << neighbors.at(v).at(n) + 1;
        n++;
        while (n < neighbors.at(v).size()){
            outFile << " " << neighbors.at(v).at(n) + 1;
            n++;
        }
    }
    outFile.close();
    return 0;
}