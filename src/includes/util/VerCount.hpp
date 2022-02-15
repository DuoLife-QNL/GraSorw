//
// Created by lihz on 2020/10/29.
//

#ifndef IOE_SORW_VERCOUNT_HPP
#define IOE_SORW_VERCOUNT_HPP

#include <string>
#include <cmath>
#include <fstream>
#include "engine/Settings.hpp"
#include "logger/logger.hpp"

size_t VerCount(std::string baseFileName){
    bool *tbl;
    vid_t from;
    vid_t to;
    size_t maxVerNum = pow(2, sizeof(vid_t) * 8);
    size_t upBound = 0;
    size_t nVertices = 0;
    tbl = new bool[maxVerNum];
    std::ifstream edgeList(baseFileName.c_str());
    if (!edgeList.good()){
        logstream(LOG_ERROR) << "Wrong edgeList file" << std::endl;
        assert(false);
    }
    while (!edgeList.eof()){
        edgeList >> from >> to;
        tbl[from] = true;
        tbl[to] = true;
        upBound += 2;
    }
    for (vid_t i = 0; i < upBound; i++){
        if (tbl[i]){
            nVertices++;
        }
    }
    return nVertices;
}

vid_t getMaxVertexId(std::string baseFileName){
    std::ifstream inFile(baseFileName);
    vid_t maxId = 0, sourceVertex, destVertex;
    std::string line;
    eid_t nReadEdge = 0;
    while (inFile.good()){
        getline(inFile, line);
        if (line.empty()){
            continue;
        }
        if (line.at(0) == '#' || line.at(0) == '%'){
            continue;
        }
        nReadEdge++;
        std::istringstream iss(line);
        iss >> sourceVertex;
        iss >> destVertex;
        if (sourceVertex > maxId)
            maxId = sourceVertex;
        if (destVertex > maxId)
            maxId = destVertex;
    }
}

#endif //IOE_SORW_VERCOUNT_HPP
