/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2020/11/21
 */

#ifndef IOE_SORW_CMPEDGECUT_HPP
#define IOE_SORW_CMPEDGECUT_HPP
#include <fstream>
#include <vector>
#include <engine/Settings.hpp>

bid_t getBlock(vid_t vertex, bid_t nBlocks, std::vector<vid_t>startVertex){
    bid_t b = 0;
    while (b < nBlocks){
        if (vertex < startVertex[b + 1])
            return b;
        b++;
    }
}

double cmpEdgeCut(const std::string& graphFileName, const std::string& brfName){
    std::vector<vid_t>startVertex;
    std::ifstream brf(brfName);
    std::string line;
    while (brf.good()){
        getline(brf, line);
        if (line.empty())
            continue;
        std::istringstream tmp(line);
        vid_t vertex;
        tmp >> vertex;
        startVertex.emplace_back(vertex);
    }
    brf.close();
    bid_t nBlocks = startVertex.size();

    std::ifstream graph(graphFileName);
    vid_t from;
    vid_t to;
    eid_t edgeIn = 0;
    eid_t edgeCut = 0;
    while (getline(graph, line)){
        if (line.empty() || line.at(0) =='#' || line.at(0) == '%')
            continue;
        std::istringstream tmp(line);
        tmp >> from >> to;
        if (getBlock(from, nBlocks, startVertex) == getBlock(to, nBlocks, startVertex)){
            edgeIn++;
        }else{
            edgeCut++;
        }
    }
    return static_cast<double>(edgeCut) / static_cast<double>(edgeIn + edgeCut);
}

#endif //IOE_SORW_CMPEDGECUT_HPP
