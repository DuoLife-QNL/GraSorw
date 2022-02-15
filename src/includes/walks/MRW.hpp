//
// Created by Ethan on 2021/1/7.
//

#ifndef IOE_SORW_MRW_HPP
#define IOE_SORW_MRW_HPP

#include "engine/Settings.hpp"
#include <vector>
#include <bits/stdc++.h>
#include "IO/VertexIO.hpp"
#include "util/RandNum.hpp"

class NodeSampler{
private:
    RandNum mainRand;

public:
    NodeSampler():mainRand(979873879865674341){};

    /* return vertex id */
    vid_t sampleOneRoot(vid_t nVertices){
        return mainRand.iRand(nVertices);
    }

    /* return index of the vector */
    vid_t sampleOneNode(const std::vector<VertexInfo> &sampleSet, vid_t sumDegree){
        vid_t degree = mainRand.iRand(sumDegree);
        vid_t sum = 0;
        for (vid_t v = 0; v < sampleSet.size(); v++){
            vid_t preSum = sum;
            sum += sampleSet.at(v).outDegree;
            if (sum < preSum)
                abort();

            if (degree < sum){
                return v;
            }
        }
        abort();
    }

    vid_t sampleNextNode(const VertexInfo& v){
        vid_t index = mainRand.iRand(v.outDegree);
        return v.csr[index];
    }
};

class MultiRW{
private:
    metrics &m;

    std::vector<VertexInfo> sampleSet;
    vid_t sumDegree = 0;
    vid_t nRoots;
    vid_t nodeBudget;
    vid_t nSampledNodes = 0;

    NodeSampler nodeSampler;

    VertexIO vertexIo;

    void sampleRootNodes(){
        while (nSampledNodes < nRoots){
            VertexInfo v(nodeSampler.sampleOneRoot(vertexIo.getNTotalVertices()));
            vertexIo.getVertexCSR(v, m);
            if (!v.outDegree)
                continue;
            v.delAfterUse = false;
            sampleSet.emplace_back(v);
            sumDegree += v.outDegree;

            vid_t sum = 0;
            for (const auto& it: sampleSet){
                sum += it.outDegree;
            }

            nSampledNodes ++;
        }
    }

public:
    MultiRW(vid_t _nRoots, vid_t _nodeBudget, int begPosFile, int csrFile, vid_t nTotalVertices, metrics &_m)
            : m(_m), vertexIo(begPosFile, csrFile, 0, nullptr, nTotalVertices){
        nRoots = _nRoots;
        nodeBudget = _nodeBudget;
    }

    void run(){
        m.start_time("sample-root-nodes");
        sampleRootNodes();
        m.stop_time("sample-root-nodes");
        while (nSampledNodes < nodeBudget){
            m.start_time("sample-one-node");
            vid_t index = nodeSampler.sampleOneNode(sampleSet, sumDegree);
            m.stop_time("sample-one-node");
            m.start_time("sample-next-node");
            vid_t nextNode = nodeSampler.sampleNextNode(sampleSet.at(index));
            m.stop_time("sample-next-node");
            VertexInfo v(nextNode);
            vertexIo.getVertexCSR(v, m);
            m.start_time("update-sample-set");
            sumDegree = sumDegree - sampleSet.at(index).outDegree + v.outDegree;
            sampleSet.at(index).delAfterUse = true;
            sampleSet.erase(sampleSet.begin() + index);
            v.delAfterUse = false;
            sampleSet.emplace_back(v);
            nSampledNodes ++;
            m.stop_time("update-sample-set");
        }
    }
};

#endif //IOE_SORW_MRW_HPP
