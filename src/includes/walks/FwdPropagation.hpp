//
// Created by Ethan on 2021/1/8.
//

#ifndef IOE_SORW_FWDPROPAGATION_HPP
#define IOE_SORW_FWDPROPAGATION_HPP

#include <utility>
#include <vector>

#include "WalkManager.hpp"
#include "IO/VertexIO.hpp"
#include "BasicIncludes.hpp"
#include "sampler/FirstOrder/AliasSampler.hpp"

class FwdPropagation: public RandomWalk{
private:
    int k;  // num of sampling level
    std::vector<int>s; // sampling nums of each level
    vid_t startVertexId = 0;
public:
    void setStartVertexId(vid_t _startVertexId) {
        FwdPropagation::startVertexId = _startVertexId;
    }

private:

    vid_t *startVertex;
    bid_t nBlocks;

    VertexIO vertexIo;
    Alias1stOrderSampler alias1StOrderSampler;

public:
    FwdPropagation(int _k, std::vector<int>_s, int begPosFile, int csrFile, vid_t _nTotalVertices, bid_t _nBlocks,
                   vid_t *_startVertex, metrics &m)
                   : vertexIo(begPosFile, csrFile, _nBlocks, _startVertex, _nTotalVertices),
                     alias1StOrderSampler(&vertexIo, m, true){
        k = _k;
        s = std::move(_s);
        startVertex = _startVertex;
        nBlocks = _nBlocks;
    }

    void initWalks(vid_t nRootVertex, WalkManager &walkManager, metrics &m) override{
        tid_t nThreads = get_option_int("nThreads", omp_get_max_threads());
        omp_set_num_threads(nThreads);
        m.start_time("initWalks");

        walkManager.currentNWalks = nRootVertex;
        for (vid_t n = 0; n < nRootVertex; n++){
            vid_t sourceVertexId = startVertexId + n;
            VertexInfo v(sourceVertexId);
            vertexIo.getVertexCSR(v, m);
            if (!v.outDegree){
                walkManager.currentNWalks -=1;
                continue;
            }
            bid_t resideBlock = vertexIo.getBlock(sourceVertexId);
            vid_t currentVertexId = sourceVertexId - startVertex[resideBlock];
            WalkDataType walk = walkManager.encode(sourceVertexId, sourceVertexId, currentVertexId, 0);
            walkManager.moveWalk(walk, resideBlock, omp_get_thread_num(), sourceVertexId, currentVertexId, 0);
            vertexIo.vertexBuffer.streamV(v);
        }
        vid_t v = 0;
        while (true){
            if (v == startVertexId)
                v += nRootVertex;
            if (v >= vertexIo.getNTotalVertices())
                break;
            VertexInfo vertex(v);
            vertexIo.getVertexCSR(vertex, m);
            if (vertex.outDegree)
                vertexIo.vertexBuffer.streamV(vertex);
            v++;
        }
        vertexIo.vertexBuffer.sort();

        m.stop_time("initWalks");

        for (bid_t block = 0; block < nBlocks; block++){
            walkManager.nWalks[block] = walkManager.nDiskWalks[block];

            for (tid_t t = 0; t < nThreads; t++){
                walkManager.nWalks[block] += walkManager.walkPool[t][block].size_w;
            }
            if (walkManager.nWalks[block]){
                walkManager.minStep[block] = 0;
            }
        }
    }

    void updateWalk(metrics &m, const WalkDataType &walk, bid_t execBlock, eid_t *beg_pos, vid_t *csr, WalkManager &walkManager) override {
        tid_t threadId = omp_get_thread_num();
        vid_t currentVertex = walkManager.getCurrentVertex(walk);
        vid_t sourceVertex = walkManager.getSourceVertex(walk);
        hid_t nHops = walkManager.getHops(walk);
        int nSample = s.at(nHops);

        vid_t absCurrentVertex = startVertex[execBlock] + currentVertex;
        VertexInfo v(absCurrentVertex);

//#pragma omp parallel for schedule(static)
        for (int i = 0; i < nSample; i++){
            vid_t destVertexId = alias1StOrderSampler.sampleDestVertex(v, m);
            if (nHops == k - 1)
                continue;
            bid_t destVertexResideBlock = vertexIo.getBlock(destVertexId);
            WalkDataType newWalk = walkManager.encode(sourceVertex, absCurrentVertex, destVertexId - startVertex[destVertexResideBlock], nHops +1);
            if (destVertexResideBlock == execBlock){
                updateWalk(m, newWalk, execBlock, nullptr, nullptr, walkManager);
            }else{
                walkManager.moveWalk(newWalk, destVertexResideBlock, threadId, absCurrentVertex, destVertexId - startVertex[destVertexResideBlock], 0);
                walkManager.isModified[destVertexResideBlock] = true;
            }
        }
    }

    void setIOInfo(bid_t _nInMemBlocks, vid_t **_csrBuf, eid_t **_begPosBuf, bid_t *_inMemIndex) override {
        vertexIo.setIOInfo(_nInMemBlocks, _csrBuf, _begPosBuf, _inMemIndex);
    }
};

#endif //IOE_SORW_FWDPROPAGATION_HPP
