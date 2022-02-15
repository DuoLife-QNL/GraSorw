/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/8/10 11:35
 */

#ifndef IOE_SORW_DEEPWALK_HPP
#define IOE_SORW_DEEPWALK_HPP
#include <string>
#include <fstream>
#include <cmath>
#include <fstream>

#include "BasicIncludes.hpp"
#include "util/toplist.hpp"
#include "util/comperror.hpp"
#include "walks/WalkManager.hpp"
#include "engine/Settings.hpp"
#include "IO/VertexIO.hpp"
#include "walks/FixedLengthFORW.hpp"

typedef unsigned VertexDataType;

class DeepWalk : public FixedLengthFORW{
public:
    int nWalkPerVertex;

private:
    Alias1stOrderSampler alias1StOrderSampler;
    std::vector<Alias1stOrderSampler> aliasSamplers;

public:

    DeepWalk(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int walkNumPerVertex,
             int begPosFile, int csrFile, vid_t nTotalVertices, metrics &m)
             : FixedLengthFORW(nBlocks, startVertex, walkLength, begPosFile, csrFile, nTotalVertices),
             alias1StOrderSampler(&vertexIO, m, true) {
        this->nWalkPerVertex = walkNumPerVertex;
        aliasSamplers.assign(100, alias1StOrderSampler);
        tid_t nThreads = omp_get_max_threads();
    }

    vid_t sampleStartEdge(metrics &m, VertexInfo &srcVertex, bool &noOutEdge) {
        if (!srcVertex.csr){
            vertexIO.getVertexCSR(srcVertex, m);
        }
        if (!srcVertex.outDegree){
            noOutEdge = true;
            return INVALID_VID;
        }else{
            vid_t destVertex = alias1StOrderSampler.sampleDestVertex(srcVertex, m);
            noOutEdge = false;
            return destVertex;
        }
    }

    void initWalks(vid_t nVertices, WalkManager &walkManager, metrics &m) {

        tid_t nThreads = get_option_int("nThreads", omp_get_max_threads());
        omp_set_num_threads(nThreads);
        m.start_time("initWalks");
        walkManager.currentNWalks = nVertices * nWalkPerVertex;
#if MULTI_THREAD
#pragma omp parallel for schedule(static)
#endif
        for (vid_t sourceVertex = 0; sourceVertex < nVertices; sourceVertex++){
            vid_t currentVertex_abs = sourceVertex;
            WalkDataType walk = WalkManager::encode(sourceVertex, INVALID_VID, currentVertex_abs, INVALID_BID, INVALID_BID, 0);
            bid_t resideBlock = vertexIO.getBlock(currentVertex_abs);
            for (wid_t w = 0; w < nWalkPerVertex; w++){
                walkManager.moveWalk(walk, resideBlock, omp_get_thread_num(), INVALID_VID, currentVertex_abs, 0, INVALID_BID);
            }
        }
        //    vertexIO.vertexBuffer.sort();
        m.stop_time("initWalks");


        std::cout << "currentNWalks before run: " << walkManager.currentNWalks << std::endl;

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

    void updateWalk(metrics &m, const WalkDataType &walk, bid_t execBlock, eid_t *beg_pos, vid_t *csr,
                         WalkManager &walkManager) {
        tid_t threadId = omp_get_thread_num();
        vid_t currentVertex = walkManager.getCurrentVertex(walk);
        vid_t previousVertex = 0;
        hid_t previousNHops = walkManager.getHops(walk);
        hid_t forwardHops = 0;
        vid_t nextVertex = 0;

//        /* change currentVertex to absolute value*/
//        currentVertex += startVertex[execBlock];
        while (currentVertex >= startVertex[execBlock] && currentVertex < startVertex[execBlock + 1] && (previousNHops + forwardHops) < walkLength){
#if ACT_VER
            vertexIO.markActivated(currentVertex);
#endif
            //        metric.at(threadId).start_time("sampling");
            nextVertex = sampleDestVertex(m, previousVertex, currentVertex, beg_pos, csr, execBlock);
            //        metric.at(threadId).stop_time("sampling");
            /* NOTE: currentVertex is absolute value here */
            if (nextVertex == INVALID_VID){
                return;
            }
            previousVertex = currentVertex;
            currentVertex = nextVertex;
            forwardHops++;
        }
        if (previousNHops + forwardHops < walkLength){
            bid_t resideBlock = getBlock(currentVertex);
            walkManager.moveWalk(walk, resideBlock, threadId, previousVertex, currentVertex,
                                 forwardHops, execBlock);
            walkManager.setMinStep(resideBlock, previousNHops + forwardHops);
            walkManager.isModified[resideBlock] = true;
        }
    }

    vid_t sampleDestVertex(metrics &m, vid_t preVertex, vid_t curVertex, eid_t *beg_pos, vid_t *csr, bid_t execBlock){
        VertexInfo curVertexInfo(curVertex, execBlock);
        return alias1StOrderSampler.sampleDestVertex(curVertexInfo, m);
    }

};

#endif //IOE_SORW_DEEPWALK_HPP
