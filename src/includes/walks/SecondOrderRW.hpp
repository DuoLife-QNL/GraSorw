//
// Created by lihz on 2020/10/22.
//

#ifndef GRAPHWALKER_SECONDORDERRW_HPP
#define GRAPHWALKER_SECONDORDERRW_HPP

#include "walks/WalkManager.hpp"
#include "engine/Settings.hpp"
#include "BasicIncludes.hpp"
#include "IO/VertexIO.hpp"
#include "walks/RandomWalk.hpp"
#include <cmath>

/**
 * @nBlocks: total number of blocks
 * @startVertex: start vertex of a block = startVertex[blockId]
 * @walkLength: for now we only consider Node2Vec, so each walk will stop
 *              with a specific walk length
 *
 * @initWalks(): put each walk on its source vertex by encoding it to 64-bit and
 *               store it into the corresponding block
 * @updateWalk(): update one walk until it reaches the boundary of the block
 * */
class SecondOrderRW: public RandomWalk {
public:
    std::vector<metrics> metric;

    bid_t nBlocks;
    vid_t *startVertex;
    hid_t walkLength;

    int csrFile;

    VertexIO vertexIO;

#if DEBUG
    uint64_t totalSteps = 0;
    std::vector<uint64_t> totalSteps_t;
#endif

    SecondOrderRW(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int begPosFile, int csrFile,
                  vid_t nTotalVertices);

    virtual vid_t sampleStartEdge(metrics &m, VertexInfo &srcVertex, bool &noOutEdge) = 0;

    void initWalks(vid_t nVertices, WalkManager &walkManager, metrics &m);
#if BI_BLOCK
    virtual void initWalks(bid_t block, WalkManager &walkManager, metrics &m) = 0;
#endif

    void updateWalk(metrics &m, const WalkDataType &walk, bid_t exec_block, eid_t *beg_pos, vid_t *csr,
                    WalkManager &walkManager);

#if BI_BLOCK
    virtual void
    updateWalk(metrics &m, WalkDataType &walk, bid_t staticBlock, bid_t dynamicBlock, WalkManager &walkManager,
               bool dynamicBlockInMem, bool &walkFinished) {
        logstream(LOG_ERROR) << "No definition of function : updateWalk!" << std::endl;
    }

    virtual void
    updateWalk(metrics &m, WalkDataType &walk, bid_t staticBlock, bid_t dynamicBlock, WalkManager &walkManager,
               bool dynamicBlockInMem, bool staticBlockInMem, bool &walkFinished) {
        logstream(LOG_ERROR) << "No definition of function : updateWalk!" << std::endl;
    }

    virtual void
    updateWalk(metrics &m, WalkDataType &walk, bid_t staticBlock, bid_t dynamicBlock, WalkManager &walkManager,
               bool &walkFinished) {
        logstream(LOG_ERROR) << "No definition of function : updateWalk!" << std::endl;
    }

    virtual void
    processWalk(metrics &m, WalkDataType &walk, bid_t staticBlock, bid_t dynamicBlock, WalkManager &walkManager, bool dynamicBlockInMem, bool staticBlockInMem);

    virtual void
    processWalk(metrics &m, WalkDataType &walk, bid_t staticBlock, bid_t dynamicBlock, WalkManager &walkManager);
#endif

    bid_t getBlock(vid_t v);

    static inline eid_t getVertexOutDegree(vid_t currentVertex, eid_t *beg_pos);

    bool allWalksFinished(WalkManager &walkManager);

    virtual double rejectAccRatio(VertexInfo &preV, VertexInfo &destV, metrics &m);

    virtual double rejectAccRatio(VertexInfo &preV, VertexInfo &curV, VertexInfo &destV, metrics &m){
        logstream(LOG_ERROR) << "Reject accept ratio is not defined" << std::endl;
        return 0;
    }

    void setIOInfo(bid_t _nInMemBlocks, vid_t **_csrBuf, eid_t **_begPosBuf, bid_t *_inMemIndex);
};

SecondOrderRW::SecondOrderRW(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int begPosFile, int csrFile,
                             vid_t nTotalVertices)
                             : vertexIO(begPosFile, csrFile, nBlocks, startVertex, nTotalVertices) {
    this->nBlocks = nBlocks;
    this->startVertex = startVertex;
    this->walkLength = walkLength;
    this->csrFile = csrFile;
#if DEBUG
    totalSteps_t.assign(100, 0);
#endif
}

void SecondOrderRW::initWalks(vid_t nVertices, WalkManager &walkManager, metrics &m) {
    logstream(LOG_ERROR) << "No definition of function : initWalks!" << std::endl;
}

void
SecondOrderRW::updateWalk(metrics &m, const WalkDataType &walk, bid_t exec_block, eid_t *beg_pos, vid_t *csr,
                          WalkManager &walkManager) {
    logstream(LOG_ERROR) << "No definition of function : updateWalk!" << std::endl;
}

bid_t SecondOrderRW::getBlock(vid_t v) {
    for (bid_t i = 0; i < nBlocks; i++){
        if (startVertex[i + 1] > v){
            return i;
        }
    }
    logstream(LOG_ERROR) << "Vertex" << v << "out of range!" << std::endl;
    abort();
}

bool SecondOrderRW::allWalksFinished(WalkManager &walkManager) {
    wid_t remainWalkNum = walkManager.currentNWalks;
    return (remainWalkNum == 0);
}


eid_t SecondOrderRW::getVertexOutDegree(vid_t currentVertex, eid_t *beg_pos) {
    return beg_pos[currentVertex + 1] - beg_pos[currentVertex];
}

double SecondOrderRW::rejectAccRatio(VertexInfo &preV, VertexInfo &destV, metrics &m) {
    logstream(LOG_ERROR) << "Reject accept ratio is not defined" << std::endl;
    return 0;
}

void SecondOrderRW::setIOInfo(bid_t _nInMemBlocks, vid_t **_csrBuf, eid_t **_begPosBuf, bid_t *_inMemIndex) {
    vertexIO.setIOInfo(_nInMemBlocks, _csrBuf, _begPosBuf, _inMemIndex);
}

#if BI_BLOCK
void SecondOrderRW::processWalk(metrics &m, WalkDataType &walk, bid_t staticBlock, bid_t dynamicBlock,
                                WalkManager &walkManager, bool dynamicBlockInMem, bool staticBlockInMem) {
    logstream(LOG_ERROR) << "Function processWalk is not defined" << std::endl;
}
void SecondOrderRW::processWalk(metrics &m, WalkDataType &walk, bid_t staticBlock, bid_t dynamicBlock,
                                WalkManager &walkManager) {

}
#endif

#endif //GRAPHWALKER_SECONDORDERRW_HPP
