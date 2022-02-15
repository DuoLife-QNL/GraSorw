/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/8/10 11:49
 */

#ifndef IOE_SORW_FIXEDLENGTHFORW_HPP
#define IOE_SORW_FIXEDLENGTHFORW_HPP
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
class FixedLengthFORW: public RandomWalk {
public:
    std::vector<metrics> metric;

    bid_t nBlocks;
    vid_t *startVertex;
    hid_t walkLength;

    int csrFile;

    VertexIO vertexIO;

    FixedLengthFORW(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int begPosFile, int csrFile,
                  vid_t nTotalVertices);

    virtual vid_t sampleStartEdge(metrics &m, VertexInfo &srcVertex, bool &noOutEdge) = 0;

    void initWalks(vid_t nVertices, WalkManager &walkManager, metrics &m);

    void updateWalk(metrics &m, const WalkDataType &walk, bid_t exec_block, eid_t *beg_pos, vid_t *csr,
                    WalkManager &walkManager);

    bid_t getBlock(vid_t v);

    static inline eid_t getVertexOutDegree(vid_t currentVertex, eid_t *beg_pos);

    bool allWalksFinished(WalkManager &walkManager);

    void setIOInfo(bid_t _nInMemBlocks, vid_t **_csrBuf, eid_t **_begPosBuf, bid_t *_inMemIndex);

    void getGraphInfo_lineGraph(vid_t &nVertices, eid_t &nEdges){
        eid_t edgeNum = 0;
        vid_t verNum = 0;
        for (vid_t v = 0; v < vertexIO.getNTotalVertices(); v++){
            VertexInfo vertexInfo(v);
            vertexIO.getVertexOutDegree(vertexInfo, AUTO);
            verNum += vertexInfo.outDegree;
            if (vertexInfo.outDegree == 1){
                continue;
            }else{
                edgeNum += (vertexInfo.outDegree * (vertexInfo.outDegree - 1)) / 2;
            }
        }
        nVertices = verNum / 2;
        nEdges = edgeNum;
    }

};

FixedLengthFORW::FixedLengthFORW(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int begPosFile, int csrFile,
                             vid_t nTotalVertices)
                             : vertexIO(begPosFile, csrFile, nBlocks, startVertex, nTotalVertices) {
    this->nBlocks = nBlocks;
    this->startVertex = startVertex;
    this->walkLength = walkLength;
    this->csrFile = csrFile;
}

void FixedLengthFORW::initWalks(vid_t nVertices, WalkManager &walkManager, metrics &m) {
    logstream(LOG_ERROR) << "No definition of function : initWalks!" << std::endl;
}

void
FixedLengthFORW::updateWalk(metrics &m, const WalkDataType &walk, bid_t exec_block, eid_t *beg_pos, vid_t *csr,
                          WalkManager &walkManager) {
    logstream(LOG_ERROR) << "No definition of function : updateWalk!" << std::endl;
}

bid_t FixedLengthFORW::getBlock(vid_t v) {
    for (bid_t i = 0; i < nBlocks; i++){
        if (startVertex[i + 1] > v){
            return i;
        }
    }
    logstream(LOG_ERROR) << "Vertex" << v << "out of range!" << std::endl;
    abort();
}

bool FixedLengthFORW::allWalksFinished(WalkManager &walkManager) {
    wid_t remainWalkNum = walkManager.currentNWalks;
    return (remainWalkNum == 0);
}


eid_t FixedLengthFORW::getVertexOutDegree(vid_t currentVertex, eid_t *beg_pos) {
    return beg_pos[currentVertex + 1] - beg_pos[currentVertex];
}


void FixedLengthFORW::setIOInfo(bid_t _nInMemBlocks, vid_t **_csrBuf, eid_t **_begPosBuf, bid_t *_inMemIndex) {
    vertexIO.setIOInfo(_nInMemBlocks, _csrBuf, _begPosBuf, _inMemIndex);
}


#endif //IOE_SORW_FIXEDLENGTHFORW_HPP
