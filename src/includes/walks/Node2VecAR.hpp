/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/4/23 17:02
 *
 * node2vec using autoregressive sampling model
 */

#ifndef IOE_SORW_NODE2VECAR_HPP
#define IOE_SORW_NODE2VECAR_HPP
#include "WalkManager.hpp"
#include "SecondOrderRW.hpp"
#include "IO/VertexIO.hpp"

class Node2VecAR : public SecondOrderRW {
protected:
    double alpha;
    double decayFactor;

    std::vector<RandNum> randNums;

public:

    int nWalkPerVertex;

    Node2VecAR(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int walkNumPerVertex, double _alpha, int begPosFile,
               int csrFile, eid_t maxOutDegree, vid_t nTotalVertices, const double &_decayFactor);

    double biasedWeight(VertexInfo &preV, VertexInfo &curV, VertexInfo &forwardV, metrics &m){
        assert(curV.outDegree);
        if (hasEdge(preV, forwardV.vertexId, m)){
            assert(preV.outDegree);
            return (1 - alpha) / curV.outDegree + alpha / preV.outDegree;
        }else{
            return (1 - alpha) / curV.outDegree;
        }
    }

    bool hasEdge(VertexInfo &srcVertex, vid_t destVertex, metrics &m);
#if BI_BLOCK
    void updateWalk(metrics &m, WalkDataType &walk, bid_t staticBlock, bid_t dynamicBlock, WalkManager &walkManager,
                    bool dynamicBlockInMem, bool staticBlockInMem, bool &walkFinished) override {
        bid_t curBlock = WalkManager::getCurBlock(walk);
        VertexInfo curVertexInfo(WalkManager::getCurrentVertex(walk) + startVertex[curBlock], curBlock, true);
        VertexInfo preVertexInfo(WalkManager::getPreviousVertex(walk), WalkManager::getPreBlock(walk), true);
        if (curBlock == dynamicBlock){
            curVertexInfo.setResideBlockInMem(dynamicBlockInMem);
            preVertexInfo.setResideBlockInMem(staticBlockInMem);
        }else{
            preVertexInfo.setResideBlockInMem(dynamicBlockInMem);
            curVertexInfo.setResideBlockInMem(staticBlockInMem);
        }
        hid_t hops = WalkManager::getHops(walk);
        vid_t nextVertex = sampleDestVertex(m, preVertexInfo, curVertexInfo);
        hops++;
#if DEBUG
        tid_t t = omp_get_thread_num();
        totalSteps_t.at(t)++;
#endif
        while (true){
            if (hops >= walkLength){
                walkFinished = true;
                return;
            }
            preVertexInfo.resetVertex(curVertexInfo.vertexId, curVertexInfo.resideBlockId, curVertexInfo.isResideBlockInMem());
            if (vertexIO.vertexInBlock(nextVertex, staticBlock)){
                curVertexInfo.resetVertex(nextVertex, staticBlock, staticBlockInMem);
            }else if(vertexIO.vertexInBlock(nextVertex, dynamicBlock)){
                curVertexInfo.resetVertex(nextVertex, dynamicBlock, dynamicBlockInMem);
            }else{
                break;
            }

            nextVertex = sampleDestVertex(m, preVertexInfo, curVertexInfo);
            hops++;
#if DEBUG
            totalSteps_t.at(t)++;
#endif
        }
        curVertexInfo.resetVertex(nextVertex, getBlock(nextVertex));
        walk = WalkManager::encode(WalkManager::getSouceVertex(walk), preVertexInfo.vertexId, curVertexInfo.vertexId - startVertex[curVertexInfo.resideBlockId], preVertexInfo.resideBlockId, curVertexInfo.resideBlockId, hops);
        //        wid_t storeBlock = preVertexInfo.vertexId < curVertexInfo.vertexId ? preVertexInfo.vertexId : curVertexInfo.vertexId;
        //        walkManager.moveWalk(walk, threadId, storeBlock);
        walkFinished = false;
    }

    void processWalk(metrics &m, WalkDataType &walk, bid_t staticBlock, bid_t dynamicBlock,
                     WalkManager &walkManager, bool dynamicBlockInMem, bool staticBlockInMem) override{
        bool walkFinished = false;

        updateWalk(m, walk, staticBlock, dynamicBlock, walkManager, dynamicBlockInMem, staticBlockInMem, walkFinished);
        if (walkFinished){
            return;
        }
        tid_t threadId = omp_get_thread_num();
        bid_t cur = WalkManager::getCurBlock(walk);
        bid_t pre = WalkManager::getPreBlock(walk);
        if (cur < staticBlock){
            walkManager.moveWalk(walk, threadId, cur);
        }else if(cur > staticBlock && cur < dynamicBlock){
            if (pre == staticBlock){
                walkManager.moveWalk(walk, threadId, staticBlock);
            }else{
                walkManager.moveWalk(walk, threadId, cur);
            }
        }else{
            assert(cur > dynamicBlock);
            if (pre == staticBlock){
#if MULTI_BUFFER
                walkManager.extendBucket[cur][threadId].push_back(walk);
#else
                walkManager.moveWalk(walk, threadId, staticBlock);
#endif
            }else{
                walkManager.moveWalk(walk, threadId, dynamicBlock);
            }
        }
    }


    virtual vid_t sampleDestVertex(metrics &m, VertexInfo &preVertexInfo, VertexInfo &curVertexInfo);
#endif

#if (PLAIN | PLAIN_BUCKET)
    void updateWalk(metrics &m, const WalkDataType &walk, bid_t execBlock, eid_t *beg_pos, vid_t *csr,
                    WalkManager &walkManager) override;
    virtual vid_t sampleDestVertex(metrics &m, vid_t preVertex, vid_t curVertex, eid_t *beg_pos, vid_t *csr, bid_t execBlock);
#endif
};

Node2VecAR::Node2VecAR(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int walkNumPerVertex, double _alpha,
                       int begPosFile, int csrFile, eid_t maxOutDegree, vid_t nTotalVertices, const double &_decayFactor)
        : SecondOrderRW(nBlocks, startVertex, walkLength, begPosFile, csrFile, nTotalVertices) {
    this->nWalkPerVertex = walkNumPerVertex;
    alpha = _alpha;
    decayFactor = _decayFactor;
    randNums.resize(100, RandNum(66532897483413));

    tid_t nThreads = omp_get_max_threads();
    for (tid_t t = 0; t < nThreads; t++){
        metrics mt("thread" + std::to_string(t));
        metric.emplace_back(mt);
    }
}

bool Node2VecAR::hasEdge(VertexInfo &srcVertex, vid_t destVertex, metrics &m) {

    if (!srcVertex.csr){
        tid_t t = omp_get_thread_num();
#if TIME_COST_INFO
        vertexIO.metric.at(t).start_time("get-csr");
#endif
#if PLAIN_BUCKET
        vertexIO.getPreVertexCSR(srcVertex, m);
#else
        vertexIO.getVertexCSR(srcVertex, m);
#endif
#if TIME_COST_INFO
        vertexIO.metric.at(t).stop_time("get-csr");
#endif
    }

//    tid_t t = omp_get_thread_num();
//    vertexIO.metric.at(t).start_time("bin-check");
    vid_t left = 0;
    vid_t right = srcVertex.outDegree;
    vid_t mid;
    while (left < right){
        mid = left + (right - left) / 2;
        if (srcVertex.csr[mid] == destVertex){
//            vertexIO.metric.at(t).stop_time("bin-check");
            return true;
        }else if (destVertex > srcVertex.csr[mid]){
            left = mid + 1;
        }else{
            right = mid;
        }
    }
//    vertexIO.metric.at(t).stop_time("bin-check");
    return false;

//    for (vid_t e = 0; e < srcVertex.outDegree; e++){
//        if (srcVertex.csr[e] == destVertex){
//            return true;
//        }
//    }
//    return false;
}

#if (PLAIN | PLAIN_BUCKET)
void Node2VecAR::updateWalk(metrics &m, const WalkDataType &walk, bid_t execBlock, eid_t *beg_pos, vid_t *csr,
                     WalkManager &walkManager) {
    tid_t threadId = omp_get_thread_num();
    vid_t currentVertex = WalkManager::getCurrentVertex(walk);
    vid_t previousVertex = WalkManager::getPreviousVertex(walk);
    hid_t previousNHops = WalkManager::getHops(walk);
    hid_t forwardHops = 0;

    /* change currentVertex to absolute value*/
    currentVertex += startVertex[execBlock];
    while (currentVertex >= startVertex[execBlock] && currentVertex < startVertex[execBlock + 1] && (previousNHops + forwardHops) < walkLength){
//        metric.at(threadId).start_time("sampling");
        vid_t nextVertex = sampleDestVertex(m, previousVertex, currentVertex, beg_pos, csr, execBlock);
//        metric.at(threadId).stop_time("sampling");
        /* NOTE: currentVertex is absolute value here */
        previousVertex = currentVertex;
        currentVertex = nextVertex;
        forwardHops++;
    }
    if (previousNHops + forwardHops < walkLength){
        bid_t resideBlock = getBlock(currentVertex);
        walkManager.moveWalk(walk, resideBlock, threadId, previousVertex, currentVertex - startVertex[resideBlock],
                             forwardHops, execBlock);
        walkManager.setMinStep(resideBlock, previousNHops + forwardHops);
        walkManager.isModified[resideBlock] = true;
    }
}

vid_t Node2VecAR::sampleDestVertex(metrics &m, vid_t preVertex, vid_t curVertex, eid_t *beg_pos, vid_t *csr, bid_t execBlock) {
    logstream(LOG_ERROR) << "sample destination vertex function not implemented" << std::endl;
    return 0;
}
#endif

#if BI_BLOCK
vid_t Node2VecAR::sampleDestVertex(metrics &m, VertexInfo &preVertexInfo, VertexInfo &curVertexInfo) {
    return 0;
}
#endif


#endif //IOE_SORW_NODE2VECAR_HPP
