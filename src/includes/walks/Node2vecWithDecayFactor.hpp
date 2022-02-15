/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/9/27 16:13
 */

#ifndef IOE_SORW_NODE2VECWITHDECAYFACTOR_HPP
#define IOE_SORW_NODE2VECWITHDECAYFACTOR_HPP

#include "WalkManager.hpp"
#include "SecondOrderRW.hpp"
#include "IO/VertexIO.hpp"

class Node2VecWithDecayFactor : public SecondOrderRW {
private:
    double p;
    double q;
    double decayFactor;
    std::vector<RandNum> randNums;

public:

    Node2VecWithDecayFactor(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, double p, double q, double decayFactor,
                            int begPosFile, int csrFile, eid_t maxOutDegree, vid_t nTotalVertices);

    double biasedWeight(VertexInfo &preV, VertexInfo &forwardV, metrics &m);

    bool hasEdge(VertexInfo &srcVertex, vid_t destVertex, metrics &m);
#if (PLAIN | PLAIN_BUCKET)
    void updateWalk(metrics &m, const WalkDataType &walk, bid_t execBlock, eid_t *beg_pos, vid_t *csr,
                    WalkManager &walkManager) override;
    virtual vid_t sampleDestVertex(metrics &m, vid_t preVertex, vid_t curVertex, eid_t *beg_pos, vid_t *csr, bid_t execBlock);
#endif
#if BI_BLOCK
    void updateWalk(metrics &m, WalkDataType &walk, bid_t staticBlock, bid_t dynamicBlock, WalkManager &walkManager,
                    bool dynamicBlockInMem, bool &walkFinished) override {
        tid_t t = omp_get_thread_num();
        if (randNums.at(t).dRand() > decayFactor){
            walkFinished = true;
            return;
        }
        bid_t curBlock = WalkManager::getCurBlock(walk);
        VertexInfo curVertexInfo(WalkManager::getCurrentVertex(walk), curBlock, true);
        VertexInfo preVertexInfo(WalkManager::getPreviousVertex(walk), WalkManager::getPreBlock(walk), true);
        if (curBlock == dynamicBlock){
            curVertexInfo.setResideBlockInMem(dynamicBlockInMem);
        }else{
            preVertexInfo.setResideBlockInMem(dynamicBlockInMem);
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
                curVertexInfo.resetVertex(nextVertex, staticBlock, true);
            }else if(vertexIO.vertexInBlock(nextVertex, dynamicBlock)){
                curVertexInfo.resetVertex(nextVertex, dynamicBlock, dynamicBlockInMem);
            }else{
                break;
            }
            if (randNums.at(t).dRand() > decayFactor){
                walkFinished = true;
#if DEBUG
                m.add("finished-walk", 1, INTEGER);
#endif
                return;
            }

            nextVertex = sampleDestVertex(m, preVertexInfo, curVertexInfo);
            hops++;
#if DEBUG
            totalSteps_t.at(t)++;
#endif
        }
        curVertexInfo.resetVertex(nextVertex, getBlock(nextVertex));
        walk = WalkManager::encode(WalkManager::getSouceVertex(walk), preVertexInfo.vertexId, curVertexInfo.vertexId, preVertexInfo.resideBlockId, curVertexInfo.resideBlockId, hops);
//        wid_t storeBlock = preVertexInfo.vertexId < curVertexInfo.vertexId ? preVertexInfo.vertexId : curVertexInfo.vertexId;
//        walkManager.moveWalk(walk, threadId, storeBlock);
        walkFinished = false;
    }

    void updateWalk(metrics &m, WalkDataType &walk, bid_t staticBlock, bid_t dynamicBlock, WalkManager &walkManager,
                    bool dynamicBlockInMem, bool staticBlockInMem, bool &walkFinished) override {
        tid_t t = omp_get_thread_num();
        if (randNums.at(t).dRand() > decayFactor){
            walkFinished = true;
            return;
        }
        bid_t curBlock = WalkManager::getCurBlock(walk);
        VertexInfo curVertexInfo(WalkManager::getCurrentVertex(walk), curBlock, true);
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
#if BLOCK_UTI | ACT_VER
        vertexIO.markActivated(curVertexInfo.vertexId);
        vertexIO.markActivated(preVertexInfo.vertexId);
#endif
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
            if (randNums.at(t).dRand() > decayFactor){
                walkFinished = true;
#if DEBUG
                m.add("finished-walk", 1, INTEGER);
#endif
                return;
            }

            nextVertex = sampleDestVertex(m, preVertexInfo, curVertexInfo);
#if BLOCK_UTI | ACT_VER
            vertexIO.markActivated(curVertexInfo.vertexId);
#endif
            hops++;
#if DEBUG
            totalSteps_t.at(t)++;
#endif
        }
        curVertexInfo.resetVertex(nextVertex, getBlock(nextVertex));
        walk = WalkManager::encode(WalkManager::getSouceVertex(walk), preVertexInfo.vertexId, curVertexInfo.vertexId, preVertexInfo.resideBlockId, curVertexInfo.resideBlockId, hops);
        //        wid_t storeBlock = preVertexInfo.vertexId < curVertexInfo.vertexId ? preVertexInfo.vertexId : curVertexInfo.vertexId;
        //        walkManager.moveWalk(walk, threadId, storeBlock);
        walkFinished = false;
    }

    void updateWalk(metrics &m, WalkDataType &walk, bid_t staticBlock, bid_t dynamicBlock, WalkManager &walkManager,
                    bool &walkFinished) override {
        bid_t curBlock = WalkManager::getCurBlock(walk);
        VertexInfo curVertexInfo(WalkManager::getCurrentVertex(walk) + startVertex[curBlock], curBlock, true);
        VertexInfo preVertexInfo(WalkManager::getPreviousVertex(walk), WalkManager::getPreBlock(walk), true);
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
            preVertexInfo.resetVertex(curVertexInfo.vertexId, curVertexInfo.resideBlockId, true);
            if (vertexIO.vertexInBlock(nextVertex, staticBlock)){
                curVertexInfo.resetVertex(nextVertex, staticBlock, true);
            }else if(vertexIO.vertexInBlock(nextVertex, dynamicBlock)){
                curVertexInfo.resetVertex(nextVertex, dynamicBlock, true);
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
                if (WalkManager::getPreviousVertex(walk) == 548142 && WalkManager::getCurrentVertex(walk) == 77859511){
                    int a = 1;
                }
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
};

Node2VecWithDecayFactor::Node2VecWithDecayFactor(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, double p,
                                                 double q, double decayFactor,
                                                 int begPosFile, int csrFile, eid_t maxOutDegree, vid_t nTotalVertices)
        : SecondOrderRW(nBlocks, startVertex, walkLength, begPosFile, csrFile, nTotalVertices) {
    this->p = p;
    this->q = q;
    this->decayFactor = decayFactor;

    tid_t nThreads = omp_get_max_threads();
    for (tid_t t = 0; t < nThreads; t++){
        metrics mt("thread" + std::to_string(t));
        metric.emplace_back(mt);
    }

    randNums.resize(100, RandNum(66532897483413));
    for (tid_t t = 0; t < nThreads; t++){
        metrics mt("thread" + std::to_string(t));
        metric.emplace_back(mt);
    }
}
#if (PLAIN | PLAIN_BUCKET)
void
Node2VecWithDecayFactor::updateWalk(metrics &m, const WalkDataType &walk, bid_t execBlock, eid_t *beg_pos, vid_t *csr,
                     WalkManager &walkManager) {
    tid_t threadId = omp_get_thread_num();
    vid_t currentVertex = walkManager.getCurrentVertex(walk);
    vid_t previousVertex = walkManager.getPreviousVertex(walk);
    hid_t previousNHops = walkManager.getHops(walk);
    hid_t forwardHops = 0;

    bool walkFinished = false;
    while (currentVertex >= startVertex[execBlock] && currentVertex < startVertex[execBlock + 1] && (previousNHops + forwardHops) < walkLength){
//        metric.at(threadId).start_time("sampling");
        if(randNums.at(threadId).dRand() > decayFactor){
            walkFinished = true;
            break;
        }
        vid_t nextVertex = sampleDestVertex(m, previousVertex, currentVertex, beg_pos, csr, execBlock);
//        metric.at(threadId).stop_time("sampling");
        /* NOTE: currentVertex is absolute value here */
        previousVertex = currentVertex;
        currentVertex = nextVertex;
        forwardHops++;
    }
    if (walkFinished){
        return;
    }
    if (previousNHops + forwardHops < walkLength){
        bid_t resideBlock = getBlock(currentVertex);
        walkManager.moveWalk(walk, resideBlock, threadId, previousVertex, currentVertex,
                             forwardHops, execBlock);
        walkManager.setMinStep(resideBlock, previousNHops + forwardHops);
        walkManager.isModified[resideBlock] = true;
    }
}
#endif

double Node2VecWithDecayFactor::biasedWeight(VertexInfo &preV, VertexInfo &forwardV, metrics &m) {
    if (preV.vertexId == forwardV.vertexId){
        return 1.0 / p;
    }
    if (hasEdge(preV, forwardV.vertexId, m)){
        return 1.0;
    }
    return 1.0 / q;
}


bool Node2VecWithDecayFactor::hasEdge(VertexInfo &srcVertex, vid_t destVertex, metrics &m) {
    if (!srcVertex.csr){
//        assert(srcVertex.resideBlockId == staticBlock_global || srcVertex.resideBlockId == dynamicBlock_global);
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
vid_t
Node2VecWithDecayFactor::sampleDestVertex(metrics &m, vid_t preVertex, vid_t curVertex, eid_t *beg_pos, vid_t *csr, bid_t execBlock) {
    logstream(LOG_ERROR) << "sample destination vertex function not implemented" << std::endl;
    return 0;
}
#endif

#if BI_BLOCK
vid_t Node2VecWithDecayFactor::sampleDestVertex(metrics &m, VertexInfo &preVertexInfo, VertexInfo &curVertexInfo) {
    return 0;
}
#endif

#endif //IOE_SORW_NODE2VECWITHDECAYFACTOR_HPP
