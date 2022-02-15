/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/4/23 17:04
 */

#ifndef IOE_SORW_NODE2VECARREJECTION_HPP
#define IOE_SORW_NODE2VECARREJECTION_HPP

#include <vector>
#include <algorithm>
#include "BasicIncludes.hpp"
#include "util/RandNum.hpp"
#include "util/toplist.hpp"
#include "sampler/FirstOrder/AliasSampler.hpp"
#include "sampler/SecondOrder/RejectSampler.hpp"
#include "walks/Node2VecAR.hpp"

class Node2VecARReject : public Node2VecAR{
public:
    Node2VecARReject(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int walkNumPerVertex, eid_t maxOutDegree,
                     double alpha, int begPosFile, int csrFile, vid_t nTotalVertices, metrics &m,
                     const double &decayFactor);

    ~Node2VecARReject();

    vid_t sampleStartEdge(metrics &m, VertexInfo &srcVertex, bool &noOutEdge) override;
#if (PLAIN | PLAIN_BUCKET)
    vid_t sampleDestVertex(metrics &m, vid_t preVertex, vid_t curVertex, eid_t *beg_pos,
                           vid_t *csr, bid_t execBlock) override;
    void initWalks(vid_t nVertices, WalkManager &walkManager, metrics &m) override;
#endif
#if BI_BLOCK
    vid_t sampleDestVertex(metrics &m, VertexInfo &preVertexInfo, VertexInfo &curVertexInfo) override {
        return soRejectSampler._sampleDestVertex(m, preVertexInfo, curVertexInfo, alias1StOrderSampler,
                                                &SecondOrderRW::rejectAccRatio,
                                                *this);
    }

    void initWalks(bid_t block, WalkManager &walkManager, metrics &m) override{

        tid_t nThreads = get_option_int("nThreads", omp_get_max_threads());
        omp_set_num_threads(nThreads);
        m.start_time("initWalks");
        vid_t startVertexId = startVertex[block];
        vid_t endVertexId = startVertex[block + 1];
#if MULTI_THREAD
#pragma omp parallel for schedule(dynamic)
#endif
        for (vid_t sourceVertex = startVertexId; sourceVertex < endVertexId; sourceVertex++){
            VertexInfo srcVertex(sourceVertex);
            bool noOutEdge = false;
            for (wid_t w = 0; w < nWalkPerVertex; w++){
                /* the absolute id of current vertex */
                vid_t currentVertex_abs = sampleStartEdge(m, srcVertex, noOutEdge);
                if (noOutEdge){
                    break;
                }
                bid_t preVertexResideBlock = block;
                bid_t resideBlock = getBlock(currentVertex_abs);
                vid_t currentVertex = currentVertex_abs - startVertex[resideBlock];
                WalkDataType walk = WalkManager::encode(sourceVertex, sourceVertex, currentVertex, preVertexResideBlock, resideBlock, 0);
                bool walkFinished = false;
                if (resideBlock == block){
                    updateWalk(m, walk, block, block, walkManager, true, true, walkFinished);
                }
                if (walkFinished){
                    continue;
                }
                preVertexResideBlock = WalkManager::getPreBlock(walk);
                resideBlock = WalkManager::getCurBlock(walk);
                bid_t toBlock = preVertexResideBlock < resideBlock ? preVertexResideBlock : resideBlock;
                walkManager.moveWalk(walk, omp_get_thread_num(), toBlock);
            }
        }
        m.stop_time("initWalks");
    }
#endif


    double rejectAccRatio(VertexInfo &preV, VertexInfo &curV, VertexInfo &destV, metrics &m) override;

private:
    double *nodeProbTbl;
    Alias1stOrderSampler alias1StOrderSampler;
    SORejectSampler soRejectSampler;
};

Node2VecARReject::Node2VecARReject(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int walkNumPerVertex,
                                   eid_t maxOutDegree,
                                   double alpha, int begPosFile, int csrFile, vid_t nTotalVertices, metrics &m,
                                   const double &decayFactor)
        : Node2VecAR(nBlocks, startVertex, walkLength, walkNumPerVertex, alpha, begPosFile, csrFile, maxOutDegree,
                     nTotalVertices, decayFactor),
          alias1StOrderSampler(&vertexIO, m, true) {
    nodeProbTbl = new double[maxOutDegree];
}

#if (PLAIN | PLAIN_BUCKET)
vid_t
Node2VecARReject::sampleDestVertex(metrics &m, vid_t preVertex, vid_t curVertex, eid_t *beg_pos,
                                 vid_t *csr, bid_t execBlock) {
    return soRejectSampler._sampleDestVertex(m, preVertex, curVertex, alias1StOrderSampler,
                                            &SecondOrderRW::rejectAccRatio,
                                            *this, execBlock);
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
void Node2VecARReject::initWalks(vid_t nVertices, WalkManager &walkManager, metrics &m) {
    std::cout << "Node2Vec using autoregressive sampling start with reject sampler" << std::endl;

    tid_t nThreads = get_option_int("nThreads", omp_get_max_threads());
    omp_set_num_threads(nThreads);
    m.start_time("initWalks");
    walkManager.currentNWalks = nVertices * nWalkPerVertex;
#pragma omp parallel for schedule(static)
    for (vid_t sourceVertex = 0; sourceVertex < nVertices; sourceVertex++){
        VertexInfo srcVertex(sourceVertex);
        bool noOutEdge = false;
        for (wid_t w = 0; w < nWalkPerVertex; w++){
            /* the absolute id of current vertex */
            vid_t currentVertex_abs = sampleStartEdge(m, srcVertex, noOutEdge);
            if (noOutEdge){
                walkManager.currentNWalks -= 1;
                continue;
            }
            bid_t preVertexResideBlock = getBlock(sourceVertex);
            bid_t resideBlock = getBlock(currentVertex_abs);
//            if (resideBlock == getBlock(sourceVertex)){
//                m.add("edge-in",1);
//            }else{
//                m.add("edge-cut", 1);
//            }
            vid_t currentVertex = currentVertex_abs - startVertex[resideBlock];
            WalkDataType walk = WalkManager::encode(sourceVertex, sourceVertex, currentVertex, preVertexResideBlock, resideBlock, 0);

            walkManager.moveWalk(walk, resideBlock, omp_get_thread_num(), sourceVertex, currentVertex, 0, preVertexResideBlock);
        }
//        if (!noOutEdge){
//            vertexIO.vertexBuffer.streamV(srcVertex);
//        }
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
#endif

Node2VecARReject::~Node2VecARReject() {
    delete [] nodeProbTbl;
}

vid_t Node2VecARReject::sampleStartEdge(metrics &m, VertexInfo &srcVertex, bool &noOutEdge) {
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

double Node2VecARReject::rejectAccRatio(VertexInfo &preV, VertexInfo &curV, VertexInfo &destV, metrics &m) {
    double bw = biasedWeight(preV, curV, destV, m);
    double result = bw / ((1 - alpha) / curV.outDegree + alpha / preV.outDegree);
    return bw / ((1 - alpha) / curV.outDegree + alpha / preV.outDegree);
}


#pragma clang diagnostic pop
#endif //IOE_SORW_NODE2VECARREJECTION_HPP
