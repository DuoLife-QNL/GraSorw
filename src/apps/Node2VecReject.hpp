//
// Created by lihz on 2020/11/4.
//

#ifndef IOE_SORW_NODE2VECREJECT_HPP
#define IOE_SORW_NODE2VECREJECT_HPP


#include <vector>
#include <algorithm>

#include "BasicIncludes.hpp"
#include "walks/Node2Vec.hpp"
#include "util/RandNum.hpp"
#include "util/toplist.hpp"
#include "sampler/FirstOrder/AliasSampler.hpp"
#include "sampler/SecondOrder/RejectSampler.hpp"

class Node2VecReject : public Node2Vec{
public:
    int nWalkPerVertex;
    Node2VecReject(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int walkNumPerVertex,
                   eid_t maxOutDegree, double p, double q, int begPosFile, int csrFile,
                   vid_t nTotalVertices, metrics &m);

    ~Node2VecReject();

    vid_t sampleStartEdge(metrics &m, VertexInfo &srcVertex, bool &noOutEdge) override;

#if BI_BLOCK
    vid_t sampleDestVertex(metrics &m, VertexInfo &preVertexInfo, VertexInfo &curVertexInfo) override {
        tid_t t = omp_get_thread_num();
        return soRejectSampler.sampleDestVertex(m, preVertexInfo, curVertexInfo, aliasSamplers.at(t),
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
            VertexInfo srcVertex(sourceVertex, block, true);
            staticBlock_global = block;
            bool noOutEdge = false;
            for (wid_t w = 0; w < nWalkPerVertex; w++){
                /* the absolute id of current vertex */
                vid_t currentVertex_abs = sampleStartEdge(m, srcVertex, noOutEdge);
                if (noOutEdge){
                    break;
                }
#if COMPUTE_EDGE_CUT
                /* compute edge-cut */
                for (eid_t e = 0; e < srcVertex.outDegree; e++){
                    if (vertexIO.vertexInBlock(srcVertex.csr[e], block)){
                        edgeIn++;
                    }else{
                        edgeCut++;
                    }
                }
#endif

                bid_t preVertexResideBlock = block;
                bid_t resideBlock = getBlock(currentVertex_abs);
                vid_t currentVertex = currentVertex_abs;
                WalkDataType walk = WalkManager::encode(sourceVertex, sourceVertex, currentVertex, preVertexResideBlock, resideBlock, 0);
                vid_t v = WalkManager::getCurrentVertex(walk);
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
#if (PLAIN | PLAIN_BUCKET)
    void initWalks(vid_t nVertices, WalkManager &walkManager, metrics &m) override;
    vid_t sampleDestVertex(metrics &m, vid_t preVertex, vid_t curVertex, eid_t *beg_pos,
                           vid_t *csr, bid_t execBlock) override;
#endif
    double rejectAccRatio(VertexInfo &preV, VertexInfo &destV, metrics &m) override ;

private:
    double *nodeProbTbl;
    Alias1stOrderSampler alias1StOrderSampler;
    std::vector<Alias1stOrderSampler> aliasSamplers;
    SORejectSampler soRejectSampler;
    double accRatioConst;

    static inline double calAccRatioConst(double p, double q);
};

Node2VecReject::Node2VecReject(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int walkNumPerVertex,
                               eid_t maxOutDegree, double p, double q, int begPosFile, int csrFile,
                               vid_t nTotalVertices, metrics &m)
        : Node2Vec(nBlocks, startVertex, walkLength, p, q, begPosFile, csrFile, maxOutDegree, nTotalVertices),
          alias1StOrderSampler(&vertexIO, m, true) {
    nWalkPerVertex = walkNumPerVertex;
    nodeProbTbl = new double[maxOutDegree];
    aliasSamplers.assign(100, alias1StOrderSampler);
    accRatioConst = calAccRatioConst(p, q);
}

Node2VecReject::~Node2VecReject() {
    delete [] nodeProbTbl;
}

#if (PLAIN | PLAIN_BUCKET)
#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
void Node2VecReject::initWalks(vid_t nVertices, WalkManager &walkManager, metrics &m) {
#if STATICCACHE
    m.start_time("initCache");
    /* 统计每个节点的度数 */
    logstream(LOG_INFO) << "initiating static buffer" << std::endl;
    std::vector<VertexInfo> vertexDegree;
    vertexDegree.resize(nVertices);
#if MULTI_THREAD
#pragma omp parallel for schedule(dynamic)
#endif
    for (vid_t v = 0; v < nVertices; v++){
        vertexDegree[v].vertexId = v;
        vertexIO.getVertexOutDegree(vertexDegree.at(v),AUTO);
    }
    std::sort(vertexDegree.begin(), vertexDegree.end(), sortByVertexDegreeDes);
    eid_t edgesPerBlock = blockSize_kb_global * 1024 / sizeof(vid_t);
    eid_t accDeg = 0;
    for (const auto &v: vertexDegree){
        accDeg += v.outDegree;
        staticBufferSize ++;
        vertexIO.vertexBuffer.needVertex(v.vertexId);
        if (accDeg > edgesPerBlock){
            break;
        }
    }
    vertexIO.vertexBuffer.setBufferSize(staticBufferSize);
    m.stop_time("initCache");
#endif

    std::cout << "Node2Vec start with reject sampler" << std::endl;

    tid_t nThreads = get_option_int("nThreads", omp_get_max_threads());
    omp_set_num_threads(nThreads);
    m.start_time("initWalks");
    walkManager.currentNWalks = (wid_t)nVertices * nWalkPerVertex;
#if MULTI_THREAD
#pragma omp parallel for schedule(dynamic)
#endif
    for (vid_t sourceVertex = 0; sourceVertex < nVertices; sourceVertex++){
        VertexInfo srcVertex(sourceVertex);
        bool noOutEdge = false;
        for (wid_t w = 0; w < nWalkPerVertex; w++){
            /* the absolute id of current vertex */
            vid_t currentVertex_abs = sampleStartEdge(m, srcVertex, noOutEdge);
            if (noOutEdge){
#pragma omp critical
                {
                    walkManager.currentNWalks --;
                }
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
#if STATICCACHE
        m.start_time("initCache");
        if (vertexIO.vertexBuffer.shouldCache(sourceVertex)) {
#pragma omp critical
            {
                vertexIO.vertexBuffer.streamV(srcVertex);
            }
        }
        m.stop_time("initCache");
#endif
    }
#if STATICCACHE
    m.start_time("initCache");
    vertexIO.vertexBuffer.sort();
    assert(vertexIO.vertexBuffer.getAddedVerticesNum() == staticBufferSize);
    vertexIO.vertexBuffer.bufferInitilized = true;
    m.stop_time("initCache");
#endif
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

vid_t
Node2VecReject::sampleDestVertex(metrics &m, vid_t preVertex, vid_t curVertex, eid_t *beg_pos,
                                 vid_t *csr, bid_t execBlock) {
    return soRejectSampler.sampleDestVertex(m, preVertex, curVertex, alias1StOrderSampler,
                                            &SecondOrderRW::rejectAccRatio,
                                            *this, execBlock);
}
#endif

vid_t Node2VecReject::sampleStartEdge(metrics &m, VertexInfo &srcVertex, bool &noOutEdge) {
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

double Node2VecReject::rejectAccRatio(VertexInfo &preV, VertexInfo &destV, metrics &m) {
//    tid_t t = omp_get_thread_num();
//    metric.at(t).start_time("compute-acc-ratio");
    double result = biasedWeight(preV, destV, m) * accRatioConst;
//    metric.at(t).stop_time("compute-acc-ratio");
    return result;
}

inline double Node2VecReject::calAccRatioConst(double p, double q) {
    double min = p;
    if (q < min)    min = q;
    if (1 < min)    min = 1;
    return min;
}


#pragma clang diagnostic pop


#endif //IOE_SORW_NODE2VECREJECT_HPP
