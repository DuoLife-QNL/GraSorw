/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/4/21 14:50
 */

#ifndef IOE_SORW_PPRONNV_HPP
#define IOE_SORW_PPRONNV_HPP

#include <vector>
#include <algorithm>
#include "BasicIncludes.hpp"
#include "util/RandNum.hpp"
#include "util/toplist.hpp"
#include "sampler/FirstOrder/AliasSampler.hpp"
#include "sampler/SecondOrder/RejectSampler.hpp"
#include "walks/Node2vecWithDecayFactor.hpp"

class PPRonNV : public Node2VecWithDecayFactor{
public:
    vid_t nStartVertices;
    wid_t nWalksPerSource;
    std::vector<vid_t> startVertices;

    PPRonNV(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int walkNumPerVertex, eid_t maxOutDegree,
            double p, double q, double decayFactor, int begPosFile, int csrFile,
            const std::string &startVerticesFileName, const vid_t &nStartVertices, vid_t nTotalVertices,
            metrics &m);

    ~PPRonNV();

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
        for (auto source: startVertices){
            if (source < startVertexId){
                continue;
            }
            if (source >= endVertexId){
                continue;
            }
            VertexInfo srcVertex(source);
            vertexIO.getVertexCSR(srcVertex, m);
#pragma omp parallel for schedule(dynamic)
            for (wid_t w = 0; w < nWalksPerSource; w++){
                bool noOutEdge = false;
                vid_t currentVertex_abs = sampleStartEdge(m, srcVertex, noOutEdge);
                if (noOutEdge){
                    continue;
                }
                bid_t preVertexResideBlock = block;
                bid_t resideBlock = getBlock(currentVertex_abs);
                vid_t currentVertex = currentVertex_abs;
                WalkDataType walk = WalkManager::encode(source, source, currentVertex, preVertexResideBlock, resideBlock, 1);
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
    double accRatioConst;
    inline double calAccRatioConst(double p, double q){
        double min = p;
        if (q < min)    min = q;
        if (1 < min)    min = 1;
        return min;
    }
};

PPRonNV::PPRonNV(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int walkNumPerVertex, eid_t maxOutDegree,
                 double p, double q, double decayFactor, int begPosFile, int csrFile,
                 const std::string &startVerticesFileName, const vid_t &nStartVertices, vid_t nTotalVertices,
                 metrics &m)
        : Node2VecWithDecayFactor(nBlocks, startVertex, walkLength, p, q, decayFactor, begPosFile, csrFile, maxOutDegree,
                                  nTotalVertices),
          alias1StOrderSampler(&vertexIO, m, true) {
    nWalksPerSource = walkNumPerVertex;
    this->nStartVertices = nStartVertices;

    std::ifstream starts(startVerticesFileName);
    vid_t source;
    for (vid_t v = 0; v < nStartVertices; v++){
        starts >> source;
        startVertices.push_back(source);
    }
    starts.close();
    nodeProbTbl = new double[maxOutDegree];
    accRatioConst = calAccRatioConst(p, q);
}

#if (PLAIN | PLAIN_BUCKET)
vid_t
PPRonNV::sampleDestVertex(metrics &m, vid_t preVertex, vid_t curVertex, eid_t *beg_pos,
                          vid_t *csr, bid_t execBlock) {
    return soRejectSampler._sampleDestVertex(m, preVertex, curVertex, alias1StOrderSampler,
                                            &SecondOrderRW::rejectAccRatio,
                                            *this, execBlock);
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
void PPRonNV::initWalks(vid_t nVertices, WalkManager &walkManager, metrics &m) {
#if STATICCACHE
    m.start_time("init-cache");
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
    m.stop_time("init-cache");
#endif

    std::cout << "Initiating walks" << std::endl;

    tid_t nThreads = get_option_int("nThreads", omp_get_max_threads());
    omp_set_num_threads(nThreads);
    m.start_time("initWalks");
    walkManager.currentNWalks = nStartVertices * nWalksPerSource;
    for (auto sourceVertex: startVertices){
        VertexInfo srcVertex(sourceVertex);
        vertexIO.getVertexCSR(srcVertex, m);
        bool noOutEdge = false;
#if MULTI_THREAD
#pragma omp parallel for schedule(dynamic)
#endif
        for (wid_t w = 0; w < nWalksPerSource; w++){
            tid_t t = omp_get_thread_num();
            /* the absolute id of current vertex */
            vid_t currentVertex_abs = sampleStartEdge(m, srcVertex, noOutEdge);
            if (noOutEdge){
#pragma omp critical
                {
                    walkManager.currentNWalks--;
                }
                continue;
            }
#if DEBUG
            totalSteps_t.at(t)++;
#endif
            bid_t preVertexResideBlock = getBlock(sourceVertex);
            bid_t resideBlock = getBlock(currentVertex_abs);
            vid_t currentVertex = currentVertex_abs;
            WalkDataType walk = WalkManager::encode(sourceVertex, sourceVertex, currentVertex, preVertexResideBlock, resideBlock, 1);

            walkManager.moveWalk(walk, resideBlock, omp_get_thread_num(), sourceVertex, currentVertex, 0, preVertexResideBlock);
        }
    }

#if STATICCACHE
    m.start_time("init-cache");
    std::cout << "Filling Static Buffer" << std::endl;
#if MULTI_THREAD
#pragma omp parallel for schedule(dynamic)
#endif
    for (vid_t sourceVertex = 0; sourceVertex < nVertices; sourceVertex++){
        VertexInfo srcVertex(sourceVertex);
        if (vertexIO.vertexBuffer.shouldCache(sourceVertex)){
#pragma omp critical
            {
                vertexIO.getVertexCSR(srcVertex, m, DISK);
                vertexIO.vertexBuffer.streamV(srcVertex);
            }
        }
    }
    vertexIO.vertexBuffer.sort();
    assert(vertexIO.vertexBuffer.getAddedVerticesNum() == staticBufferSize);
    vertexIO.vertexBuffer.bufferInitilized = true;
    m.stop_time("init-cache");
#endif
    m.stop_time("initWalks");


    std::cout << "currentNWalks before run: " << walkManager.currentNWalks << std::endl;

    for (bid_t block = 0; block < nBlocks; block++){
        walkManager.nWalks[block] = walkManager.nDiskWalks[block];

        for (tid_t t = 0; t < nThreads; t++){
            walkManager.nWalks[block] += walkManager.walkPool[t][block].size_w;
        }
        if (walkManager.nWalks[block]){
            walkManager.minStep[block] = 1;
        }
    }
}
#endif

PPRonNV::~PPRonNV() {
    delete [] nodeProbTbl;
}

vid_t PPRonNV::sampleStartEdge(metrics &m, VertexInfo &srcVertex, bool &noOutEdge) {
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

double PPRonNV::rejectAccRatio(VertexInfo &preV, VertexInfo &curV, VertexInfo &destV, metrics &m) {
    double result = biasedWeight(preV, destV, m) * accRatioConst;
    return result;
}


#pragma clang diagnostic pop
#endif //IOE_SORW_PPRONNV_HPP
