//
// Created by lihz on 2020/10/28.
//
#ifndef IOE_SORW_NODE2VECONLINE_HPP
#define IOE_SORW_NODE2VECONLINE_HPP

#include <vector>

#include "BasicIncludes.hpp"
#include "walks/Node2Vec.hpp"
#include "util/RandNum.hpp"
#include "util/toplist.hpp"

class Node2vecOnline : public Node2Vec{
public:
    Node2vecOnline(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int walkNumPerVertex,
                   eid_t maxOutDegree, double p, double q, int begPosFile, int csrFile,
                   vid_t nTotalVertices);

    ~Node2vecOnline();

    vid_t sampleStartEdge(metrics &m, VertexInfo &srcVertex, bool &noOutEdge) override;

    vid_t sampleDestVertex(metrics &m, vid_t preVertex, vid_t curVertex, eid_t *beg_pos,
                           vid_t *csr, bid_t execBlock) override;

    void initWalks(vid_t nVertices, WalkManager &soWalkManager, metrics &m) override;


private:
    double *nodeProbTbl;
    RandNum mainRand;
};

vid_t
Node2vecOnline::sampleDestVertex(metrics &m, vid_t preVertex, vid_t curVertex, eid_t *beg_pos,
                                 vid_t *csr, bid_t execBlock) {

    eid_t startEdge = beg_pos[curVertex - startVertex[execBlock]];
    eid_t startDestVertexOffset = startEdge - beg_pos[0];
//    nodeProbTbl[0] = biasedWeight(m, preVertex, csr[startDestVertexOffset]);
    eid_t outDegree =  getVertexOutDegree(curVertex - startVertex[execBlock], beg_pos);
//    m.start_time("compute-biased-weight");
    for (int i = 1; i < outDegree; i++){
//        nodeProbTbl[i] = nodeProbTbl[i - 1] + biasedWeight(m, preVertex, csr[startDestVertexOffset + i]);
    }
//    m.stop_time("compute-biased-weight");
    double probRand = mainRand.dRand() * nodeProbTbl[outDegree - 1];

    for (int i = 0; i < outDegree - 1; i++){
        if (probRand < nodeProbTbl[i]){
            return csr[startDestVertexOffset + i];
        }
    }
    return csr[startDestVertexOffset + outDegree - 1];
}

Node2vecOnline::Node2vecOnline(bid_t nBlocks, vid_t *startVertex, hid_t walkLength, int walkNumPerVertex,
                               eid_t maxOutDegree, double p, double q, int begPosFile, int csrFile,
                               vid_t nTotalVertices)
        : Node2Vec(nBlocks, startVertex, walkLength, walkNumPerVertex, p, q, begPosFile, csrFile, maxOutDegree, nTotalVertices), mainRand(9798798452435) {
    nodeProbTbl = new double[maxOutDegree];
}

Node2vecOnline::~Node2vecOnline() {
    delete [] nodeProbTbl;
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
void Node2vecOnline::initWalks(vid_t nVertices, WalkManager &soWalkManager, metrics &m) {
    std::cout << "Node2Vec start with online sampler" << std::endl;

    tid_t nThreads = get_option_int("nThreads", omp_get_max_threads());
    soWalkManager.currentNWalks = nVertices * nWalkPerVertex;
    omp_set_num_threads(nThreads);
#pragma omp parallel for schedule(static)
    for (vid_t sourceVertex = 0; sourceVertex < nVertices; sourceVertex++){
        for (wid_t w = 0; w < nWalkPerVertex; w++){
            bool noOutEdge = false;
            /* the absolute id of current vertex */
            vid_t currentVertex_abs = sampleStartEdge(m, sourceVertex, noOutEdge);
            if (noOutEdge){
                soWalkManager.currentNWalks -= 1;
                continue;
            }
            bid_t resideBlock = getBlock(currentVertex_abs);
            vid_t currentVertex = currentVertex_abs - startVertex[resideBlock];
            WalkDataType walk = soWalkManager.encode(sourceVertex, sourceVertex, currentVertex, 0, 0);
            soWalkManager.moveWalk(walk, resideBlock, omp_get_thread_num(), sourceVertex, currentVertex, 0, 0);
        }
    }

    for (bid_t block = 0; block < nBlocks; block++){
        soWalkManager.nWalks[block] = soWalkManager.nDiskWalks[block];

        for (tid_t t = 0; t < nThreads; t++){
            soWalkManager.nWalks[block] += soWalkManager.walkPool[t][block].size_w;
        }
        if (soWalkManager.nWalks[block]){
            soWalkManager.minStep[block] = 0;
        }
    }
}

vid_t Node2vecOnline::sampleStartEdge(metrics &m, VertexInfo &srcVertex, bool &noOutEdge) {
    vertexIO.getVertexCSR(srcVertex, m);
    if (srcVertex.outDegree == 0){
        noOutEdge = true;
        return 0;
    }
    vid_t secondVertex = srcVertex.csr[mainRand.iRand(outDegree)];
    noOutEdge = false;
    return secondVertex;
}
#pragma clang diagnostic pop

#endif //IOE_SORW_NODE2VECONLINE_HPP