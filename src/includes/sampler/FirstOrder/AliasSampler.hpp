/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2020/11/5
 */

#ifndef IOE_SORW_ALIASSAMPLER_HPP
#define IOE_SORW_ALIASSAMPLER_HPP
#include "engine/Settings.hpp"
#include "IO/VertexIO.hpp"
#include "metrics/metrics.hpp"
#include "util/RandNum.hpp"

class Alias1stOrderSampler{
private:
    VertexIO *pVertexIO;

    vid_t **aliasTbl;
    float **probTbl;

    std::vector<RandNum> mainRand;
//    RandNum mainRand;

    bool unWeightedGraph;

public:
    Alias1stOrderSampler(VertexIO *_pVertexIO, metrics &m, bool _unWeightedGraph = false);

    vid_t sampleDestVertex(VertexInfo &v, metrics &m);

};

Alias1stOrderSampler::Alias1stOrderSampler(VertexIO *_pVertexIO, metrics &m, bool _unWeightedGraph){
    unWeightedGraph = _unWeightedGraph;
    pVertexIO = _pVertexIO;
    vid_t nTotalVertices = _pVertexIO->getNTotalVertices();
    mainRand.resize(100, RandNum(987865565678));
}

vid_t Alias1stOrderSampler::sampleDestVertex(VertexInfo &v, metrics &m) {
#if UNWEIGHTED_GRAPH
        if (!v.csr) {
//#if METRIC_GET_CSR
//            if (v.resideBlockId == dynamicBlock_global){
//                auto getCSRStart = std::chrono::system_clock::now();
//                pVertexIO->getVertexCSR(v, m);
//                auto getCSREnd = std::chrono::system_clock::now();
//                getCSR_dynamicBlock_global_ms += static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(getCSREnd - getCSRStart).count());
//            }else{
//                assert(v.resideBlockId == staticBlock_global);
//                pVertexIO->getVertexCSR(v, m);
//            }
//#else
#if TIME_COST_INFO
            tid_t t = omp_get_thread_num();
            pVertexIO->metric.at(t).start_time("get-csr");
#endif
            pVertexIO->getVertexCSR(v, m);
#if TIME_COST_INFO
            pVertexIO->metric.at(t).stop_time("get-csr");
#endif
//#endif
        }
        if (!v.outDegree) {
            return INVALID_VID;
//            abort();
        }
        tid_t t = omp_get_thread_num();
        return v.csr[int(floor(mainRand.at(t).dRand() * v.outDegree))];
#endif
}

#endif //IOE_SORW_ALIASSAMPLER_HPP
