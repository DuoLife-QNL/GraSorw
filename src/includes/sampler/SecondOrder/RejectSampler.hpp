/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2020/11/6
 */

#ifndef IOE_SORW_REJECTSAMPLER_HPP
#define IOE_SORW_REJECTSAMPLER_HPP
#include "metrics/metrics.hpp"
#include "engine/Settings.hpp"
#include "IO/VertexIO.hpp"
#include "sampler/FirstOrder/AliasSampler.hpp"
#include "logger/logger.hpp"
#include "util/RandNum.hpp"
#include "walks/SecondOrderRW.hpp"

class SORejectSampler{
public:
    SORejectSampler();
    vid_t sampleDestVertex(metrics &m, vid_t previousVertex, vid_t currentVertex, Alias1stOrderSampler &aliasSampler,
                           double (SecondOrderRW::*accRatio)(VertexInfo &, VertexInfo &, metrics &), SecondOrderRW &SORW,
                           const bid_t &execBlock);
    vid_t _sampleDestVertex(metrics &m, vid_t previousVertex, vid_t currentVertex, Alias1stOrderSampler &aliasSampler,
                           double (SecondOrderRW::*accRatio)(VertexInfo &, VertexInfo &, VertexInfo &, metrics &), SecondOrderRW &SORW,
                           const bid_t &execBlock){
        VertexInfo curVertexInfo(currentVertex, execBlock, true);
        VertexInfo destVertexInfo;
        VertexInfo preVertexInfo(previousVertex);
        preVertexInfo.delAfterUse = false;
//    tid_t t = omp_get_thread_num();
        do {
//        SORW.metric.at(t).start_time("1-orderAliasSampler");
            destVertexInfo.vertexId = aliasSampler.sampleDestVertex(curVertexInfo, m);
//        SORW.metric.at(t).stop_time("1-orderAliasSampler");
            if (!curVertexInfo.outDegree){
                logstream(LOG_ERROR) << "current vertex has no out edge" << std::endl;
                abort();
            }
        }while (!accept((SORW.*accRatio)(preVertexInfo, curVertexInfo, destVertexInfo, m)));

        return destVertexInfo.vertexId;
    }

    vid_t sampleDestVertex(metrics &m, VertexInfo &preVertexInfo, VertexInfo &curVertexInfo, Alias1stOrderSampler &aliasSampler,
                           double (SecondOrderRW::*accRatio)(VertexInfo &, VertexInfo &, metrics &), SecondOrderRW &SORW){
        VertexInfo destVertexInfo;
        do {
            destVertexInfo.vertexId = aliasSampler.sampleDestVertex(curVertexInfo, m);
        }while (!accept((SORW.*accRatio)(preVertexInfo, destVertexInfo, m)));

        return destVertexInfo.vertexId;
    }

    vid_t _sampleDestVertex(metrics &m, VertexInfo &preVertexInfo, VertexInfo &curVertexInfo, Alias1stOrderSampler &aliasSampler,
                           double (SecondOrderRW::*accRatio)(VertexInfo &, VertexInfo &, VertexInfo &, metrics &), SecondOrderRW &SORW){
        VertexInfo destVertexInfo;
        do {
            destVertexInfo.resetVertex(aliasSampler.sampleDestVertex(curVertexInfo, m));
        }while (!accept((SORW.*accRatio)(preVertexInfo, curVertexInfo, destVertexInfo, m)));

        return destVertexInfo.vertexId;
    }

private:
    std::vector<RandNum> randNum;

    bool accept(double accRatio);
};

SORejectSampler::SORejectSampler() {
    randNum.resize(100, RandNum(9898676785859));
}

vid_t SORejectSampler::sampleDestVertex(metrics &m, vid_t previousVertex, vid_t currentVertex,
                                        Alias1stOrderSampler &aliasSampler,
                                        double (SecondOrderRW::*accRatio)(VertexInfo &, VertexInfo &, metrics &), SecondOrderRW &SORW,
                                        const bid_t &execBlock) {
    VertexInfo curVertexInfo(currentVertex, execBlock, true);
    VertexInfo destVertexInfo;
    VertexInfo preVertexInfo(previousVertex);
//    tid_t t = omp_get_thread_num();
    do {
//        SORW.metric.at(t).start_time("1-orderAliasSampler");
         destVertexInfo.vertexId = aliasSampler.sampleDestVertex(curVertexInfo, m);
//        SORW.metric.at(t).stop_time("1-orderAliasSampler");
        if (!curVertexInfo.outDegree){
            logstream(LOG_ERROR) << "current vertex has no out edge" << std::endl;
            abort();
        }
    }while (!accept((SORW.*accRatio)(preVertexInfo, destVertexInfo, m)));

    return destVertexInfo.vertexId;
}

bool SORejectSampler::accept(double accRatio) {
    tid_t t = omp_get_thread_num();
    return randNum.at(t).dRand() < accRatio;
}

#endif //IOE_SORW_REJECTSAMPLER_HPP
