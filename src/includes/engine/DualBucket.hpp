/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/3/26 10:48
 */

#ifndef IOE_SORW_DUALBUCKET_HPP
#define IOE_SORW_DUALBUCKET_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <fcntl.h>
#include <unistd.h>
#include <cassert>
#include <omp.h>
#include <vector>
#include <map>
#include <ctime>
#include <sys/mman.h>
#include <asm/mman.h>
#include <numeric>

#include "api/filename.hpp"
#include "api/io.hpp"
#include "metrics/metrics.hpp"
#include "api/pthread_tools.hpp"
#include "walks/SecondOrderRW.hpp"
#include "Settings.hpp"
#include "logger/logger.hpp"
#include "walks/RandomWalk.hpp"
#include "util/Timer.hpp"

#if BI_BLOCK
#if DEBUG
class MemoryComponent{
public:
    double graphfile = 0;
    double currentWalkArray = 0;
    double bucket = 0;
    double walkPool = 0;
};
#endif

class DualBucketEngine {
public:
    std::string base_filename;
    // unsigned membudget_mb;
    unsigned long long blocksize_kb;
    bid_t nblocks;
    vid_t nvertices;
    tid_t exec_threads;
    vid_t *blocks;
    timeval start;

    /* Ｉn memory blocks */
    bid_t nmblocks; //number of in memory blocks
    vid_t **csrbuf;
    eid_t **beg_posbuf;
    bid_t cmblocks; //current number of in memory blocks
    bid_t *inMemIndex;
    int beg_posf, csrf;

    bid_t staticBlockIndex = 0;
    bid_t dynamicBlockIndex = 1;

    /* State */
    bid_t exec_block;

    /* Metrics */
    metrics &m;
    double memUsage;
    double maxMemUsage = 0;
    long long memCount = 0;
#if DEBUG
    std::vector<double> maxMemUsage_staticBlock;
    MemoryComponent memoryComponent_b0;
    std::vector<MemoryComponent> blockMemComp;
#endif

    WalkManager *walk_manager;

    void print_config();

    double runtime();

public:

    /**
     * Initialize GraphChi engine
     * @param base_filename prefix of the graph files
     * @param nblocks number of shards
     * @param selective_scheduling if true, uses selective scheduling
     */
    DualBucketEngine(const std::string &_base_filename, unsigned long long int _blocksize_kb,
                     bid_t _nblocks, bid_t _nmblocks, metrics &_m);

    virtual ~DualBucketEngine();

    void load_block_range(std::string base_filename, unsigned long long blocksize_kb, vid_t *&blocks,
                          bool allowfail = false);

    void loadSubGraph(bid_t p, eid_t *&beg_pos, vid_t *&csr, vid_t *nverts, eid_t *nedges);

    void loadSubGraph(bid_t p, vid_t *nverts, eid_t *nedges, bid_t doNotUnloadBlockId) {
        if (inMemIndex[p] == nmblocks) {//the block is not in memory
            bid_t swapin;
            if (cmblocks < nmblocks) {
                swapin = cmblocks++;
            } else {
                bid_t minmwb = chooseSwapOutExcept(doNotUnloadBlockId);
                swapin = unLoadBlock_id(minmwb);
                // munmap(beg_posbuf[swapin], sizeof(eid_t)*(blocks[minmwb+1] - blocks[minmwb] + 1));
            }
            loadSubGraph(p, beg_posbuf[swapin], csrbuf[swapin], nverts, nedges);
            inMemIndex[p] = swapin;
        }
    }

    void loadBlock(bid_t p, eid_t *&beg_pos, vid_t *&csr) {

        // m.start_time("__g_loadSubGraph_malloc_begpos");
        /* read beg_pos file */

        vid_t nverts = blocks[p + 1] - blocks[p];
        beg_pos = (eid_t *) malloc((nverts + 1) * sizeof(eid_t));
        // m.stop_time("__g_loadSubGraph_malloc_begpos");
        // beg_pos=(eid_t *)mmap(NULL,(size_t)(*nverts+1)*sizeof(eid_t),
        //         PROT_READ | PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, 0, 0);
        // if(beg_pos == MAP_FAILED)
        // {
        //     printf("%ld\n",(size_t)(*nverts+1)*sizeof(eid_t));
        //     perror("beg_pos alloc mmap");
        //     exit(-1);
        // }
//        m.start_time("z__g_loadSubGraph_read_begpos");
        preada(beg_posf, beg_pos, (size_t) (nverts + 1) * sizeof(eid_t), (size_t) blocks[p] * sizeof(eid_t));
//        m.stop_time("z__g_loadSubGraph_read_begpos");
        /* read csr file */
        m.start_time("z__g_loadSubGraph_realloc_csr");
        eid_t nedges = beg_pos[nverts] - beg_pos[0];
        if (nedges * sizeof(vid_t) > blocksize_kb * 1024) {
            csr = (vid_t *) realloc(csr, (nedges) * sizeof(vid_t));
        }
        m.stop_time("z__g_loadSubGraph_realloc_csr");
        m.start_time("z__g_loadSubGraph_read_csr");
        preada(csrf, csr, (nedges) * sizeof(vid_t), beg_pos[0] * sizeof(vid_t));
        m.stop_time("z__g_loadSubGraph_read_csr");

        m.stop_time("g_loadSubGraph");
    }

    void swapBlock(bid_t swapInBlockId, bid_t swapOutBlockId = INVALID_BID) {
        m.start_time("5_LoadBlocks");
        bid_t swapOutBlockIndex;
        if (swapInBlockId == swapOutBlockId){
            return;
        }
        if (cmblocks < nmblocks) {
            swapOutBlockIndex = cmblocks++;
        } else {
            assert(swapOutBlockId != INVALID_BID);
            swapOutBlockIndex = inMemIndex[swapOutBlockId];
            inMemIndex[swapOutBlockId] = nmblocks;
            unLoadBlock_index(swapOutBlockIndex);
        }
        loadBlock(swapInBlockId, beg_posbuf[swapOutBlockIndex], csrbuf[swapOutBlockIndex]);
        inMemIndex[swapInBlockId] = swapOutBlockIndex;
        m.stop_time("5_LoadBlocks");
    }

    void unLoadAllBlock() {
        for (bid_t b = 0; b < nblocks; b++) {
            inMemIndex[b] = nmblocks;
        }
        for (bid_t b = 0; b < nmblocks; b++) {
            unLoadBlock_index(b);
        }
        cmblocks = 0;
    }

    bid_t unLoadBlock_id(bid_t p) {
        bid_t unLoadIndex = inMemIndex[p];
        inMemIndex[p] = nmblocks;
        assert(unLoadIndex < nmblocks);
        if (beg_posbuf[unLoadIndex] != nullptr) free(beg_posbuf[unLoadIndex]);
        return unLoadIndex;
    }

    void unLoadBlock_index(bid_t blockIndex) {
        assert(blockIndex < nmblocks);
        if (beg_posbuf[blockIndex] != nullptr) {
            free(beg_posbuf[blockIndex]);
        }
    }

    void loadPreVertexSubGraph(bid_t p, eid_t *&beg_pos, vid_t *&csr, vid_t *nverts, eid_t *nedges);

    void findSubGraph(bid_t p, eid_t *&beg_pos, vid_t *&csr, vid_t *nverts, eid_t *nedges);

    bid_t swapOut();

    bid_t chooseSwapOutExcept(bid_t doNotUnloadBlockId) {
        m.start_time("z_g_swapOut");
        wid_t minmw = 0xffffffff;
        bid_t minmwb = 0;
        for (bid_t b = 0; b < nblocks; b++) {
            if (b == doNotUnloadBlockId) {
                continue;
            }
            if (inMemIndex[b] < nmblocks && walk_manager->nWalks[b] < minmw) {
                minmw = walk_manager->nWalks[b];
                minmwb = b;
            }
        }
        m.start_time("z_g_swapOut");
        return minmwb;
    }

    vid_t num_vertices() {
        return blocks[nblocks];
    }

    void run(vid_t nRootVertices, SecondOrderRW &userprogram, float prob, float th);

#if FULLY_LOAD

    void run_fullyLoadTest(vid_t nRootVertices, SecondOrderRW &userprogram, float prob, float th);

#elif NO_LOAD_TEST
    void run_noLoadTest(vid_t nRootVertices, SecondOrderRW &userprogram, float prob, float th);
#elif ONDEMAND_LOAD
    void run_onDemandLoad(vid_t nRootVertices, SecondOrderRW &userprogram, float prob, float th);
#endif

#define PRE_TRAIN_NO_LOAD 0
#define PRE_TRAIN_FULL_LOAD 1

    void preTrain(SecondOrderRW &userprogram){
        std::vector<eid_t> nEdges;
        nEdges.assign(nblocks, 0);
        vid_t tmp;
        for (bid_t b = 0; b < nblocks; b++){
            userprogram.vertexIO.getBlockSize(b, tmp, nEdges.at(b));
        }
        auto baseRandNum = RandNum(754388217423);
        std::vector<RandNum> randNums;
        /* randnum for each thread */
        for (int i = 0; i < 100; i++){
            auto randNum = RandNum(3245523322 + baseRandNum.iRand(664384234));
            randNums.push_back(randNum);
        }

        /* 设置输出文件 */
        std::ofstream preTrainLog_fullLoad("/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/pre-train-full-load.csv");
        preTrainLog_fullLoad << "blockId,ratio,ncsr,time" << std::endl;
        std::ofstream preTrainLog_noLoad("/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/pre-train-no-load.csv");
        preTrainLog_noLoad << "blockId,ratio,ncsr,time" << std::endl;

        /* 设定mini-test的节点数量及数量变化幅度 */
        const auto startRatio = 0.01;
        const auto ratioStep = 0.01;
        const auto stopRatio = 1;
        Timer loadBlockTimer, csrTimer;

        /* full load test */
        const int forwardLoadTime = 1;
        const int reverseLoadTime = 2;
        const int totalLoadTime = forwardLoadTime + reverseLoadTime + 1;
        typedef struct loadCSRLog{
            double ratio;
            vid_t nverts;
            double time;
        }loadCSRLog;
        std::vector<std::vector<double>> loadBlockTimeLog;
        std::vector<std::vector<loadCSRLog>> csrTime_fullLoad;
        loadBlockTimeLog.resize(nblocks);
        csrTime_fullLoad.resize(nblocks);
        /* 正序加载测试 */
        for (int i = 0; i < forwardLoadTime; i++){
            bid_t preBlockId = INVALID_BID;
            for (bid_t b = 0; b < nblocks; b++){
                loadBlockTimer.start();
                swapBlock(b, preBlockId);
                loadBlockTimer.stop();
                loadBlockTimeLog.at(b).push_back(loadBlockTimer.duration_s());
                preBlockId = b;
            }
            unLoadAllBlock();
        }
        /* 逆序加载测试 */
        for (int i = 0; i < reverseLoadTime; i++){
            bid_t preBlockId = INVALID_BID;
            for (bid_t b = nblocks - 1; b >= 0 && b < nblocks; b--){
                loadBlockTimer.start();
                swapBlock(b, preBlockId);
                loadBlockTimer.stop();
                loadBlockTimeLog.at(b).push_back(loadBlockTimer.duration_s());
                preBlockId = b;
            }
            unLoadAllBlock();
        }
        /* 计算内存中获取节点信息总时间，本次采用顺序加载，加载时间同样作为取平均的依据 */
        bid_t preBlockId = INVALID_BID;
        for (bid_t b = 0; b < nblocks; b++){
            loadBlockTimer.start();
            swapBlock(b, preBlockId);
            userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
            loadBlockTimer.stop();
            loadBlockTimeLog.at(b).push_back(loadBlockTimer.duration_s());
            auto loadBlockTime = std::accumulate(loadBlockTimeLog.at(b).begin(), loadBlockTimeLog.at(b).end(), 0.0) / loadBlockTimeLog.at(b).size();
            preBlockId = b;
            auto ratio = startRatio;
            vid_t nverts = blocks[b + 1] - blocks[b];
            while (ratio < stopRatio){
                auto nSampleV = static_cast<vid_t>(ratio * nverts);
                csrTimer.start();
#pragma omp parallel for schedule(dynamic)
                for (vid_t v = 0; v < nSampleV; v++){
                    tid_t t = omp_get_thread_num();
                    vid_t nodeId = blocks[0] + randNums[t].iRand(blocks[1] - blocks[0]);
                    VertexInfo node(nodeId, 0, true);
                    userprogram.vertexIO.getVertexCSR(node, m);
                }
                csrTimer.stop();
                loadCSRLog temp = {.ratio = ratio, .nverts = nSampleV, .time = csrTimer.duration_s() + loadBlockTime};
                csrTime_fullLoad.at(b).push_back(temp);
                ratio += ratioStep;
            }
        }
        unLoadAllBlock();
        /* 将Log输出到文件 */
        for (bid_t b = 0; b < nblocks; b++){
            for (const auto & log: csrTime_fullLoad.at(b)){
                preTrainLog_fullLoad << b << ',' << log.ratio << ',' << log.nverts << ',' << log.time << std::endl;
            }
        }
        preTrainLog_fullLoad.close();

        /* no-load测试 */
        std::vector<std::vector<loadCSRLog>> csrTime_noLoad;
        csrTime_noLoad.resize(nblocks);
        auto ratio = startRatio;
        while (ratio < stopRatio){
            for (bid_t b = 0; b < nblocks; b++){
                auto nSampleV = static_cast<vid_t>(ratio * (blocks[b + 1] - blocks[b]));
                csrTimer.start();
#pragma omp parallel for schedule(dynamic)
                for (vid_t v = 0; v < nSampleV; v++){
                    tid_t t = omp_get_thread_num();
                    vid_t nodeId = blocks[0] + randNums[t].iRand(blocks[1] - blocks[0]);
                    VertexInfo node(nodeId, 0, false);
                    userprogram.vertexIO.getVertexCSR(node, m);
                }
                csrTimer.stop();
                loadCSRLog temp = {.ratio = ratio, .nverts = nSampleV, .time = csrTimer.duration_s()};
                csrTime_noLoad.at(b).push_back(temp);
            }
            ratio += ratioStep;
        }
        for (bid_t b = 0; b < nblocks; b++){
            for (const auto & log: csrTime_noLoad.at(b)){
                preTrainLog_noLoad << b << ',' << log.ratio << ',' << log.nverts << ',' << log.time << std::endl;
            }
        }
        preTrainLog_noLoad.close();

//#if PRE_TRAIN_FULL_LOAD
//        swapBlock(0);
//        userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
//#endif
//        auto nverts = blocks[1] - blocks[0];
//        auto increment = 0.001;
//        auto ratio = increment;
//#if PRE_TRAIN_NO_LOAD
//
//#elif PRE_TRAIN_FULL_LOAD
//
//#endif
//        Timer noLoadTimer;
//        while (ratio < 10){
//            vid_t nSampledV = ratio * nverts;
//            noLoadTimer.start();
//#pragma omp parallel for schedule(dynamic)
//            for (vid_t v = 0; v < nSampledV; v++){
//                tid_t t = omp_get_thread_num();
//                vid_t nodeId = blocks[0] + randNums[t].iRand(blocks[1] - blocks[0]);
//                VertexInfo node(nodeId, 0, false);
//                userprogram.vertexIO.getVertexCSR(node, m);
//            }
//            noLoadTimer.stop();
//#if PRE_TRAIN_FULL_LOAD
//            preTrainLog_fullLoad << nSampledV << ',' << noLoadTimer.duration_s() + loadBlockTimer.duration_s() << std::endl;
//#elif PRE_TRAIN_NO_LOAD
//            preTrainLog_noLoad << nSampledV << ',' << noLoadTimer.duration_s() << std::endl;
//#endif
//            ratio += increment;
//        }
//        ratio = 0.6;
//        vid_t nSampledV = ratio * nverts;
//        noLoadTimer.start();
//        #pragma omp parallel for schedule(dynamic)
//        for (vid_t v = 0; v < nSampledV; v++){
//            tid_t t = omp_get_thread_num();
//            vid_t nodeId = blocks[0] + randNums[t].iRand(blocks[1] - blocks[0]);
//            VertexInfo node(nodeId, 0, false);
//            userprogram.vertexIO.getVertexCSR(node, m);
//        }
//        noLoadTimer.stop();
//        //            preTrainLog_fullLoad << nSampledV << ',' << noLoadTimer.duration_s() + loadBlockTimer.duration_s() << std::endl;
//        preTrainLog_noLoad << nSampledV << ',' << noLoadTimer.duration_s() << std::endl;
    }

    void preTrain2(SecondOrderRW &userprogram){
//        /* 采样多少个节点来描述一条射线 */
//        const int N = 10;
        /* 统计子图信息 */
        std::vector<eid_t> nEdges;
        nEdges.assign(nblocks, 0);
        vid_t tmp;
        for (bid_t b = 0; b < nblocks; b++){
            userprogram.vertexIO.getBlockSize(b, tmp, nEdges.at(b));
        }
        auto baseRandNum = RandNum(754388217423);

        const auto startRatio = 0.01;
        const auto ratioStep = 0.01;
        const auto stopRatio = 1;
        Timer csrTimer;

        const bid_t dynamicBlockId = 16;
        bid_t preBlockId = INVALID_BID;
        swapBlock(dynamicBlockId, preBlockId);
        userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);

        /* 将所有节点进行分类，domains[b]中的节点均存在属于子图b的邻居 */
        std::vector<std::vector<vid_t>> domains;
        domains.resize(nblocks);
        for (vid_t v = blocks[dynamicBlockId]; v < blocks[dynamicBlockId +1]; v++){
            VertexInfo vertexInfo(v, dynamicBlockId, true);
            userprogram.vertexIO.getVertexCSR(vertexInfo, m);

            /* 将v放入应有的domain中*/
            vid_t neighborIndex = 0;
            bid_t staticBlockId = 0;
            while (neighborIndex < vertexInfo.outDegree && staticBlockId < dynamicBlockId){
                vid_t neighbor = vertexInfo.csr[neighborIndex];
                if (userprogram.vertexIO.vertexInBlock(neighbor, staticBlockId)){
                    domains.at(staticBlockId).push_back(v);
                    neighborIndex++;
                    while (userprogram.vertexIO.vertexInBlock(vertexInfo.csr[neighborIndex], staticBlockId)){
                        neighborIndex++;
                    }
                    staticBlockId++;
                }else{
                    staticBlockId++;
                }
            }
        }

        /* 在不同的domain中采样 */

    }

    /* use array instead of vector to store bucket */
    void run_array(vid_t nRootVertices, SecondOrderRW &userprogram, float prob, float th) {
        // srand((unsigned)time(NULL));
        m.start_time("2_StartWalks");
        bid_t preInitBlock = 0;
        for (bid_t b = 0; b < nblocks; b++) {
            swapBlock(b, preInitBlock);
            userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
            userprogram.initWalks(b, *walk_manager, m);
#if INFO
            logstream(LOG_INFO) << "initiated block" << b << std::endl;
#endif
            preInitBlock = b;
        }
        walk_manager->updateWalkNum();
        std::cout << "currentNWalks before run: " << walk_manager->currentNWalks << std::endl;
        m.stop_time("2_StartWalks");

        gettimeofday(&start, NULL);
        m.start_time("1_ProcessTime");

        vid_t nverts;
        eid_t nedges;
        int blockcount = 0;
        bid_t staticBlock, dynamicBlock;
        bid_t preStaticBlock = 0, preDynamicBlock = 1;
        bid_t lostBlockId = INVALID_BID;
        unLoadAllBlock();
        while (!userprogram.allWalksFinished(*walk_manager)) {
            for (staticBlock = 0; staticBlock < nblocks - 1; staticBlock++) {
                blockcount++;
                if (!walk_manager->nWalks[staticBlock]) {
                    continue;
                }
                if (inMemIndex[staticBlock] < nmblocks) {
                    lostBlockId = preStaticBlock;
                } else {
                    swapBlock(staticBlock, preStaticBlock);
                }
                wid_t nWalks = walk_manager->getCurrentWalks(staticBlock);
                walk_manager->clearRecoredWalkNum(staticBlock);
//             if(blockcount % (nBlocks/100+1)==1)
//             if(blockcount % (1024*1024*1024/nedges+1) == 1)
                {
                    logstream(LOG_DEBUG) << runtime() << "s : staticBlockCount: " << blockcount << std::endl;
                    logstream(LOG_INFO) << "currentNWalks = " << walk_manager->currentNWalks << ", nwalks["
                                        << staticBlock << "] = " << nWalks << std::endl;
                }
                walk_manager->collectBucket(staticBlock);
                for (dynamicBlock = staticBlock + 1; dynamicBlock < nblocks; dynamicBlock++) {
                    if (!walk_manager->preBucket[dynamicBlock].size() &&
                        !walk_manager->curBucket[dynamicBlock].size()) {
                        continue;
                    }
                    if (preDynamicBlock == staticBlock) {
                        assert(lostBlockId != INVALID_BID);
                        swapBlock(dynamicBlock, lostBlockId);
                        lostBlockId = INVALID_BID;
                    } else {
                        swapBlock(dynamicBlock, preDynamicBlock);
                    }
                    userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
                    m.start_time("4_UpdateWalk");
#pragma omp parallel for schedule(dynamic)
                    for (wid_t w = 0; w < walk_manager->preBucket[dynamicBlock].size(); w++) {
                        userprogram.processWalk(m, walk_manager->preBucket[dynamicBlock][w], staticBlock, dynamicBlock,
                                                *walk_manager);
                    }
#pragma omp parallel for schedule(dynamic)
                    for (wid_t w = 0; w < walk_manager->curBucket[dynamicBlock].size(); w++) {
                        userprogram.processWalk(m, walk_manager->curBucket[dynamicBlock][w], staticBlock, dynamicBlock,
                                                *walk_manager);
                    }
                    preDynamicBlock = dynamicBlock;
                    m.stop_time("4_UpdateWalk");
                }
#if INFO | DEBUG
                double usage = get_memory();
                if (usage > maxMemUsage) {
                    maxMemUsage = usage;
                }
#endif
#if DEBUG
                if (usage > maxMemUsage_staticBlock.at(staticBlock)){
                maxMemUsage_staticBlock.at(staticBlock) = usage;
//                if (staticBlock == 0){
//                    /* two block currently in memory is staticBlock and preDynamicBlock */
//                    vid_t nVertsStaticBlock = blocks[staticBlock + 1] - blocks[staticBlock];
//                    vid_t nVertsPreDynamicBlock = blocks[preDynamicBlock + 1] - blocks[preDynamicBlock];
//                    bid_t staticBlockInMemIndex = 0;
//                    bid_t preDynamicBlockInMemIndex = 1;
//                    /* FIXME: 这里实际上不是很准确，如果涉及到lostBlock，这里将统计错误，但考虑到子图大小差别不大，这里误差应该较小 */
//                    if (inMemIndex[preDynamicBlock] == 0){
//                        preDynamicBlockInMemIndex = 0;
//                        staticBlockInMemIndex = 1;
//                    }
//                    eid_t nEdgesStaticBlock = beg_posbuf[staticBlockInMemIndex][nVertsStaticBlock] - beg_posbuf[staticBlockInMemIndex][0];
//                    eid_t nEdgesPreDynamicBlock = beg_posbuf[preDynamicBlockInMemIndex][nVertsPreDynamicBlock] - beg_posbuf[preDynamicBlockInMemIndex][0];
//                    size_t csrByte = (nEdgesPreDynamicBlock + nEdgesStaticBlock) * sizeof(vid_t);
//                    size_t begPosByte = (nVertsStaticBlock + nVertsPreDynamicBlock) * sizeof(vid_t);
//                    memoryComponent_b0.graphfile = static_cast<double>(csrByte + begPosByte) / 1024 / 1024;
//                    size_t currentWalkByte = walk_manager->curBlockNWalks * sizeof(WalkDataType);
//                    memoryComponent_b0.currentWalkArray = static_cast<double>(currentWalkByte) / 1024 / 1024;
//                    size_t bucketsByte = 0;
//                    for (bid_t b = 0; b < nblocks; b++){
//                        /* vec.size() and vec.capacity differ a lot */
//                        bucketsByte += walk_manager->preBucket[b].capacity() * sizeof(WalkDataType);
//                        bucketsByte += walk_manager->curBucket[b].capacity() * sizeof(WalkDataType);
//                    }
//                    memoryComponent_b0.bucket = static_cast<double>(bucketsByte) / 1024 / 1024;
//                    size_t walkPoolSize_b = 0;
//                    for (tid_t t = 0; t < 16; t++){
//                        for (bid_t b = 0; b < nblocks; b++){
//                            walkPoolSize_b += walk_manager->walkPool[t][b].size_w * sizeof(WalkDataType);
//                        }
//                    }
//                    memoryComponent_b0.walkPool = static_cast<double>(walkPoolSize_b) / 1024 / 1024;
//                }
                /* two block currently in memory is staticBlock and preDynamicBlock */
                vid_t nVertsStaticBlock = blocks[staticBlock + 1] - blocks[staticBlock];
                vid_t nVertsPreDynamicBlock = blocks[preDynamicBlock + 1] - blocks[preDynamicBlock];
                bid_t staticBlockInMemIndex = 0;
                bid_t preDynamicBlockInMemIndex = 1;
                /* FIXME: 这里实际上不是很准确，如果涉及到lostBlock，这里将统计错误，但考虑到子图大小差别不大，这里误差应该较小 */
                if (inMemIndex[preDynamicBlock] == 0){
                    preDynamicBlockInMemIndex = 0;
                    staticBlockInMemIndex = 1;
                }
                eid_t nEdgesStaticBlock = beg_posbuf[staticBlockInMemIndex][nVertsStaticBlock] - beg_posbuf[staticBlockInMemIndex][0];
                eid_t nEdgesPreDynamicBlock = beg_posbuf[preDynamicBlockInMemIndex][nVertsPreDynamicBlock] - beg_posbuf[preDynamicBlockInMemIndex][0];
                size_t csrByte = (nEdgesPreDynamicBlock + nEdgesStaticBlock) * sizeof(vid_t);
                size_t begPosByte = (nVertsStaticBlock + nVertsPreDynamicBlock) * sizeof(vid_t);
                blockMemComp.at(staticBlock).graphfile = static_cast<double>(csrByte + begPosByte) / 1024 / 1024;
                size_t currentWalkByte = walk_manager->curBlockNWalks * sizeof(WalkDataType);
                blockMemComp.at(staticBlock).currentWalkArray = static_cast<double>(currentWalkByte) / 1024 / 1024;
                size_t bucketsByte = 0;
                for (bid_t b = 0; b < nblocks; b++){
                    /* vec.size() and vec.capacity differ a lot */
                    bucketsByte += walk_manager->preBucket[b].size() * sizeof(WalkDataType);
                    bucketsByte += walk_manager->curBucket[b].size() * sizeof(WalkDataType);
                }
                blockMemComp.at(staticBlock).bucket = static_cast<double>(bucketsByte) / 1024 / 1024;
                size_t walkPoolSize_b = 0;
                for (tid_t t = 0; t < 16; t++){
                    for (bid_t b = 0; b < nblocks; b++){
                        walkPoolSize_b += walk_manager->walkPool[t][b].size_w * sizeof(WalkDataType);
                    }
                }
                blockMemComp.at(staticBlock).walkPool = static_cast<double>(walkPoolSize_b) / 1024 / 1024;
            }
#endif
#if INFO | DEBUG
                memUsage = ((memUsage * memCount) + usage) / (memCount + 1);
                memCount++;
#endif
                walk_manager->updateWalkNum();
                walk_manager->clearBucket();
                preStaticBlock = staticBlock;
            }
        } // For block loop
#if DEBUG
        for (bid_t b = 0; b < nblocks - 1; b++){
        m.set("maxMem_SBlock=" + std::to_string(b), maxMemUsage_staticBlock.at(b));
    }
    std::ofstream memCompLog("log/dual-bucket/memComp.csv");
    for (bid_t b = 0; b < nblocks - 1; b++){
        memCompLog << blockMemComp.at(b).graphfile << "," << blockMemComp.at(b).bucket << ","
                   << blockMemComp.at(b).currentWalkArray << "," << blockMemComp.at(b).walkPool
                   << std::endl;
    }
    uint64_t totalSteps = 0;
    for (tid_t t = 0; t < 16; t++){
        totalSteps += userprogram.totalSteps_t.at(t);
    }
    m.set("total-steps", totalSteps);
#endif
        m.set("avg-mem-usage", memUsage);
        m.set("max-mem-usage", maxMemUsage);
        m.stop_time("1_ProcessTime");
    }
};

DualBucketEngine::DualBucketEngine(const std::string &_base_filename, unsigned long long int _blocksize_kb,
                                   bid_t _nblocks, bid_t _nmblocks, metrics &_m)
        : base_filename(_base_filename), blocksize_kb(_blocksize_kb), nblocks(_nblocks), nmblocks(_nmblocks), m(_m) {
    // membudget_mb = get_option_int("membudget_mb", 1024);
    exec_threads = get_option_int("nThreads", omp_get_max_threads());
    omp_set_num_threads(exec_threads);
    load_block_range(base_filename, blocksize_kb, blocks);
    logstream(LOG_INFO) << "block_range loaded!" << std::endl;
    nvertices = num_vertices();
    walk_manager = new WalkManager(m, nblocks, exec_threads, base_filename);
    logstream(LOG_INFO) << "walk_manager created!" << std::endl;
    if (nmblocks < 2) {
        if (nblocks > 1) {
            logstream(LOG_ERROR) << "In memory block number should at least be 2" << std::endl;
            abort();
        }
    }
    csrbuf = (vid_t **) malloc(nmblocks * sizeof(vid_t *));
    for (bid_t b = 0; b < nmblocks; b++) {
        csrbuf[b] = (vid_t *) malloc(blocksize_kb * 1024);
        // csrbuf[b] = (vid_t *)mmap(NULL, blocksize_kb*1024,
        //         PROT_READ | PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS
        //         //| MAP_HUGETLB | MAP_HUGE_2MB
        //         , 0, 0);
        // if(csrbuf[b] == MAP_FAILED){
        //     printf("%lld\n",b*blocksize_kb*1024);
        //     perror("csrbuf alloc mmap");
        //     exit(-1);
        // }
    }
    logstream(LOG_INFO) << "csrbuf malloced!" << std::endl;
    beg_posbuf = (eid_t **) malloc(nmblocks * sizeof(eid_t *));
    for (bid_t b = 0; b < nmblocks; b++) beg_posbuf[b] = nullptr;
    inMemIndex = (bid_t *) malloc(nblocks * sizeof(bid_t));
    for (bid_t b = 0; b < nblocks; b++) inMemIndex[b] = nmblocks;
    cmblocks = 0;

    // m.start_time("__g_loadSubGraph_filename");
    std::string invlname = fidname(base_filename, 0); //only 1 file
    std::string beg_posname = invlname + ".beg_pos";
    std::string csrname = invlname + ".csr";
    // m.stop_time("__g_loadSubGraph_filename");
    // m.start_time("__g_loadSubGraph_open_begpos");
    beg_posf = open(beg_posname.c_str(), O_RDONLY | O_CREAT, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
    // m.stop_time("__g_loadSubGraph_open_begpos");
    // m.start_time("__g_loadSubGraph_open_csr");
    csrf = open(csrname.c_str(), O_RDONLY | O_CREAT, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
    // m.stop_time("__g_loadSubGraph_open_csr");
    // m.start_time("__g_loadSubGraph_if_open_success");
    if (csrf < 0 || beg_posf < 0) {
        logstream(LOG_FATAL) << "Could not load :" << csrname << " or " << beg_posname << ", error: " << strerror(errno)
                             << std::endl;
    }
    assert(csrf > 0 && beg_posf > 0);
    // m.stop_time("__g_loadSubGraph_if_open_success");

    _m.set("file", _base_filename);
    _m.set("engine", "default");
    _m.set("nBlocks", (size_t) nblocks);

    print_config();
#if DEBUG
    maxMemUsage_staticBlock.assign(_nblocks, 0);
    MemoryComponent memoryComponent;
    blockMemComp.assign(_nblocks, memoryComponent);
#endif
}

DualBucketEngine::~DualBucketEngine() {
    delete walk_manager;

    if (inMemIndex != NULL) free(inMemIndex);
    if (blocks != NULL) free(blocks);

    for (bid_t b = 0; b < cmblocks; b++) {
        if (beg_posbuf[b] != NULL) free(beg_posbuf[b]);
        if (csrbuf[b] != NULL) free(csrbuf[b]);
        // munmap(csrbuf[b], blocksize_kb*1024);
    }
    if (beg_posbuf != NULL) free(beg_posbuf);
    if (csrbuf != NULL) free(csrbuf);

    close(beg_posf);
    close(csrf);
}

void
DualBucketEngine::load_block_range(std::string base_filename, unsigned long long int blocksize_kb, vid_t *&blocks,
                                   bool allowfail) {
    std::string blockrangefile = blockrangename(base_filename, blocksize_kb);
    std::ifstream brf(blockrangefile.c_str());

    if (!brf.good()) {
        logstream(LOG_ERROR) << "Could not load block range file: " << blockrangefile << std::endl;
    }
    assert(brf.good());

    /* block中存的实际上是vertex id */
    blocks = (vid_t *) malloc((nblocks + 1) * sizeof(vid_t));
    vid_t en;
    for (bid_t i = 0; i < nblocks + 1; i++) {
        assert(!brf.eof());
        brf >> en;
        blocks[i] = en;
    }
    for (bid_t i = nblocks - 1; i < nblocks; i++) {
        logstream(LOG_INFO) << "last shard: " << blocks[i] << " - " << blocks[i + 1] << std::endl;
    }
    brf.close();
}

void DualBucketEngine::loadSubGraph(bid_t p, eid_t *&beg_pos, vid_t *&csr, vid_t *nverts, eid_t *nedges) {
    m.start_time("g_loadSubGraph");


    // m.start_time("__g_loadSubGraph_malloc_begpos");
    /* read beg_pos file */
    *nverts = blocks[p + 1] - blocks[p];
    beg_pos = (eid_t *) malloc((*nverts + 1) * sizeof(eid_t));
    // m.stop_time("__g_loadSubGraph_malloc_begpos");
    // beg_pos=(eid_t *)mmap(NULL,(size_t)(*nverts+1)*sizeof(eid_t),
    //         PROT_READ | PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, 0, 0);
    // if(beg_pos == MAP_FAILED)
    // {
    //     printf("%ld\n",(size_t)(*nverts+1)*sizeof(eid_t));
    //     perror("beg_pos alloc mmap");
    //     exit(-1);
    // }
//    m.start_time("z__g_loadSubGraph_read_begpos");
    preada(beg_posf, beg_pos, (size_t) (*nverts + 1) * sizeof(eid_t), (size_t) blocks[p] * sizeof(eid_t));
//    m.stop_time("z__g_loadSubGraph_read_begpos");
    /* read csr file */
    m.start_time("z__g_loadSubGraph_realloc_csr");
    *nedges = beg_pos[*nverts] - beg_pos[0];
    if (*nedges * sizeof(vid_t) > blocksize_kb * 1024) {
        csr = (vid_t *) realloc(csr, (*nedges) * sizeof(vid_t));
    }
    m.stop_time("z__g_loadSubGraph_realloc_csr");
    m.start_time("z__g_loadSubGraph_read_csr");
    preada(csrf, csr, (*nedges) * sizeof(vid_t), beg_pos[0] * sizeof(vid_t));
    m.stop_time("z__g_loadSubGraph_read_csr");

    m.stop_time("g_loadSubGraph");
}

void DualBucketEngine::loadPreVertexSubGraph(bid_t p, eid_t *&beg_pos, vid_t *&csr, vid_t *nverts, eid_t *nedges) {
    m.start_time("loadPreSubGraph");


    // m.start_time("__g_loadSubGraph_malloc_begpos");
    /* read beg_pos file */
    *nverts = blocks[p + 1] - blocks[p];
    beg_pos = (eid_t *) malloc((*nverts + 1) * sizeof(eid_t));
    // m.stop_time("__g_loadSubGraph_malloc_begpos");
    // beg_pos=(eid_t *)mmap(NULL,(size_t)(*nverts+1)*sizeof(eid_t),
    //         PROT_READ | PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, 0, 0);
    // if(beg_pos == MAP_FAILED)
    // {
    //     printf("%ld\n",(size_t)(*nverts+1)*sizeof(eid_t));
    //     perror("beg_pos alloc mmap");
    //     exit(-1);
    // }
    m.start_time("loadPreSubGraph_read_begpos");
    preada(beg_posf, beg_pos, (size_t) (*nverts + 1) * sizeof(eid_t), (size_t) blocks[p] * sizeof(eid_t));
    m.stop_time("loadPreSubGraph_read_begpos");
    /* read csr file */
    *nedges = beg_pos[*nverts] - beg_pos[0];
    csr = (vid_t *) malloc((*nedges) * sizeof(vid_t));
    m.start_time("loadPreSubGraph_read_csr");
    preada(csrf, csr, (*nedges) * sizeof(vid_t), beg_pos[0] * sizeof(vid_t));
    m.stop_time("loadPreSubGraph_read_csr");

    m.stop_time("loadPreSubGraph");
}

void DualBucketEngine::findSubGraph(bid_t p, eid_t *&beg_pos, vid_t *&csr, vid_t *nverts, eid_t *nedges) {
    m.start_time("2_findSubGraph");
    if (inMemIndex[p] == nmblocks) {//the block is not in memory
        bid_t swapin;
        if (cmblocks < nmblocks) {
            swapin = cmblocks++;
        } else {
            bid_t minmwb = swapOut();
            swapin = inMemIndex[minmwb];
            inMemIndex[minmwb] = nmblocks;
            assert(swapin < nmblocks);
            if (beg_posbuf[swapin] != NULL) free(beg_posbuf[swapin]);
            // munmap(beg_posbuf[swapin], sizeof(eid_t)*(blocks[minmwb+1] - blocks[minmwb] + 1));
        }
        loadSubGraph(p, beg_posbuf[swapin], csrbuf[swapin], nverts, nedges);
        inMemIndex[p] = swapin;
    } else {
    }
    beg_pos = beg_posbuf[inMemIndex[p]];
    csr = csrbuf[inMemIndex[p]];
    m.stop_time("2_findSubGraph");
}


bid_t DualBucketEngine::swapOut() {
    m.start_time("z_g_swapOut");
    wid_t minmw = 0xffffffff;
    bid_t minmwb = 0;
    for (bid_t b = 0; b < nblocks; b++) {
        if (inMemIndex[b] < nmblocks && walk_manager->nWalks[b] < minmw) {
            minmw = walk_manager->nWalks[b];
            minmwb = b;
        }
    }
    m.start_time("z_g_swapOut");
    return minmwb;
}

#if FULLY_LOAD
void DualBucketEngine::run_fullyLoadTest(vid_t nRootVertices, SecondOrderRW &userprogram, float prob, float th) {
    // srand((unsigned)time(NULL));
    std::vector<eid_t> nEdges;
    nEdges.assign(nblocks, 0);
    vid_t tmp;
    for (bid_t b = 0; b < nblocks; b++){
        userprogram.vertexIO.getBlockSize(b, tmp, nEdges.at(b));
    }
    m.start_time("2_StartWalks");

    std::ofstream blockLogFile_dynamicBlocks;
    std::ofstream blockLogFile_staticBlocks;
    /* 以追加的方式写log文件 */
    blockLogFile_dynamicBlocks.open(loadTestOutPutFileName_dynamicBlocks);
    blockLogFile_dynamicBlocks << "dynamic-block-id,static-block-id,nverts,nedges,nwalks,load-block-time,exec-block-time" << std::endl;

    blockLogFile_staticBlocks.open(loadTestOutPutFileName_staticBlocks);
    if (!blockLogFile_staticBlocks.is_open()){
        abort();
    }
    blockLogFile_staticBlocks << "static-block-id,nverts,nedges,nwalks,load-block-time,exec-block-time" << std::endl;

    /* initiate walks */
    bid_t preInitBlock = INVALID_BID;
    for (bid_t b = 0; b < nblocks; b++) {
        swapBlock(b, preInitBlock);
        userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
        userprogram.initWalks(b, *walk_manager, m);
#if INFO
        logstream(LOG_INFO) << "initiated block" << b << std::endl;
#endif
        preInitBlock = b;
    }
    walk_manager->updateWalkNum();
    std::cout << "currentNWalks before run: " << walk_manager->currentNWalks << std::endl;
    m.stop_time("2_StartWalks");

    gettimeofday(&start, NULL);
    m.start_time("1_ProcessTime");

    vid_t nverts;
    eid_t nedges;
    int blockcount = 0;
    bid_t staticBlock, dynamicBlock;
    bid_t preStaticBlock = INVALID_BID, preDynamicBlock = INVALID_BID;
    bid_t lostBlockId = INVALID_BID;
    unLoadAllBlock();
    if (nblocks == 1) {
        eid_t *beg_pos;
        vid_t *csr;
        loadBlock(0, beg_pos, csr);
        wid_t nWalks = walk_manager->getCurrentWalks(staticBlock);
        {
            logstream(LOG_DEBUG) << runtime() << "s : staticBlockCount: " << blockcount << std::endl;
            logstream(LOG_INFO) << "currentNWalks = " << walk_manager->currentNWalks << ", nwalks[" << staticBlock
                                << "] = " << nWalks << std::endl;
        }
#pragma for parallel (dynamic)
        for (wid_t w = 0; w < nWalks; w++) {
            userprogram.updateWalk(m, walk_manager->currentWalks[w], 0, beg_pos, csr, *walk_manager);
        }
    }
    while (!userprogram.allWalksFinished(*walk_manager)) {
        for (staticBlock = 0; staticBlock < nblocks - 1; staticBlock++) {
            blockcount++;
            if (!walk_manager->nWalks[staticBlock]) {
                continue;
            }
            Timer loadTimer_static;
            loadTimer_static.start();
            if (inMemIndex[staticBlock] < nmblocks) {
                lostBlockId = preStaticBlock;
            } else {
                swapBlock(staticBlock, preStaticBlock);
            }
            loadTimer_static.stop();
            wid_t nWalks = walk_manager->getCurrentWalks(staticBlock);
            walk_manager->clearRecoredWalkNum(staticBlock);
//             if(blockcount % (nBlocks/100+1)==1)
//             if(blockcount % (1024*1024*1024/nedges+1) == 1)
            {
                logstream(LOG_DEBUG) << runtime() << "s : staticBlockCount: " << blockcount << std::endl;
                logstream(LOG_INFO) << "currentNWalks = " << walk_manager->currentNWalks << ", nwalks[" << staticBlock
                                    << "] = " << nWalks << std::endl;
            }
            walk_manager->collectBucket(staticBlock);
            Timer executeTimer_static;
            executeTimer_static.start();
            for (dynamicBlock = staticBlock + 1; dynamicBlock < nblocks; dynamicBlock++) {
#if MULTI_BUFFER
                m.start_time("extend-buffer");
                walk_manager->exendBucket(dynamicBlock);
                walk_manager->clearExtendBucket(dynamicBlock);
                m.stop_time("extend-buffer");
#endif
                if (!walk_manager->bucket[dynamicBlock].size()) {
                    continue;
                }
                wid_t totalWalks = walk_manager->bucket[dynamicBlock].size();
                auto clockStart = std::chrono::system_clock::now();
                if (preDynamicBlock == staticBlock) {
                    assert(lostBlockId != INVALID_BID);
                    swapBlock(dynamicBlock, lostBlockId);
                    lostBlockId = INVALID_BID;
                } else {
                    swapBlock(dynamicBlock, preDynamicBlock);
                }
                auto swapBlockEnd = std::chrono::system_clock::now();
                userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
                m.start_time("4_UpdateWalk");
#if MULTI_THREAD
#pragma omp parallel for schedule(dynamic)
#endif
                for (wid_t w = 0; w < walk_manager->bucket[dynamicBlock].size(); w++) {
                    userprogram.processWalk(m, walk_manager->bucket[dynamicBlock][w], staticBlock, dynamicBlock,
                                            *walk_manager, true, true);
                }
                preDynamicBlock = dynamicBlock;
                m.stop_time("4_UpdateWalk");
                auto execBucketEnd = std::chrono::system_clock::now();
                double swapBlock_s = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
                        swapBlockEnd - clockStart).count())
                                * std::chrono::microseconds::period::num
                                / std::chrono::microseconds::period::den;
                double execBucket_s = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
                        execBucketEnd - clockStart).count())
                                      * std::chrono::microseconds::period::num
                                      / std::chrono::microseconds::period::den;
                blockLogFile_dynamicBlocks << dynamicBlock << ","
                             << staticBlock << ","
                             << blocks[dynamicBlock + 1] - blocks[dynamicBlock] << ","
                             << nEdges.at(dynamicBlock) << ","
                             << totalWalks << ","
                             << swapBlock_s << ","
                             << execBucket_s << std::endl;
            }
            executeTimer_static.stop();
            blockLogFile_staticBlocks << staticBlock << ","
                                      << blocks[staticBlock + 1] - blocks[staticBlock] << ","
                                      << nEdges.at(staticBlock) << ","
                                      << nWalks << ","
                                      << loadTimer_static.duration_s() << ","
                                      << executeTimer_static.duration_s() << std::endl;
#if INFO | DEBUG
            double usage = get_memory();
            if (usage > maxMemUsage) {
                maxMemUsage = usage;
            }
#endif
#if DEBUG
            if (usage > maxMemUsage_staticBlock.at(staticBlock)){
                maxMemUsage_staticBlock.at(staticBlock) = usage;
//                if (staticBlock == 0){
//                    /* two block currently in memory is staticBlock and preDynamicBlock */
//                    vid_t nVertsStaticBlock = blocks[staticBlock + 1] - blocks[staticBlock];
//                    vid_t nVertsPreDynamicBlock = blocks[preDynamicBlock + 1] - blocks[preDynamicBlock];
//                    bid_t staticBlockInMemIndex = 0;
//                    bid_t preDynamicBlockInMemIndex = 1;
//                    /* FIXME: 这里实际上不是很准确，如果涉及到lostBlock，这里将统计错误，但考虑到子图大小差别不大，这里误差应该较小 */
//                    if (inMemIndex[preDynamicBlock] == 0){
//                        preDynamicBlockInMemIndex = 0;
//                        staticBlockInMemIndex = 1;
//                    }
//                    eid_t nEdgesStaticBlock = beg_posbuf[staticBlockInMemIndex][nVertsStaticBlock] - beg_posbuf[staticBlockInMemIndex][0];
//                    eid_t nEdgesPreDynamicBlock = beg_posbuf[preDynamicBlockInMemIndex][nVertsPreDynamicBlock] - beg_posbuf[preDynamicBlockInMemIndex][0];
//                    size_t csrByte = (nEdgesPreDynamicBlock + nEdgesStaticBlock) * sizeof(vid_t);
//                    size_t begPosByte = (nVertsStaticBlock + nVertsPreDynamicBlock) * sizeof(vid_t);
//                    memoryComponent_b0.graphfile = static_cast<double>(csrByte + begPosByte) / 1024 / 1024;
//                    size_t currentWalkByte = walk_manager->curBlockNWalks * sizeof(WalkDataType);
//                    memoryComponent_b0.currentWalkArray = static_cast<double>(currentWalkByte) / 1024 / 1024;
//                    size_t bucketsByte = 0;
//                    for (bid_t b = 0; b < nblocks; b++){
//                        /* vec.size() and vec.capacity differ a lot */
//                        bucketsByte += walk_manager->preBucket[b].capacity() * sizeof(WalkDataType);
//                        bucketsByte += walk_manager->curBucket[b].capacity() * sizeof(WalkDataType);
//                    }
//                    memoryComponent_b0.bucket = static_cast<double>(bucketsByte) / 1024 / 1024;
//                    size_t walkPoolSize_b = 0;
//                    for (tid_t t = 0; t < 16; t++){
//                        for (bid_t b = 0; b < nblocks; b++){
//                            walkPoolSize_b += walk_manager->walkPool[t][b].size_w * sizeof(WalkDataType);
//                        }
//                    }
//                    memoryComponent_b0.walkPool = static_cast<double>(walkPoolSize_b) / 1024 / 1024;
//                }
                /* two block currently in memory is staticBlock and preDynamicBlock */
                vid_t nVertsStaticBlock = blocks[staticBlock + 1] - blocks[staticBlock];
                vid_t nVertsPreDynamicBlock = blocks[preDynamicBlock + 1] - blocks[preDynamicBlock];
                bid_t staticBlockInMemIndex = 0;
                bid_t preDynamicBlockInMemIndex = 1;
                /* FIXME: 这里实际上不是很准确，如果涉及到lostBlock，这里将统计错误，但考虑到子图大小差别不大，这里误差应该较小 */
                if (inMemIndex[preDynamicBlock] == 0){
                    preDynamicBlockInMemIndex = 0;
                    staticBlockInMemIndex = 1;
                }
                eid_t nEdgesStaticBlock = beg_posbuf[staticBlockInMemIndex][nVertsStaticBlock] - beg_posbuf[staticBlockInMemIndex][0];
                eid_t nEdgesPreDynamicBlock = beg_posbuf[preDynamicBlockInMemIndex][nVertsPreDynamicBlock] - beg_posbuf[preDynamicBlockInMemIndex][0];
                size_t csrByte = (nEdgesPreDynamicBlock + nEdgesStaticBlock) * sizeof(vid_t);
                size_t begPosByte = (nVertsStaticBlock + nVertsPreDynamicBlock) * sizeof(vid_t);
                blockMemComp.at(staticBlock).graphfile = static_cast<double>(csrByte + begPosByte) / 1024 / 1024;
                size_t currentWalkByte = walk_manager->curBlockNWalks * sizeof(WalkDataType);
                blockMemComp.at(staticBlock).currentWalkArray = static_cast<double>(currentWalkByte) / 1024 / 1024;
                size_t bucketsByte = 0;
                for (bid_t b = 0; b < nblocks; b++){
                    /* vec.size() and vec.capacity differ a lot */
                    bucketsByte += walk_manager->preBucket[b].size() * sizeof(WalkDataType);
                    bucketsByte += walk_manager->curBucket[b].size() * sizeof(WalkDataType);
                }
                blockMemComp.at(staticBlock).bucket = static_cast<double>(bucketsByte) / 1024 / 1024;
                size_t walkPoolSize_b = 0;
                for (tid_t t = 0; t < 16; t++){
                    for (bid_t b = 0; b < nblocks; b++){
                        walkPoolSize_b += walk_manager->walkPool[t][b].size_w * sizeof(WalkDataType);
                    }
                }
                blockMemComp.at(staticBlock).walkPool = static_cast<double>(walkPoolSize_b) / 1024 / 1024;
            }
#endif
#if INFO | DEBUG
            memUsage = ((memUsage * memCount) + usage) / (memCount + 1);
            memCount++;
#endif
            walk_manager->updateWalkNum();
            walk_manager->clearBucket();
            preStaticBlock = staticBlock;
        }
    } // For block loop
    blockLogFile_dynamicBlocks.close();
#if DEBUG
    for (bid_t b = 0; b < nblocks - 1; b++){
        m.set("maxMem_SBlock=" + std::to_string(b), maxMemUsage_staticBlock.at(b));
    }
    std::ofstream memCompLog("log/dual-bucket/memComp.csv");
    for (bid_t b = 0; b < nblocks - 1; b++){
        memCompLog << blockMemComp.at(b).graphfile << "," << blockMemComp.at(b).bucket << ","
                   << blockMemComp.at(b).currentWalkArray << "," << blockMemComp.at(b).walkPool
                   << std::endl;
    }
    uint64_t totalSteps = 0;
    for (tid_t t = 0; t < 16; t++){
        totalSteps += userprogram.totalSteps_t.at(t);
    }
    m.set("total-steps", totalSteps);
#endif
    m.set("avg-mem-usage", memUsage);
    m.set("max-mem-usage", maxMemUsage);
    m.stop_time("1_ProcessTime");
}
#endif

#if NO_LOAD_TEST
void DualBucketEngine::run_noLoadTest(vid_t nRootVertices, SecondOrderRW &userprogram, float prob, float th) {
    std::vector<std::vector<double>> ths;
    ths.resize(nblocks);
    std::ifstream ths_txt(thsFile_dynamic);
    for (bid_t dynamicBlock = 1; dynamicBlock < nblocks; dynamicBlock++){
        ths.at(dynamicBlock).resize(dynamicBlock);
        for (bid_t staticBlock = 0; staticBlock < dynamicBlock; staticBlock++){
            double temp;
            ths_txt >> temp;
            ths.at(dynamicBlock).at(staticBlock) = temp;
        }
    }
    ths_txt.close();
    // srand((unsigned)time(NULL));
    std::vector<eid_t> nEdges;
    nEdges.assign(nblocks, 0);
    vid_t tmp;
    for (bid_t b = 0; b < nblocks; b++){
        userprogram.vertexIO.getBlockSize(b, tmp, nEdges.at(b));
    }
    m.start_time("2_StartWalks");
    std::ofstream blockLogFile;
//    /* 以追加的方式写文件 */
//    blockLogFile.open(LOAD_TRADEOFF_FILE, std::ios_base::app);
    /* 以覆盖方式写文件 */
    blockLogFile.open(loadTestOutPutFileName_dynamicBlocks);
    blockLogFile << "static-block-id,dynamic-block-id,nverts,nedges,nwalks,exec-bucket-time" << std::endl;
    bid_t preInitBlock = 0;
    for (bid_t b = 0; b < nblocks; b++) {
        swapBlock(b, preInitBlock);
        userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
        userprogram.initWalks(b, *walk_manager, m);
#if INFO
        logstream(LOG_INFO) << "initiated block" << b << std::endl;
#endif
        preInitBlock = b;
    }
    walk_manager->updateWalkNum();
    std::cout << "currentNWalks before run: " << walk_manager->currentNWalks << std::endl;
    m.stop_time("2_StartWalks");

    gettimeofday(&start, NULL);
    m.start_time("1_ProcessTime");

    vid_t nverts;
    eid_t nedges;
    int blockcount = 0;
    bid_t staticBlock, dynamicBlock;
    bid_t preStaticBlock = 0, preDynamicBlock = 1;
    bid_t lostBlockId = INVALID_BID;
    unLoadAllBlock();
    if (nblocks == 1) {
        eid_t *beg_pos;
        vid_t *csr;
        loadBlock(0, beg_pos, csr);
        wid_t nWalks = walk_manager->getCurrentWalks(staticBlock);
        {
            logstream(LOG_DEBUG) << runtime() << "s : staticBlockCount: " << blockcount << std::endl;
            logstream(LOG_INFO) << "currentNWalks = " << walk_manager->currentNWalks << ", nwalks[" << staticBlock
                                << "] = " << nWalks << std::endl;
        }
#pragma for parallel (dynamic)
        for (wid_t w = 0; w < nWalks; w++) {
            userprogram.updateWalk(m, walk_manager->currentWalks[w], 0, beg_pos, csr, *walk_manager);
        }
    }
    while (!userprogram.allWalksFinished(*walk_manager)) {
        for (staticBlock = 0; staticBlock < nblocks - 1; staticBlock++) {
            blockcount++;
            if (!walk_manager->nWalks[staticBlock]) {
                continue;
            }
            if (inMemIndex[staticBlock] < nmblocks) {
                lostBlockId = preStaticBlock;
            } else {
                swapBlock(staticBlock, preStaticBlock);
                if (userprogram.vertexIO.blockInMem(preStaticBlock)){
                    lostBlockId = preStaticBlock;
                }
//                std::cout << "static block swap in: " << staticBlock << ", swap out: " << preStaticBlock <<std::endl;
            }
#if OUTPUT_NO_LOAD_DATA
            staticBlock_global = staticBlock;
#endif
            wid_t nWalks = walk_manager->getCurrentWalks(staticBlock);
            walk_manager->clearRecoredWalkNum(staticBlock);
//             if(blockcount % (nBlocks/100+1)==1)
//             if(blockcount % (1024*1024*1024/nedges+1) == 1)
            {
                logstream(LOG_DEBUG) << runtime() << "s : staticBlockCount: " << blockcount << std::endl;
                logstream(LOG_INFO) << "currentNWalks = " << walk_manager->currentNWalks << ", nwalks[" << staticBlock
                                    << "] = " << nWalks << std::endl;
            }
            walk_manager->collectBucket(staticBlock);
            for (dynamicBlock = staticBlock + 1; dynamicBlock < nblocks; dynamicBlock++) {
                if (!walk_manager->bucket[dynamicBlock].size()) {
                    continue;
                }
#if OUTPUT_NO_LOAD_DATA
                dynamicBlock_global = dynamicBlock;
                getCSR_dynamicBlock_global_ms = 0;
#endif
                wid_t totalWalks = walk_manager->bucket[dynamicBlock].size();
                bool loadDynamicBlock = true;
#if USE_TRAINED_THS
                if (ths.at(dynamicBlock).at(staticBlock) > 0 && totalWalks < (blocks[dynamicBlock + 1] - blocks[dynamicBlock]) * ths.at(dynamicBlock).at(staticBlock) && !userprogram.vertexIO.blockInMem(dynamicBlock)){
#elif USE_PRE_DEFINE_TH
                if (totalWalks < (blocks[dynamicBlock + 1] - blocks[dynamicBlock]) * NO_LOAD_TH && !userprogram.vertexIO.blockInMem(dynamicBlock)){
#endif
                    loadDynamicBlock = false;
                }else{
                    if (preDynamicBlock == staticBlock) {
                        assert(lostBlockId != INVALID_BID);
                        swapBlock(dynamicBlock, lostBlockId);
                        lostBlockId = INVALID_BID;
                    } else {
                        if (lostBlockId != INVALID_BID && !userprogram.vertexIO.blockInMem(preDynamicBlock) && userprogram.vertexIO.blockInMem(lostBlockId)){
                            swapBlock(dynamicBlock, lostBlockId);
                            lostBlockId = INVALID_BID;
                        }else{
                            swapBlock(dynamicBlock, preDynamicBlock);
                        }
//                        std::cout << "dynamic block swap in: " << dynamicBlock << ", swap out: " << preDynamicBlock << std::endl;
                    }
                    userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
                }
#if OUTPUT_NO_LOAD_DATA
                auto execBucketStart = std::chrono::system_clock::now();
#endif
                m.start_time("4_UpdateWalk");
#pragma omp parallel for schedule(dynamic)
                for (wid_t w = 0; w < walk_manager->bucket[dynamicBlock].size(); w++) {
                    userprogram.processWalk(m, walk_manager->bucket[dynamicBlock][w], staticBlock, dynamicBlock,
                                            *walk_manager, loadDynamicBlock);
                }
                if (loadDynamicBlock){
                    preDynamicBlock = dynamicBlock;
                }
                m.stop_time("4_UpdateWalk");
#if OUTPUT_NO_LOAD_DATA
                auto execBucketEnd = std::chrono::system_clock::now();
                double execBucket_s = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
                        execBucketEnd - execBucketStart).count())
                                      * std::chrono::microseconds::period::num
                                      / std::chrono::microseconds::period::den;
                if (!loadDynamicBlock){
                    blockLogFile << staticBlock << ","
                                 << dynamicBlock << ","
                                 << blocks[dynamicBlock + 1] - blocks[dynamicBlock] << ","
                                 << nEdges.at(dynamicBlock) << ","
                                 << walk_manager->bucket[dynamicBlock].size() << ","
                                 << execBucket_s << std::endl;
                }
#endif
            }
#if INFO | DEBUG
            double usage = get_memory();
            if (usage > maxMemUsage) {
                maxMemUsage = usage;
            }
#endif
#if DEBUG
            if (usage > maxMemUsage_staticBlock.at(staticBlock)){
                maxMemUsage_staticBlock.at(staticBlock) = usage;
//                if (staticBlock == 0){
//                    /* two block currently in memory is staticBlock and preDynamicBlock */
//                    vid_t nVertsStaticBlock = blocks[staticBlock + 1] - blocks[staticBlock];
//                    vid_t nVertsPreDynamicBlock = blocks[preDynamicBlock + 1] - blocks[preDynamicBlock];
//                    bid_t staticBlockInMemIndex = 0;
//                    bid_t preDynamicBlockInMemIndex = 1;
//                    /* FIXME: 这里实际上不是很准确，如果涉及到lostBlock，这里将统计错误，但考虑到子图大小差别不大，这里误差应该较小 */
//                    if (inMemIndex[preDynamicBlock] == 0){
//                        preDynamicBlockInMemIndex = 0;
//                        staticBlockInMemIndex = 1;
//                    }
//                    eid_t nEdgesStaticBlock = beg_posbuf[staticBlockInMemIndex][nVertsStaticBlock] - beg_posbuf[staticBlockInMemIndex][0];
//                    eid_t nEdgesPreDynamicBlock = beg_posbuf[preDynamicBlockInMemIndex][nVertsPreDynamicBlock] - beg_posbuf[preDynamicBlockInMemIndex][0];
//                    size_t csrByte = (nEdgesPreDynamicBlock + nEdgesStaticBlock) * sizeof(vid_t);
//                    size_t begPosByte = (nVertsStaticBlock + nVertsPreDynamicBlock) * sizeof(vid_t);
//                    memoryComponent_b0.graphfile = static_cast<double>(csrByte + begPosByte) / 1024 / 1024;
//                    size_t currentWalkByte = walk_manager->curBlockNWalks * sizeof(WalkDataType);
//                    memoryComponent_b0.currentWalkArray = static_cast<double>(currentWalkByte) / 1024 / 1024;
//                    size_t bucketsByte = 0;
//                    for (bid_t b = 0; b < nblocks; b++){
//                        /* vec.size() and vec.capacity differ a lot */
//                        bucketsByte += walk_manager->preBucket[b].capacity() * sizeof(WalkDataType);
//                        bucketsByte += walk_manager->curBucket[b].capacity() * sizeof(WalkDataType);
//                    }
//                    memoryComponent_b0.bucket = static_cast<double>(bucketsByte) / 1024 / 1024;
//                    size_t walkPoolSize_b = 0;
//                    for (tid_t t = 0; t < 16; t++){
//                        for (bid_t b = 0; b < nblocks; b++){
//                            walkPoolSize_b += walk_manager->walkPool[t][b].size_w * sizeof(WalkDataType);
//                        }
//                    }
//                    memoryComponent_b0.walkPool = static_cast<double>(walkPoolSize_b) / 1024 / 1024;
//                }
                /* two block currently in memory is staticBlock and preDynamicBlock */
                vid_t nVertsStaticBlock = blocks[staticBlock + 1] - blocks[staticBlock];
                vid_t nVertsPreDynamicBlock = blocks[preDynamicBlock + 1] - blocks[preDynamicBlock];
                bid_t staticBlockInMemIndex = 0;
                bid_t preDynamicBlockInMemIndex = 1;
                /* FIXME: 这里实际上不是很准确，如果涉及到lostBlock，这里将统计错误，但考虑到子图大小差别不大，这里误差应该较小 */
                if (inMemIndex[preDynamicBlock] == 0){
                    preDynamicBlockInMemIndex = 0;
                    staticBlockInMemIndex = 1;
                }
                eid_t nEdgesStaticBlock = beg_posbuf[staticBlockInMemIndex][nVertsStaticBlock] - beg_posbuf[staticBlockInMemIndex][0];
                eid_t nEdgesPreDynamicBlock = beg_posbuf[preDynamicBlockInMemIndex][nVertsPreDynamicBlock] - beg_posbuf[preDynamicBlockInMemIndex][0];
                size_t csrByte = (nEdgesPreDynamicBlock + nEdgesStaticBlock) * sizeof(vid_t);
                size_t begPosByte = (nVertsStaticBlock + nVertsPreDynamicBlock) * sizeof(vid_t);
                blockMemComp.at(staticBlock).graphfile = static_cast<double>(csrByte + begPosByte) / 1024 / 1024;
                size_t currentWalkByte = walk_manager->curBlockNWalks * sizeof(WalkDataType);
                blockMemComp.at(staticBlock).currentWalkArray = static_cast<double>(currentWalkByte) / 1024 / 1024;
                size_t bucketsByte = 0;
                for (bid_t b = 0; b < nblocks; b++){
                    /* vec.size() and vec.capacity differ a lot */
                    bucketsByte += walk_manager->preBucket[b].size() * sizeof(WalkDataType);
                    bucketsByte += walk_manager->curBucket[b].size() * sizeof(WalkDataType);
                }
                blockMemComp.at(staticBlock).bucket = static_cast<double>(bucketsByte) / 1024 / 1024;
                size_t walkPoolSize_b = 0;
                for (tid_t t = 0; t < 16; t++){
                    for (bid_t b = 0; b < nblocks; b++){
                        walkPoolSize_b += walk_manager->walkPool[t][b].size_w * sizeof(WalkDataType);
                    }
                }
                blockMemComp.at(staticBlock).walkPool = static_cast<double>(walkPoolSize_b) / 1024 / 1024;
            }
#endif
#if INFO | DEBUG
            memUsage = ((memUsage * memCount) + usage) / (memCount + 1);
            memCount++;
#endif
            walk_manager->updateWalkNum();
            walk_manager->clearBucket();
            preStaticBlock = staticBlock;
        }
    } // For block loop
    blockLogFile.close();
#if DEBUG
    for (bid_t b = 0; b < nblocks - 1; b++){
        m.set("maxMem_SBlock=" + std::to_string(b), maxMemUsage_staticBlock.at(b));
    }
    std::ofstream memCompLog("log/dual-bucket/memComp.csv");
    for (bid_t b = 0; b < nblocks - 1; b++){
        memCompLog << blockMemComp.at(b).graphfile << "," << blockMemComp.at(b).bucket << ","
                   << blockMemComp.at(b).currentWalkArray << "," << blockMemComp.at(b).walkPool
                   << std::endl;
    }
    uint64_t totalSteps = 0;
    for (tid_t t = 0; t < 16; t++){
        totalSteps += userprogram.totalSteps_t.at(t);
    }
    m.set("total-steps", totalSteps);
#endif
    m.set("avg-mem-usage", memUsage);
    m.set("max-mem-usage", maxMemUsage);
    m.stop_time("1_ProcessTime");
}
#endif

#if ONDEMAND_LOAD
void DualBucketEngine::run_onDemandLoad(vid_t nRootVertices, SecondOrderRW &userprogram, float prob, float th) {
        std::vector<double> ths_dynamic;
        ths_dynamic.resize(nblocks);
        std::ifstream ths_txt_dynamic(thsFile_dynamic);
        if (!ths_txt_dynamic.is_open()){
            assert(false && "Ths file error");
        }
        /* dynamic block读进来的第一个数代表block 0的 th，但由于block 0不会称为dynamic block，故这个th没有作用 */
        for (bid_t dynamicBlock = 0; dynamicBlock < nblocks; dynamicBlock++){
                double temp;
                ths_txt_dynamic >> temp;
                ths_dynamic.at(dynamicBlock) = temp;
        }
        ths_txt_dynamic.close();

        std::vector<double> ths_static;
        ths_static.resize(nblocks);
        std::ifstream ths_txt_static(thsFile_static);
        if (!ths_txt_static.is_open()){
            assert(false && "Ths file (static) error");
        }
        for (bid_t staticBlock = 0; staticBlock < nblocks; staticBlock++){
            double temp;
            ths_txt_static >> temp;
            ths_static.at(staticBlock) = temp;
        }
        ths_txt_static.close();

        // srand((unsigned)time(NULL));
        /* store total edge number in each block */
        std::vector<eid_t> nEdges;
        nEdges.assign(nblocks, 0);
        vid_t tmp;
        for (bid_t b = 0; b < nblocks; b++){
            userprogram.vertexIO.getBlockSize(b, tmp, nEdges.at(b));
        }
        m.start_time("2_StartWalks");
        std::ofstream blockLogFile_dynamicBlocks;
        std::ofstream blockLogFile_staticBlocks;
        //    /* 以追加的方式写文件 */
        //    blockLogFile_dynamicBlocks.open(LOAD_TRADEOFF_FILE, std::ios_base::app);
        /* 以覆盖方式写文件 */
#if OUTPUT_ON_DEMAND_DATA
        blockLogFile_dynamicBlocks.open(loadTestOutPutFileName_dynamicBlocks);
        if (!blockLogFile_dynamicBlocks.is_open()){
            abort();
        }
        blockLogFile_staticBlocks.open(loadTestOutPutFileName_staticBlocks);
        if (!blockLogFile_staticBlocks.is_open()){
            abort();
        }
#if NO_LOAD_TH
        blockLogFile_dynamicBlocks << "static-block-id,dynamic-block-id,nverts,nedges,nwalks,map-get-time,exec-bucket-time" << std::endl;
#else
        blockLogFile_dynamicBlocks << "static-block-id,dynamic-block-id,nverts,nedges,nwalks,load-block-time,exec-block-time" << std::endl;
        blockLogFile_staticBlocks << "static-block-id,nverts,nedges,nwalks,load-block-time,exec-block-time" << std::endl;
#endif
#endif
        bid_t preInitBlock = INVALID_BID;
        for (bid_t b = 0; b < nblocks; b++) {
            swapBlock(b, preInitBlock);
            userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
            userprogram.initWalks(b, *walk_manager, m);
#if INFO
            logstream(LOG_INFO) << "initiated block" << b << std::endl;
#endif
            preInitBlock = b;
        }
#if COMPUTE_EDGE_CUT
        std::cout << "Edge-Cut: " << edgeCut << std::endl;
        std::cout << "Edge-In: " << edgeIn << std::endl;
        std::cout << "Edge-Cut Percentage: " << static_cast<double>(edgeCut) / (static_cast<double>(edgeIn) + static_cast<double>(edgeCut)) << std::endl;
        return;
#endif
        walk_manager->updateWalkNum();
        std::cout << "currentNWalks before run: " << walk_manager->currentNWalks << std::endl;
        m.stop_time("2_StartWalks");

        gettimeofday(&start, NULL);
        m.start_time("1_ProcessTime");

        vid_t nverts;
        eid_t nedges;
        int blockcount = 0;
        bid_t staticBlock, dynamicBlock;
        bid_t preStaticBlock = INVALID_BID, preDynamicBlock = INVALID_BID;
        bid_t lostBlockId = INVALID_BID;
        unLoadAllBlock();
        if (nblocks == 1) {
            eid_t *beg_pos;
            vid_t *csr;
            loadBlock(0, beg_pos, csr);
            wid_t nWalks = walk_manager->getCurrentWalks(staticBlock);
            {
                logstream(LOG_DEBUG) << runtime() << "s : staticBlockCount: " << blockcount << std::endl;
                logstream(LOG_INFO) << "currentNWalks = " << walk_manager->currentNWalks << ", nwalks[" << staticBlock
                << "] = " << nWalks << std::endl;
            }
#pragma for parallel (dynamic)
for (wid_t w = 0; w < nWalks; w++) {
                userprogram.updateWalk(m, walk_manager->currentWalks[w], 0, beg_pos, csr, *walk_manager);
            }
        }

#if BLOCK_UTI
        std::ofstream blockUtiFile(blockUtiFileName);
        if (!blockUtiFile.is_open()){
            abort();
        }
        blockUtiFile << "IO Num," << "Block Id," << "Block Type," << "Used Edges," << "Total Load Edges" << std::endl;
#endif
#if ACT_VER
        std::ofstream actVerFile(actVerFileName);
        if (!actVerFile.is_open()){
            abort();
        }
        actVerFile << "IO Num," << "Block Id," << "Block Type," << "Activated Vertices," << "Total Vertices" << std::endl;
#endif
        /* main while start */
#if BLOCK_UTI
        long blockIONum = 0;
#endif
#if ACT_VER
        long blockIONum_ver = 0;
#endif
        while (!userprogram.allWalksFinished(*walk_manager)) {
            for (staticBlock = 0; staticBlock < nblocks - 1; staticBlock++) {
                blockcount++;
                Timer loadTimer_staticBlock;
                if (!walk_manager->nWalks[staticBlock]) {
                    continue;
                }
#if SWAP_BLOCK_DEBUG
                std::cout << "static-block = " << staticBlock << std::endl;
#endif
                wid_t nWalks = walk_manager->getCurrentWalks(staticBlock);
                bool loadStaticBlock = true;
#if ONDEMAND_LOAD & ONDEMAND_LOAD_STATIC
                if (inMemIndex[staticBlock] < nmblocks) {
                    lostBlockId = preStaticBlock;
#if SWAP_BLOCK_DEBUG
                    std::cout << "lost-block = pre-static-block = " << lostBlockId << std::endl;
#endif
                }else if ((staticBlock > 0 && ths_static.at(staticBlock) > 0 && nWalks < (blocks[staticBlock + 1] - blocks[staticBlock]) * ths_static.at(staticBlock))){
                    loadStaticBlock = false;
                    m.start_time("map-get");
#if OUTPUT_ON_DEMAND_DATA
                    loadTimer_staticBlock.start();
#endif
                    walk_manager->mapBlockVertex(staticBlock, userprogram.vertexIO);
                    userprogram.vertexIO.getVertexOnDemand(staticBlock, false);
#if OUTPUT_ON_DEMAND_DATA
                    loadTimer_staticBlock.stop();
#endif
                    m.stop_time("map-get");
                }else if (lostBlockId != INVALID_BID){
                    swapBlock(staticBlock, lostBlockId);
                    lostBlockId = INVALID_BID;
#if SWAP_BLOCK_DEBUG
                    std::cout << "swap out lost-block-id: " << lostBlockId << ", swap in static-block: " << staticBlock << std::endl;
                    std::cout << "lost-block = invalid" << std::endl;
#endif
                }else{
                    swapBlock(staticBlock, preStaticBlock);
#if SWAP_BLOCK_DEBUG
                    std::cout << "swap out pre-static-block: " << preStaticBlock << ", swap in static-block: " << staticBlock << std::endl;
#endif
                    if (userprogram.vertexIO.blockInMem(preStaticBlock)){
                        lostBlockId = preStaticBlock;
#if SWAP_BLOCK_DEBUG
                        std::cout << "pos2: lost-block = pre-static-block = " << lostBlockId << std::endl;
#endif
                    }
                    //                std::cout << "static block swap in: " << staticBlock << ", swap out: " << preStaticBlock <<std::endl;
                }
#else
                if (inMemIndex[staticBlock] < nmblocks) {
                    lostBlockId = preStaticBlock;
                } else {
                    swapBlock(staticBlock, preStaticBlock);
                    if (userprogram.vertexIO.blockInMem(preStaticBlock)){
                        lostBlockId = preStaticBlock;
                    }
                    //                std::cout << "static block swap in: " << staticBlock << ", swap out: " << preStaticBlock <<std::endl;
                }
#endif

#if OUTPUT_ON_DEMAND_DATA
                staticBlock_global = staticBlock;
#endif

                walk_manager->clearRecoredWalkNum(staticBlock);
                //             if(blockcount % (nBlocks/100+1)==1)
                //             if(blockcount % (1024*1024*1024/nedges+1) == 1)
                {
                    logstream(LOG_DEBUG) << runtime() << "s : staticBlockCount: " << blockcount << std::endl;
                    logstream(LOG_INFO) << "currentNWalks = " << walk_manager->currentNWalks << ", nwalks[" << staticBlock
                    << "] = " << nWalks << std::endl;
                }
                walk_manager->collectBucket(staticBlock);
#if OUTPUT_ON_DEMAND_DATA
                Timer executeTimer_staticBlock;
                executeTimer_staticBlock.start();
#endif
                for (dynamicBlock = staticBlock + 1; dynamicBlock < nblocks; dynamicBlock++) {

#if MULTI_BUFFER
                    m.start_time("extend-buffer");
                    walk_manager->exendBucket(dynamicBlock);
                    walk_manager->clearExtendBucket(dynamicBlock);
                    m.stop_time("extend-buffer");
#endif
                    if (!walk_manager->bucket[dynamicBlock].size()) {
                        continue;
                    }
#if SWAP_BLOCK_DEBUG
                    std::cout << "dynamic-block = " << dynamicBlock << std::endl;
#endif
#if OUTPUT_ON_DEMAND_DATA
                    dynamicBlock_global = dynamicBlock;
                    getCSR_dynamicBlock_global_ms = 0;
#endif
                    wid_t totalWalks = walk_manager->bucket[dynamicBlock].size();
                    bool loadDynamicBlock = true;
                    Timer mapAndGetTimer;
                    Timer loadBlockTimer;
#if USE_TRAINED_THS
                    if (ths_dynamic.at(dynamicBlock) > 0 && totalWalks < (blocks[dynamicBlock + 1] - blocks[dynamicBlock]) * ths_dynamic.at(dynamicBlock) && !userprogram.vertexIO.blockInMem(dynamicBlock)){
#elif USE_PRE_DEFINE_TH
                        if (totalWalks < (blocks[dynamicBlock + 1] - blocks[dynamicBlock]) * NO_LOAD_TH && !userprogram.vertexIO.blockInMem(dynamicBlock)){
#endif
                            loadDynamicBlock = false;
                            m.start_time("map-get");
                            mapAndGetTimer.start();
                            walk_manager->mapBucketVertex(dynamicBlock, userprogram.vertexIO);
                            userprogram.vertexIO.getVertexOnDemand(dynamicBlock);
                            mapAndGetTimer.stop();
                            m.stop_time("map-get");
                            logstream(LOG_DEBUG) << "OL-block " << dynamicBlock <<std::endl;
                        }else{
                            loadBlockTimer.start();
                            if (preDynamicBlock == staticBlock) {
                                assert(lostBlockId != INVALID_BID);
                                swapBlock(dynamicBlock, lostBlockId);
#if SWAP_BLOCK_DEBUG
                                std::cout << "pos1: swap out lost-block: " << lostBlockId << ", swap in dynamic-block: " << dynamicBlock << std::endl;
#endif
                                lostBlockId = INVALID_BID;
#if SWAP_BLOCK_DEBUG
                                std::cout << "pos1: lost-block = invalid" << std::endl;
#endif
                            } else {
                                if (lostBlockId != INVALID_BID && !userprogram.vertexIO.blockInMem(preDynamicBlock) && userprogram.vertexIO.blockInMem(lostBlockId)){
                                    swapBlock(dynamicBlock, lostBlockId);
#if SWAP_BLOCK_DEBUG
                                    std::cout << "pos2: swap out lost-block: " << lostBlockId << ", swap in dynamic-block: " << dynamicBlock << std::endl;
#endif
                                    lostBlockId = INVALID_BID;
#if SWAP_BLOCK_DEBUG
                                    std::cout << "pos2: lost-block = invalid" << std::endl;
#endif
                                }else{
                                    swapBlock(dynamicBlock, preDynamicBlock);
#if SWAP_BLOCK_DEBUG
                                    std::cout << "swap out pre-dynamic-block: " << preDynamicBlock << ", swap in dynamic-block: " << dynamicBlock << std::endl;
#endif
                                }
                                //                        std::cout << "dynamic block swap in: " << dynamicBlock << ", swap out: " << preDynamicBlock << std::endl;
                            }
                            userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
                            loadBlockTimer.stop();
                        }
#if OUTPUT_ON_DEMAND_DATA
                        Timer execBucketTimer;
                        execBucketTimer.start();
#endif
                        m.start_time("4_UpdateWalk");
#if MULTI_THREAD
#pragma omp parallel for schedule(dynamic)
#endif
                        for (wid_t w = 0; w < walk_manager->bucket[dynamicBlock].size(); w++) {
                            if (w == 18579507 && staticBlock == 2 && dynamicBlock == 12) {
                                WalkDataType walk = walk_manager->bucket[dynamicBlock][w];
                                vid_t preVId = WalkManager::getPreviousVertex(walk);
                                wid_t preBId = WalkManager::getPreBlock(walk);
                                int a = 1;
                            }
                            userprogram.processWalk(m, walk_manager->bucket[dynamicBlock][w], staticBlock, dynamicBlock,
                                                    *walk_manager, loadDynamicBlock, loadStaticBlock);
                        }
                        if (loadDynamicBlock){
                            preDynamicBlock = dynamicBlock;
#if SWAP_BLOCK_DEBUG
                            std::cout << "pre-dynamic-block = dynamic-block = " << preDynamicBlock << std::endl;
#endif
#if BLOCK_UTI
                            blockUtiFile << blockIONum << ","
                                         << dynamicBlock << ","
                                         << "dynamic-block" << ","
                                         << userprogram.vertexIO.getUsedEdgesNum(dynamicBlock) << ","
                                         << nEdges.at(dynamicBlock) << std::endl;
                            userprogram.vertexIO.clearBlockMap(dynamicBlock);
                            blockIONum++;
#endif
                        }else{
                            userprogram.vertexIO.clearVertexMapAndCSR(dynamicBlock);
#if BLOCK_UTI
                            blockUtiFile << blockIONum << ","
                                         << dynamicBlock << ","
                                         << "dynamic-block" << ","
                                         << userprogram.vertexIO.getUsedEdgesNum(dynamicBlock) << ","
                                         << userprogram.vertexIO.getUsedEdgesNum(dynamicBlock) << std::endl;
                            userprogram.vertexIO.clearBlockMap(dynamicBlock);
                            blockIONum++;
#endif
                        }
                        m.stop_time("4_UpdateWalk");
#if OUTPUT_ON_DEMAND_DATA
                        execBucketTimer.stop();
//#if NO_LOAD_TH
                        if (!loadDynamicBlock){
                            blockLogFile_dynamicBlocks << staticBlock << ","
                            << dynamicBlock << ","
                            << blocks[dynamicBlock + 1] - blocks[dynamicBlock] << ","
                            << nEdges.at(dynamicBlock) << ","
                            << walk_manager->bucket[dynamicBlock].size() << ","
                            << mapAndGetTimer.duration_s() << ","
                            << execBucketTimer.duration_s() << std::endl;
                        }
//#else
//                        blockLogFile_dynamicBlocks << staticBlock << ","
//                        << dynamicBlock << ","
//                        << blocks[dynamicBlock + 1] - blocks[dynamicBlock] << ","
//                        << nEdges.at(dynamicBlock) << ","
//                        << walk_manager->bucket[dynamicBlock].size() << ","
//                        << loadBlockTimer.duration_s() << ","
//                        << execBucketTimer.duration_s() << std::endl;
//#endif
#endif
                    }
#if OUTPUT_ON_DEMAND_DATA
                    executeTimer_staticBlock.stop();
                    if (!loadStaticBlock){
                        blockLogFile_staticBlocks << staticBlock << ","
                                                   << blocks[staticBlock + 1] - blocks[staticBlock] << ","
                                                   << nEdges.at(staticBlock) << ","
                                                   << nWalks << ","
                                                   << loadTimer_staticBlock.duration_s() << ","
                                                   << executeTimer_staticBlock.duration_s() << std::endl;
                    }
#endif
#if INFO | DEBUG
double usage = get_memory();
                    if (usage > maxMemUsage) {
                        maxMemUsage = usage;
                    }
#endif
#if DEBUG
                    if (usage > maxMemUsage_staticBlock.at(staticBlock)){
                        maxMemUsage_staticBlock.at(staticBlock) = usage;
                        //                if (staticBlock == 0){
                        //                    /* two block currently in memory is staticBlock and preDynamicBlock */
                        //                    vid_t nVertsStaticBlock = blocks[staticBlock + 1] - blocks[staticBlock];
                        //                    vid_t nVertsPreDynamicBlock = blocks[preDynamicBlock + 1] - blocks[preDynamicBlock];
                        //                    bid_t staticBlockInMemIndex = 0;
                        //                    bid_t preDynamicBlockInMemIndex = 1;
                        //                    /* FIXME: 这里实际上不是很准确，如果涉及到lostBlock，这里将统计错误，但考虑到子图大小差别不大，这里误差应该较小 */
                        //                    if (inMemIndex[preDynamicBlock] == 0){
                        //                        preDynamicBlockInMemIndex = 0;
                        //                        staticBlockInMemIndex = 1;
                        //                    }
                        //                    eid_t nEdgesStaticBlock = beg_posbuf[staticBlockInMemIndex][nVertsStaticBlock] - beg_posbuf[staticBlockInMemIndex][0];
                        //                    eid_t nEdgesPreDynamicBlock = beg_posbuf[preDynamicBlockInMemIndex][nVertsPreDynamicBlock] - beg_posbuf[preDynamicBlockInMemIndex][0];
                        //                    size_t csrByte = (nEdgesPreDynamicBlock + nEdgesStaticBlock) * sizeof(vid_t);
                        //                    size_t begPosByte = (nVertsStaticBlock + nVertsPreDynamicBlock) * sizeof(vid_t);
                        //                    memoryComponent_b0.graphfile = static_cast<double>(csrByte + begPosByte) / 1024 / 1024;
                        //                    size_t currentWalkByte = walk_manager->curBlockNWalks * sizeof(WalkDataType);
                        //                    memoryComponent_b0.currentWalkArray = static_cast<double>(currentWalkByte) / 1024 / 1024;
                        //                    size_t bucketsByte = 0;
                        //                    for (bid_t b = 0; b < nblocks; b++){
                        //                        /* vec.size() and vec.capacity differ a lot */
                        //                        bucketsByte += walk_manager->preBucket[b].capacity() * sizeof(WalkDataType);
                        //                        bucketsByte += walk_manager->curBucket[b].capacity() * sizeof(WalkDataType);
                        //                    }
                        //                    memoryComponent_b0.bucket = static_cast<double>(bucketsByte) / 1024 / 1024;
                        //                    size_t walkPoolSize_b = 0;
                        //                    for (tid_t t = 0; t < 16; t++){
                        //                        for (bid_t b = 0; b < nblocks; b++){
                        //                            walkPoolSize_b += walk_manager->walkPool[t][b].size_w * sizeof(WalkDataType);
                        //                        }
                        //                    }
                        //                    memoryComponent_b0.walkPool = static_cast<double>(walkPoolSize_b) / 1024 / 1024;
                        //                }
                        /* two block currently in memory is staticBlock and preDynamicBlock */
                        vid_t nVertsStaticBlock = blocks[staticBlock + 1] - blocks[staticBlock];
                        vid_t nVertsPreDynamicBlock = blocks[preDynamicBlock + 1] - blocks[preDynamicBlock];
                        bid_t staticBlockInMemIndex = 0;
                        bid_t preDynamicBlockInMemIndex = 1;
                        /* FIXME: 这里实际上不是很准确，如果涉及到lostBlock，这里将统计错误，但考虑到子图大小差别不大，这里误差应该较小 */
                        if (inMemIndex[preDynamicBlock] == 0){
                            preDynamicBlockInMemIndex = 0;
                            staticBlockInMemIndex = 1;
                        }
                        eid_t nEdgesStaticBlock = beg_posbuf[staticBlockInMemIndex][nVertsStaticBlock] - beg_posbuf[staticBlockInMemIndex][0];
                        eid_t nEdgesPreDynamicBlock = beg_posbuf[preDynamicBlockInMemIndex][nVertsPreDynamicBlock] - beg_posbuf[preDynamicBlockInMemIndex][0];
                        size_t csrByte = (nEdgesPreDynamicBlock + nEdgesStaticBlock) * sizeof(vid_t);
                        size_t begPosByte = (nVertsStaticBlock + nVertsPreDynamicBlock) * sizeof(vid_t);
                        blockMemComp.at(staticBlock).graphfile = static_cast<double>(csrByte + begPosByte) / 1024 / 1024;
                        size_t currentWalkByte = walk_manager->curBlockNWalks * sizeof(WalkDataType);
                        blockMemComp.at(staticBlock).currentWalkArray = static_cast<double>(currentWalkByte) / 1024 / 1024;
                        size_t bucketsByte = 0;
                        for (bid_t b = 0; b < nblocks; b++){
                            /* vec.size() and vec.capacity differ a lot */
                            bucketsByte += walk_manager->preBucket[b].size() * sizeof(WalkDataType);
                            bucketsByte += walk_manager->curBucket[b].size() * sizeof(WalkDataType);
                        }
                        blockMemComp.at(staticBlock).bucket = static_cast<double>(bucketsByte) / 1024 / 1024;
                        size_t walkPoolSize_b = 0;
                        for (tid_t t = 0; t < 16; t++){
                            for (bid_t b = 0; b < nblocks; b++){
                                walkPoolSize_b += walk_manager->walkPool[t][b].size_w * sizeof(WalkDataType);
                            }
                        }
                        blockMemComp.at(staticBlock).walkPool = static_cast<double>(walkPoolSize_b) / 1024 / 1024;
                    }
#endif
#if INFO | DEBUG
                    memUsage = ((memUsage * memCount) + usage) / (memCount + 1);
                    memCount++;
#endif
                    walk_manager->updateWalkNum();
                    walk_manager->clearBucket();
                    if (loadStaticBlock){
                        preStaticBlock = staticBlock;
#if SWAP_BLOCK_DEBUG
                        std::cout << "pre-static-block = static-block = " << preStaticBlock << std::endl;
#endif
#if BLOCK_UTI
                        blockUtiFile << blockIONum << ","
                                     << staticBlock << ","
                                     << "static-block" << ","
                                     << userprogram.vertexIO.getUsedEdgesNum(staticBlock) << ","
                                     << nEdges.at(staticBlock) << std::endl;
                        userprogram.vertexIO.clearBlockMap(staticBlock);
                        blockIONum++;
#endif
#if ACT_VER
                        actVerFile << blockIONum_ver << ","
                                   << staticBlock << ","
                                   << "static-block" << ","
                                   << userprogram.vertexIO.getActivatedVerticesNum(staticBlock) << ","
                                   << userprogram.vertexIO.getStartVertexId(staticBlock + 1) - userprogram.vertexIO.getStartVertexId(staticBlock) << std::endl;
                        userprogram.vertexIO.clearBlockMap(staticBlock);
                        blockIONum_ver++;
#endif

                    }else{
                        userprogram.vertexIO.clearVertexMapAndCSR(staticBlock, false);
#if BLOCK_UTI
                        blockUtiFile << blockIONum << ","
                                     << staticBlock << ","
                                     << "static-block" << ","
                                     << userprogram.vertexIO.getUsedEdgesNum(staticBlock) << ","
                                     << userprogram.vertexIO.getUsedEdgesNum(staticBlock) << std::endl;
                        userprogram.vertexIO.clearBlockMap(staticBlock);
                        blockIONum++;
#endif
                    }
                }
            } // For block loop
            blockLogFile_dynamicBlocks.close();
#if DEBUG
            for (bid_t b = 0; b < nblocks - 1; b++){
                m.set("maxMem_SBlock=" + std::to_string(b), maxMemUsage_staticBlock.at(b));
            }
            std::ofstream memCompLog("log/dual-bucket/memComp.csv");
            for (bid_t b = 0; b < nblocks - 1; b++){
                memCompLog << blockMemComp.at(b).graphfile << "," << blockMemComp.at(b).bucket << ","
                << blockMemComp.at(b).currentWalkArray << "," << blockMemComp.at(b).walkPool
                << std::endl;
            }
            uint64_t totalSteps = 0;
            for (tid_t t = 0; t < 16; t++){
                totalSteps += userprogram.totalSteps_t.at(t);
            }
            m.set("total-steps", totalSteps);
#endif
            m.set("avg-mem-usage", memUsage);
            m.set("max-mem-usage", maxMemUsage);
            m.stop_time("1_ProcessTime");
}
#endif

        void DualBucketEngine::run(vid_t nRootVertices, SecondOrderRW &userprogram, float prob, float th) {
#if FULLY_LOAD
            run_fullyLoadTest(nRootVertices, userprogram, prob, th);
#elif NO_LOAD_TEST
            run_noLoadTest(nRootVertices, userprogram, prob, th);
#elif ONDEMAND_LOAD
            run_onDemandLoad(nRootVertices, userprogram, prob, th);
#else
            // srand((unsigned)time(NULL));
            m.start_time("2_StartWalks");
            bid_t preInitBlock = INVALID_BID;
            for (bid_t b = 0; b < nblocks; b++) {
                swapBlock(b, preInitBlock);
                userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
                userprogram.initWalks(b, *walk_manager, m);
#if INFO
                logstream(LOG_INFO) << "initiated block" << b << std::endl;
#endif
                preInitBlock = b;
            }
            walk_manager->updateWalkNum();
            std::cout << "currentNWalks before run: " << walk_manager->currentNWalks << std::endl;
            m.stop_time("2_StartWalks");

            gettimeofday(&start, NULL);
            m.start_time("1_ProcessTime");

            vid_t nverts;
            eid_t nedges;
            int blockcount = 0;
            bid_t staticBlock, dynamicBlock;
            bid_t preStaticBlock = INVALID_BID, preDynamicBlock = INVALID_BID;
            bid_t lostBlockId = INVALID_BID;
            unLoadAllBlock();
            if (nblocks == 1) {
                //        eid_t *beg_pos;
                //        vid_t *csr;
                //        loadBlock(0, beg_pos, csr);
                //        wid_t nWalks = walk_manager->getCurrentWalks(staticBlock);
                //        {
                //            logstream(LOG_DEBUG) << runtime() << "s : staticBlockCount: " << blockcount << std::endl;
                //            logstream(LOG_INFO) << "currentNWalks = " << walk_manager->currentNWalks << ", nwalks[" << staticBlock
                //                                << "] = " << nWalks << std::endl;
                //        }
                //#pragma omp parallel for schedule(dynamic)
                //        for (wid_t w = 0; w < nWalks; w++) {
                //            userprogram.updateWalk(m, walk_manager->currentWalks[w], 0, beg_pos, csr, *walk_manager);
                //        }
                return;
            }
            while (!userprogram.allWalksFinished(*walk_manager)) {
                for (staticBlock = 0; staticBlock < nblocks - 1; staticBlock++) {
                    blockcount++;
                    if (!walk_manager->nWalks[staticBlock]) {
                        continue;
                    }
                    if (inMemIndex[staticBlock] < nmblocks) {
                        lostBlockId = preStaticBlock;
                    } else {
                        swapBlock(staticBlock, preStaticBlock);
                    }
                    wid_t nWalks = walk_manager->getCurrentWalks(staticBlock);
                    walk_manager->clearRecoredWalkNum(staticBlock);
                    //             if(blockcount % (nBlocks/100+1)==1)
                    //             if(blockcount % (1024*1024*1024/nedges+1) == 1)
                    {
                        logstream(LOG_DEBUG) << runtime() << "s : staticBlockCount: " << blockcount << std::endl;
                        logstream(LOG_INFO) << "currentNWalks = " << walk_manager->currentNWalks << ", nwalks[" << staticBlock
                        << "] = " << nWalks << std::endl;
                    }
                    walk_manager->collectBucket(staticBlock);
                    for (dynamicBlock = staticBlock + 1; dynamicBlock < nblocks; dynamicBlock++) {
#if MULTI_BUFFER
                        m.start_time("extend-buffer");
                        walk_manager->exendBucket(dynamicBlock);
                        walk_manager->clearExtendBucket(dynamicBlock);
                        m.stop_time("extend-buffer");
#endif
                        if (!walk_manager->bucket[dynamicBlock].size()) {
                            continue;
                        }
                        if (preDynamicBlock == staticBlock) {
                            assert(lostBlockId != INVALID_BID);
                            swapBlock(dynamicBlock, lostBlockId);
                            lostBlockId = INVALID_BID;
                        } else {
                            swapBlock(dynamicBlock, preDynamicBlock);
                        }
                        userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
                        m.start_time("4_UpdateWalk");

#pragma omp parallel for schedule(dynamic)
                        for (wid_t w = 0; w < walk_manager->bucket[dynamicBlock].size(); w++) {
                            userprogram.processWalk(m, walk_manager->bucket[dynamicBlock][w], staticBlock, dynamicBlock,
                                                    *walk_manager, true, true);
                        }
                        preDynamicBlock = dynamicBlock;
                        m.stop_time("4_UpdateWalk");
                    }
#if INFO | DEBUG
                    double usage = get_memory();
                    if (usage > maxMemUsage) {
                        maxMemUsage = usage;
                    }
#endif
#if DEBUG
                    if (usage > maxMemUsage_staticBlock.at(staticBlock)){
                        maxMemUsage_staticBlock.at(staticBlock) = usage;
                        //                if (staticBlock == 0){
                        //                    /* two block currently in memory is staticBlock and preDynamicBlock */
                        //                    vid_t nVertsStaticBlock = blocks[staticBlock + 1] - blocks[staticBlock];
                        //                    vid_t nVertsPreDynamicBlock = blocks[preDynamicBlock + 1] - blocks[preDynamicBlock];
                        //                    bid_t staticBlockInMemIndex = 0;
                        //                    bid_t preDynamicBlockInMemIndex = 1;
                        //                    /* FIXME: 这里实际上不是很准确，如果涉及到lostBlock，这里将统计错误，但考虑到子图大小差别不大，这里误差应该较小 */
                        //                    if (inMemIndex[preDynamicBlock] == 0){
                        //                        preDynamicBlockInMemIndex = 0;
                        //                        staticBlockInMemIndex = 1;
                        //                    }
                        //                    eid_t nEdgesStaticBlock = beg_posbuf[staticBlockInMemIndex][nVertsStaticBlock] - beg_posbuf[staticBlockInMemIndex][0];
                        //                    eid_t nEdgesPreDynamicBlock = beg_posbuf[preDynamicBlockInMemIndex][nVertsPreDynamicBlock] - beg_posbuf[preDynamicBlockInMemIndex][0];
                        //                    size_t csrByte = (nEdgesPreDynamicBlock + nEdgesStaticBlock) * sizeof(vid_t);
                        //                    size_t begPosByte = (nVertsStaticBlock + nVertsPreDynamicBlock) * sizeof(vid_t);
                        //                    memoryComponent_b0.graphfile = static_cast<double>(csrByte + begPosByte) / 1024 / 1024;
                        //                    size_t currentWalkByte = walk_manager->curBlockNWalks * sizeof(WalkDataType);
                        //                    memoryComponent_b0.currentWalkArray = static_cast<double>(currentWalkByte) / 1024 / 1024;
                        //                    size_t bucketsByte = 0;
                        //                    for (bid_t b = 0; b < nblocks; b++){
                        //                        /* vec.size() and vec.capacity differ a lot */
                        //                        bucketsByte += walk_manager->preBucket[b].capacity() * sizeof(WalkDataType);
                        //                        bucketsByte += walk_manager->curBucket[b].capacity() * sizeof(WalkDataType);
                        //                    }
                        //                    memoryComponent_b0.bucket = static_cast<double>(bucketsByte) / 1024 / 1024;
                        //                    size_t walkPoolSize_b = 0;
                        //                    for (tid_t t = 0; t < 16; t++){
                        //                        for (bid_t b = 0; b < nblocks; b++){
                        //                            walkPoolSize_b += walk_manager->walkPool[t][b].size_w * sizeof(WalkDataType);
                        //                        }
                        //                    }
                        //                    memoryComponent_b0.walkPool = static_cast<double>(walkPoolSize_b) / 1024 / 1024;
                        //                }
                        /* two block currently in memory is staticBlock and preDynamicBlock */
                        vid_t nVertsStaticBlock = blocks[staticBlock + 1] - blocks[staticBlock];
                        vid_t nVertsPreDynamicBlock = blocks[preDynamicBlock + 1] - blocks[preDynamicBlock];
                        bid_t staticBlockInMemIndex = 0;
                        bid_t preDynamicBlockInMemIndex = 1;
                        /* FIXME: 这里实际上不是很准确，如果涉及到lostBlock，这里将统计错误，但考虑到子图大小差别不大，这里误差应该较小 */
                        if (inMemIndex[preDynamicBlock] == 0){
                            preDynamicBlockInMemIndex = 0;
                            staticBlockInMemIndex = 1;
                        }
                        eid_t nEdgesStaticBlock = beg_posbuf[staticBlockInMemIndex][nVertsStaticBlock] - beg_posbuf[staticBlockInMemIndex][0];
                        eid_t nEdgesPreDynamicBlock = beg_posbuf[preDynamicBlockInMemIndex][nVertsPreDynamicBlock] - beg_posbuf[preDynamicBlockInMemIndex][0];
                        size_t csrByte = (nEdgesPreDynamicBlock + nEdgesStaticBlock) * sizeof(vid_t);
                        size_t begPosByte = (nVertsStaticBlock + nVertsPreDynamicBlock) * sizeof(vid_t);
                        blockMemComp.at(staticBlock).graphfile = static_cast<double>(csrByte + begPosByte) / 1024 / 1024;
                        size_t currentWalkByte = walk_manager->curBlockNWalks * sizeof(WalkDataType);
                        blockMemComp.at(staticBlock).currentWalkArray = static_cast<double>(currentWalkByte) / 1024 / 1024;
                        size_t bucketsByte = 0;
                        for (bid_t b = 0; b < nblocks; b++){
                            /* vec.size() and vec.capacity differ a lot */
                            bucketsByte += walk_manager->preBucket[b].size() * sizeof(WalkDataType);
                            bucketsByte += walk_manager->curBucket[b].size() * sizeof(WalkDataType);
                        }
                        blockMemComp.at(staticBlock).bucket = static_cast<double>(bucketsByte) / 1024 / 1024;
                        size_t walkPoolSize_b = 0;
                        for (tid_t t = 0; t < 16; t++){
                            for (bid_t b = 0; b < nblocks; b++){
                                walkPoolSize_b += walk_manager->walkPool[t][b].size_w * sizeof(WalkDataType);
                            }
                        }
                        blockMemComp.at(staticBlock).walkPool = static_cast<double>(walkPoolSize_b) / 1024 / 1024;
                    }
#endif
#if INFO | DEBUG
                    memUsage = ((memUsage * memCount) + usage) / (memCount + 1);
                    memCount++;
#endif
                    walk_manager->updateWalkNum();
                    walk_manager->clearBucket();
                    preStaticBlock = staticBlock;
                }
            } // For block loop
#if DEBUG
for (bid_t b = 0; b < nblocks - 1; b++){
                m.set("maxMem_SBlock=" + std::to_string(b), maxMemUsage_staticBlock.at(b));
            }
            std::ofstream memCompLog("log/dual-bucket/memComp.csv");
            for (bid_t b = 0; b < nblocks - 1; b++){
                memCompLog << blockMemComp.at(b).graphfile << "," << blockMemComp.at(b).bucket << ","
                << blockMemComp.at(b).currentWalkArray << "," << blockMemComp.at(b).walkPool
                << std::endl;
            }
            uint64_t totalSteps = 0;
            for (tid_t t = 0; t < 16; t++){
                totalSteps += userprogram.totalSteps_t.at(t);
            }
            m.set("total-steps", totalSteps);
#endif
            m.set("avg-mem-usage", memUsage);
            m.set("max-mem-usage", maxMemUsage);
            m.stop_time("1_ProcessTime");
#endif
    }

void DualBucketEngine::print_config() {
    logstream(LOG_INFO) << "Engine configuration: " << std::endl;
    logstream(LOG_INFO) << " exec_threads = " << (int) exec_threads << std::endl;
    logstream(LOG_INFO) << " blocksize_kb = " << blocksize_kb << "kb" << std::endl;
    logstream(LOG_INFO) << " number of total blocks = " << nblocks << std::endl;
    logstream(LOG_INFO) << " number of in-memory blocks = " << nmblocks << std::endl;
}

double DualBucketEngine::runtime() {
    timeval end;
    gettimeofday(&end, NULL);
    return end.tv_sec - start.tv_sec + ((double) (end.tv_usec - start.tv_usec)) / 1.0E6;
}

#endif
#endif //IOE_SORW_DUALBUCKET_HPP
