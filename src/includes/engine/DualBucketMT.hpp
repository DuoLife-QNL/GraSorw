/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/4/29 17:15
 * Multi-thread disk based version of collecting bucket
 */

#ifndef IOE_SORW_DUALBUCKETMT_HPP
#define IOE_SORW_DUALBUCKETMT_HPP

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

#include "api/filename.hpp"
#include "api/io.hpp"
#include "metrics/metrics.hpp"
#include "api/pthread_tools.hpp"
#include "walks/SecondOrderRW.hpp"
#include "Settings.hpp"
#include "logger/logger.hpp"
#include "walks/RandomWalk.hpp"

#if BI_BLOCK
#if 0
class DualBucketEngine_MT {
public:
    std::string base_filename;
    // unsigned membudget_mb;
    unsigned long long blocksize_kb;
    bid_t nblocks;
    vid_t nvertices;
    tid_t exec_threads;
    vid_t* blocks;
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
    DualBucketEngine_MT(const std::string& _base_filename, unsigned long long int _blocksize_kb,
                     bid_t _nblocks, bid_t _nmblocks, metrics &_m);

    virtual ~DualBucketEngine_MT();

    void load_block_range(std::string base_filename, unsigned long long blocksize_kb, vid_t * &blocks, bool allowfail=false);

    void loadSubGraph(bid_t p, eid_t * &beg_pos, vid_t * &csr, vid_t *nverts, eid_t *nedges);

    void loadSubGraph(bid_t p, vid_t *nverts, eid_t *nedges, bid_t doNotUnloadBlockId){
        if(inMemIndex[p] == nmblocks){//the block is not in memory
            bid_t swapin;
            if(cmblocks < nmblocks){
                swapin = cmblocks++;
            }else{
                bid_t minmwb = chooseSwapOutExcept(doNotUnloadBlockId);
                swapin = unLoadBlock_id(minmwb);
                // munmap(beg_posbuf[swapin], sizeof(eid_t)*(blocks[minmwb+1] - blocks[minmwb] + 1));
            }
            loadSubGraph(p, beg_posbuf[swapin], csrbuf[swapin], nverts, nedges);
            inMemIndex[p] = swapin;
        }
    }

    void loadBlock(bid_t p, eid_t* &beg_pos, vid_t* &csr){

        // m.start_time("__g_loadSubGraph_malloc_begpos");
        /* read beg_pos file */

        vid_t nverts = blocks[p+1] - blocks[p];
        beg_pos = (eid_t*) malloc((nverts+1)*sizeof(eid_t));
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
        preada(beg_posf, beg_pos, (size_t)(nverts+1)*sizeof(eid_t), (size_t)blocks[p]*sizeof(eid_t));
//        m.stop_time("z__g_loadSubGraph_read_begpos");
        /* read csr file */
        m.start_time("z__g_loadSubGraph_realloc_csr");
        eid_t nedges = beg_pos[nverts] - beg_pos[0];
        if(nedges*sizeof(vid_t) > blocksize_kb*1024){
            csr = (vid_t*)realloc(csr, (nedges)*sizeof(vid_t) );
        }
        m.stop_time("z__g_loadSubGraph_realloc_csr");
        m.start_time("z__g_loadSubGraph_read_csr");
        preada(csrf, csr, (nedges)*sizeof(vid_t), beg_pos[0]*sizeof(vid_t));
        m.stop_time("z__g_loadSubGraph_read_csr");

        m.stop_time("g_loadSubGraph");
    }

    void swapBlock(bid_t swapInBlockId, bid_t swapOutBlockId) {
        m.start_time("5_LoadBlocks");
        bid_t swapOutBlockIndex = inMemIndex[swapOutBlockId];
        if(cmblocks < nmblocks) {
            swapOutBlockIndex = cmblocks++;
        }else{
            inMemIndex[swapOutBlockId] = nmblocks;
            unLoadBlock_index(swapOutBlockIndex);
        }
        loadBlock(swapInBlockId, beg_posbuf[swapOutBlockIndex], csrbuf[swapOutBlockIndex]);
        inMemIndex[swapInBlockId] = swapOutBlockIndex;
        m.stop_time("5_LoadBlocks");
    }

    void unLoadAllBlock(){
        for (bid_t b = 0; b < nblocks; b++){
            inMemIndex[b] = nmblocks;
        }
        for (bid_t b = 0; b < nmblocks; b++){
            unLoadBlock_index(b);
        }
        cmblocks = 0;
    }

    bid_t unLoadBlock_id(bid_t p){
        bid_t unLoadIndex = inMemIndex[p];
        inMemIndex[p] = nmblocks;
        assert(unLoadIndex < nmblocks);
        if(beg_posbuf[unLoadIndex] != nullptr) free(beg_posbuf[unLoadIndex]);
        return unLoadIndex;
    }

    void unLoadBlock_index(bid_t blockIndex){
        assert(blockIndex < nmblocks);
        if(beg_posbuf[blockIndex] != nullptr){
            free(beg_posbuf[blockIndex]);
        }
    }

    void loadPreVertexSubGraph(bid_t p, eid_t * &beg_pos, vid_t * &csr, vid_t *nverts, eid_t *nedges);

    void findSubGraph(bid_t p, eid_t * &beg_pos, vid_t * &csr, vid_t *nverts, eid_t *nedges);

    bid_t swapOut();

    bid_t chooseSwapOutExcept(bid_t doNotUnloadBlockId) {
        m.start_time("z_g_swapOut");
        wid_t minmw = 0xffffffff;
        bid_t minmwb = 0;
        for(bid_t b = 0; b < nblocks; b++){
            if (b == doNotUnloadBlockId){
                continue;
            }
            if(inMemIndex[b]<nmblocks && walk_manager->nWalks[b] < minmw){
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
    void exec_updates(SecondOrderRW &userprogram, wid_t nwalks, eid_t *&beg_pos, vid_t *&csr,
                      const vid_t &nBlockVertices, float th);

    void processWalk(SecondOrderRW &userprogram, WalkDataType &walk, bid_t staticBlock, bid_t dynamicBlock){
        bool walkFinished = false;

        userprogram.updateWalk(m, walk, staticBlock, dynamicBlock, *walk_manager, walkFinished);
        if (walkFinished){
            return;
        }
        tid_t threadId = omp_get_thread_num();
        bid_t cur = WalkManager::getCurBlock(walk);
        bid_t pre = WalkManager::getPreBlock(walk);
        if (cur < staticBlock){
            walk_manager->moveWalk(walk, threadId, cur);
        }else if(cur > staticBlock && cur < dynamicBlock){
            if (pre == staticBlock){
                walk_manager->moveWalk(walk, threadId, staticBlock);
            }else{
                walk_manager->moveWalk(walk, threadId, cur);
            }
        }else{
            assert(cur > dynamicBlock);
            if (pre == staticBlock){
#pragma omp critical
                {
                    walk_manager->addWalk2BucketFile(walk, cur, CUR);
                    walk_manager->curBucketsSize[cur]++;
                }
            }else{
                walk_manager->moveWalk(walk, threadId, dynamicBlock);
            }
        }
    }

    void run(vid_t nRootVertices, SecondOrderRW &userprogram, float prob, float th);

    void checkMemoryUsage(){
#if INFO | DEBUG
        double usage = get_memory();
        if (usage > maxMemUsage){
            maxMemUsage = usage;
        }
        memUsage = ((memUsage * memCount) + usage) / (memCount + 1);
        memCount ++;
#endif
    }
};

DualBucketEngine_MT::DualBucketEngine_MT(const std::string& _base_filename, unsigned long long int _blocksize_kb,
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
    if (nmblocks < 2){
        logstream(LOG_ERROR) << "In memory block number should at least be 2" << std::endl;
        abort();
    }
    csrbuf = (vid_t**)malloc(nmblocks*sizeof(vid_t*));
    for(bid_t b = 0; b < nmblocks; b++){
        csrbuf[b] = (vid_t*)malloc(blocksize_kb*1024);
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
    beg_posbuf = (eid_t**)malloc(nmblocks*sizeof(eid_t*));
    for (bid_t b = 0; b < nmblocks; b++) beg_posbuf[b] = nullptr;
    inMemIndex = (bid_t*)malloc(nblocks*sizeof(bid_t));
    for(bid_t b = 0; b < nblocks; b++)  inMemIndex[b] = nmblocks;
    cmblocks = 0;

    // m.start_time("__g_loadSubGraph_filename");
    std::string invlname = fidname( base_filename, 0 ); //only 1 file
    std::string beg_posname = invlname + ".beg_pos";
    std::string csrname = invlname + ".csr";
    // m.stop_time("__g_loadSubGraph_filename");
    // m.start_time("__g_loadSubGraph_open_begpos");
    beg_posf = open(beg_posname.c_str(),O_RDONLY | O_CREAT, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
    // m.stop_time("__g_loadSubGraph_open_begpos");
    // m.start_time("__g_loadSubGraph_open_csr");
    csrf = open(csrname.c_str(),O_RDONLY | O_CREAT, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
    // m.stop_time("__g_loadSubGraph_open_csr");
    // m.start_time("__g_loadSubGraph_if_open_success");
    if (csrf < 0 || beg_posf < 0) {
        logstream(LOG_FATAL) << "Could not load :" << csrname << " or " << beg_posname << ", error: " << strerror(errno) << std::endl;
    }
    assert(csrf > 0 && beg_posf > 0);
    // m.stop_time("__g_loadSubGraph_if_open_success");

    _m.set("file", _base_filename);
    _m.set("engine", "default");
    _m.set("nBlocks", (size_t)nblocks);

    print_config();
#if DEBUG
    maxMemUsage_staticBlock.assign(_nblocks, 0);
    MemoryComponent memoryComponent;
    blockMemComp.assign(_nblocks, memoryComponent);
#endif
}

DualBucketEngine_MT::~DualBucketEngine_MT() {
    delete walk_manager;

    if(inMemIndex != NULL) free(inMemIndex);
    if(blocks != NULL) free(blocks);

    for(bid_t b = 0; b < cmblocks; b++){
        if(beg_posbuf[b] != NULL)   free(beg_posbuf[b]);
        if(csrbuf[b] != NULL)   free(csrbuf[b]);
        // munmap(csrbuf[b], blocksize_kb*1024);
    }
    if(beg_posbuf != NULL) free(beg_posbuf);
    if(csrbuf != NULL) free(csrbuf);

    close(beg_posf);
    close(csrf);
}

void
DualBucketEngine_MT::load_block_range(std::string base_filename, unsigned long long int blocksize_kb, vid_t *&blocks,
                                   bool allowfail) {
    std::string blockrangefile = blockrangename(base_filename, blocksize_kb);
    std::ifstream brf(blockrangefile.c_str());

    if (!brf.good()) {
        logstream(LOG_ERROR) << "Could not load block range file: " << blockrangefile << std::endl;
    }
    assert(brf.good());

    /* block中存的实际上是vertex id */
    blocks = (vid_t*)malloc((nblocks+1)*sizeof(vid_t));
    vid_t en;
    for(bid_t i=0; i < nblocks+1; i++) {
        assert(!brf.eof());
        brf >> en;
        blocks[i] = en;
    }
    for(bid_t i=nblocks-1; i < nblocks; i++) {
        logstream(LOG_INFO) << "last shard: " << blocks[i] << " - " << blocks[i+1] << std::endl;
    }
    brf.close();
}

void DualBucketEngine_MT::loadSubGraph(bid_t p, eid_t *&beg_pos, vid_t *&csr, vid_t *nverts, eid_t *nedges) {
    m.start_time("g_loadSubGraph");


    // m.start_time("__g_loadSubGraph_malloc_begpos");
    /* read beg_pos file */
    *nverts = blocks[p+1] - blocks[p];
    beg_pos = (eid_t*) malloc((*nverts+1)*sizeof(eid_t));
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
    preada(beg_posf, beg_pos, (size_t)(*nverts+1)*sizeof(eid_t), (size_t)blocks[p]*sizeof(eid_t));
//    m.stop_time("z__g_loadSubGraph_read_begpos");
    /* read csr file */
    m.start_time("z__g_loadSubGraph_realloc_csr");
    *nedges = beg_pos[*nverts] - beg_pos[0];
    if(*nedges*sizeof(vid_t) > blocksize_kb*1024){
        csr = (vid_t*)realloc(csr, (*nedges)*sizeof(vid_t) );
    }
    m.stop_time("z__g_loadSubGraph_realloc_csr");
    m.start_time("z__g_loadSubGraph_read_csr");
    preada(csrf, csr, (*nedges)*sizeof(vid_t), beg_pos[0]*sizeof(vid_t));
    m.stop_time("z__g_loadSubGraph_read_csr");

    m.stop_time("g_loadSubGraph");
}

void DualBucketEngine_MT::loadPreVertexSubGraph(bid_t p, eid_t *&beg_pos, vid_t *&csr, vid_t *nverts, eid_t *nedges) {
    m.start_time("loadPreSubGraph");


    // m.start_time("__g_loadSubGraph_malloc_begpos");
    /* read beg_pos file */
    *nverts = blocks[p+1] - blocks[p];
    beg_pos = (eid_t*) malloc((*nverts+1)*sizeof(eid_t));
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
    preada(beg_posf, beg_pos, (size_t)(*nverts+1)*sizeof(eid_t), (size_t)blocks[p]*sizeof(eid_t));
    m.stop_time("loadPreSubGraph_read_begpos");
    /* read csr file */
    *nedges = beg_pos[*nverts] - beg_pos[0];
    csr = (vid_t*)malloc((*nedges)*sizeof(vid_t));
    m.start_time("loadPreSubGraph_read_csr");
    preada(csrf, csr, (*nedges)*sizeof(vid_t), beg_pos[0]*sizeof(vid_t));
    m.stop_time("loadPreSubGraph_read_csr");

    m.stop_time("loadPreSubGraph");
}

void DualBucketEngine_MT::findSubGraph(bid_t p, eid_t *&beg_pos, vid_t *&csr, vid_t *nverts, eid_t *nedges) {
    m.start_time("2_findSubGraph");
    if(inMemIndex[p] == nmblocks){//the block is not in memory
        bid_t swapin;
        if(cmblocks < nmblocks){
            swapin = cmblocks++;
        }else{
            bid_t minmwb = swapOut();
            swapin = inMemIndex[minmwb];
            inMemIndex[minmwb] = nmblocks;
            assert(swapin < nmblocks);
            if(beg_posbuf[swapin] != NULL) free(beg_posbuf[swapin]);
            // munmap(beg_posbuf[swapin], sizeof(eid_t)*(blocks[minmwb+1] - blocks[minmwb] + 1));
        }
        loadSubGraph(p, beg_posbuf[swapin], csrbuf[swapin], nverts, nedges);
        inMemIndex[p] = swapin;
    }else{
    }
    beg_pos = beg_posbuf[ inMemIndex[p] ];
    csr = csrbuf[ inMemIndex[p] ];
    m.stop_time("2_findSubGraph");
}


bid_t DualBucketEngine_MT::swapOut() {
    m.start_time("z_g_swapOut");
    wid_t minmw = 0xffffffff;
    bid_t minmwb = 0;
    for(bid_t b = 0; b < nblocks; b++){
        if(inMemIndex[b]<nmblocks && walk_manager->nWalks[b] < minmw){
            minmw = walk_manager->nWalks[b];
            minmwb = b;
        }
    }
    m.start_time("z_g_swapOut");
    return minmwb;
}

void DualBucketEngine_MT::exec_updates(SecondOrderRW &userprogram, wid_t nwalks, eid_t *&beg_pos, vid_t *&csr,
                                    const vid_t &nBlockVertices, float th) { //, VertexDataType* vertex_value){
    // unsigned count = walk_manager->readblockWalks(exec_block);
    if(nwalks < 100) omp_set_num_threads(1);
    userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
    userprogram.vertexIO.setPreBlockInfo(nullptr, nullptr, INVALID_BID);
    if (nwalks > th * nBlockVertices){
        m.add("use-pre-bucket-sort", 1);
        for (bid_t p = 0; p < walk_manager->nBlocks; p++){
            walk_manager->currentWalkBucket[p].clear();
        }
//        m.start_time("3_CollectBucket");
        for (wid_t w = 0; w < nwalks; w++){
            WalkDataType walk = walk_manager->currentWalks[w];
            bid_t preResideBlock = walk_manager->getPreBlock(walk);
            walk_manager->currentWalkBucket[preResideBlock].push_back(walk);
        }
//        m.stop_time("3_CollectBucket");
        for (bid_t p = 0; p < walk_manager->nBlocks; p++){
            if (true || walk_manager->currentWalkBucket[p].size() < 0.45 * (blocks[p + 1] - blocks[p])){
                m.start_time("exec_updates_2od_no-load");
#pragma omp parallel for schedule(static)
                for(wid_t i = 0; i < walk_manager->currentWalkBucket[p].size(); i++ ){
                    userprogram.updateWalk(m, walk_manager->currentWalkBucket[p][i], exec_block, beg_pos, csr, *walk_manager);//, vertex_value);
                }
                m.stop_time("exec_updates_2od_no-load");
                continue;
            }
            bool deletePreVertexBlock = false;
            vid_t nverts = 0;
            vid_t *preCsr = nullptr;
            eid_t nedges = 0;
            eid_t *preBegPos = nullptr;
            if (userprogram.vertexIO.blockInMem(p)){ // preVertex Block is currently in memory
                deletePreVertexBlock = false;
                userprogram.vertexIO.getPreBlockInfoFromMem(p);
            }else{
                deletePreVertexBlock = true;
                m.add("load-pre-block", 1);
                loadPreVertexSubGraph(p, preBegPos, preCsr, &nverts, &nedges);
                userprogram.vertexIO.setPreBlockInfo(preCsr, preBegPos, p);
            }
            m.start_time("exec_updates_2od_fully-load");
#pragma omp parallel for schedule(dynamic)
            for (wid_t i = 0; i < walk_manager->currentWalkBucket[p].size(); i++){
                userprogram.updateWalk(m, walk_manager->currentWalkBucket[p][i], exec_block, beg_pos, csr, *walk_manager);
            }

            if (deletePreVertexBlock){
                free(preBegPos);
                free(preCsr);
            }
            m.stop_time("exec_updates_2od_fully-load");
        }
    }else{
        m.start_time("5_exec_updates");
#pragma omp parallel for schedule(static)
        for(wid_t i = 0; i < nwalks; i++ ){
            WalkDataType walk = walk_manager->currentWalks[i];
            userprogram.updateWalk(m, walk, exec_block, beg_pos, csr, *walk_manager);//, vertex_value);
        }
        m.stop_time("5_exec_updates");
    }


    // walk_manager->writeblockWalks(exec_block);
}

void DualBucketEngine_MT::run(vid_t nRootVertices, SecondOrderRW &userprogram, float prob, float th) {
    // srand((unsigned)time(NULL));
    m.start_time("2_StartWalks");
    bid_t preInitBlock = 0;
    for (bid_t b = 0; b < nblocks; b++){
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
    while(!userprogram.allWalksFinished(*walk_manager) ){
        for (staticBlock = 0; staticBlock < nblocks - 1; staticBlock ++){
            blockcount++;
            if (!walk_manager->nWalks[staticBlock]){
                continue;
            }
            if (inMemIndex[staticBlock] < nmblocks){
                lostBlockId = preStaticBlock;
            }else{
                swapBlock(staticBlock, preStaticBlock);
            }
            wid_t nWalks = walk_manager->nWalks[staticBlock];
//             if(blockcount % (nBlocks/100+1)==1)
//             if(blockcount % (1024*1024*1024/nedges+1) == 1)
            {
                logstream(LOG_DEBUG) << runtime() << "s : staticBlockCount: " << blockcount << std::endl;
                logstream(LOG_INFO) << "currentNWalks = " << walk_manager->currentNWalks << ", nwalks[" << staticBlock << "] = " << nWalks << std::endl;
            }
            walk_manager->collectBuckets2Disk(staticBlock);
            walk_manager->clearRecoredWalkNum(staticBlock);
            for (dynamicBlock = staticBlock + 1; dynamicBlock < nblocks; dynamicBlock++){
                wid_t preBucketSize = walk_manager->getBucketSize(dynamicBlock, PRE);
                wid_t curBucketSize = walk_manager->getBucketSize(dynamicBlock, CUR);
                if (!preBucketSize && !curBucketSize){
                    continue;
                }
                if (preDynamicBlock == staticBlock){
                    assert(lostBlockId != INVALID_BID);
                    swapBlock(dynamicBlock, lostBlockId);
                    lostBlockId = INVALID_BID;
                }else{
                    swapBlock(dynamicBlock, preDynamicBlock);
                }
                userprogram.setIOInfo(nmblocks, csrbuf, beg_posbuf, inMemIndex);
                m.start_time("4_UpdateWalk");

                auto preBucketWalks = new WalkDataType [preBucketSize];
                walk_manager->readBucket(preBucketWalks, dynamicBlock, PRE);
//#pragma omp parallel for schedule(dynamic)
                for (wid_t w = 0; w < preBucketSize; w++){
                    processWalk(userprogram, preBucketWalks[w], staticBlock, dynamicBlock);
                }
                checkMemoryUsage();
                delete[] preBucketWalks;
                auto curBucketWalks = new WalkDataType [curBucketSize];
                curBucketSize = walk_manager->getBucketSize(dynamicBlock, CUR);
                walk_manager->readBucket(curBucketWalks, dynamicBlock, CUR);
//#pragma omp parallel for schedule(dynamic)
                for (wid_t w = 0; w < curBucketSize; w++){
                    processWalk(userprogram, curBucketWalks[w], staticBlock, dynamicBlock);
                }
                checkMemoryUsage();
                delete[] curBucketWalks;
                preDynamicBlock = dynamicBlock;
                m.stop_time("4_UpdateWalk");
            }
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

void DualBucketEngine_MT::print_config() {
    logstream(LOG_INFO) << "Engine configuration: " << std::endl;
    logstream(LOG_INFO) << " exec_threads = " << (int)exec_threads << std::endl;
    logstream(LOG_INFO) << " blocksize_kb = " << blocksize_kb << "kb" << std::endl;
    logstream(LOG_INFO) << " number of total blocks = " << nblocks << std::endl;
    logstream(LOG_INFO) << " number of in-memory blocks = " << nmblocks << std::endl;
}

double DualBucketEngine_MT::runtime() {
    timeval end;
    gettimeofday(&end, NULL);
    return end.tv_sec-start.tv_sec+ ((double)(end.tv_usec-start.tv_usec))/1.0E6;
}
#endif
#endif
#endif //IOE_SORW_DUALBUCKETMT_HPP
