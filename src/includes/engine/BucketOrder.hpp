/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/4/23 14:09
 */

#ifndef IOE_SORW_BUCKETORDER_HPP
#define IOE_SORW_BUCKETORDER_HPP
#if PLAIN_BUCKET
class BucketEngine {
public:
    std::string base_filename;
    // unsigned membudget_mb;
    unsigned long long blocksize_kb;
    bid_t nblocks;
    vid_t nvertices;
    tid_t exec_threads;
    vid_t* blocks{};
    timeval start{};

    /* Ｉn memory blocks */
    bid_t nmblocks; //number of in memory blocks
    vid_t **csrbuf;
    eid_t **beg_posbuf;
    bid_t cmblocks; //current number of in memory blocks
    bid_t *inMemIndex;
    int beg_posf, csrf;

    /* State */
    bid_t exec_block{};

    /* Metrics */
    metrics &m;
    double memUsage{};
    long long memCount = 0;

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
    BucketEngine(std::string _base_filename, unsigned long long int _blocksize_kb,
                      bid_t _nblocks, bid_t _nmblocks, metrics &_m);

    virtual ~BucketEngine();

    void load_block_range(std::string base_filename, unsigned long long blocksize_kb, vid_t * &blocks, bool allowfail=false);

    void loadSubGraph(bid_t p, eid_t * &beg_pos, vid_t * &csr, vid_t *nverts, eid_t *nedges);

    void loadPreVertexSubGraph(bid_t p, eid_t * &beg_pos, vid_t * &csr, vid_t *nverts, eid_t *nedges);

    void findSubGraph(bid_t p, eid_t * &beg_pos, vid_t * &csr, vid_t *nverts, eid_t *nedges);

    bid_t swapOut();

    virtual size_t num_vertices() {
        return blocks[nblocks];
    }
    void exec_updates(SecondOrderRW &userprogram, wid_t nwalks, eid_t *&beg_pos, vid_t *&csr,
                      const vid_t &nBlockVertices, float th);
    void run(vid_t nRootVertices, SecondOrderRW &userprogram, float prob, float th);
};

BucketEngine::BucketEngine(std::string _base_filename, unsigned long long int _blocksize_kb,
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
}

BucketEngine::~BucketEngine() {
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
BucketEngine::load_block_range(std::string base_filename, unsigned long long int blocksize_kb, vid_t *&blocks,
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

void BucketEngine::loadSubGraph(bid_t p, eid_t *&beg_pos, vid_t *&csr, vid_t *nverts, eid_t *nedges) {
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

void BucketEngine::loadPreVertexSubGraph(bid_t p, eid_t *&beg_pos, vid_t *&csr, vid_t *nverts, eid_t *nedges) {
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

void BucketEngine::findSubGraph(bid_t p, eid_t *&beg_pos, vid_t *&csr, vid_t *nverts, eid_t *nedges) {
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


bid_t BucketEngine::swapOut() {
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

void BucketEngine::exec_updates(SecondOrderRW &userprogram, wid_t nwalks, eid_t *&beg_pos, vid_t *&csr,
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
            bid_t preResideBlock = WalkManager::getPreBlock(walk);
            walk_manager->currentWalkBucket[preResideBlock].push_back(walk);
        }
//        m.stop_time("3_CollectBucket");
        for (bid_t p = 0; p < walk_manager->nBlocks; p++){
//            if (walk_manager->currentWalkBucket[p].size() < 0.45 * (blocks[p + 1] - blocks[p])){
//                m.start_time("exec_updates_2od_no-load");
//#pragma omp parallel for schedule(static)
//                for(wid_t i = 0; i < walk_manager->currentWalkBucket[p].size(); i++ ){
//                    userprogram.updateWalk(m, walk_manager->currentWalkBucket[p][i], exec_block, beg_pos, csr, *walk_manager);//, vertex_value);
//                }
//                m.stop_time("exec_updates_2od_no-load");
//                continue;
//            }
            if (walk_manager->currentWalkBucket[p].empty()){
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
        m.start_time("exec_updates_no_bucket");
#pragma omp parallel for schedule(dynamic)
        for(wid_t i = 0; i < nwalks; i++ ){
            WalkDataType walk = walk_manager->currentWalks[i];
            userprogram.updateWalk(m, walk, exec_block, beg_pos, csr, *walk_manager);//, vertex_value);
        }
        m.stop_time("exec_updates_no_bucket");
    }


    // walk_manager->writeblockWalks(exec_block);
}

void BucketEngine::run(vid_t nRootVertices, SecondOrderRW &userprogram, float prob, float th) {
    // srand((unsigned)time(NULL));
    m.start_time("2_StartWalks");
    userprogram.initWalks(nRootVertices, *walk_manager, m);
    m.stop_time("2_StartWalks");

    gettimeofday(&start, NULL);
    m.start_time("1_ProcessTime");

    vid_t nverts, *csr;
    eid_t nedges, *beg_pos;
    /*loadOnDemand -- block loop */
    int blockcount = 0;
    while(!userprogram.allWalksFinished(*walk_manager) ){
        blockcount++;
        m.start_time("1_chooseBlock");
        exec_block = walk_manager->chooseBlock(prob);
        m.stop_time("1_chooseBlock");
        findSubGraph(exec_block, beg_pos, csr, &nverts, &nedges);

        /*load walks info*/
        // walk_manager->loadWalkPool(exec_block);
        wid_t nwalks;
        nwalks = walk_manager->getCurrentWalks(exec_block);

        // if(blockcount % (nBlocks/100+1)==1)
//        if(blockcount % (1024*1024*1024/nedges+1) == 1)
        {
            logstream(LOG_DEBUG) << runtime() << "s : blockcount: " << blockcount << std::endl;
            logstream(LOG_INFO) << "nverts = " << nverts << ", nedges = " << nedges << std::endl;
            logstream(LOG_INFO) << "currentNWalks = " << walk_manager->currentNWalks << ", nwalks[" << exec_block << "] = " << nwalks << std::endl;
        }

        double usage = get_memory();
        memUsage = ((memUsage * memCount) + usage) / (memCount + 1);
        memCount ++;

        exec_updates(userprogram, nwalks, beg_pos, csr, nverts, th);
        walk_manager->updateWalkNum(exec_block);
        // userprogram.compUtilization(beg_pos[nverts] - beg_pos[0]);

    } // For block loop
#if DEBUG
    uint64_t totalSteps = 0;
    for (tid_t t = 0; t < 16; t++){
        totalSteps += userprogram.totalSteps_t.at(t);
    }
    m.set("total-steps", totalSteps);
#endif
    m.set("avg-mem-usage", memUsage);
    m.stop_time("1_ProcessTime");
}

void BucketEngine::print_config() {
    logstream(LOG_INFO) << "Engine configuration: " << std::endl;
    logstream(LOG_INFO) << " exec_threads = " << (int)exec_threads << std::endl;
    logstream(LOG_INFO) << " blocksize_kb = " << blocksize_kb << "kb" << std::endl;
    logstream(LOG_INFO) << " number of total blocks = " << nblocks << std::endl;
    logstream(LOG_INFO) << " number of in-memory blocks = " << nmblocks << std::endl;
}

double BucketEngine::runtime() {
    timeval end;
    gettimeofday(&end, NULL);
    return end.tv_sec-start.tv_sec+ ((double)(end.tv_usec-start.tv_usec))/1.0E6;
}

#endif
#endif //IOE_SORW_BUCKETORDER_HPP
