#include <iostream>
#include "BasicIncludes.hpp"
#include "util/toplist.hpp"
#include "api/MetricReport.hpp"
//#include "Node2VecOnline.hpp"
#include "Node2VecReject.hpp"
#include "Node2VecARRejection.hpp"
#include "PPRonNV.hpp"
#include "LinkSCANReject.hpp"
#include "preprocess/cmpEdgeCut.hpp"
#include "DeepWalk.hpp"
#include "engine/SOGraphWalker.hpp"
#include "engine/GraphWalker.hpp"
#include "engine/DualBucket.hpp"
#include "engine/DualBucketMT.hpp"
#include "engine/BucketOrder.hpp"
#include "engine/PreTrain.hpp"

metrics runProgram(){
    metrics m("SORW");

    std::string fileName = get_option_string("file", "/home/lihz/Code/GraphWalker/Dataset/test_mini.txt");
    std::string startVerticesFileName = get_option_string("starts-vertices-file", "conf/PPRStarts.cnf");
    loadTestOutPutFileName_dynamicBlocks = get_option_string("load-test-output-file-dynamic", LOAD_TRADEOFF_FILE);
    loadTestOutPutFileName_staticBlocks = get_option_string("load-test-output-file-static", LOAD_TRADEOFF_FILE);
    blockUtiFileName = get_option_string("block-uti-file", BLOCK_UTI_FILE_NAME);
#if (NO_LOAD_TEST | ONDEMAND_LOAD) & USE_TRAINED_THS
    thsFile_dynamic = get_option_string("ths-file-dynamic", THS_FILE);
    thsFile_static = get_option_string("ths-file-static", THS_FILE);
#endif
    int nWalksPerVertex = get_option_int("walks-per-vertex", 1);
    int nStartVertices = get_option_int("num-start-vertices", 100);
    double decayFactor = get_option_float("decay-factor", 0.85);
    hid_t walkLength = get_option_int("walk-length", 20);
    unsigned long long blocksize_kb = get_option_long("blocksize_kb", 0); // Size of block, represented in KB
    float probChooseMinSteps = get_option_float("prob", 0);
    bid_t nmblocks = get_option_int("nmblocks", 2); // number of in-memory blocks
    float p = get_option_float("p", 1);
    float q = get_option_float("q", 1);
    float alpha = get_option_float("alpha", 0.2);
    float th = get_option_float("bucket-th", 0);
    std::cout << "bucket-th = " << th << std::endl;
    eid_t maxOutDegree = 0;

    m.set("walks-per-vertex", nWalksPerVertex);
    m.set("walkLength", (int) walkLength);
    m.set("p", p);
    m.set("q", q);
    m.set("alpha", alpha);

    PAR_METHOD parMethod = DEFAULT;
    /* graph partition command */
    std::string grpParCmd = get_option_string("graph-partition-method", "default");
    std::string parFileName = get_option_string("graph-partition-file", "ERROR");
    int preParNBlocks = get_option_int("pre-partition-block-nums", 0);
    if (grpParCmd == "metis"){
        parMethod = METIS;
    }

    m.set("walks-per-vertex", nWalksPerVertex);
    m.set("walkLength", (int) walkLength);
    m.set("p", p);
    m.set("q", q);

    bid_t nBlocks;
    vid_t nVertices;
    nVertices = get_option_long("num-vertices", 0);
    if (!nVertices){
        nVertices = get_num_vertices(fileName);
    }
    if (parMethod == METIS){
        /* use graph partition from metis */
        convert_if_notexists(fileName, blocksize_kb, maxOutDegree, METIS, parFileName, preParNBlocks);
        fileName.append(".MetisReOrder");
        nBlocks = preParNBlocks;
        if(nmblocks > preParNBlocks) nmblocks = nBlocks;

    }else{
        vid_t nwalks = nVertices * nWalksPerVertex;
        if(blocksize_kb == 0){
            int dom = 0;
            while(nwalks){
                dom++;
                nwalks /= 10;
            }
            blocksize_kb = pow(2, 10 + dom);
        }
        if(nmblocks == 0){
            nmblocks = MEM_BUDGET / blocksize_kb;
        }
        nBlocks = convert_if_notexists(fileName, blocksize_kb, maxOutDegree, DEFAULT);
        blockSize_kb_global = blocksize_kb;
        if(nmblocks > nBlocks) nmblocks = nBlocks;
    }
    logstream(LOG_DEBUG) << "nBlocks nmblocks : " << nBlocks << " " << nmblocks << std::endl;

#if BI_BLOCK
    DualBucketEngine engine(fileName, blocksize_kb, nBlocks, nmblocks, m);
//    DualBucketEngine_MT engine(fileName, blocksize_kb, nBlocks, nmblocks, m);
#elif PLAIN_BUCKET
    BucketEngine engine(fileName, blocksize_kb, nBlocks, nmblocks, m);
#elif PLAIN
    SOGraphWalkerEngine engine(fileName, blocksize_kb, nBlocks, nmblocks, m);
#endif
#if FIRST_ORDER_ENGINE
    GraphWalkerEngine engine(fileName, blocksize_kb, nBlocks, nmblocks, m);
#endif
#if DEEPWALK
    DeepWalk program(nBlocks, engine.blocks, walkLength, nWalksPerVertex, engine.beg_posf, engine.csrf, nVertices, m);
#elif RWNV
    Node2VecReject program(nBlocks, engine.blocks, walkLength, nWalksPerVertex, maxOutDegree, p, q, engine.beg_posf,
                           engine.csrf, nVertices, m);
#elif PRNV
    PPRonNV program(nBlocks, engine.blocks, walkLength, nWalksPerVertex, maxOutDegree, p, q, decayFactor,
                                 engine.beg_posf, engine.csrf, startVerticesFileName, nStartVertices,
                                 nVertices, m);
#endif
    engine.run(nVertices, program, probChooseMinSteps, th);
    metrics_report(m);
#if TIME_COST_INFO
    tid_t nThreads = engine.exec_threads;
    double totTime = 0;
    std::string metricItemName = "get-csr";
    for (tid_t t = 0; t < nThreads; t++){
        auto time = program.vertexIO.metric.at(t).get(metricItemName).value;
        std::cout << "thread" + std::to_string(t) + " " + metricItemName + ": " << time
                  << "(" << program.vertexIO.metric.at(t).get(metricItemName).count << ")" << std::endl;
        totTime += time;
    }
    std::cout << "avg " + metricItemName + " time: " << totTime / nThreads << std::endl;
#endif

    return m;
}

const std::vector<std::string>itemName{
        "1_ProcessTime",
        "2_StartWalks",
        "g_loadSubGraph",
        "z_w_readWalksfromDisk",
        "4_writeWalks2Disk_",
        "5_exec_updates"
};

const std::vector<std::string>displayName{
        "total time",
        "start walks time",
        "block IO",
        "walk pool I",
        "walk pool O",
        "exec update"
};

class Unit{
public:
    double time;
    size_t count;
    Unit(double _time, size_t _count){
        time = _time;
        count = _count;
    }
};

class OneExec{
public:
    std::vector<Unit>result;
    tid_t nThreads;

    explicit OneExec(metrics &m, tid_t _nThreads){
        nThreads = _nThreads;
        for (const auto& s:itemName){
            double time = 0;
            size_t count = 0;
            time = m.get(s).value;
            count = m.get(s).count;
            result.emplace_back(time, count);
        }
    }
};

class MetricResult{
private:
    int nThreads;
public:
    std::vector<OneExec>results;

    explicit MetricResult(tid_t _nThreads){
        nThreads = _nThreads;
    }

    void addResult(metrics &m){
        results.emplace_back(OneExec(m, nThreads));
    }

    OneExec getResult(int round){
        return results.at(round);
    }

    Unit getAvg(int itemId){
        double time = 0;
        double count = 0;
        int round = results.size();
        for (const auto& r: results){
            time += r.result.at(itemId).time / round;
            count += static_cast<double>(r.result.at(itemId).count) / round;
        }
        return {time, static_cast<size_t>(count)};
    }

    void show(){
        for (int i = 0; i < itemName.size(); i++){
            Unit itemAvg(getAvg(i));
            std::cout << displayName.at(i) << ": " << itemAvg.time << "  " << itemAvg.count << std::endl;
        }
    }

};

int main(int argc, const char ** argv){
    set_argc(argc, argv);
    tid_t nthreads = get_option_int("execthreads", omp_get_max_threads());
    MetricResult metricResult(nthreads);
    int execRound = 1;
    for (int round = 0; round < execRound; round++){
        omp_set_num_threads(nthreads);
        metrics m = runProgram();
        metricResult.addResult(m);
    }
    metricResult.show();

    for (int round = 0; round < execRound; round++){
        double runTime = metricResult.getResult(round).result.at(0).time;
        std::cout << "round " << round << " run time: " << runTime << std::endl;
    }

    return 0;


}
