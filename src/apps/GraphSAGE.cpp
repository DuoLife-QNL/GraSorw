//
// Created by Ethan on 2021/1/8.
//

#include "walks/FwdPropagation.hpp"
#include "util/toplist.hpp"
#include "api/MetricReport.hpp"

int main(int argc, const char **argv){
    set_argc(argc, argv);
    metrics m("GraphSAGE");

    std::string fileName = get_option_string("file", "/home/lihz/Code/GraphWalker/Dataset/test_mini.txt");
    unsigned long long blocksize_kb = get_option_long("blocksize_kb", 0); // Size of block, represented in KB
    bid_t nmblocks = get_option_int("nmblocks", 0); // number of in-memory blocks
    eid_t maxOutDegree = 0;

    int k = 4;
    std::vector<int> s{ 50, 20, 10, 5 };

    vid_t nVertices = get_num_vertices(fileName);
    bid_t nblocks = convert_if_notexists(fileName, blocksize_kb, maxOutDegree, DEFAULT);
    if(nmblocks > nblocks) nmblocks = nblocks;

    logstream(LOG_DEBUG) << "nBlocks nmblocks : " << nblocks << " " << nmblocks << std::endl;

    GraphWalkerEngine engine(fileName, blocksize_kb, nblocks, nmblocks, m);

    FwdPropagation program(k, s, engine.beg_posf, engine.csrf, nVertices, nblocks, engine.blocks, m);
    program.setStartVertexId(0);

    engine.run(nVertices, program, 0, 0);

    metrics_report(m);
    return 0;
}