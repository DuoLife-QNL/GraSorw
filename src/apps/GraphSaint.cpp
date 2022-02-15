//
// Created by Ethan on 2021/1/7.
//

#include <string>

#include "walks/MRW.hpp"
#include "BasicIncludes.hpp"
#include "api/MetricReport.hpp"
#include "util/toplist.hpp"

int main(int argc, const char **argv){
    set_argc(argc, argv);
    metrics m("GraphSaint");

    std::string fileName = get_option_string("file", "/home/lihz/Code/GraphWalker/Dataset/test_mini.txt");
    vid_t nRoots = get_option_int("root-num", 100);
    vid_t nodeBudget = get_option_int("node-budget", 200);
    vid_t nVertices = get_num_vertices(fileName);

    std::string invlname = fidname(fileName, 0);
    std::string beg_posname = invlname + ".beg_pos";
    std::string csrname = invlname + ".csr";
    int beg_posf = open(beg_posname.c_str(),O_RDONLY | O_CREAT, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
    int csrf = open(csrname.c_str(),O_RDONLY | O_CREAT, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);

    MultiRW program(nRoots, nodeBudget, beg_posf, csrf, nVertices, m);

    m.start_time("run");
    program.run();
    m.stop_time("run");

    metrics_report(m);

    return 0;

}