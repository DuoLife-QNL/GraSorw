/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/4/22 10:40
 */
#include <iostream>
#include <string>
#include <set>
#include <fstream>
#include <unistd.h>
#include "engine/Settings.hpp"
#include "util/RandNum.hpp"

const vid_t nStarts = 100;

int main(int argc, char *argv[]){
    int opt;
    vid_t nVertices;
    opt = getopt (argc, argv, "n:");
    if(opt > 0){
        nVertices = atoi(optarg);
    }else{
        std::cout << "Use -n to declare graph vertex num!" << std::endl;
        return 0;
    }
    RandNum randNum(time(0));
    std::set<vid_t> starts;
    do {
        starts.insert(randNum.iRand(nVertices));
    } while (starts.size() < nStarts);
    std::ofstream startsFile("conf/PPRStarts/PPRStarts.txt");
    for (auto v: starts){
        startsFile << v << std::endl;
    }
    startsFile.close();
}