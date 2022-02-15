/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/3/7 15:58
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

int main(){
    std::string edgeListFileName = "/mnt/home/idmg/lhz/CLionProjects/IOE-SORW/Dataset/soc-LiveJournal1.txt.unDir";
    std::string rabbitTmpFileName = "/mnt/home/idmg/lhz/CLionProjects/IOE-SORW/Dataset/soc-LiveJournal1.txt.unDir.rabbitOrder-tmp";
    std::string outPutFileName = "/mnt/home/idmg/lhz/CLionProjects/IOE-SORW/Dataset/soc-LiveJournal1.txt.unDir.rabbitOrder";

    std::ifstream rabbitOrder(rabbitTmpFileName);
    std::string line;
    std::vector<uint32_t> newId;
    uint32_t tmp;
    while (rabbitOrder.good()){
        getline(rabbitOrder, line);
        if (line.empty()){
            continue;
        }
        if (line.at(0) == '#' || line.at(0) == '%'){
            continue;
        }
        istringstream iss(line);
        iss >> tmp;
        newId.push_back(tmp);
    }
    rabbitOrder.close();

    std::ifstream oldEdgeList(edgeListFileName);
    std::ofstream newEdgeList(outPutFileName);
    uint32_t source, dest;
    while (oldEdgeList.good()){
        getline(oldEdgeList, line);
        if (line.empty()){
            continue;
        }
        if (line.at(0) == '#' || line.at(0) == '%'){
            continue;
        }
        istringstream iss(line);
        iss >> source >> dest;
        newEdgeList << newId.at(source) << " " << newId.at(dest) << std::endl;
    }
    oldEdgeList.close();
    newEdgeList.close();

}