/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/3/7 16:25
 */

#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <sstream>

int main(){
    std::string oriFileName = "/mnt/home/idmg/lhz/CLionProjects/IOE-SORW/Dataset/soc-LiveJournal1.txt.unDir.rabbitOrder";
    std::string sortFileName = oriFileName + ".sort";
    uint32_t maxVertexId = 4847571;

    std::ifstream oriFile(oriFileName);
    std::string line;
    std::vector<std::list<uint32_t>> edgeList;
    edgeList.resize(maxVertexId + 1);
    uint32_t source, dest;
    while (oriFile.good()){
        getline(oriFile, line);
        if (line.empty()){
            continue;
        }
        if (line.at(0) == '#' || line.at(0) == '%'){
            continue;
        }
        std::istringstream iss(line);
        iss >> source >> dest;
        edgeList.at(source).push_back(dest);
    }
    oriFile.close();

    std::ofstream sortFile(sortFileName);
    uint32_t src = 0;
    for (const auto &nodeOutList: edgeList){
        for (const auto &d: nodeOutList){
            sortFile << src << " " << d << std::endl;
        }
        src++;
    }
    sortFile.close();
}