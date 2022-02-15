//
// Created by lihz on 2020/10/29.
//

#include <cstdio>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <unistd.h>
#include <sstream>
#include <map>
#include <set>

#include "engine/Settings.hpp"

using namespace std;

vid_t *nodeid;
vid_t *degrees;
vid_t *nodeindex;
vid_t **neighbors;

void stupidConversion(int argc, char *argv[]){
    char baseFileName[100];
    char unDirFileName[100];
    vid_t n;
    eid_t e;
    vector<vid_t> xs, ys;

    strcpy(baseFileName, argv[1]);
    strcpy(unDirFileName, baseFileName);
    strcat(unDirFileName, ".unDir");
    ofstream unDirFile(unDirFileName);
    FILE *fp = fopen(baseFileName, "r");

    /**
     * load edges from disk to memory
     * @x: source vertex
     * @y: destination vertex
     * @maxi: the max node id in base edgeList file
     */
    vid_t maxi = 0, x, y;
    e = 0;

    char s[1024];
    while(fgets(s, 1024, fp) != NULL) {
        if (s[0] == '#') continue; // Comment
        if (s[0] == '%') continue; // Comment
        char *t1, *t2;
        t1 = strtok(s, "\t, ");
        t2 = strtok(NULL, "\t, ");
        if (t1 == NULL || t2 == NULL ) {
            cout << "Input file is not in right format. "
                 << "Expecting <from> <to>. "
                 << "Current line: " << s << endl;
            assert(false);
        }
        x = atoi(t1);
        y = atoi(t2);
        if (x == y){
            continue;
        }
        e++;
        if (x > maxi) maxi = x;
        if (y > maxi) maxi = y;
        xs.push_back(x);
        ys.push_back(y);
    }
    std::cout << "scanned edges: " << e << std::endl;

    fclose(fp);

    /**
     * the map table mapping old id to new id
     * nodeid[oldid] = newid
     */
    nodeid = static_cast<vid_t *>(malloc((maxi + 1) * sizeof(vid_t)));
    for (vid_t i = 0; i <= maxi; i++)
        nodeid[i] = -1;
    n = 0;
    /**
     * assign nodes with unique nodeid from 0 to e-1 (edge nums - 1)
     * give id to nodes increasingly following the order when they first appear
     * @n: max id + 1
     */
    for (eid_t j = 0; j < e; j++) {
        if (nodeid[xs[j]] == -1) {
            nodeid[xs[j]] = n;
            ++n;
        }
        if (nodeid[ys[j]] == -1) {
            nodeid[ys[j]] = n;
            ++n;
        }
    }

    std::cout << "max new id " << n - 1 << std::endl;

    degrees = static_cast<vid_t *>(malloc(n * sizeof(vid_t)));
    nodeindex = static_cast<vid_t *>(malloc(n * sizeof(vid_t)));
    neighbors = static_cast<vid_t **>(malloc(n * sizeof(vid_t *)));
    memset(degrees, 0, n * sizeof(vid_t));
    memset(nodeindex, 0, n * sizeof(vid_t));

    /* The degree of the node includes in-degree and out-degree */
    /* FIXME: 这里如果原本就有双向边，会导致输出两遍同一条边，使输出图非简单无向图（两节点间仅有一边）*/
    for (eid_t j = 0; j < e; j++) {
        degrees[nodeid[xs[j]]]++;
        degrees[nodeid[ys[j]]]++;
    }

    for (vid_t i = 0; i < n; i++)
        neighbors[i] = static_cast<vid_t *>(malloc(degrees[i] * sizeof(vid_t)));


    for (eid_t j = 0; j < e; j++) {
        neighbors[nodeid[xs[j]]][nodeindex[nodeid[xs[j]]]++] = nodeid[ys[j]];
        neighbors[nodeid[ys[j]]][nodeindex[nodeid[ys[j]]]++] = nodeid[xs[j]];
    }
    for (vid_t i = 0; i < n; i++) {
        nodeindex[i] = 0;
        if (degrees[i] == 0) {
            continue;
        }
        /* sort node's neighbors by nodeid */
        sort(neighbors[i], neighbors[i] + degrees[i]);
        /**
         * each node has at least one neighbor, according to the aforementioned * definition of neighbor
         */
        nodeindex[i] = 1;
        for (vid_t j = 1; j < degrees[i]; j++)
            if (neighbors[i][j] != neighbors[i][j - 1])
                nodeindex[i]++;
    }

    vid_t preFrom, preTo;
    unDirFile << 0 << " " << neighbors[0][0] << std::endl;
    preFrom = 0;
    preTo = neighbors[0][0];
    eid_t totalNEdge = 1;
    for (vid_t from = 0; from < n; from++){
        for (vid_t toIndex = 0; toIndex < degrees[from]; toIndex++){
            vid_t to = neighbors[from][toIndex];
            if (from != preFrom || to != preTo){
                unDirFile << from << " " << to << endl;
                totalNEdge++;
            }
            if (toIndex == 0){
                preFrom = from;
            }
            preTo = to;
        }
    }
    std::cout << "total vertex number: " << n << std::endl;
    std::cout << "total edge number: " << totalNEdge << std::endl;
    unDirFile.close();
}

void convertWithSet(int argc, char *argv[]){
    int opt;
    string baseFileName;
    opt = getopt (argc, argv, "f:");
    if(opt){
        baseFileName = optarg;
    }else{
        cout << "Use -f to declare input file name!" << endl;
        return;
    }
    string outFileName = baseFileName + ".unDir";
    ifstream inFile(baseFileName);
    vid_t maxOriId = 0, sourceVertex, destVertex;
    string line;
//    vector<vid_t> oriSrcV, oriDestV;
    eid_t nReadEdge = 0;
    while (inFile.good()){
        getline(inFile, line);
        if (line.empty()){
            continue;
        }
        if (line.at(0) == '#' || line.at(0) == '%'){
            continue;
        }
        istringstream iss(line);
        iss >> sourceVertex;
        iss >> destVertex;
        if (sourceVertex == destVertex){
            continue;
        }
        nReadEdge++;
        if (sourceVertex > maxOriId)
            maxOriId = sourceVertex;
        if (destVertex > maxOriId)
            maxOriId = destVertex;
//        oriSrcV.push_back(sourceVertex);
//        oriDestV.push_back(destVertex);
        if (nReadEdge % 1000000 == 0){
            std::cout << "read edge " << nReadEdge << std::endl;
        }
    }

    vector<vid_t> newId;
    newId.assign(maxOriId + 1, INVALID_VID);
    vid_t oriSrcId, oriDestId;
    vid_t newSrcId, newDestId;
    vid_t assignNewId = 0;

    vector<set<vid_t>> unDirGraph;
    unDirGraph.resize(maxOriId + 1);

    inFile.clear();
    inFile.seekg(0);
    eid_t e = 0;
    while (e < nReadEdge){
        getline(inFile, line);
        if (line.empty()){
            continue;
        }
        if (line.at(0) == '#' || line.at(0) == '%'){
            continue;
        }
        istringstream iss(line);
        iss >> oriSrcId;
        iss >> oriDestId;
        if (sourceVertex == destVertex){
            continue;
        }
        e++;

        if (newId.at(oriSrcId) == INVALID_VID){
            newId.at(oriSrcId) = assignNewId;
            assignNewId++;
        }
        if (newId.at(oriDestId) == INVALID_VID){
            newId.at(oriDestId) = assignNewId;
            assignNewId++;
        }

        newSrcId = newId.at(oriSrcId);
        newDestId = newId.at(oriDestId);
        unDirGraph.at(newSrcId).insert(newDestId);
        unDirGraph.at(newDestId).insert(newSrcId);
        if (e % 1000000 == 0){
            std::cout << "insert edge " << e << std::endl;
        }
    }
    std::cout << "max assigned new id = " << assignNewId - 1 << std::endl;

    std::ofstream unDirFile(outFileName);
    for (vid_t src = 0; src < unDirGraph.size(); src++){
        for (auto dest: unDirGraph.at(src)){
            unDirFile << src << " " << dest << std::endl;
        }
        if (src % 10000 == 0){
            std::cout << "output src: " << src << std::endl;
        }
    }
    unDirFile.close();
}

int main (int argc, char *argv[])
{
    stupidConversion(argc, argv);

    return 0;
}