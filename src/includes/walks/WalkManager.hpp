#ifndef DEF_GRAPHWALKER_WALK
#define DEF_GRAPHWALKER_WALK

#include <iostream>
#include <cstdio>
#include <sstream>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <string>
#include <queue>
#include <chrono>
#include <thread>

#include "metrics/metrics.hpp"
#include "api/filename.hpp"
#include "api/io.hpp"
#include "engine/WalkBuffer.hpp"
#include "preprocess/conversions.hpp"
#include "util/CircularQueue.hpp"
#include "metrics/FileSize.hpp"
#include "metrics/Colors.hpp"
#include "IO/VertexIO.hpp"


class WalkManager{
protected:
    std::string baseFile;
    metrics &m;

private:
    CircularQueue walkQueue;
public:
    tid_t nThreads;
    bid_t nBlocks;
    wid_t *nWalks;
    wid_t *nDiskWalks;
    hid_t* minStep;
    WalkBuffer **walkPool;

    WalkDataType *currentWalks = nullptr;
    std::vector<WalkDataType> *currentWalkBucket;

    std::vector<WalkDataType> *bucket;
    std::vector<WalkDataType> *curBucket;
    std::vector<WalkDataType> *preBucket;

    std::vector<WalkDataType> **extendBucket;

    wid_t currentNWalks = 0;
    wid_t curBlockNWalks;
    bool *isModified;

    int *preBucketFiles;
    int *curBucketFiles;
    wid_t *preBucketsSize;
    wid_t *curBucketsSize;

    WalkManager(metrics &m, bid_t nBlocks, tid_t nThreads, std::string baseFile);

    static WalkDataType encode(vid_t sourceId, vid_t preId, vid_t curId, bid_t preBlock, bid_t curBlock, hid_t hop){
        /* & MAX: in case that no such entry (item) in walk representation */
        WalkDataType walk = (((WalkDataType)sourceId & MAX_SRC_VERTEX) << SRC_VERTEX_S)
                            | (((WalkDataType)preId & MAX_PRE_VERTEX) << PRE_VERTEX_S)
                            | (((WalkDataType)curId & MAX_CUR_VERTEX) << CUR_VERTEX_S)
                            | (((WalkDataType)preBlock & MAX_PRE_BLOCK) << PRE_BLOCK_S)
                            | (((WalkDataType)curBlock & MAX_CUR_VERTEX) << CUR_BLOCK_S)
                            | ((WalkDataType)hop & MAX_HOP);
        return walk;
    }


    static vid_t getSouceVertex(const WalkDataType &walk){
        return (vid_t)(walk >> SRC_VERTEX_S);
    }


    static vid_t getPreviousVertex(const WalkDataType &walk){
        return (vid_t)(walk >> PRE_VERTEX_S);
    }


    static vid_t getCurrentVertex(const WalkDataType &walk){
        return (vid_t)(walk >> CUR_VERTEX_S) & MAX_CUR_VERTEX;
    }

    static bid_t getPreBlock(const WalkDataType &walk){
        return (bid_t)(walk >> PRE_BLOCK_S) & MAX_PRE_BLOCK;
    }

    static bid_t getCurBlock(const WalkDataType &walk){
        return (bid_t)(walk >> CUR_BLOCK_S) & MAX_CUR_BLOCK;
    }

    static hid_t getHops(const WalkDataType &walk){
        return (hid_t)(walk >> HOP_S) & MAX_HOP;
    }


	wid_t getCurrentWalks(bid_t p);

    void clearRecoredWalkNum(bid_t block){
        nWalks[block] = 0;
    }

    /* nAppendHops: num of newly walked steps*/
    WalkDataType reEncode(WalkDataType walk, vid_t previousVertex, vid_t currentVertex, hid_t nAppendHops,
                          bid_t preResideBlock);

    void moveWalk(WalkDataType walk, bid_t toBlock, tid_t thread, vid_t previousVertex, vid_t currentVertex,
                  hid_t nAppendHops, bid_t preResideBlock);

    void moveWalk(const WalkDataType walk, const tid_t thread, const bid_t toBlock){
        if(walkPool[thread][toBlock].size_w == WALK_BUFFER_SIZE){
            writeWalks2Disk(thread, toBlock);
        }
        assert(walkPool[thread][toBlock].size_w < WALK_BUFFER_SIZE);
        walkPool[thread][toBlock].push_back(walk);
        isModified[toBlock] = true;
    }

    void writeWalks2Disk(tid_t thread, bid_t block);

    wid_t getNCurrentWalks(bid_t block);

    void readWalksFromDisk(bid_t block);

    void updateWalkNum(bid_t block);

    void updateWalkNum(){
        wid_t newNWalks = 0;
        for (bid_t b = 0; b < nBlocks; b++){
            if (isModified[b]){
                isModified[b] = false;
                wid_t newBlockWalkNum = 0;
                newBlockWalkNum = nDiskWalks[b];
                for(tid_t t = 0; t < nThreads; t++){
                    newBlockWalkNum += walkPool[t][b].size_w;
                }
                nWalks[b] = newBlockWalkNum;
                newNWalks += newBlockWalkNum;
            }else{
                newNWalks += nWalks[b];
            }
        }
        currentNWalks = newNWalks;
        free(currentWalks);
        currentWalks = nullptr;
    }

    void collectBucket(bid_t staticBlock){
        m.add("total-collected-walks", curBlockNWalks);
        m.start_time("3_CollectBucket");
        for (wid_t w = 0; w < curBlockNWalks; w++){
            WalkDataType walk = currentWalks[w];
            bid_t preBlock = WalkManager::getPreBlock(walk);
            if (preBlock == staticBlock){
                bucket[WalkManager::getCurBlock(walk)].push_back(walk);
            }else{
                bucket[preBlock].push_back(walk);
            }
        }
        m.stop_time("3_CollectBucket");
    }

    void mapBucketVertex(bid_t dynamicBlock, VertexIO &vertexIo){
        m.start_time("map-vertex");
        for (const auto &walk: bucket[dynamicBlock]){
            bid_t preBlock = WalkManager::getPreBlock(walk);
            vid_t v;
            if (preBlock == dynamicBlock){
                v = WalkManager::getPreviousVertex(walk);
            }else{
                v = WalkManager::getCurrentVertex(walk);
            }
            vertexIo.needVertex(v);
        }
        m.stop_time("map-vertex");
    }

    void mapBlockVertex(bid_t staticBlock, VertexIO &vertexIo){
        m.start_time("map-vertex");
        for (wid_t w = 0; w < curBlockNWalks; w++){
            WalkDataType walk = currentWalks[w];
            bid_t preBlock = WalkManager::getPreBlock(walk);
            vid_t v;
            if (preBlock == staticBlock){
                v = WalkManager::getPreviousVertex(walk);
            }else{
                v = WalkManager::getCurrentVertex(walk);
            }
            vertexIo.needVertex(v);
        }
        m.stop_time("map-vertex");
    }

    void clearBucket(){
        for (bid_t b = 0; b < nBlocks; b++){
            bucket[b].clear();
        }
    }

    void exendBucket(bid_t dynamicBlock){
        for (tid_t t = 0; t < nThreads; t++){
            bucket[dynamicBlock].insert(bucket[dynamicBlock].end(), extendBucket[dynamicBlock][t].begin(), extendBucket[dynamicBlock][t].end());
        }
    }

    void clearExtendBucket(bid_t dynamicBlock){
        for (tid_t t = 0; t < nThreads; t++){
            extendBucket[dynamicBlock][t].clear();
        }
    }

    void setMinStep(bid_t p, hid_t hop );

    bid_t blockWithMaxWalks();

    bid_t blockWithMinStep();

    bid_t chooseBlock(float prob);
#if BI_BLOCK
    void collectBuckets2Disk(bid_t staticBlock){
//        m.start_time("3_CollectBucket");
        walkQueue.ClearQueue();
        wid_t inCount = 0;
        wid_t outCount = 0;
        bool allWalksRead = false;
        for (bid_t b = 0; b < nBlocks; b++){
            ftruncate(preBucketFiles[b], 0);
            ftruncate(curBucketFiles[b], 0);
            preBucketsSize[b] = 0;
            curBucketsSize[b] = 0;
        }
#pragma omp parallel
        {
            tid_t t = omp_get_thread_num();
#pragma omp single
            {
                if (omp_get_num_threads() < 2){
                    assert(false);
                }
            };
#pragma omp barrier
            if (t == 0){
                for (tid_t tid = 0; tid < nThreads; tid++){
                    if (walkPool[tid][staticBlock].size_w > 0){
                        for(wid_t w = 0; w < walkPool[tid][staticBlock].size_w; w++){
                            while (true){
                                if (!walkQueue.QueueFull()){
                                    WalkDataType &walk = walkPool[tid][staticBlock][w];
                                    walkQueue.Enqueue(walk);
//                                    uint64_t high = walk >> 64;
//                                    std::cout << "enqueued walk" << inCount << "high: " << std::hex << high << std::endl;
//                                    uint64_t low = (uint64_t)walk;
//                                    std::cout << "enqueued walk" << inCount << " low: " << std::hex << low << std::endl;
                                    std::cout << RED << "walk id " << inCount << " enqueued " << std::endl;
                                    std::this_thread::sleep_for(std::chrono::milliseconds(20));
                                    break;
                                }
                            }
                            inCount++;
//                            std::cout << "added to queue: " << inCount << std::endl;
                        }
                        walkPool[tid][staticBlock].size_w = 0;
//                        std::cout << "walk pool added" << std::endl;
                    }
                    if (nDiskWalks[staticBlock] > 0){
                        std::string walksFile = walksname(baseFile, staticBlock);
                        int f = open(walksFile.c_str(), O_RDWR, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
                        assert(f > 0);
                        WalkDataType *diskWalks;
                        diskWalks = (WalkDataType *)mmap(nullptr, nDiskWalks[staticBlock] * sizeof(WalkDataType), PROT_READ, MAP_PRIVATE, f, 0);
                        for (wid_t w = 0; w < nDiskWalks[staticBlock]; w++){
                            while (true){
                                if (!walkQueue.QueueFull()){
                                    walkQueue.Enqueue(diskWalks[w]);
                                    std::cout << RED << "walk id " << inCount << " enqueued " << std::endl;
                                    std::this_thread::sleep_for(std::chrono::milliseconds(20));
                                    break;
                                }
                            }
                            inCount++;
                        }
                        munmap(diskWalks, nDiskWalks[staticBlock] * sizeof(WalkDataType));
                        ftruncate(f, 0);
                        close(f);
                        unlink(walksFile.c_str());
                        nDiskWalks[staticBlock] = 0;
                        allWalksRead = true;
                    }
                }
            }else if (t == 1){
                while (!allWalksRead || !walkQueue.QueueEmpty()){
                    WalkDataType walk;
                    bool dequeued = false;
#pragma omp critical (Queue)
                    {
                        if (!walkQueue.QueueEmpty()){
                            walkQueue.DeQueue(walk);
                            outCount++;
                            dequeued = true;
                        }
                    }
                    if (dequeued){
                        bid_t preBlock = getPreBlock(walk);
                        bid_t curBlock = getCurBlock(walk);
//                        std::cout << "dequeued walk " << outCount - 1 << " pre block = " << preBlock << std::endl;
//                        std::cout << "dequeued walk " << outCount - 1 << " cur block = " << curBlock << std::endl;
                        if (preBlock == staticBlock){
//                        bid_t curBlockId = getCurBlock(walk);

                            uint64_t high = walk >> 64;
//                            std::cout << "dequeued walk " << outCount - 1 << ", high: " << std::hex << high << std::endl;
//                            uint64_t low = (uint64_t)walk;
//                            std::cout << "dequeued walk " << outCount - 1 << ", low: " << std::hex << low << std::endl;
//                            std::cout << "adding one walk to curBucket " << curBlock << std::endl;
                            addWalk2BucketFile(walk, curBlock, CUR);
                            std::string bucketFileName = bucketName(baseFile, curBlock, CUR);
                            int bucketFileSize = file_size(bucketFileName.c_str());
                            std::cout << GREEN << "walk id " << outCount - 1 << " dequeued" << std::endl
                                      << "    " << "to CUR bucket: " << curBlock << std::endl
                                      << "    " << "bucket file size: " << bucketFileSize << std::endl;
#pragma omp critical (curBucketsSize)
                            {
                                curBucketsSize[curBlock]++;
                            };
                        }else{
                            addWalk2BucketFile(walk, preBlock, PRE);
                            std::string bucketFileName = bucketName(baseFile, preBlock, PRE);
                            int bucketFileSize = file_size(bucketFileName.c_str());
                            std::cout << GREEN << "walk id " << outCount - 1 << " dequeued" << std::endl
                                      << "    " << "to PRE bucket: " << preBlock << std::endl
                                      << "    " << "bucket file size: " << bucketFileSize << std::endl;
#pragma omp critical (preBucketsSize)
                            {
                                preBucketsSize[preBlock]++;
                            };
                        }
                    }

                }
            }
        };
        assert(inCount == outCount);
        assert(inCount == nWalks[staticBlock]);
        std::vector<int>precur = {PRE, CUR};
        for (auto flag: precur){
            for (bid_t b = 0; b < nBlocks; b ++){
                int f = getBucketFile(b, flag);
                std::string fileName = bucketName(baseFile, b, flag);
                int fileSize = file_size(fileName.c_str());
                int bucketSize = getBucketSize(b, flag);
                assert(fileSize == bucketSize * sizeof(WalkDataType));
            }
        }
//        m.stop_time("3_CollectBucket");
    }

    void readBucket(WalkDataType *walks, bid_t dynamicBlock, int preOrCurFlag){
        std::string bucketFile = bucketName(baseFile, dynamicBlock, preOrCurFlag);
//        int f = open(bucketFile.c_str(), O_RDWR, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
//        assert(f > 0);
        int f = getBucketFile(dynamicBlock, preOrCurFlag);
        int fileSize = file_size(bucketFile.c_str());
        wid_t bucketSize = getBucketSize(dynamicBlock, preOrCurFlag);
        logstream(DEBUG) << "file: " << bucketFile << " size: " << fileSize;
        assert(fileSize == bucketSize * sizeof(WalkDataType));
        preada(f, walks, bucketSize * sizeof(WalkDataType));
        /* clear file */
        ftruncate(f,0);
        close(f);
        /* remove the walk file*/
        unlink(bucketFile.c_str());
    }

    wid_t getBucketSize(bid_t dynamicBlock, int preOrCurFlag){
        if (preOrCurFlag == PRE){
            return preBucketsSize[dynamicBlock];
        }else if (preOrCurFlag == CUR){
            return curBucketsSize[dynamicBlock];
        }else{
            assert(false);
        }
    }

    void addWalk2BucketFile(const WalkDataType &walk, bid_t dynamicBlock, int preOrCurFlag){
        int bucketFile = getBucketFile(dynamicBlock, preOrCurFlag);
        pwritea(bucketFile, &walk, sizeof(WalkDataType));
    }

    int getBucketFile(bid_t dynamicBlock, int preOrCurFlag) const{
        if (preOrCurFlag == PRE){
            return preBucketFiles[dynamicBlock];
        }else if (preOrCurFlag == CUR){
            return curBucketFiles[dynamicBlock];
        }else{
            assert(false);
        }
    }
#endif
};

WalkManager::WalkManager(metrics &m, bid_t nBlocks, tid_t nThreads, std::string baseFile)
        :baseFile(baseFile), nBlocks(nBlocks), nThreads(nThreads), m(m){
    walkPool = new WalkBuffer*[nThreads];
    for(tid_t i = 0; i < nThreads; i++)
        walkPool[i] = new WalkBuffer[nBlocks];

    nWalks = (wid_t*)malloc(nBlocks*sizeof(wid_t));
    nDiskWalks = (wid_t*)malloc(nBlocks*sizeof(wid_t));
    minStep = (hid_t*)malloc(nBlocks*sizeof(hid_t));
    memset(nWalks, 0, nBlocks*sizeof(wid_t));
    memset(nDiskWalks, 0, nBlocks*sizeof(wid_t));
    memset(minStep, 0xffff, nBlocks*sizeof(hid_t));
    currentNWalks = 0;

    rm_dir((baseFile+"_GraSorw/walks/").c_str());
    mkdir((baseFile+"_GraSorw/walks/").c_str(), 0777);

    isModified = (bool*)malloc(nBlocks*sizeof(bool));
    memset(isModified, false, nBlocks*sizeof(bool));
    currentWalkBucket = new std::vector<WalkDataType>[nBlocks];
    bucket = new std::vector<WalkDataType>[nBlocks];
    curBucket = new std::vector<WalkDataType>[nBlocks];
    preBucket = new std::vector<WalkDataType>[nBlocks];
    extendBucket = new std::vector<WalkDataType> *[nBlocks];
    for (bid_t b = 0; b < nBlocks; b++){
        extendBucket[b] = new std::vector<WalkDataType>[nThreads];
    }
    walkQueue.InitQueue();
    preBucketFiles = new int[nBlocks];
    curBucketFiles = new int[nBlocks];
    preBucketsSize = new wid_t [nBlocks];
    curBucketsSize = new wid_t [nBlocks];
    for (bid_t b = 0; b < nBlocks; b++){
        std::string prebucketFile = bucketName(baseFile, b, PRE);
        preBucketFiles[b] = open(prebucketFile.c_str(), O_WRONLY | O_CREAT | O_APPEND, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
        std::string curBucketFile = bucketName(baseFile, b, PRE);
        curBucketFiles[b] = open(curBucketFile.c_str(), O_WRONLY | O_CREAT | O_APPEND, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);

        preBucketsSize[b] = 0;
        curBucketsSize[b] = 0;
    }
}

WalkDataType
WalkManager::reEncode(WalkDataType walk, vid_t previousVertex, vid_t currentVertex, hid_t nAppendHops,
                      bid_t preResideBlock) {
    walk += nAppendHops;
    vid_t source = getSouceVertex(walk);
    return encode(source, previousVertex, currentVertex, preResideBlock, 0, getHops(walk));
}

void WalkManager::moveWalk(WalkDataType walk, bid_t toBlock, tid_t thread, vid_t previousVertex, vid_t currentVertex,
                           hid_t nAppendHops, bid_t preResideBlock) {
    if(walkPool[thread][toBlock].size_w == WALK_BUFFER_SIZE){
        writeWalks2Disk(thread, toBlock);
    }
    assert(walkPool[thread][toBlock].size_w < WALK_BUFFER_SIZE);
    walk = reEncode(walk, previousVertex, currentVertex, nAppendHops, preResideBlock);
    walkPool[thread][toBlock].push_back(walk);
}

void WalkManager::writeWalks2Disk(tid_t thread, bid_t block) {
//    m.start_time("4_writeWalks2Disk");
    std::string walksFile = walksname(baseFile, block );
    int f = open(walksFile.c_str(), O_WRONLY | O_CREAT | O_APPEND, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
    pwritea(f, &walkPool[thread][block][0], walkPool[thread][block].size_w * sizeof(WalkDataType) );
    nDiskWalks[block] += walkPool[thread][block].size_w;
    walkPool[thread][block].size_w = 0;
    close(f);
//    m.stop_time("4_writeWalks2Disk");
}

wid_t WalkManager::getNCurrentWalks(bid_t block) {
//    m.start_time("3_getcurrentwalks");
    currentWalks = (WalkDataType*)malloc(nWalks[block] * sizeof(WalkDataType));
    if(nDiskWalks[block] > 0){
        readWalksFromDisk(block);
    }
    wid_t count = nDiskWalks[block];
    for(tid_t t = 0; t < nThreads; t++){
        if(walkPool[t][block].size_w > 0){
            for(wid_t w = 0; w < walkPool[t][block].size_w; w++)
                currentWalks[count + w] = walkPool[t][block][w];
            count += walkPool[t][block].size_w;
            /* TODO: why set size_w to 0 here? */
            walkPool[t][block].size_w = 0;
        }
    }
    if (count != nWalks[block]) {
        logstream(LOG_DEBUG) << "read walks count = " << count << ", recorded nWalks[p] = " << nWalks[block] << ", nDiskWalks[block]" << nDiskWalks[block] << std::endl;
    }
    nDiskWalks[block] = 0;
//    m.stop_time("3_getcurrentwalks");
    return count;
}

void WalkManager::readWalksFromDisk(bid_t block) {
    m.start_time("z_w_readWalksfromDisk");

    std::string walksFile = walksname(baseFile, block);
    int f = open(walksFile.c_str(), O_RDWR, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
    if (f < 0) {
        logstream(LOG_FATAL) << "Could not load :" << walksFile << " error: " << strerror(errno) << std::endl;
    }
    assert(f > 0);
    /* read from file*/
    preada(f, &currentWalks[0], nDiskWalks[block] * sizeof(WalkDataType), 0);
    /* clear file */
    ftruncate(f,0);
    close(f);
    /* remove the walk file*/
    unlink(walksFile.c_str());

    m.stop_time("z_w_readWalksfromDisk");
}

void WalkManager::updateWalkNum(bid_t block) {
    m.start_time("6_updateWalkNum");
    wid_t forwardWalks = 0;
    for(bid_t b = 0; b < nBlocks; b++){
        if(isModified[b] && b != block){
            isModified[b] = false;
            wid_t newWalkNum = 0;
            newWalkNum = nDiskWalks[b];
            for(tid_t t = 0; t < nThreads; t++){
                newWalkNum += walkPool[t][b].size_w;
            }
            if(newWalkNum < nWalks[b]){
                logstream(LOG_DEBUG) << " b = " << b << ", newWalkNum = " << newWalkNum << ", nWalks[b] = " << nWalks[b] << std::endl;
                assert(false);
            }
            forwardWalks += newWalkNum - nWalks[b];
            nWalks[b] = newWalkNum;
        }
    }


    m.start_time("z_w_clear_curwalks");
    currentNWalks += forwardWalks;
    currentNWalks -= nWalks[block];
    nWalks[block] = 0;
    minStep[block] = 0xffff;
    free(currentWalks);
    currentWalks = NULL;
    m.stop_time("z_w_clear_curwalks");

    m.stop_time("6_updateWalkNum");

//    wid_t newNWalks = 0;
//    for (bid_t b = 0; b < nBlocks; b++){
//        if (isModified[b]){
//            isModified[b] = false;
//            wid_t newBlockWalkNum = 0;
//            newBlockWalkNum = nDiskWalks[b];
//            for(tid_t t = 0; t < nThreads; t++){
//                newBlockWalkNum += walkPool[t][b].size_w;
//            }
//            nWalks[b] = newBlockWalkNum;
//            newNWalks += newBlockWalkNum;
//        }else{
//            newNWalks += nWalks[b];
//        }
//    }
//    currentNWalks = newNWalks;
//    minStep[block] = 0xffff;
//    free(currentWalks);
//    currentWalks = nullptr;
}



bid_t WalkManager::blockWithMaxWalks() {
    wid_t maxw = 0, maxp = 0;
    for(bid_t p = 0; p < nBlocks; p++) {
        if(maxw < nWalks[p] ){
            maxw = nWalks[p];
            maxp = p;
        }
    }
    return maxp;
}

bid_t WalkManager::blockWithMinStep() {
    hid_t mins = 0xffff, minp = 0;
    for(bid_t p = 0; p < nBlocks; p++) {
        if(mins > minStep[p] ){
            mins = minStep[p];
            minp = p;
        }
    }
    if(nWalks[minp] > 0)
        return minp;
    return blockWithMaxWalks();
}

bid_t WalkManager::chooseBlock(float prob) {
    // return blockWithMaxWeight();//////////////
    float cc = ((float)rand())/RAND_MAX;
    if( cc < prob ){
        return blockWithMinStep();
    }
    return blockWithMaxWalks();
}

void WalkManager::setMinStep(bid_t p, hid_t hop) {
    if(minStep[p] > hop)
    {
#pragma omp critical
        {
            minStep[p] = hop;
        }
    }
}

wid_t WalkManager::getCurrentWalks(bid_t p) {
//	m.start_time("3_getcurrentwalks");
	currentWalks = (WalkDataType*)malloc(nWalks[p] * sizeof(WalkDataType));
	if(nDiskWalks[p] > 0){
		readWalksFromDisk(p);
	}
	wid_t count = nDiskWalks[p];
	for(tid_t t = 0; t < nThreads; t++){
		if(walkPool[t][p].size_w > 0){
			for(wid_t w = 0; w < walkPool[t][p].size_w; w++)
				currentWalks[count + w] = walkPool[t][p][w];
			count += walkPool[t][p].size_w;
			walkPool[t][p].size_w = 0;
		}
	}
	if (count != nWalks[p]) {
		logstream(LOG_DEBUG) << "read walks count = " << count << ", recorded nWalks[p] = " << nWalks[p] << ", disk nWalks[p]" << nDiskWalks[p] << std::endl;
	}
	nDiskWalks[p] = 0;
//	m.stop_time("3_getcurrentwalks");
	curBlockNWalks = count;
	return count;
}


#endif