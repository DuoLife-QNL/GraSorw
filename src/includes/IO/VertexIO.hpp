//
// Created by lihz on 2020/11/4.
//

#ifndef IOE_SORW_VERTEXIO_HPP
#define IOE_SORW_VERTEXIO_HPP

#define TWOLEVELCACHE 0
#define METRICS 0
#define BUFFER_SIZE 300

typedef enum GET_CSR_METHOD{
    DISK,
    MEMORY,
    ON_DEMAND_CACHE,
    AUTO
}GET_CSR_METHOD;

#include <map>
#include "engine/Settings.hpp"
#include "metrics/metrics.hpp"
#include "api/io.hpp"
#include "IO/VertexInfo.hpp"

#include "IO/VertexBuffer.hpp"

class CacheNode{
public:
    CacheNode(vid_t _vertex, eid_t _outDegree, vid_t *_csr);
    ~CacheNode();
    vid_t vertex;
    eid_t outDegree;
    vid_t *csr;
    CacheNode *pre, *next;
};

CacheNode::CacheNode(vid_t _vertex, eid_t _outDegree, vid_t *_csr) {
    vertex = _vertex;
    outDegree = _outDegree;
    csr = _csr;
    pre = nullptr;
    next = nullptr;
}

CacheNode::~CacheNode() {
    delete [] csr;
}


class CSRBuffer{
private:
    vid_t bufferSize = 0;
    CacheNode *head = nullptr;
    CacheNode *tail = nullptr;
    std::map<vid_t, CacheNode *> mp;

    void remove(CacheNode *node);
    void setHead(CacheNode *node);

public:
    CSRBuffer();
    explicit CSRBuffer(vid_t _bufferSize);
    void setBufferSize(vid_t _bufferSize);
    vid_t *get(vid_t vertex, eid_t &outDegree);
    void add(vid_t vertex, vid_t *csr, vid_t outDegree);
};

CSRBuffer::CSRBuffer() = default;

CSRBuffer::CSRBuffer(vid_t _bufferSize) {
    bufferSize = _bufferSize;
}

void CSRBuffer::setBufferSize(vid_t _bufferSize) {
    bufferSize = _bufferSize;
}

void CSRBuffer::remove(CacheNode *node) {
    if (node->pre != nullptr){
        node->pre->next = node->next;
    }else{
        head = node->next;
    }
    if (node->next != nullptr){
        node->next->pre = node->pre;
    }else{
        tail = node->pre;
    }
}

void CSRBuffer::setHead(CacheNode *node) {
    node->next = head;
    node->pre = nullptr;
    if (head != nullptr){
        head->pre = node;
    }
    head = node;
    if (tail == nullptr){
        tail = head;
    }
}

vid_t *CSRBuffer::get(vid_t vertex, eid_t &outDegree) {
    auto it = mp.find(vertex);
    if (it != mp.end()){
        CacheNode *node = it->second;
        remove(node);
        setHead(node);
        outDegree = node->outDegree;
        return node->csr;
    }else{
        return nullptr;
    }
}

void CSRBuffer::add(vid_t vertex, vid_t *csr, vid_t outDegree) {
    auto *newNode = new CacheNode(vertex, outDegree, csr);
    if (mp.size() >= bufferSize){
        auto it = mp.find(tail->vertex);
        delete it->second;
        remove(tail);
        mp.erase(it);
    }
    setHead(newNode);
    mp[vertex] = newNode;
}

#if false
class VertexIO{
private:
    int begPosFile, csrFile;
    bid_t nBlocks;
    vid_t *startVertex;
    vid_t nTotalVertex;
    CSRBuffer csrBuffer;

    /**
     * If vertexIO has not stored CSR of any vertex, or blocks in memory
     * have been changed, flag should be set to false, meaning that the
     * outDegree and *csr here is not reliable and we have to get vertex
     * CSR from disk or block memory. After that, we set flag to true.
     */
    bool flag = false;

    /* buffer for getVertexCSR */
    vid_t currentVertex = 0;
    vid_t *csr = nullptr;
#if !TWOLEVELCACHE
    vid_t *constCSR = nullptr;
#endif
    vid_t outDegree = 0;

    /* find csr segment directly from memory */
    bid_t nInMemBlocks = 0;
    vid_t **csrBuf = nullptr;
    eid_t **begPosBuf = nullptr;
    bid_t *inMemIndex = nullptr;

public:
    VertexIO(int _begPosFile, int _csrFile, bid_t _nBlocks, vid_t *_startVertex, vid_t _nTotalVertices, vid_t maxOutDegree);

    void setIOInfo(bid_t _nInMemBlocks, vid_t **_csrBuf, eid_t **_begPosBuf, bid_t *_inMemIndex);

    vid_t *getVertexCSR(vid_t vertex, vid_t &_outDegree, metrics &m, bool useBlockInMem, bool preV);

    vid_t getOutDegree(vid_t vertex) const;

    vid_t getNTotalVertices() const;

    bid_t getBlock(vid_t v);

    bool blockInMem(vid_t vertex);

};

VertexIO::VertexIO(int _begPosFile, int _csrFile, bid_t _nBlocks, vid_t *_startVertex, vid_t _nTotalVertices, vid_t maxOutDegree)
        : csrBuffer(BUFFER_SIZE) {
    begPosFile = _begPosFile;
    csrFile = _csrFile;
    nBlocks = _nBlocks;
    startVertex = _startVertex;
    nTotalVertex = _nTotalVertices;
#if !TWOLEVELCACHE
    csr = new vid_t [maxOutDegree];
    constCSR = csr;
#endif
}

#if TWOLEVELCACHE
vid_t * VertexIO::getVertexCSR(vid_t vertex, vid_t &_outDegree, metrics &m,
                               bool useBlockInMem = true,
                               bool preV = false) {

    bid_t resideBlock = getBlock(vertex);
    /* hit first level cache */
    if (flag && vertex == currentVertex){
#if METRICS
        m.add("hit0", 1);
#endif
        _outDegree = outDegree;
        return csr;
    }

    currentVertex = vertex;
    if (!useBlockInMem || !inMemIndex ||  inMemIndex[resideBlock] == nInMemBlocks){
        /* block not in memory, or we don't use blocks in memory */
//        if (preV){
//            m.add("getPreDisk", 1);
//        }
        eid_t *begPos;
        /* check Buffer */
#if METRICS
        m.start_time("getVertexCSR_Cache");
#endif
        csr = csrBuffer.get(currentVertex, outDegree);
#if METRICS
        m.stop_time("getVertexCSR_Cache");
#endif
        if (csr == nullptr){
            /**
             * eid_t[0]: start edge ID of vertex "vertex"
             * eid_t[1]: start edge ID of vertex (vertex + 1)
             * so eid_[1] - eid[0] = out degree of vertex "vertex"
             */
            begPos = new eid_t[2];
#if METRICS
            m.start_time("getVertexCSR_Disk");
#endif
            preada(begPosFile, begPos, 2 * sizeof(eid_t), vertex * sizeof(eid_t));
            outDegree = begPos[1] - begPos[0];
            auto *vertexCSR = new vid_t[outDegree];
            preada(csrFile, vertexCSR, outDegree * sizeof(vid_t), begPos[0] * sizeof(vid_t));
            delete [] begPos;
            csr = vertexCSR;
#if METRICS
            m.stop_time("getVertexCSR_Disk");
#endif

//            m.start_time("maintainCache");
            csrBuffer.add(currentVertex, vertexCSR, outDegree);
//            m.stop_time("maintainCache");
        }else{
//            m.add("hit", 1);
        }
    }else{  //get the csr from block info in memory
        if (preV){
#if METRICS
            m.add("getPreMemory", 1);
#endif
        }
        assert(csr);
//        m.start_time("getVertexCSR_MemoryBlock");
        vid_t vertexOffset = vertex - startVertex[resideBlock];
        bid_t blockOffset = inMemIndex[resideBlock];
        eid_t *blockBegPos = begPosBuf[blockOffset];
        vid_t *blockCSR = csrBuf[blockOffset];
        eid_t verOutOffset = blockBegPos[vertexOffset] - blockBegPos[0];
        csr = blockCSR + verOutOffset;
        outDegree = blockBegPos[vertexOffset + 1] - blockBegPos[vertexOffset];
//        m.stop_time("getVertexCSR_MemoryBlock");
    }
    flag = true;
    _outDegree = outDegree;
    return csr;
}

#else
vid_t * VertexIO::getVertexCSR(vid_t vertex, vid_t &_outDegree, metrics &m,
                               bool useBlockInMem = true,
                               bool preV = false) {
    /* hit first level cache */
    if (flag && vertex == currentVertex){
//        m.add("hit0", 1);
        _outDegree = outDegree;
        return csr;
    }

    bid_t resideBlock = getBlock(vertex);

    currentVertex = vertex;
    if (!useBlockInMem || !inMemIndex ||  inMemIndex[resideBlock] == nInMemBlocks){
        /* block not in memory, or we don't use blocks in memory */
        if (preV){
#if METRICS
            m.add("getPreDisk", 1);
#endif
        }
        eid_t *begPos;
        /**
         * eid_t[0]: start edge ID of vertex "vertex"
         * eid_t[1]: start edge ID of vertex (vertex + 1)
         * so eid_[1] - eid[0] = out degree of vertex "vertex"
         */
        begPos = new eid_t[2];
#if METRICS
            m.start_time("getVertexCSR_Disk");
#endif
        preada(begPosFile, begPos, 2 * sizeof(eid_t), vertex * sizeof(eid_t));
        outDegree = begPos[1] - begPos[0];
        csr = constCSR;
        preada(csrFile, csr, outDegree * sizeof(vid_t), begPos[0] * sizeof(vid_t));
        delete [] begPos;
#if METRICS
            m.stop_time("getVertexCSR_Disk");
#endif
    }else{  //get the csr from block info in memory
        if (preV){
#if METRICS
            m.add("getPreMemory", 1);
#endif
        }
        assert(csr);
#if METRICS
        m.start_time("getVertexCSR_MemoryBlock");
#endif
        vid_t vertexOffset = vertex - startVertex[resideBlock];
        bid_t blockOffset = inMemIndex[resideBlock];
        eid_t *blockBegPos = begPosBuf[blockOffset];
        vid_t *blockCSR = csrBuf[blockOffset];
        eid_t verOutOffset = blockBegPos[vertexOffset] - blockBegPos[0];
        csr = blockCSR + verOutOffset;
        outDegree = blockBegPos[vertexOffset + 1] - blockBegPos[vertexOffset];
#if METRICS
        m.stop_time("getVertexCSR_MemoryBlock");
#endif
    }
    flag = true;
    _outDegree = outDegree;
    return csr;
}
#endif

void VertexIO::setIOInfo(bid_t _nInMemBlocks, vid_t **_csrBuf, eid_t **_begPosBuf, bid_t *_inMemIndex) {
    this->nInMemBlocks = _nInMemBlocks;
    this->csrBuf = _csrBuf;
    this->begPosBuf = _begPosBuf;
    this->inMemIndex = _inMemIndex;
    flag = false;
}

vid_t VertexIO::getOutDegree(vid_t vertex) const {
    if (flag && vertex == currentVertex){
        return outDegree;
    }else{
        eid_t *begPos;
        /**
         * eid_t[0]: start edge ID of vertex "vertex"
         * eid_t[1]: start edge ID of vertex (vertex + 1)
         * so eid_[1] - eid[0] = out degree of vertex "vertex"
         */
        begPos = new eid_t[2];
        preada(begPosFile, begPos, 2 * sizeof(eid_t), vertex * sizeof(eid_t));
        vid_t result = begPos[1] - begPos[0];
        delete []begPos;
        return result;
    }
}

vid_t VertexIO::getNTotalVertices() const {
    return nTotalVertex;
}

bid_t VertexIO::getBlock(vid_t v) {
    for (bid_t i = 0; i < nBlocks; i++){
        if (startVertex[i + 1] > v){
            return i;
        }
    }
    logstream(LOG_ERROR) << "Vertex" << v << "out of range!" << std::endl;
    abort();
}

bool VertexIO::blockInMem(vid_t vertex) {
    if (!inMemIndex ||  inMemIndex[getBlock(vertex)] == nInMemBlocks){
        return false;
    }else{
        return true;
    }
}

#else

class VertexIO{
private:
    int begPosFile, csrFile;
    bid_t nBlocks;
    vid_t *startVertex;
    vid_t nTotalVertex;
#if CACHE
    CSRBuffer csrBuffer;
#endif
    /* find csr segment directly from memory */
    bid_t nInMemBlocks = 0;
    vid_t **csrBuf = nullptr;
    eid_t **begPosBuf = nullptr;
    bid_t *inMemIndex = nullptr;

    vid_t *preBlockCSR = nullptr;
    eid_t *preBlockBegPos = nullptr;
    bid_t preBlockId = INVALID_BID;

    /* used for on-demand loading */
    bool *vertexMap;
        /* dynamic block on-demand loading */
    VertexInfo **onDemandVertexInfo_dynamic;
            /* 标记按需加载CSR链表已清空，防止内存泄漏 */
    bool onDemandCSRCleared_dynamic = true;
    bid_t onDemandLoadBlockId_dynamic = INVALID_BID;
        /* static block on-demand loading*/
    VertexInfo **onDemandVertexInfo_static;
    bool onDemandCSRCleared_static = true;
    bid_t onDemandLoadBlockId_static = INVALID_BID;

#if BLOCK_UTI |ACT_VER
    /*
     * vertex map for block utilization statistics
     * used to mark which vertices are used
     */
    bool *vertexMap_bu;
#endif

public:
    std::vector<metrics> metric;

    VertexBuffer vertexBuffer;

    VertexIO(int _begPosFile, int _csrFile, bid_t _nBlocks, vid_t *_startVertex, vid_t _nTotalVertices);

    void setIOInfo(bid_t _nInMemBlocks, vid_t **_csrBuf, eid_t **_begPosBuf, bid_t *_inMemIndex);

    void getVertexCSR(VertexInfo &v, metrics &m, GET_CSR_METHOD getCsrMethod = AUTO);

    void getVertexOutDegree(VertexInfo &v, GET_CSR_METHOD getCsrMethod);

    void getPreVertexCSR(VertexInfo &preV, metrics &m);

    vid_t getNTotalVertices() const;

    bid_t getBlock(vid_t v);

    bool blockInMem(VertexInfo &v);

    bool blockInMem(bid_t &blockId);

    void setPreBlockInfo(vid_t *_preBlockCSR, eid_t *_preBlockBegPos, bid_t _preBlockId);

    void getPreBlockInfoFromMem(bid_t _preBlockId);

#if BLOCK_UTI | ACT_VER
    void clearBlockMap(bid_t blockId){
        memset(vertexMap_bu + getStartVertexId(blockId), false, getBlockVertexNum(blockId) * sizeof(false));
    }

    void markActivated(vid_t vertexId){
        vertexMap_bu[vertexId] = true;
    }

    eid_t getUsedEdgesNum(bid_t blockId){
        eid_t usedEdgeNum = 0;
        vid_t startVertexId = getStartVertexId(blockId);
        vid_t endVertexId = getStartVertexId(blockId + 1);
        for (vid_t vertexId = startVertexId; vertexId < endVertexId; vertexId++){
            if (!vertexMap_bu[vertexId]){
                continue;
            }
            VertexInfo v(vertexId);
            getVertexOutDegree(v, AUTO);
            usedEdgeNum += v.outDegree;
        }
        return usedEdgeNum;
    }

    vid_t getActivatedVerticesNum(bid_t blockId){
        vid_t activedVertices = 0;
        vid_t startVertexId = getStartVertexId(blockId);
        vid_t endVertexId = getStartVertexId(blockId + 1);
        for (vid_t vertexId = startVertexId; vertexId < endVertexId; vertexId++){
            if (!vertexMap_bu[vertexId]){
                continue;
            }
            activedVertices++;
        }
        return activedVertices;
    }
#endif


    bool vertexInBlock(vid_t vertexId, bid_t blockId){
        return vertexId >= startVertex[blockId] && vertexId < startVertex[blockId + 1];
    }

    void getBlockSize(bid_t blockId, vid_t &nVertices, eid_t &nEdges){
        auto begPos_s = new eid_t;
        auto begPos_e = new eid_t;
        vid_t vertexId_s = startVertex[blockId];
        vid_t vertexId_e = startVertex[blockId + 1];
        preada(begPosFile, begPos_s, sizeof(eid_t), vertexId_s * sizeof(eid_t));
        preada(begPosFile, begPos_e, sizeof(eid_t), vertexId_e * sizeof(eid_t));
        nVertices = vertexId_e - vertexId_s;
        nEdges = *begPos_e - *begPos_s;
    }

    void needVertex(const vid_t &v){
        vertexMap[v] = true;
    }

    /* @dynamicBlock to mark whether it is on-demand loading the dynamic block or the static block */
    void getVertexOnDemand(const bid_t &b, bool dynamicBlock = true){
        if ((dynamicBlock && !onDemandCSRCleared_dynamic) || (!dynamicBlock && !onDemandCSRCleared_static)){
            assert(false && "Clear on demand CSR table first!");
        }
        if (dynamicBlock){
            vid_t vs = startVertex[b];
            vid_t ve = startVertex[b + 1];
            vid_t blockVertexNum = ve - vs;
            onDemandVertexInfo_dynamic = new VertexInfo *[blockVertexNum];
            /* v: vertex offset to vs */
            for (vid_t v = 0; v < blockVertexNum; v++){
                onDemandVertexInfo_dynamic[v] = nullptr;
                if (vertexMap[v + vs]){
                    auto *pVertexInfo = new VertexInfo(v + vs, b, false);
                    getVertexCSR(*pVertexInfo, metric[0], DISK);
                    onDemandVertexInfo_dynamic[v] = pVertexInfo;
                }
            }
            onDemandCSRCleared_dynamic = false;
            onDemandLoadBlockId_dynamic = b;
        }else{
            vid_t vs = startVertex[b];
            vid_t ve = startVertex[b + 1];
            vid_t blockVertexNum = ve - vs;
            onDemandVertexInfo_static = new VertexInfo *[blockVertexNum];
            /* v: vertex offset to vs */
            for (vid_t v = 0; v < blockVertexNum; v++){
                onDemandVertexInfo_static[v] = nullptr;
                if (vertexMap[v + vs]){
                    auto *pVertexInfo = new VertexInfo(v + vs, b, false);
                    getVertexCSR(*pVertexInfo, metric[0], DISK);
                    onDemandVertexInfo_static[v] = pVertexInfo;
                }
            }
            onDemandCSRCleared_static = false;
            onDemandLoadBlockId_static = b;
        }

    }
    /* @dynamicBlock to mark whether it is on-demand loading the dynamic block or the static block */
    void clearVertexMapAndCSR(const bid_t &b, bool dynamicBlock = true){
        if (dynamicBlock){
            vid_t vs = startVertex[onDemandLoadBlockId_dynamic];
            vid_t ve = startVertex[onDemandLoadBlockId_dynamic + 1];
            vid_t blockVertexNum = ve - vs;
            memset(vertexMap + vs, false, blockVertexNum * sizeof(bool));
            for (vid_t v = 0; v < blockVertexNum; v++){
                if (onDemandVertexInfo_dynamic[v]){
                    delete onDemandVertexInfo_dynamic[v];
                }
            }
            delete onDemandVertexInfo_dynamic;
            onDemandCSRCleared_dynamic = true;
            onDemandLoadBlockId_dynamic = INVALID_BID;
        }else{
            vid_t vs = startVertex[onDemandLoadBlockId_static];
            vid_t ve = startVertex[onDemandLoadBlockId_static + 1];
            vid_t blockVertexNum = ve - vs;
            memset(vertexMap + vs, false, blockVertexNum * sizeof(bool));
            for (vid_t v = 0; v < blockVertexNum; v++){
                if (onDemandVertexInfo_static[v]){
                    delete onDemandVertexInfo_static[v];
                }
            }
            delete onDemandVertexInfo_static;
            onDemandCSRCleared_static = true;
            onDemandLoadBlockId_static = INVALID_BID;
        }
    }

    vid_t getStartVertexId(bid_t blockId){
        return startVertex[blockId];
    }

    vid_t getBlockVertexNum(bid_t blockId){
        return startVertex[blockId + 1] - startVertex[blockId];
    }
};

VertexIO::VertexIO(int _begPosFile, int _csrFile, bid_t _nBlocks, vid_t *_startVertex, vid_t _nTotalVertices): vertexBuffer(_nTotalVertices) {
    begPosFile = _begPosFile;
    csrFile = _csrFile;
    nBlocks = _nBlocks;
    startVertex = _startVertex;
    nTotalVertex = _nTotalVertices;

    tid_t nThreads = omp_get_max_threads();
    for (tid_t t = 0; t < nThreads; t++){
        metrics mt("thread" + std::to_string(t));
        metric.emplace_back(mt);
    }
    vertexMap = new bool[nTotalVertex];
#if CACHE
    csrBuffer.setBufferSize(blockSize_kb_global * 1024 / sizeof(vid_t));
#endif
#if BLOCK_UTI | ACT_VER
    vertexMap_bu = new bool [getNTotalVertices()];
    memset(vertexMap_bu, false, getNTotalVertices() * sizeof(bool));
#endif
}

void VertexIO::setIOInfo(bid_t _nInMemBlocks, vid_t **_csrBuf, eid_t **_begPosBuf, bid_t *_inMemIndex) {
    this->nInMemBlocks = _nInMemBlocks;
    this->csrBuf = _csrBuf;
    this->begPosBuf = _begPosBuf;
    this->inMemIndex = _inMemIndex;
}

vid_t VertexIO::getNTotalVertices() const {
    return nTotalVertex;
}

bid_t VertexIO::getBlock(vid_t v) {
    for (bid_t i = 0; i < nBlocks; i++){
        if (startVertex[i + 1] > v){
            return i;
        }
    }
    logstream(LOG_ERROR) << "Vertex" << v << "out of range!" << std::endl;
    abort();
}

bool VertexIO::blockInMem(VertexInfo &v) {
    if (v.resideBlockId == INVALID_BID){
        v.resideBlockId = getBlock(v.vertexId);
    }
    if (!inMemIndex || inMemIndex[v.resideBlockId] == nInMemBlocks){
        v.setResideBlockInMem(false);
        return false;
    }else{
        v.setResideBlockInMem(true);
        return true;
    }
}
#if CACHE
void VertexIO::getVertexCSR(VertexInfo &v, metrics &m, GET_CSR_METHOD getCsrMethod) {
    /* 判断CSR是否存在的逻辑放到外部，只有需要CSR的时候才会获取 */
    if (v.csr){
        abort();
    }
    vid_t *csr;
    eid_t outDegree;
    if (v.vertexId == INVALID_VID){
        abort();
    }

    if ((v.isResideBlockInMem() || blockInMem(v))) {  //get the csr from block info in memory
        bid_t resideBlock = v.resideBlockId;
        vid_t vertexOffset = v.vertexId - startVertex[resideBlock];
        bid_t blockOffset = inMemIndex[resideBlock];
        eid_t *blockBegPos = begPosBuf[blockOffset];
        vid_t *blockCSR = csrBuf[blockOffset];
        eid_t verOutOffset = blockBegPos[vertexOffset] - blockBegPos[0];
        csr = blockCSR + verOutOffset;
        outDegree = blockBegPos[vertexOffset + 1] - blockBegPos[vertexOffset];
        v.delAfterUse = false;

    }else{
        /* block not in memory */
        csr = csrBuffer.get(v.vertexId, outDegree);
        if (csr == nullptr){
            /* vertex info not in cache */

            eid_t *begPos;
            /**
             * eid_t[0]: start edge ID of vertex "vertex"
             * eid_t[1]: start edge ID of vertex (vertex + 1)
             * so eid_[1] - eid[0] = out degree of vertex "vertex"
             */
            begPos = new eid_t[2];
            preada(begPosFile, begPos, 2 * sizeof(eid_t), v.vertexId * sizeof(eid_t));
            outDegree = begPos[1] - begPos[0];
            csr = new vid_t[outDegree];
            preada(csrFile, csr, outDegree * sizeof(vid_t), begPos[0] * sizeof(vid_t));
            delete[] begPos;
            csrBuffer.add(v.vertexId, csr, outDegree);
        }
        v.delAfterUse = false;
    }

    v.csr = csr;
    v.outDegree = outDegree;

}
#else
void VertexIO::getVertexCSR(VertexInfo &v, metrics &m, GET_CSR_METHOD getCsrMethod) {
    if (v.vertexId == 548142){
        int a = 1;
    }
    if (v.csr){
        abort();
    }
    vid_t *csr;
    eid_t outDegree;
    if (v.vertexId == INVALID_VID){
        abort();
    }

    if ((v.isResideBlockInMem() || blockInMem(v)) && (getCsrMethod == AUTO || getCsrMethod == MEMORY)) {  //get the csr from block info in memory
#if METRICS
        m.start_time("getVertexCSR_MemoryBlock");
#endif
        bid_t resideBlock = v.resideBlockId;
        vid_t vertexOffset = v.vertexId - startVertex[resideBlock];
        bid_t blockOffset = inMemIndex[resideBlock];
        eid_t *blockBegPos = begPosBuf[blockOffset];
        vid_t *blockCSR = csrBuf[blockOffset];
        eid_t verOutOffset = blockBegPos[vertexOffset] - blockBegPos[0];
        csr = blockCSR + verOutOffset;
        outDegree = blockBegPos[vertexOffset + 1] - blockBegPos[vertexOffset];
        v.delAfterUse = false;
#if METRICS
        m.stop_time("getVertexCSR_MemoryBlock");
#endif
    }else if(getCsrMethod == AUTO &&
                (
                 /* vertex is in dynamic block and has been on-demand loaded into memory */
                 (
                  !onDemandCSRCleared_dynamic &&
                  v.resideBlockId == onDemandLoadBlockId_dynamic &&
                  onDemandVertexInfo_dynamic[v.vertexId - startVertex[onDemandLoadBlockId_dynamic]]
                 )
                 ||
                 /* vertex is in static block and has been on-demand loaded into memory */
                 (
                  !onDemandCSRCleared_static &&
                  v.resideBlockId == onDemandLoadBlockId_static &&
                  onDemandVertexInfo_static[v.vertexId - startVertex[onDemandLoadBlockId_static]]
                 )
                )
            ){
        if (v.resideBlockId == onDemandLoadBlockId_dynamic){
            // get from on demand CSR cache (dynamic block)
            vid_t vertexOffset = v.vertexId - startVertex[onDemandLoadBlockId_dynamic];
            assert(v.vertexId == onDemandVertexInfo_dynamic[vertexOffset]->vertexId);
            csr = onDemandVertexInfo_dynamic[vertexOffset]->csr;
            outDegree = onDemandVertexInfo_dynamic[vertexOffset]->outDegree;
            v.delAfterUse = false;
        }else{
            vid_t vertexOffset = v.vertexId - startVertex[onDemandLoadBlockId_static];
            assert(v.vertexId == onDemandVertexInfo_static[vertexOffset]->vertexId);
            csr = onDemandVertexInfo_static[vertexOffset]->csr;
            outDegree = onDemandVertexInfo_static[vertexOffset]->outDegree;
            v.delAfterUse = false;
        }

    }else{
        /* block not in memory */

#if PLAIN & STATICCACHE
        #if METRICS
        m.start_time("hit-buffer");
        if (vertexBuffer.getVertex(v)){
            m.stop_time("hit-buffer");
            return;
        }
        m.stop_time("hit-buffer");
#else
        if (vertexBuffer.getVertex(v)){
//            tid_t t = omp_get_thread_num();
//            if (t == 1){
//                m.add("get-csr-cache-th-1", 1);
//            }else if (t == 2){
//                m.add("get-csr-cache-th-2", 1);
//            }
//            return;
        }
#endif
#endif

        eid_t *begPos;
        /**
         * eid_t[0]: start edge ID of vertex "vertex"
         * eid_t[1]: start edge ID of vertex (vertex + 1)
         * so eid_[1] - eid[0] = out degree of vertex "vertex"
         */
        begPos = new eid_t[2];
#if METRICS
        m.start_time("getVertexCSR_Disk");
#endif
        preada(begPosFile, begPos, 2 * sizeof(eid_t), v.vertexId * sizeof(eid_t));
        outDegree = begPos[1] - begPos[0];
        csr = new vid_t[outDegree];
        preada(csrFile, csr, outDegree * sizeof(vid_t), begPos[0] * sizeof(vid_t));
        delete[] begPos;
        v.delAfterUse = true;
#if METRICS
        m.stop_time("getVertexCSR_Disk");
#endif
//        tid_t t = omp_get_thread_num();
//        if (t == 1){
//            m.add("get-csr-disk-th-1", 1);
//        }else if (t == 2){
//            m.add("get-csr-disk-th-2", 1);
//        }
    }

    v.csr = csr;
    v.outDegree = outDegree;

}

/* only fill the outDegree in VertexInfo, the csr in VertexInfo remains NULL */
void VertexIO::getVertexOutDegree(VertexInfo &v, GET_CSR_METHOD getCsrMethod){
    if (v.csr){
        abort();
    }
    eid_t outDegree;
    if (v.vertexId == INVALID_VID){
        abort();
    }

    if ((v.isResideBlockInMem() || blockInMem(v)) && (getCsrMethod == AUTO || getCsrMethod == MEMORY)) {  //get the csr from block info in memory
        bid_t resideBlock = v.resideBlockId;
        vid_t vertexOffset = v.vertexId - startVertex[resideBlock];
        bid_t blockOffset = inMemIndex[resideBlock];
        eid_t *blockBegPos = begPosBuf[blockOffset];
        eid_t verOutOffset = blockBegPos[vertexOffset] - blockBegPos[0];
        outDegree = blockBegPos[vertexOffset + 1] - blockBegPos[vertexOffset];
    }else if(getCsrMethod == AUTO &&
             (
                     /* vertex is in dynamic block and has been on-demand loaded into memory */
                     (
                             !onDemandCSRCleared_dynamic &&
                             v.resideBlockId == onDemandLoadBlockId_dynamic &&
                             onDemandVertexInfo_dynamic[v.vertexId - startVertex[onDemandLoadBlockId_dynamic]]
                     )
                     ||
                     /* vertex is in static block and has been on-demand loaded into memory */
                     (
                             !onDemandCSRCleared_static &&
                             v.resideBlockId == onDemandLoadBlockId_static &&
                             onDemandVertexInfo_static[v.vertexId - startVertex[onDemandLoadBlockId_static]]
                     )
             )
            ) {
        if (v.resideBlockId == onDemandLoadBlockId_dynamic) {
            // get from on demand CSR cache (dynamic block)
            vid_t vertexOffset = v.vertexId - startVertex[onDemandLoadBlockId_dynamic];
            assert(v.vertexId == onDemandVertexInfo_dynamic[vertexOffset]->vertexId);
            outDegree = onDemandVertexInfo_dynamic[vertexOffset]->outDegree;
            v.delAfterUse = false;
        } else {
            vid_t vertexOffset = v.vertexId - startVertex[onDemandLoadBlockId_static];
            assert(v.vertexId == onDemandVertexInfo_static[vertexOffset]->vertexId);
            outDegree = onDemandVertexInfo_static[vertexOffset]->outDegree;
            v.delAfterUse = false;
        }
    }else{
        /* block not in memory */

        eid_t *begPos;
        /**
         * eid_t[0]: start edge ID of vertex "vertex"
         * eid_t[1]: start edge ID of vertex (vertex + 1)
         * so eid_[1] - eid[0] = out degree of vertex "vertex"
         */
        begPos = new eid_t[2];
        preada(begPosFile, begPos, 2 * sizeof(eid_t), v.vertexId * sizeof(eid_t));
        outDegree = begPos[1] - begPos[0];
        delete[] begPos;
        v.delAfterUse = true;
    }
    v.outDegree = outDegree;
}
#endif

void VertexIO::getPreVertexCSR(VertexInfo &preV, metrics &m) {
    tid_t t = omp_get_thread_num();
    if (preBlockId != INVALID_BID && preV.vertexId >= startVertex[preBlockId] && preV.vertexId < startVertex[preBlockId + 1]){
//        metric.at(t).add("get-pre-csr-pre-block", 1);
        eid_t VertexOffset = preV.vertexId - startVertex[preBlockId];
        vid_t csrOffset = preBlockBegPos[VertexOffset] - preBlockBegPos[0];
        preV.csr = preBlockCSR + csrOffset;
        preV.outDegree = preBlockBegPos[VertexOffset + 1] - preBlockBegPos[VertexOffset];
        preV.delAfterUse = false;
    }else{
//        metric.at(t).add("get-pre-csr-general", 1);
        getVertexCSR(preV, m);
    }
}

void VertexIO::setPreBlockInfo(vid_t *_preBlockCSR, eid_t *_preBlockBegPos, bid_t _preBlockId) {
    preBlockCSR = _preBlockCSR;
    preBlockBegPos = _preBlockBegPos;
    preBlockId = _preBlockId;
}

bool VertexIO::blockInMem(bid_t &blockId) {
    if (blockId == INVALID_BID || !inMemIndex || inMemIndex[blockId] == nInMemBlocks){
        return false;
    }else{
        return true;
    }
}

void VertexIO::getPreBlockInfoFromMem(bid_t _preBlockId) {
    assert(blockInMem(_preBlockId));
    preBlockId = _preBlockId;
    bid_t blockOffset = inMemIndex[preBlockId];
    preBlockBegPos = begPosBuf[blockOffset];
    preBlockCSR = csrBuf[blockOffset];
}

#endif
#endif //IOE_SORW_VERTEXIO_HPP
