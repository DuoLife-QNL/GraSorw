//
// Created by Ethan on 2020/12/28.
//

#ifndef IOE_SORW_VERTEXINFO_H
#define IOE_SORW_VERTEXINFO_H


class VertexInfo{
public:
    /**
     * TODO: should set vertexId to private in the future
     * so that only reset() is able to change the vertex info
     */
    vid_t vertexId;
    vid_t *csr;
    vid_t outDegree;
    bid_t resideBlockId= INVALID_BID;
private:
    /* 当不知道是否在内存中时，缺省填false */
    bool resideBlockInMem = false;
public:
    bool isResideBlockInMem() const {
        return resideBlockInMem;
    }

public:
    void setResideBlockInMem(const bool& _resideBlockInMem) {
        resideBlockInMem = _resideBlockInMem;
        if (resideBlockInMem){
            delAfterUse = false;
        }else{
            delAfterUse = true;
        }
    }

    bool delAfterUse = true;

    VertexInfo(){
        vertexId = INVALID_VID;
        csr = nullptr;
        outDegree = 0;
    };

    explicit VertexInfo(vid_t _vertexId){
        vertexId = _vertexId;
        csr = nullptr;
        outDegree = 0;
    }

    VertexInfo(vid_t _vertexId, bid_t _resideBlockId){
        vertexId = _vertexId;
        csr = nullptr;
        outDegree = 0;
        resideBlockId = _resideBlockId;
    }

    VertexInfo(vid_t _vertexId, bid_t _resideBlockId, bool _resideBlockInMem){
        vertexId = _vertexId;
        csr = nullptr;
        outDegree = 0;
        resideBlockId = _resideBlockId;
        resideBlockInMem = _resideBlockInMem;
        if (resideBlockInMem){
            delAfterUse = false;
        }
    }

    void resetVertex(vid_t v){
        if (csr){
            if (delAfterUse){
                delete csr;
            }
        }
        csr = nullptr;
        vertexId = v;
        outDegree = 0;
        resideBlockId = INVALID_BID;
        resideBlockInMem = false;
        delAfterUse = false;
    }

    void resetVertex(vid_t _vertexId, bid_t _resideBlockId){
        if (delAfterUse){
            delete csr;
        }
        vertexId = _vertexId;
        outDegree = 0;
        resideBlockId = _resideBlockId;
        resideBlockInMem = false;
        delAfterUse = false;
    }

    void resetVertex(vid_t _vertexId, bid_t _resideBlockId, bool _resideBlockInMem){
        if (delAfterUse){
            delete csr;
        }
        csr = nullptr;
        vertexId = _vertexId;
        outDegree = 0;
        resideBlockId = _resideBlockId;
        resideBlockInMem = _resideBlockInMem;
        delAfterUse = false;
        if (!resideBlockInMem){
            delAfterUse = true;
        }
    }

    virtual ~VertexInfo() {
        if (delAfterUse){
            delete csr;
        }
    };
};

bool sortByVertexDegreeDes(VertexInfo v1, VertexInfo v2){
    return v1.outDegree > v2.outDegree;
}


#endif //IOE_SORW_VERTEXINFO_H
