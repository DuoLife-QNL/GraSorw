//
// Created by Ethan on 2020/12/28.
//

#ifndef IOE_SORW_VERTEXBUFFER_HPP
#define IOE_SORW_VERTEXBUFFER_HPP

#include "IO/VertexInfo.hpp"
#include "engine/Settings.hpp"
#include <algorithm>
vid_t staticBufferSize = 0;

class VertexBuffer{
private:
    std::vector<bool> vertexMap;
    std::vector<VertexInfo> buffer;
    vid_t minDegree = 0xffffffff;
    size_t minDegPos = 0;
    vid_t addedVertices = 0;

    void updateMin(){
        minDegree = 0xffffffff;
        vid_t thisOutDegree;
        for (size_t i = 0; i < buffer.size(); i++){
            thisOutDegree = buffer.at(i).outDegree;
            if (thisOutDegree < minDegree){
                minDegree = thisOutDegree;
                minDegPos = i;
            }
        }
    }

    void insert(const VertexInfo& v){
        buffer.at(minDegPos).delAfterUse = true;
        buffer.at(minDegPos) = v;
        updateMin();
    }

    static bool sortByVertexId(const VertexInfo& u, const VertexInfo& v){
        return u.vertexId < v.vertexId;
    }

public:
    bool bufferInitilized = false;
    vid_t getAddedVerticesNum(){
        return addedVertices;
    }
    VertexBuffer(vid_t totalVertices){
        vertexMap.assign(totalVertices, false);
    }

    void needVertex(vid_t vertexId){
        vertexMap.at(vertexId) = true;
    }

    bool shouldCache(vid_t vertexId){
        return vertexMap.at(vertexId);
    }

    void setBufferSize(vid_t bufferSize){
        buffer.resize(bufferSize);
    }

    void streamV(VertexInfo& v){
        v.delAfterUse = false;
        buffer.at(addedVertices) = v;
        addedVertices ++;
//        if (buffer.size() < staticBufferSize){
//            v.delAfterUse = false;
//            buffer.push_back(v);
//            if (v.outDegree < minDegree){
//                minDegree = v.outDegree;
//                minDegPos = buffer.size() - 1;
//            }
//        }else{
//            if (v.outDegree > minDegree){
//                v.delAfterUse = false;
//                insert(v);
//            }
//        }
    }

    void sort(){
        std::sort(buffer.begin(), buffer.end(), sortByVertexId);
    }

    bool getVertex(VertexInfo &v){
        if (!vertexMap.at(v.vertexId) || !bufferInitilized){
            return false;
        }
        v.delAfterUse =false;
        vid_t vertexId = v.vertexId;
        size_t mid, left = 0;
        size_t right = buffer.size();
        vid_t midVertexId;
        while (left < right){
            mid = left +  (right - left) / 2;
            midVertexId = buffer.at(mid).vertexId;
            if (vertexId > midVertexId){
                left = mid + 1;
            }else if (vertexId < midVertexId){
                right = mid;
            }else{
                v.csr = buffer.at(mid).csr;
                v.outDegree = buffer.at(mid).outDegree;
                return true;
            }
        }
    }
};

#endif //IOE_SORW_VERTEXBUFFER_HPP
