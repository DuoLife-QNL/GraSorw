#ifndef DEF_WALK_BUFFER
#define DEF_WALK_BUFFER


#include <cstring>
#include "Settings.hpp"
#include <cstdlib>
#include <cassert>

/**
 * second order walk buffer
 */
class WalkBuffer{
public:
    bool malloced;
    wid_t size_w;
    WalkDataType *walks;

    WalkBuffer(){
        size_w = 0;
        malloced = false;
    }

    ~WalkBuffer(){
        if(size_w > 0 && walks != NULL){
            free(walks);
        }
    }

    WalkDataType & operator[] (int i){
        return walks[i];
    }

    void push_back(WalkDataType w){
        if(!malloced){
            assert(size_w == 0);
            walks = (WalkDataType*)malloc(WALK_BUFFER_SIZE * sizeof(WalkDataType));
            malloced = true;
        }
        walks[size_w++] = w;
    }
};

#endif