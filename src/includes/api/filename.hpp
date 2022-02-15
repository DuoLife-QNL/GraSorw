#ifndef GraSorw_FILENAMES_DEF
#define GraSorw_FILENAMES_DEF

#include <fstream>
#include <fcntl.h>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>

#include "engine/Settings.hpp"
#include "logger/logger.hpp"

#define PRE 01
#define CUR 02

static std::string fidname( std::string basefilename, bid_t p ){
    std::stringstream ss;
    ss << basefilename;
    ss << "_GraSorw/graphinfo/file";
    ss << "_" << p;
    return ss.str();
}

static std::string walksname( std::string basefilename, bid_t p ){
    std::stringstream ss;
    ss << basefilename;
    ss << "_GraSorw/walks/pool";
    ss << "_" << p << ".walks";
    return ss.str();
}

static std::string bucketName(std::string baseFileName, bid_t dynamicBlock, int preOrCurFlag){
    std::stringstream ss;
    ss << baseFileName;
    if (preOrCurFlag == PRE){
        ss << "_GraSorw/walks/preBucket";
    }else if (preOrCurFlag == CUR){
        ss << "_GraSorw/walks/curBucket";
    }else{
        assert(false);
    }
    ss << "_" << dynamicBlock << ".walks";
    return ss.str();
}

static std::string filerangename(std::string basefilename, uint16_t filesize_GB){
    std::stringstream ss;
    ss << basefilename;
    ss << "_GraSorw/filesize_" << filesize_GB << "GB.filerange";
    return ss.str();
}

static std::string maxOutDegreeName(std::string basefilename){
    std::stringstream ss;
    ss << basefilename;
    ss <<"_GraSorw/maxOutDegree.dat";
    return ss.str();
}

static std::string blockrangename(std::string basefilename, unsigned long long blocksize_KB){
    std::stringstream ss;
    ss << basefilename;
    ss << "_GraSorw/blocksize_" << blocksize_KB << "KB.blockrange";
    return ss.str();
}

static std::string nverticesname(std::string basefilename) {
    std::stringstream ss;
    ss << basefilename;
    ss << "_GraSorw/N.nvertices"; 
    return ss.str();
}

/**
 * Configuration file name
 */
static std::string configname() {
    char * chi_root = getenv("GRAPHCHI_ROOT");
    if (chi_root != NULL) {
        return std::string(chi_root) + "/conf/graphchi.cnf";
    } else {
        return "conf/GraSorw.cnf";
    }
}

/**
 * Configuration file name - local version which can
 * override the version in the version control.
 */
static std::string configlocalname() {
    char * chi_root = getenv("GRAPHCHI_ROOT");
    if (chi_root != NULL) {
        return std::string(chi_root) + "/conf/GraSorw.local.cnf";
    } else {
        return "conf/GraSorw.local.cnf";
    }
}

#endif