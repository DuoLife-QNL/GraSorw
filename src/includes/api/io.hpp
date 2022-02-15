
#ifndef DEF_IOUTIL_HPP
#define DEF_IOUTIL_HPP

#include <unistd.h>
#include <cassert>
#include <cstdlib>
#include <cerrno>
#include <zlib.h>
#include <fcntl.h>
#include "logger/logger.hpp"


template <typename T>
void preada(int f, T * tbuf, size_t nbytes, size_t off = 0) {
    size_t nread = 0;
    T * buf = (T*)tbuf;
    while(nread<nbytes) {
        ssize_t a = pread(f, buf, nbytes - nread, off + nread);
        if (a == (-1)) {
            logstream(LOG_INFO) << "Error, could not read: " << strerror(errno) << "; file-desc: " << f << std::endl;
            logstream(LOG_INFO) << "Pread arguments: " << f << " tbuf: " << tbuf << " nbytes: " << nbytes << " off: " << off << std::endl;
            assert(a != (-1));
        }
        assert(a>0);
        buf += a/sizeof(T);
        nread += a;
    }
    assert(nread <= nbytes);
}

template <typename T>
size_t readfull(int f, T ** buf) {
    off_t sz = lseek(f, 0, SEEK_END);
    lseek(f, 0, SEEK_SET);
    *buf = (T*)malloc(sz);
    preada(f, *buf, sz, 0);
    return sz;
}

template <typename T>
void pwritea(int f, T * tbuf, size_t nbytes, size_t off = 0) {
    size_t nwritten = 0;
    T * buf = (T*)tbuf;
    while(nwritten<nbytes) {
        size_t a = pwrite(f, buf, nbytes-nwritten, off+nwritten);
        assert(a>0);
        if (a == size_t(-1)) {
            logstream(LOG_ERROR) << "Could not write " << (nbytes-nwritten) << " bytes!" << " error:" <<  strerror(errno) << std::endl;
            assert(false);
        }
        buf += a;
        nwritten += a;
    }  
}

template <typename T>
void writefile(std::string fname, T * buf, T * &bufptr){
    int f = open(fname.c_str(), O_WRONLY | O_CREAT, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
    if (f < 0) {
        logstream(LOG_ERROR) << "Could not open " << fname << " error: " << strerror(errno) << std::endl;
    }
    assert(f >= 0);
    pwritea( f, buf, bufptr - buf );
    close(f);
    bufptr = buf;
}

template <typename T>
void appendfile(std::string fname, T * buf, T * &bufptr){
    int f = open(fname.c_str(), O_WRONLY | O_CREAT | O_APPEND, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
    if (f < 0) {
        logstream(LOG_ERROR) << "Could not open " << fname << " error: " << strerror(errno) << std::endl;
    }
    assert(f >= 0);
    pwritea( f, buf, bufptr - buf );
    close(f);
}

#endif