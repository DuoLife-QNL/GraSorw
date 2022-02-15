
#ifndef GRAPHWALKER_CONVERSIONS_DEF
#define GRAPHWALKER_CONVERSIONS_DEF

#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

#include <fstream>
#include <iostream>

#include "engine/Settings.hpp"
#include "logger/logger.hpp"
#include "api/filename.hpp"
#include "api/io.hpp"
#include "ReOrder.hpp"

typedef enum PAR_METHOD{
    DEFAULT,
    METIS
} PAR_METHOD;

    long long max_value(long long a, long long b){
        return (a > b ? a : b);
    }
    
    long long min_value(long long a, long long b){
        return (a < b ? a : b);
    }

    //for convert_to_csr
    eid_t max_nedges; //max num of edges in an csr file
    bid_t fid, fnum;
    std::vector<vid_t> files;
    vid_t fstv; //start vertex of current file
    eid_t fstp; //start position in csr of current file

    vid_t bstv; //start vertex of current buffer
    eid_t bstp; //start position in csr of current buffer

    vid_t curvertex;
    eid_t cpos; //position in the total csr

    //for compute_block

    int rm_dir(std::string dir_full_path){    
        DIR* dirp = opendir(dir_full_path.c_str());    
        if(!dirp){
            return -1;
        }
        struct dirent *dir;
        struct stat st;
        while((dir = readdir(dirp)) != NULL){
            if(strcmp(dir->d_name,".") == 0 || strcmp(dir->d_name,"..") == 0){
                continue;
            }    
            std::string sub_path = dir_full_path + '/' + dir->d_name;
            if(lstat(sub_path.c_str(),&st) == -1){
                // Log("rm_dir:lstat ",sub_path," error");
                logstream(LOG_WARNING) << "rm_dir:lstat " << sub_path << " error" << std::endl;
                continue;
            }    
            if(S_ISDIR(st.st_mode)){
                if(rm_dir(sub_path) == -1){ // 如果是目录文件，递归删除
                    closedir(dirp);
                    return -1;
                }
                rmdir(sub_path.c_str());
            }
            else if(S_ISREG(st.st_mode)) {
                unlink(sub_path.c_str());     // 如果是普通文件，则unlink
            }
            else{
                // Log("rm_dir:st_mode ",sub_path," error");
                logstream(LOG_WARNING) << "rm_dir:st_mode " << sub_path << " error" << std::endl;
                continue;
            }
        }
        if(rmdir(dir_full_path.c_str()) == -1){//delete dir itself.
            closedir(dirp);
            return -1;
        }
        closedir(dirp);
        return 0;
    }

    /**
     * Returns the number of shards if a file has been already
     * sharded or 0 if not found.
     */
    static bid_t find_filerange(std::string base_filename, uint16_t filesize_GB) {
        std::string filerangefile = filerangename(base_filename, filesize_GB);
        FILE *tryf = fopen(filerangefile.c_str(), "r");
        if (tryf != NULL) { // Found!
            bid_t nshards = 0;
            while(!feof(tryf)){
                char flag = fgetc(tryf);
                if(flag == '\n')
                nshards++;
            }
            fclose(tryf);
            return nshards-1;
        }
        // Not found!
        logstream(LOG_WARNING) << "Could not find shards with filesize = " << filesize_GB << "GB." << std::endl;
        return 0;
    }

    void readMaxOutDegree(std::string base_file, eid_t &maxOutDegree){
        std::string maxOutDegreeFile = maxOutDegreeName(base_file);
        std::ifstream file(maxOutDegreeFile.c_str());
        file >> maxOutDegree;
    }

    /**
     * Returns the number of blocks if a file has been already
     * 0 if not found.
     */
    static bid_t find_blockrange(std::string base_filename, unsigned long long blocksize_kb) {
        std::string blockrangefile = blockrangename(base_filename, blocksize_kb);
        FILE *tryf = fopen(blockrangefile.c_str(), "r");
        if (tryf != NULL) { // Found!
            bid_t nblocks = 0;
            while(!feof(tryf)){
                char flag = fgetc(tryf);
                if(flag == '\n')
                nblocks++;
            }
            fclose(tryf);
            return nblocks-1;
        }
        // Not found!
        logstream(LOG_WARNING) << "Could not find blocks with blocksize = " << blocksize_kb << "KB." << std::endl;
        return 0;
    }

    void writeFileRange(std::string filename, uint16_t filesize_GB){
        /*write csr file range*/
        std::string filerangefile = filerangename(filename, filesize_GB);
        std::ofstream frf(filerangefile.c_str());      
        for( bid_t p = 0; p < files.size(); p++ ){
            frf << files[p] << std::endl;
        }
        frf.close();

        /*write nvertices*/
        std::string nverticesfile = nverticesname(filename);
        std::ofstream nvf(nverticesfile.c_str());
        nvf << files[fnum] << std::endl;
        nvf.close();
    }

    void writeMaxOutDegree(std::string filename, vid_t maxOutDegree){
        std::string maxOutDegreeFile = maxOutDegreeName(filename);
        std::ofstream frf(maxOutDegreeFile.c_str());
        frf << maxOutDegree << std::endl;
        frf.close();
    }

    void bwritezero( char * beg_pos, char * &beg_posptr, eid_t count ){
        while( count-- ){
            *((eid_t*)beg_posptr) = cpos;
            beg_posptr += sizeof(eid_t);
        }
    }

    void bwrite(char * beg_pos, char * &beg_posptr, char * csr, char * &csrptr, eid_t outd, std::vector<vid_t> outv, std::string filename ){
        cpos += outd;
        *((eid_t*)beg_posptr) = cpos;
        beg_posptr += sizeof(eid_t);
        for( eid_t i = 0; i < outd; i++ ){
            *((vid_t*)csrptr) = outv[i];
            csrptr += sizeof(vid_t);
        }
    }

    void flushInvl(std::string filename, char * csr, char * &csrptr, char * beg_pos, char * &beg_posptr){
        if( cpos - fstp >= max_nedges){ //new file needed
            logstream(LOG_INFO) << "FILE_" << fid << " : [ " << fstv << " , " << curvertex-1 << " ]" << std::endl;
            fstv = curvertex;
            files.push_back(fstv);
            fstp = cpos;
            fid++;
        }
        std::string fidfile = fidname(filename,fid);
        std::string csrname = fidfile + ".csr";
        std::string beg_posname = fidfile + ".beg_pos";
        appendfile(csrname, csr, csrptr);
        appendfile(beg_posname, beg_pos, beg_posptr);
        logstream(LOG_INFO) << "Buffer_" << " : [ " << bstv << " , " << curvertex-1 << " ] finished, with csr position : [" << bstp << ", " << cpos << ")" << std::endl;

        csrptr = csr;
        beg_posptr = beg_pos;
        bstv = curvertex;
        bstp = cpos;
    }

    bid_t convert_to_csr(std::string filename, uint16_t filesize_GB, eid_t &maxOutDegree){

        FILE * inf = fopen(filename.c_str(), "r");
        if (inf == NULL) {
            logstream(LOG_FATAL) << "Could not load :" << filename << " error: " << strerror(errno) << std::endl;
        }
        assert(inf != NULL);

        rm_dir((filename+"_GraSorw/").c_str());
        mkdir((filename+"_GraSorw/").c_str(), 0777);
        mkdir((filename+"_GraSorw/graphinfo/").c_str(), 0777);

        max_nedges = (eid_t)filesize_GB * 1024 * 1024 * 1024 / sizeof(vid_t); //max number of (vertices+edges) of a shard
        logstream(LOG_INFO) << "Begin convert_to_csr, max_nedges in an csr file = " << max_nedges << std::endl;
        logstream(LOG_INFO) << "VERT_SIZE in a beg_pos buffer = " << VERT_SIZE << "EDGE_SIZE in an csr buffer = " << EDGE_SIZE << std::endl;
        
        char * csr = (char*) malloc(EDGE_SIZE*sizeof(vid_t));
        char * csrptr = csr;
        char * beg_pos = (char*) malloc(VERT_SIZE*sizeof(eid_t));
        char * beg_posptr = beg_pos;
               
        logstream(LOG_INFO) << "Reading in edge list format!" << std::endl;

        fid = 0;
        fstv = 0;
        fstp = 0;
        files.push_back(fstv);

        bstv = 0;
        bstp = 0;
        *((eid_t*)beg_posptr) = 0;
        beg_posptr += sizeof(eid_t);

        curvertex = 0;
        cpos = 0;
        vid_t max_vert = 0;
        eid_t outd = 0;
        eid_t maxoutd = 0;
        std::vector<vid_t> outv;

        char s[1024];
        while(fgets(s, 1024, inf) != NULL) {
            if (s[0] == '#') continue; // Comment
            if (s[0] == '%') continue; // Comment
            
            char *t1, *t2;
            t1 = strtok(s, "\t, ");
            t2 = strtok(NULL, "\t, ");
            if (t1 == NULL || t2 == NULL ) {
                logstream(LOG_ERROR) << "Input file is not in right format. "
                << "Expecting <from> <to>. "
                << "Current line: " << s << "\n";
                assert(false);
            }
            vid_t from = atoi(t1);
            vid_t to = atoi(t2);
            if( from == to ) continue;
            max_vert = max_value(max_vert, from);
            max_vert = max_value(max_vert, to);
            if( from == curvertex ){
                outv.push_back(to);
                outd++;
                if (outd > maxoutd){
                    maxoutd = outd;
                }
            }else{  //a new vertex
                if( cpos - bstp + outd >= EDGE_SIZE || curvertex - bstv + 1 >= VERT_SIZE ){
                    flushInvl(filename, csr, csrptr, beg_pos, beg_posptr);
                    if( outd > EDGE_SIZE){
                        logstream(LOG_ERROR) << "Too small memory capacity with EDGE_SIZE = " << EDGE_SIZE << " to support larger ourdegree of vert " << curvertex << ", with outdegree = " << outd << std::endl;
                        assert(false);
                    }
                }
                bwrite(beg_pos, beg_posptr, csr, csrptr, outd, outv, filename); //write a vertex to buffer
                if( from - curvertex > 1 ){ //there are verts with zero out-links
                    vid_t remianzero = from-curvertex-1;
                    vid_t remainsize = VERT_SIZE - (curvertex - bstv);
                    while(remianzero > remainsize){
                        bwritezero( beg_pos, beg_posptr, remainsize ); 
                        flushInvl(filename, csr, csrptr, beg_pos, beg_posptr);
                        logstream(LOG_DEBUG) << remianzero << " , remainsize =  " << remainsize << std::endl;
                        remianzero -= remainsize;
                        remainsize = VERT_SIZE;
                    }
                    bwritezero( beg_pos, beg_posptr, remianzero ); 
                }
                curvertex = from;
                outd = 1;
                outv.clear();
                outv.push_back(to);
            }
        }
        fclose(inf);
        bwrite(beg_pos, beg_posptr, csr, csrptr, outd, outv, filename);//write the last vertex to buffer

        if(max_vert > curvertex){
            logstream(LOG_INFO) << "need bwritezero, as max_vert = " << max_vert << ", curvertex = " << curvertex << std::endl;
            bwritezero( beg_pos, beg_posptr, max_vert - curvertex );
        }       
        
        flushInvl(filename, csr, csrptr, beg_pos, beg_posptr);

        files.push_back(max_vert+1);
        fnum = fid+1;
        logstream(LOG_INFO) << "Partitioned csr file number : " << fnum << std::endl;
        writeFileRange(filename,filesize_GB);
        maxOutDegree = maxoutd;
        writeMaxOutDegree(filename, maxoutd);


        //output beg_pos information
        logstream(LOG_INFO) << "nverts = " << max_vert+1 << ", " << "nedges(fnum=1)) = " << cpos << std::endl;

        if(csr!=NULL) free(csr);
        if(beg_pos!=NULL) free(beg_pos);


        return fnum;
    }

    bid_t compute_block(std::string filename, unsigned long long blocksize_kb){
        eid_t mneb = (eid_t)blocksize_kb * 1024 / sizeof(vid_t); // max number of edges in a block
        logstream(LOG_INFO) << "Begin compute_block with blocksize = " << blocksize_kb << "KB, max number of edges in a block = " << mneb << std::endl;
        
        bid_t blockid = 0;
        std::vector<vid_t> blocks;
        vid_t stvb= 0; //start vertex of current block
        eid_t bgstvb= 0; //beg_pos of the start vertex of current block
        blocks.push_back(stvb);

        eid_t * beg_pos = (eid_t*) malloc(VERT_SIZE*sizeof(eid_t));
        std::string beg_posname = fidname(filename,fid) + ".beg_pos";
        int beg_posf = open(beg_posname.c_str(),O_RDONLY | O_CREAT, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
        if (beg_posf < 0) {
            logstream(LOG_FATAL) << "Could not load :" << beg_posname << ", error: " << strerror(errno) << std::endl;
        }
        assert(beg_posf > 0);
        vid_t ttv = lseek(beg_posf, 0, SEEK_END) / sizeof(eid_t); //total vertices number in begpos file
        vid_t rv = VERT_SIZE; //read vertices number in begpos file
        vid_t nread = 0;
        while( nread < ttv ){
            if(ttv - nread < VERT_SIZE) rv = ttv - nread;
            preada( beg_posf, beg_pos, (size_t)rv*sizeof(eid_t), (size_t)nread*sizeof(eid_t) );
            logstream(LOG_DEBUG) << "nread = " << nread << ", beg_pos[0] = " << beg_pos[0] << std::endl;
            for(vid_t v = 0; v < rv; v++){
                if(beg_pos[v]-bgstvb > mneb){
                    logstream(LOG_INFO) << "Block_" << blockid << " : [" << stvb << ", " << nread+v-1 << ")" << std::endl;
                    if(beg_pos[v]-beg_pos[v-1] > mneb){
                        logstream(LOG_WARNING) << "Too small blocksize with max num of edges of a block = " << mneb << " to support larger ourdegree of vert " << nread+v << ", with outdegree = " << beg_pos[v]-beg_pos[v-1] << std::endl;
                        logstream(LOG_WARNING) << "v = " << v << ", beg_pos[v-1] = " << beg_pos[v-1] << ", beg_pos[v] = " << beg_pos[v] << std::endl;
                        blockid+=2;
                        stvb = nread+v-1;
                        blocks.push_back(stvb);
                        stvb = nread+v;
                        blocks.push_back(stvb);
                        bgstvb = beg_pos[v];
                    }else{
                        blockid++;
                        stvb = nread+v-1;
                        blocks.push_back(stvb);
                        bgstvb = beg_pos[v-1];
                    }
                }
            }
            nread += rv;
        }
        close(beg_posf);
        blocks.push_back(ttv-1);
        blockid++;

        /*write block range*/
        std::string blockrangefile = blockrangename(filename, blocksize_kb);
        std::ofstream brf(blockrangefile.c_str());      
        for( bid_t p = 0; p < blocks.size(); p++ ){
            brf << blocks[p] << std::endl;
        }
        brf.close();

        return blockid;
    }

    /**
     * Converts graph from an edge list format. Input may contain
     * value for the edges. Self-edges are ignored.
     */

    /* parameters which have a default value should be given if not using default partition method here */
    bid_t convert_if_notexists(std::string basefilename, unsigned long long blocksize_kb, eid_t &maxOutDegree,
                               PAR_METHOD parMethod, std::string graParFileName = "", bid_t nBlocks = 0) {
        if (parMethod == METIS){
            std::string reOrderFileName = basefilename + ".MetisReOrder";
            bid_t nshards = find_filerange(reOrderFileName, FILE_SIZE);
            if(nshards > 0) {
                logstream(LOG_INFO) << "Found preprocessed files for " << basefilename << ", shardsize = " << FILE_SIZE << "GB, num file=" << nshards << std::endl;
                readMaxOutDegree(basefilename, maxOutDegree);
                logstream(LOG_INFO) << "Found max out degree: " << maxOutDegree << std::endl;
            }else{
                ReOrder reOrder(basefilename, graParFileName, nBlocks);
                reOrder.genReOrderFile(reOrderFileName);
                nshards = convert_to_csr(reOrderFileName, FILE_SIZE, maxOutDegree);
                logstream(LOG_INFO) << "Successfully finished sharding for " << reOrderFileName << std::endl;
                logstream(LOG_INFO) << "Created " << nshards << " shards." << std::endl;
                /* write block range file */
                std::ofstream brf(blockrangename(reOrderFileName, blocksize_kb));
                for (bid_t b = 0; b < nBlocks; b++){
                    brf << reOrder.getStartVertex(b) << std::endl;
                }
                brf << reOrder.getNVertex();
                brf.close();
            }
            return nBlocks;
        }else{
            bid_t nshards = find_filerange(basefilename, FILE_SIZE);
            /* Check if input file is already sharded */
            if(nshards > 0) {
                logstream(LOG_INFO) << "Found preprocessed files for " << basefilename << ", shardsize = " << FILE_SIZE << "GB, num file=" << nshards << std::endl;
                readMaxOutDegree(basefilename, maxOutDegree);
                logstream(LOG_INFO) << "Found max out degree: " << maxOutDegree << std::endl;
                //return nshards;
            }else{
                logstream(LOG_INFO) << "Did not find preprocessed shards for " << basefilename  << std::endl;
                logstream(LOG_INFO) << "Will try create them now..." << std::endl;

                nshards = convert_to_csr(basefilename, FILE_SIZE, maxOutDegree);

                logstream(LOG_INFO) << "Successfully finished sharding for " << basefilename << std::endl;
                logstream(LOG_INFO) << "Created " << nshards << " shards." << std::endl;
            }

            assert(blocksize_kb > 0);
            bid_t nblocks = find_blockrange(basefilename, blocksize_kb);
            if(nblocks > 0) {
                logstream(LOG_INFO) << "Found computed blocks for " << basefilename << ", blocksize = " << blocksize_kb << "KB, num blocks=" << nblocks << std::endl;
                //return nshards;
            }else{
                logstream(LOG_INFO) << "Did not find computed blocks for " << basefilename  << std::endl;
                logstream(LOG_INFO) << "Will try compute the blcok range now..." << std::endl;

                nblocks = compute_block(basefilename, blocksize_kb);

                logstream(LOG_INFO) << "Successfully finished compute_block for " << basefilename << std::endl;
                logstream(LOG_INFO) << "computed " << nblocks << " blocks." << std::endl;
            }
            return nblocks;
        }
    }


#endif

