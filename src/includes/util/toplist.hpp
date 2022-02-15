
/**
 * @file
 * @author  Aapo Kyrola <akyrola@cs.cmu.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Aapo Kyrola, Guy Blelloch, Carlos Guestrin / Carnegie Mellon University]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 
 *
 * @section DESCRIPTION
 *
 * Tools for listing the TOP K values from a verte data file.
 */

#ifndef DEF_GRAPHCHI_TOPLIST
#define DEF_GRAPHCHI_TOPLIST

#include <vector>
#include <algorithm>
#include <errno.h>
#include <assert.h>

#include "util/merge.hpp"
#include "util/qsort.hpp"
#include "api/filename.hpp"
#include "util/VerCount.hpp"
  
    template <typename VertexDataType>
    struct vertex_value {
        vid_t vertex;
        VertexDataType value;
        vertex_value() {}
        vertex_value(vid_t v, VertexDataType x) : vertex(v), value(x) {}
    };

    static VARIABLE_IS_NOT_USED std::string filename_vertex_data(std::string basefilename) {
        std::stringstream ss;
        ss << basefilename;
        ss << "_GraphWalker/" << sizeof(VertexDataType) << "B.vvalue";
        return ss.str();
    }
    
    template <typename VertexDataType>
    bool vertex_value_greater(const vertex_value<VertexDataType> &a, const vertex_value<VertexDataType> &b) {
        return a.value > b.value;
    }

    static VARIABLE_IS_NOT_USED size_t get_num_vertices(std::string basefilename) {
        std::string numv_filename = nverticesname(basefilename);
        std::ifstream vfileF(numv_filename.c_str());
        size_t n;
        if (!vfileF.good()) {
            logstream(LOG_WARNING) << "Could not find file " << numv_filename << std::endl;
            logstream(LOG_INFO) << "Pre-counting num of vertices." << std::endl;
            n = getMaxVertexId(basefilename);
            return n;
        }
        vfileF >> n;
        vfileF.close();
        return n;
    }
     
    /**
      * Reads the vertex data file and returns top N values.
      * Vertex value type must be given as a template parameter.
      * This method has been implemented in a manner to consume very little
      * memory, i.e the whole file is not loaded into memory (unless ntop = nvertices).
      * @param basefilename name of the graph
      * @param ntop number of top values to return (if ntop is smaller than the total number of vertices, returns all in sorted order)
      * @param from first vertex to include (default, 0)
      * @param to last vertex to include (default, all)
      * @return a vector of top ntop values  
     */
    template <typename VertexDataType>
    std::vector<vertex_value<VertexDataType> > get_top_vertices(std::string basefilename, int ntop, vid_t from=0, vid_t to=0) {
        typedef vertex_value<VertexDataType> vv_t;
        
        /* Initialize striped IO manager */
        metrics m("toplist");
        // stripedio * iomgr = new stripedio(m);
        
        /* Initialize the vertex-data reader */
        vid_t readwindow = 1024 * 1024;
        size_t numvertices = get_num_vertices(basefilename);
        VertexDataType * vertexdata = (VertexDataType*)calloc(readwindow, sizeof(VertexDataType));
        VertexDataType * vertexdataptr = vertexdata;
        // vertex_data_store<VertexDataType> * vertexdata =
        // new vertex_data_store<VertexDataType>(basefilename, numvertices, iomgr);
           
        if ((size_t)ntop > numvertices) {
            ntop = (int)numvertices;
        }
                
        /* Initialize buffer */
        vv_t * buffer_idxs = (vv_t*) calloc(readwindow, sizeof(vv_t));
        vv_t * topbuf = (vv_t*) calloc(ntop, sizeof(vv_t));
        vv_t * mergearr = (vv_t*) calloc(ntop * 2, sizeof(vv_t));
        
        /* Iterate the vertex values and maintain the top-list */
        size_t idx = 0;
        vid_t st = 0;
        vid_t en = numvertices - 1;

        int count = 0;
        while(st <= numvertices - 1) {
            en = st + readwindow - 1;
            if (en >= numvertices - 1) en = numvertices - 1;
            int nt = en - st + 1;
            
            /* Load the vertex values */
            std::string vertex_value_file = filename_vertex_data(basefilename);
           int f = open(vertex_value_file.c_str(), O_RDONLY | O_CREAT, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
            if (f < 0) {
                logstream(LOG_ERROR) << "Could not open " << vertex_value_file << " error: " << strerror(errno) << std::endl;
            }
            assert(f >= 0);
            preada(f, vertexdata, nt*sizeof(VertexDataType), st*sizeof(VertexDataType));
            close(f);
            // vertexdata->load(st, en);
            
            int k = 0;
            VertexDataType minima = VertexDataType();
            if (count > 0) {
                minima = topbuf[ntop - 1].value; // Minimum value that should be even considered
            }
            vertexdataptr = vertexdata;
            for(int j=0; j < nt; j++) {
                VertexDataType val = *vertexdataptr;
                vertexdataptr++;

                if (count == 0 || (val > minima)) {
                    buffer_idxs[k] = vv_t((vid_t)idx + from, val);
                    k++;
                }
                idx++;
            }
            nt = k; /* How many were actually included */
            
            /* Sort buffer-idxs */
            quickSort(buffer_idxs, nt, vertex_value_greater<VertexDataType>);

            /* Merge the top with the current top */
            if (count == 0) {
                /* Nothing to merge, just copy */
                memcpy(topbuf, buffer_idxs, ntop * sizeof(vv_t));
            } else {
                // void merge(ET* S1, int l1, ET* S2, int l2, ET* R, F f) {
                merge<vv_t>(topbuf, ntop, buffer_idxs, std::min(ntop, nt), mergearr,  vertex_value_greater<VertexDataType>);
                memcpy(topbuf, mergearr, ntop * sizeof(vv_t));
            }
            
            count++;
            st += readwindow;
        }
                   
        /* Return */
        std::vector< vv_t > ret;
        for(int i=0; i < ntop; i++) {
            ret.push_back(topbuf[i]);
        }
        free(buffer_idxs);
        free(mergearr);
        free(topbuf);
        free(vertexdata);
        // delete iomgr;

        return ret;
    }

#endif

