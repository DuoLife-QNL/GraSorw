//
// Created by lihz on 2020/10/28.
//

#ifndef IOE_SORW_SETTINGS_HPP
#define IOE_SORW_SETTINGS_HPP

#include <vector>
#include <cstdint>

#define DEBUG 0
#define INFO 1

/* The bi-block execution engine */
#define BI_BLOCK 0
/* The plain bucket execution engine */
#define PLAIN_BUCKET 0
/* The plain execution engine, including SOGW and SGSC */
#define PLAIN 0
#if PLAIN
/* SGSC if set to 1, SOGW if set to 0 */
#define STATICCACHE 0
#endif
/* The first-order random walk  engine */
#define FIRST_ORDER_ENGINE 1

#if FIRST_ORDER_ENGINE
/* Whether output the time cost metric, for example, the time of massive light vertex IO */
#define TIME_COST_INFO 1
/* The iteration-based method */
#define ITERATION_BASED 1
#if ITERATION_BASED
/* If 0, use the alphabet algorithm */
#define SKIP_EMPTY 0
#endif
/* State-aware block scheduling mechanism used in GraphWalker */
#define STATE_AWARE 0
#endif

#define RWNV 1
#define PRNV 0
#define DEEPWALK 0

/* Learning-based block loading settings. All set to 0 if you do not use learning-based block loading */
#define FULLY_LOAD 0
#if FULLY_LOAD
#define OUTPUT_FULLY_DATA 1
#endif
#define ONDEMAND_LOAD 0
#if ONDEMAND_LOAD
#define OUTPUT_ON_DEMAND_DATA 0
#endif



#define	RAND_MAX	2147483647
#define	FILE_SIZE	1024 // GB
#define	VERT_SIZE	64 * 1024 * 1024 // 64M vertices in beg_pos buffer in preprocess
#define	EDGE_SIZE	256 * 1024 * 1024 // 256M edges in csr buffer in preprocess
#define	WALK_BUFFER_SIZE	1024 * 1024 // most 1024 * 1024 walks in a in-memory walk buffer
#define	MEM_BUDGET	44 * 1024 * 1024 // for 64GB memory machine
// #define	MEM_BUDGET	4 * 1024 * 1024 // for 8GB memory machine

typedef unsigned __int128 uint128_t;

typedef uint32_t vid_t;
typedef uint64_t eid_t;
typedef uint64_t wid_t; //type of id of walks
typedef uint16_t bid_t; //type of id of blocks
typedef uint16_t hid_t; //type of id of hops
typedef uint8_t tid_t; //type of id of threads
typedef unsigned VertexDataType;
/* second order walk data type */
typedef uint128_t WalkDataType;

#define INVALID_BID UINT16_MAX
#define INVALID_VID UINT32_MAX

#define UNWEIGHTED_GRAPH 1

#define MULTI_THREAD 1

#define CACHE 0

#if BI_BLOCK
#define HOP_S 0
#define CUR_BLOCK_S 10
#define PRE_BLOCK_S 21
#define CUR_VERTEX_S 32
#define PRE_VERTEX_S 64
#define SRC_VERTEX_S 96
#define MAX_HOP 0x3ff
#define MAX_CUR_BLOCK 0x7ff
#define MAX_PRE_BLOCK 0x7ff
#define MAX_CUR_VERTEX 0xffffffff
#define MAX_PRE_VERTEX 0xffffffff
#define MAX_SRC_VERTEX 0xffffffff
#elif (PLAIN | PLAIN_BUCKET)
#define HOP_S 0
#define CUR_BLOCK_S 16
#define PRE_BLOCK_S 16
#define CUR_VERTEX_S 32
#define PRE_VERTEX_S 64
#define SRC_VERTEX_S 96
#define MAX_HOP 0xffff
#define MAX_CUR_BLOCK 0
#define MAX_PRE_BLOCK 0xffff
#define MAX_CUR_VERTEX 0xffffffff
#define MAX_PRE_VERTEX 0xffffffff
#define MAX_SRC_VERTEX 0xffffffff
#elif FIRST_ORDER_ENGINE
#define HOP_S 0
#define CUR_BLOCK_S 32
#define PRE_BLOCK_S 32
#define CUR_VERTEX_S 32
#define PRE_VERTEX_S 64
#define SRC_VERTEX_S 64
#define MAX_HOP 0xffffffff
#define MAX_CUR_BLOCK 0
#define MAX_PRE_BLOCK 0
#define MAX_CUR_VERTEX 0xffffffff
#define MAX_PRE_VERTEX 0
#define MAX_SRC_VERTEX 0xffffffff
#endif

#if BI_BLOCK | FIRST_ORDER_ENGINE
bid_t staticBlock_global;
bid_t dynamicBlock_global;
double getCSR_dynamicBlock_global_ms;

#define MULTI_BUFFER 0

#if FULLY_LOAD
#define OUTPUT_FULLY_DATA 1
#endif

#define NO_LOAD_TEST 0

#if ONDEMAND_LOAD
// 是否static_block也用on-demand-load, only for second-order random walk
#define ONDEMAND_LOAD_STATIC 0
#define SWAP_BLOCK_DEBUG 0
// COMPUTE_EDGE_CUT only, no walk processings
#define COMPUTE_EDGE_CUT 0
#endif

// 控制对于获取csr时间的统计（统计从dynamic block中获得csr的开销）
// 真正要统计的时候这里要区分线程再取平均，暂时不用这种方式
#define METRIC_GET_CSR 1

#if NO_LOAD_TEST | ONDEMAND_LOAD
#define USE_PRE_DEFINE_TH 0
#define USE_TRAINED_THS 1
#define OUTPUT_NO_LOAD_DATA 0

#if USE_PRE_DEFINE_TH
#define NO_LOAD_TH 1
#elif USE_TRAINED_THS
#define THS_FILE "/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-all-0.txt"
#endif
#endif
#endif

std::string thsFile_dynamic;
std::string thsFile_static;
std::string loadTestOutPutFileName_dynamicBlocks;
std::string loadTestOutPutFileName_staticBlocks;
#define LOAD_TRADEOFF_FILE " "
#if FULLY_LOAD
#define LOAD_TRADEOFF_FILE "/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/full-load.csv"
#elif NO_LOAD_TEST
#define LOAD_TRADEOFF_FILE "/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/no-load.csv"
#elif ONDEMAND_LOAD
#endif

long long blockSize_kb_global = 0;

long long edgeIn = 0;
long long edgeCut = 0;

// STATISTIC
#define BLOCK_UTI 0
#define BLOCK_UTI_FILE_NAME "/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/block-uti/default.csv"
#define ACT_VER 0
#define ACT_VER_FILENAME "/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/actived-vertex/default.csv"
std::string blockUtiFileName;
std::string actVerFileName = ACT_VER_FILENAME;

#endif //IOE_SORW_SETTINGS_HPP
