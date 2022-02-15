import os
os.chdir('/home/hongzheng/Codes/CLionProjects/IOE-SORW')
LJ = 0
TW = 1
UK = 2
FR = 3
fileNames = ["LJ", "TW", "UK", "FR"]
graphPaths = [
    '/home/hongzheng/data/soc-LiveJournal1.txt.unDir',
    '/home/hongzheng/data/twitter_rv.net.unDir',
    '/home/hongzheng/data/uk_edge.txt.unDir',
    '/home/hongzheng/data/com-friendster.txt.unDir'
]

metisPartFiles = [
    None,
    '/home/hongzheng/data/twitter_rv.net.unDir.toMetis.part.18',
    '/home/hongzheng/data/uk_edge.txt.unDir.toMetis.part.25',
    None
]

thsFiles_metis = [
    None,
    '/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/twitter-metis-OL-nopair-ths.txt',
    '/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/uk-metis-OL-nopair-ths.txt',
    None
]

metisBlockNums = [None, 18, 25, None]

blockSizes = [
    '20000',
    '524288',
    '1048576',
    '524288'
]

blockNums = ['17', '18', '25', '27']
SOGW_ARNV = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/entire-performance/SOGW-ARNV'
SOGW_NVNV = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/entire-performance/SOGW-NVNV'
SOLID_ARNV = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/entire-performance/SOLID-ARNV'
SOLID_NVNV = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/entire-performance/SOLID-NVNV'
programs = [SOGW_ARNV, SOGW_NVNV, SOLID_ARNV, SOLID_NVNV]
programNames = ['SOGW_ARNV', 'SOGW_NVNV', 'SOLID_ARNV', 'SOLID_NVNV']
program1 = {
    'AR': '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/io_reducing/1_SOGW/SOGW-ARNV',
    'NV': '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/io_reducing/1_SOGW/SOGW-NVNV'
}

program2 = {
    'AR': '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/io_reducing/2_simple-bucket-based/BB-ARNV',
    'NV': '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/io_reducing/2_simple-bucket-based/BB-NVNV'
}

program3 = {
    'AR': '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/io_reducing/3_2-and-DP-MPC/BB-DP-MPC-ARNV',
    'NV': '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/io_reducing/3_2-and-DP-MPC/BB-DP-MPC-NVNV'
}

program4 = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/io_reducing/4_3-and-metis/BB-DP-MPC-NVNV'
program5 = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/io_reducing/5_3-and-OLSB/BB-DP-MPC-OLSB-NVNV'
program6 = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/io_reducing/6_3-and-metis-OLSB/BB-DP-MPC-OLSB-NVNV'
program7 = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/io_reducing/7_3-and-MBE/BB-DP-MPC-MBE-NVNV'
program8 = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/io_reducing/8_entire-SOLID/SOLID-NVNV'
exprPrograms = [program1, program2, program3, program4, program5, program6, program7, program8]
nmblocks_default = '2'
nThreads_default = '72'
walkLength_default = '80'
walksPerVertex_default = '10'
#node2vec

p_default= '1'
q_default= '1'
#autoregressive

alpha_default = '0.2'
ths_all_0 = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-all-0.txt'
ths_twitter = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-on-demand-nopair.txt'
ths_uk = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-uk.txt'
ths_friendster= '/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-friendster.txt'
thsFiles = [
    ths_all_0,
    ths_twitter,
    ths_uk,
    ths_friendster
]

logRootPath = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/paper/IO_reducing'
def makeBash(program, graphFile, args, logPath):
    bash = (
            program +
            " --file=" + graphFile +
            " --blocksize_kb=" + args['blockSize'] +
            " --nmblocks=" + args['nmblocks'] +
            " --nThreads=" + args['nThreads'] +
            " --walk-length=" + args['walkLength'] +
            " --walks-per-vertex=" + args['walksPerVertex'] +
            " --p=" + args['p'] +
            " --q=" + args['q'] +
            " --alpha=" + args['alpha'] +
            " --graph-partition-method=" + args['partitionMethod'] +
            " --graph-partition-file=" + args['partitionFile'] +
            " --pre-partition-block-nums=" + args['prePartitionBlockNums'] +
            " --load-test-output-file-static=" + args['outPutFileStatic'] + 
            " --load-test-output-file-dynamic=" + args['outPutFileDynamic'] +
            " --ths-file-static=" + args['thsFileStatic'] +
            " --ths-file-dynamic=" + args['thsFileDynamic'] +
            " > " + logPath
    )
    return bash

def setArg(
    blockSize,
    *, 
    nmblocks = nmblocks_default, 
    nThreads = nThreads_default, 
    walkLength = walkLength_default,
    walksPerVertex = walksPerVertex_default,
    p = p_default,
    q = q_default,
    alpha = alpha_default,
    partitionMethod = 'default',
    partitionFile = 'ERROR',
    prePartitionBlockNums = '0',
    outPutFileStatic = None,
    outPutFileDynamic = None,
    thsFileStatic = None,
    thsFileDynamic = None
):
    arg = {
        'blockSize': str(blockSize),
        'nmblocks': str(nmblocks),
        'nThreads': str(nThreads),
        'walkLength': str(walkLength),
        'walksPerVertex': str(walksPerVertex),
        'p': str(p),
        'q': str(q),
        'alpha': str(alpha),
        'partitionMethod' : partitionMethod,
        'partitionFile' : partitionFile,
        'prePartitionBlockNums' : str(prePartitionBlockNums),
        'outPutFileStatic' : str(outPutFileStatic),
        'outPutFileDynamic' : str(outPutFileDynamic),
        'thsFileStatic' : str(thsFileStatic),
        'thsFileDynamic' : str(thsFileDynamic)
    }
    return arg

def makeLogPath(exprId, program, graphFile, round):
    return logRootPath + '/Expr' + str(exprId) + '-' + graphFile + '-' + program + '-' + str(round) + '.log'

for exprId, program in zip(range(1, 9), exprPrograms):
    if exprId == 1:
        # prog = program['AR']
        # graphFile = graphPaths[LJ]
        # args = setArg(blockSizes[LJ])
        # for round in range(1, 4):
        #     logPath = makeLogPath(exprId, 'AR', 'LJ', round)
        #     os.system(makeBash(prog, graphFile, args, logPath))
        # prog = program['NV']
        # for blockSize, graphFile, graphName in zip(blockSizes, graphPaths, fileNames):
        #     if graphFile != graphPaths[LJ]:
        #         args = setArg(blockSize, walkLength=5, walksPerVertex=1)
        #     else:
        #         args = setArg(blockSizes[LJ])
        #     for round in range(1, 4):
        #         logPath = makeLogPath(exprId, 'NV', graphName, round)
        #         os.system(makeBash(prog, graphFile, args, logPath))
        continue
    elif exprId == 2 or exprId == 3:
        # prog = program['AR']
        # graphFile = graphPaths[LJ]
        # args = setArg(blockSizes[LJ])
        # for round in range(1, 4):
        #     logPath = makeLogPath(exprId, 'AR', 'LJ', round)
        #     os.system(makeBash(prog, graphFile, args, logPath))
        # prog = program['NV']
        # for blockSize, graphFile, graphName in zip(blockSizes, graphPaths, fileNames):
        #     args = setArg(blockSize)
        #     for round in range(1, 4):
        #         logPath = makeLogPath(exprId, 'NV', graphName, round)
        #         os.system(makeBash(prog, graphFile, args, logPath))
        continue
    elif exprId == 4:
        # startIndex = TW
        # endIndex = UK + 1
        # for blockSize, graphFile, graphName, metisPartFile, blockNum in zip(blockSizes[startIndex: endIndex], graphPaths[startIndex: endIndex], fileNames[startIndex: endIndex], metisPartFiles[startIndex: endIndex], metisBlockNums[startIndex: endIndex]):
        #     args = setArg(blockSize, partitionMethod='metis', partitionFile=metisPartFile, prePartitionBlockNums=blockNum)
        #     for round in range(1, 4):
        #         logPath = makeLogPath(exprId, 'NV', graphName, round)
        #         os.system(makeBash(program, graphFile, args, logPath))
        continue
    elif exprId == 5:
        # startIndex = TW
        # endIndex = FR + 1
        # for blockSize, graphFile, graphName, ths in zip(blockSizes[startIndex: endIndex], graphPaths[startIndex: endIndex], fileNames[startIndex: endIndex], thsFiles[startIndex: endIndex]):
        #     args = setArg(blockSize, thsFile=ths)
        #     for round in range(1, 4):
        #         logPath = makeLogPath(exprId, 'NV', graphName, round)
        #         os.system(makeBash(program, graphFile, args, logPath))
        continue
    elif exprId == 6:
        # startIndex = TW
        # endIndex = UK + 1
        # for blockSize, graphFile, graphName, metisPartFile, blockNum, ths in zip(blockSizes[startIndex: endIndex], graphPaths[startIndex: endIndex], fileNames[startIndex: endIndex], metisPartFiles[startIndex: endIndex], metisBlockNums[startIndex: endIndex], thsFiles_metis[startIndex: endIndex]):
        #     args = setArg(blockSize, partitionMethod='metis', partitionFile=metisPartFile, prePartitionBlockNums=blockNum, thsFile=ths)
        #     for round in range(1, 4):
        #         logPath = makeLogPath(exprId, 'NV', graphName, round)
        #         os.system(makeBash(program, graphFile, args, logPath))
        continue
    elif exprId == 7:
        for blockSize, graphFile, graphName in zip(blockSizes, graphPaths, fileNames):
            if (graphName == fileNames[LJ] or graphName == fileNames[FR]):
                args = setArg(blockSize)
                for round in range(1, 4):
                    logPath = makeLogPath(exprId, 'NV', graphName, round)
                    os.system(makeBash(program, graphFile, args, logPath))
    elif exprId == 8:
        for blockSize, graphFile, graphName, metisPartFile, blockNum, thsMetis, ths in zip(blockSizes, graphPaths, fileNames, metisPartFiles, metisBlockNums, thsFiles_metis, thsFiles):
            # if (graphName == fileNames[LJ]):
            #     args = setArg(blockSize)
            #     for round in range(1, 4):
            #         logPath = makeLogPath(exprId, 'NV', graphName, round)
            #         os.system(makeBash(program, graphFile, args, logPath))
            # elif graphName == fileNames[FR]:
            #     args = setArg(blockSize, thsFile=ths)
            #     for round in range(1, 4):
            #         logPath = makeLogPath(exprId, 'NV', graphName, round)
            #         os.system(makeBash(program, graphFile, args, logPath))
            # else:
            #     args = setArg(blockSize, partitionMethod='metis', partitionFile=metisPartFile, prePartitionBlockNums=blockNum, thsFile=thsMetis)
                # for round in range(1, 4):
            #         logPath = makeLogPath(exprId, 'NV', graphName, round)
            #         os.system(makeBash(program, graphFile, args, logPath))
            # if (graphName == fileNames[TW] or graphName == fileNames[UK]):
            #     args = setArg(blockSize, thsFile=ths)
            #     for round in range(1, 4):
            #         logPath = makeLogPath(exprId, 'NV', graphName, round)
            #         os.system(makeBash(program, graphFile, args, logPath))
            continue