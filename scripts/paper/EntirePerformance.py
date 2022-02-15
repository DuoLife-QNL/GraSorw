import os

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

blockSizes = [
    '20000',
    '524288',
    '1048576',
    '524288'
]

SOGW_ARNV = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/entire-performance/SOGW-ARNV'
SOGW_NVNV = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/entire-performance/SOGW-NVNV'
SOLID_ARNV = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/entire-performance/SOLID-ARNV'
SOLID_NVNV = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/entire-performance/SOLID-NVNV'

programs = [SOGW_ARNV, SOGW_NVNV, SOLID_ARNV, SOLID_NVNV]
programNames = ['SOGW_ARNV', 'SOGW_NVNV', 'SOLID_ARNV', 'SOLID_NVNV']

nmblocks = '2'
nThreads = '72'
walkLength = '80'
walksPerVertex = '10'

#node2vec
p= '1'
q= '1'

#autoregressive
alpha = '0.2'

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

logRootPath = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/paper/entire_performance'

os.chdir('/home/hongzheng/Codes/CLionProjects/IOE-SORW')
def makeBash(program, graphFile, blockSize, nmblocks, nThreads, walkLength, walksPerVertex, p, q, alpha, thsFilePath, logPath):
    bash = (
            program +
            " --file=" + graphFile +
            " --blocksize_kb=" + blockSize +
            " --nmblocks=" + nmblocks +
            " --nThreads=" + nThreads +
            " --walk-length=" + walkLength +
            " --walks-per-vertex=" + walksPerVertex +
            " --p=" + p +
            " --q=" + q +
            " --alpha=" + alpha +
            " --ths-file=" + thsFilePath +
            " > " + logPath
    )
    return bash

for fileName, graphPath, blockSize, thsFile in zip(fileNames, graphPaths, blockSizes, thsFiles):
    if fileName == fileNames[LJ]:
        continue
    for program, programName in zip(programs, programNames):
        if ((program == SOGW_ARNV or program == SOLID_ARNV) and fileName != fileNames[LJ]):
            continue
        if (program == SOGW_ARNV or program == SOGW_NVNV):
            ths = ths_all_0
            if graphPath == graphPaths[LJ]:
                walkLength = '80'
                walksPerVertex = '10'
            else:
                walkLength = '5'
                walksPerVertex = '1'
        else:
            ths = thsFile
            walkLength = '80'
            walksPerVertex = '10'

        for round in range(1, 4):
            logPath = logRootPath + '/' + fileName + '-' + programName + '-' + str(round) + '.log'
            bash = makeBash(program, graphPath, blockSize, nmblocks, nThreads, walkLength, walksPerVertex, p, q, alpha, ths, logPath)
            os.system(bash)