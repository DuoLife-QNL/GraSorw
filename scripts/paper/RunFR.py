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

logRootPath = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/friendster/2021-9-26/'
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

def makeLogPath(program, graphFile, round):
    return logRootPath + graphFile + '-' + program + '-' + str(round) + '.log'

program_FL_OUTPUT = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/dual-bucket/trade-off/map-get/Individually-trained/FL_OUTPUT'
program_OL_OUTPUT = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/dual-bucket/trade-off/map-get/Individually-trained/OL_OUTPUT'

ths_all_1 = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-all-1.txt'

outputLog_FL_Dynamic = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/friendster/2021-9-23/block-processing-logs/FL-Dynamic.csv'
outputLog_FL_Static = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/friendster/2021-9-23/block-processing-logs/FL-Static.csv'
outputLog_OL_Dynamic = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/friendster/2021-9-23/block-processing-logs/OL-TH-1-Dynamic.csv'
outputLog_OL_Static = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/friendster/2021-9-23/block-processing-logs/OL-TH-1-Static.csv'

# get log of full load
# args = setArg(blockSizes[FR], outPutFileStatic=outputLog_FL_Static, outPutFileDynamic=outputLog_FL_Dynamic)
# logPath = makeLogPath('FL-GetLog', 'FR')
# os.system(makeBash(program_FL_OUTPUT, graphPaths[FR], args, logPath))

# get log of on-demand load with th set to 1
# args = setArg(blockSizes[FR], thsFileDynamic=ths_all_1, thsFileStatic=ths_all_1, outPutFileStatic=outputLog_OL_Static, outPutFileDynamic=outputLog_OL_Dynamic)
# logPath = makeLogPath('OL-GetLog', 'FR')
# os.system(makeBash(program_OL_OUTPUT, graphPaths[FR], args, logPath))

program_new_static_th = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/dual-bucket/trade-off/map-get/static-log2/OL'
# Learning-based block loading
for round in range(1, 4):
    args = setArg(
        blockSizes[FR], 
        thsFileDynamic='/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/FR/ths-friendster-old-d.txt',
        thsFileStatic='/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/FR/ths-friendster-old-s.txt'
        )
    logPath = makeLogPath('LBBL', 'FR', round)
    os.system(makeBash(program_OL_OUTPUT, graphPaths[FR], args, logPath))
for round in range(1, 4):
    args = setArg(
        blockSizes[FR], 
        thsFileDynamic=ths_all_0,
        thsFileStatic=ths_all_0
        )
    logPath = makeLogPath('FL', 'FR', round)
    os.system(makeBash(program_OL_OUTPUT, graphPaths[FR], args, logPath))