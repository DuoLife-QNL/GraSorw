import os
# graph settings
ths_all_0 = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-all-0.txt'

class graph(object):
    ths_d = ths_all_0
    ths_s = ths_all_0

    def __init__(self, name, graphPath, nVertices, blockSize, nBlocks):
        self.name = name
        self.graphPath = graphPath
        self.blockSize = blockSize
        self.nBlocks = nBlocks
        self.nVertices = nVertices

LJ = graph('LJ', '/home/hongzheng/data/soc-LiveJournal1.txt.unDir', 4846609, 20000, 17)
TW = graph('TW', '/home/hongzheng/data/twitter_rv.net.unDir', 41652230, 524288, 18)
FR = graph('FR', '/home/hongzheng/data/com-friendster.txt.unDir', 65608366, 524288, 27)
UK = graph('UK', '/home/hongzheng/data/uk_edge.txt.unDir', 105153952, 1048576, 25)
# TODO: Complete graph info here
Kron29 = graph('KR', '/mnt/data/hongzheng/data/Kron29.txt.unDir', 276924598, 0, 0)

TW.ths_d = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-on-demand-nopair.txt'
FR.ths_d = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-friendster.txt'
UK.ths_d = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-uk.txt'

# program settings
program = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/paper/Cache/Static-Cache'
nmblocks_default = '2'
nThreads_default = '72'

## node2vec
p_default= '1'
q_default= '1'
walkLength_default = '80'
walksPerVertex_default = '10'

## autoregressive
alpha_default = '0.2'

# running specific settins
logRootPath = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/paper/Static-Cache' + '/'
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
            " --num-vertices=" +args['nVertices'] +
            " > " + logPath
    )
    return bash

def setArg(
    graph,
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
    thsFileDynamic = None,
):
    arg = {
        'blockSize': str(graph.blockSize),
        'nVertices' : str(graph.nVertices),
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

def makeLogPath(program, graphFile, round, *, prefix = 'GraSorw'):
    return logRootPath + prefix + '-' + graphFile + '-' + program + '-' + str(round) + '.log'

def run(program, graphFile, args, logPath):
    os.system(makeBash(program, graphFile, args, logPath))

g = LJ
args = setArg(g)
logPath = makeLogPath('NV', g.name, 1)
run(program, g.graphPath, args, logPath)

for g in [TW, FR, UK]:
    args = setArg(g, walkLength=10, walksPerVertex=1)
    logPath = makeLogPath('NV', g.name, 1)
    run(program, g.graphPath, args, logPath)