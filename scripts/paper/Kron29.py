from BasicSettings import *

logRootPath = '/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/paper/entire_performance' + '/'
g = Kron29

# grasorw without bucket-extending or OL
# program = program_BIB
# args = setArg(g, walkLength=30, walksPerVertex=3)
# logPath = makeLogPath(logRootPath, program.name, g.name, 1)
# run(program, g.graphPath, args, logPath)

# sogw and static-cache
for program in [program_SOGW_NVNV]:
    args = setArg(g, walkLength=5, walksPerVertex=1)
    logPath = makeLogPath(logRootPath, program.name, g.name, 1)
    run(program, g.graphPath, args, logPath)
