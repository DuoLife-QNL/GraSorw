program_Node2vec_soGraphWalker=/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/plain/release-build/Plain_Node2Vec_TCI
program_DeepWalk_GraphWalker=/home/hongzheng/Codes/CLionProjects/GraphWalker/build/release-build/deepwalk

graph_LJ=/home/hongzheng/data/soc-LiveJournal1.txt.unDir
graph_Twitter=/home/hongzheng/data/twitter_rv.net.unDir

# node2vec on lj
#for i in 1 2 3
#do
#  ${program_Node2vec_soGraphWalker} \
#  --file=${graph_LJ} \
#  --blocksize_kb=20000 \
#  --nmblocks=1 \
#  --nThreads=72 \
#  --walk-length=80 \
#  --walks-per-vertex=10 \
#  --p=1 \
#  --q=1 \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/paper/massive_io_percentage/LJ-node2vec-${i}.log
#done
#
# node2vec on twitter
for i in 1 2 3
do
  ${program_Node2vec_soGraphWalker} \
  --file=${graph_Twitter} \
  --blocksize_kb=524288 \
  --nmblocks=1 \
  --nThreads=72 \
  --walk-length=5 \
  --walks-per-vertex=1 \
  --p=1 \
  --q=1 \
  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/paper/massive_io_percentage/Twitter-node2vec-${i}.log
done

# deepwalk on lj
#for i in 1 2 3
#do
#  ${program_DeepWalk_GraphWalker} \
#  --file=${graph_LJ} \
#  --N=4846609 \
#  --R=10 \
#  --L=80 \
#  --blocksize_kb=20000 \
#  --nmblocks=1 \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/paper/massive_io_percentage/LJ-deepwalk-${i}.log
#done

# deepwalk on twitter
#for i in 1 2 3
#do
#  ${program_DeepWalk_GraphWalker} \
#  --file=${graph_Twitter} \
#  --N=41652230 \
#  --R=10 \
#  --L=80 \
#  --blocksize_kb=524288 \
#  --nmblocks=1 \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/paper/massive_io_percentage/Twitter-deepwalk-${i}.log
#done
