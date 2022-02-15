program[0]=/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/dual-bucket/DualBucket_Node2Vec_FL
program[1]=/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/dual-bucket/DualBucket_Node2Vec_NL_Test
program[2]=/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/dual-bucket/DualBucket_Node2Vec_Pairwise_OUTPUT_Twitter_18

program_name[0]=full-load
program_name[1]=no-load
program_name[2]=pairwise

graph[0]=/home/hongzheng/data/twitter_rv.net.unDir.dorder.part.18
graph[1]=/home/hongzheng/data/twitter_rv.net.unDir
graph[2]=/home/hongzheng/data/com-friendster.txt.unDir
graph[3]=/home/hongzheng/data/twitter_rv.net.unDir

graph_partition_file=/home/hongzheng/data/twitter_rv.net.unDir.toMetis.part.18

graph_name[0]=dorder
graph_name[1]=metis
graph_name[2]=friendster
graph_name[3]=twitter

for p in 2
do
  ${program[p]} \
  --file=${graph[3]} \
  --num-vertices=41652230 \
  --nmblocks=2 \
  --nThreads=36 \
  --walk-length=80 \
  --walks-per-vertex=10 \
  --load-test-output-file=/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/${program_name[p]}-${graph_name[3]}.csv
done

#for p in 0 1
#do
#  ${program[p]} \
#  --file=${graph[1]} \
#  --num-vertices=41652230 \
#  --nmblocks=2 \
#  --nThreads=36 \
#  --walk-length=80 \
#  --walks-per-vertex=10 \
#  --load-test-output-file=/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/${program_name[p]}-${graph_name[1]}.csv \
#  --graph-partition-method=metis \
#  --graph-partition-file=${graph_partition_file} \
#  --pre-partition-block-nums=18
#done
#
#for p in 0 1
#do
#  ${program[p]} \
#  --file=${graph[2]} \
#  --num-vertices=65608366 \
#  --nmblocks=2 \
#  --nThreads=36 \
#  --walk-length=80 \
#  --walks-per-vertex=10 \
#  --load-test-output-file=/home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/${program_name[p]}-${graph_name[2]}.csv
#done
