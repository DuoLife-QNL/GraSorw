program=/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/dual-bucket/trade-off/map-get/O3/release/DualBucket_Node2Vec_OL_Pairwise
program_OLSB=/home/hongzheng/Codes/CLionProjects/IOE-SORW/build/dual-bucket/trade-off/map-get/O3/release/DualBucket_Node2Vec_OLSB_Pairwise

ths_all_0=/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-all-0.txt
ths_no_pair=/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-on-demand-nopair.txt
ths_pairwise=/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-on-demand-pairwise.txt
ths_force=/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-all-force-0.01.txt
ths_light_task_20_5=/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-twitter-light-task-20-5.txt
ths_uk=/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-uk.txt
ths_friendster=/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/ths-friendster.txt
ths_tw_mts_nopair=/home/hongzheng/Codes/CLionProjects/IOE-SORW/conf/twitter-metis-OL-nopair-ths.txt


# light twitter
for i in 1 2 3
do
  ${program_OLSB} \
  --file=/home/hongzheng/data/twitter_rv.net.unDir \
  --blocksize_kb=524288 \
  --nmblocks=2 \
  --nThreads=72 \
  --walk-length=20 \
  --walks-per-vertex=5 \
  --p=1 \
  --q=1 \
  --ths-file=${ths_no_pair} \
  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/twitter/trade-off/on-demand-load/2021-8-11/light-stdth-OLSB-${i}.log
done

for i in 1 2 3
do
  ${program} \
  --file=/home/hongzheng/data/twitter_rv.net.unDir \
  --blocksize_kb=524288 \
  --nmblocks=2 \
  --nThreads=72 \
  --walk-length=20 \
  --walks-per-vertex=5 \
  --p=1 \
  --q=1 \
  --ths-file=${ths_no_pair} \
  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/twitter/trade-off/on-demand-load/2021-8-11/light-stdth-OL-${i}.log
done

for i in 1 2 3
do
  ${program} \
  --file=/home/hongzheng/data/twitter_rv.net.unDir \
  --blocksize_kb=524288 \
  --nmblocks=2 \
  --nThreads=72 \
  --walk-length=20 \
  --walks-per-vertex=5 \
  --p=1 \
  --q=1 \
  --ths-file=${ths_all_0} \
  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/twitter/trade-off/on-demand-load/2021-8-11/light-OL-th-0-${i}.log
done
#
## twitter
#for i in 1 2 3
#do
#  ${program_OLSB} \
#  --file=/home/hongzheng/data/twitter_rv.net.unDir \
#  --nmblocks=2 \
#  --nThreads=72 \
#  --walk-length=80 \
#  --walks-per-vertex=10 \
#  --p=1 \
#  --q=1 \
#  --ths-file=${ths_no_pair} \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/twitter/trade-off/on-demand-load/2021-8-10/OLSB-${i}.log
#done
#
#for i in 1 2 3
#do
#  ${program} \
#  --file=/home/hongzheng/data/twitter_rv.net.unDir \
#  --nmblocks=2 \
#  --nThreads=72 \
#  --walk-length=80 \
#  --walks-per-vertex=10 \
#  --p=1 \
#  --q=1 \
#  --ths-file=${ths_no_pair} \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/twitter/trade-off/on-demand-load/2021-8-10/OL-${i}.log
#done
#
#for i in 1 2 3
#do
#  ${program} \
#  --file=/home/hongzheng/data/twitter_rv.net.unDir \
#  --nmblocks=2 \
#  --nThreads=72 \
#  --walk-length=80 \
#  --walks-per-vertex=10 \
#  --p=1 \
#  --q=1 \
#  --ths-file=${ths_all_0} \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/twitter/trade-off/on-demand-load/2021-8-10/OL-th-0-${i}.log
#done

#uk
#for i in 1 2 3
#do
#  ${program_OLSB} \
#  --file=/home/hongzheng/data/uk_edge.txt.unDir \
#  --blocksize_kb=1048576 \
#  --nmblocks=2 \
#  --nThreads=72 \
#  --walk-length=80 \
#  --walks-per-vertex=10 \
#  --p=1 \
#  --q=1 \
#  --ths-file=${ths_uk} \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/uk/2021-8-10/OLSB-${i}.log
#done
#
#for i in 1 2 3
#do
#  ${program} \
#  --file=/home/hongzheng/data/uk_edge.txt.unDir \
#  --blocksize_kb=1048576 \
#  --nmblocks=2 \
#  --nThreads=72 \
#  --walk-length=80 \
#  --walks-per-vertex=10 \
#  --p=1 \
#  --q=1 \
#  --ths-file=${ths_uk} \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/uk/2021-8-10/OL-${i}.log
#done
#
#for i in 1 2 3
#do
#  ${program} \
#  --file=/home/hongzheng/data/uk_edge.txt.unDir \
#  --blocksize_kb=1048576 \
#  --nmblocks=2 \
#  --nThreads=72 \
#  --walk-length=80 \
#  --walks-per-vertex=10 \
#  --p=1 \
#  --q=1 \
#  --ths-file=${ths_all_0} \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/uk/2021-8-10/OL-th-0-${i}.log
#done


#for i in 1 2 3
#do
#  ${program} \
#  --file=/home/hongzheng/data/twitter_rv.net.unDir \
#  --blocksize_kb=524288 \
#  --nmblocks=2 \
#  --nThreads=72 \
#  --walk-length=80 \
#  --walks-per-vertex=10 \
#  --p=1 \
#  --q=1 \
#  --ths-file=${ths_pairwise} \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/twitter/trade-off/on-demand-load/2021-7-19/ths-pairwise-${i}.log
#done
#
#for i in 1 2 3
#do
#  ${program} \
#  --file=/home/hongzheng/data/twitter_rv.net.unDir \
#  --blocksize_kb=524288 \
#  --nmblocks=2 \
#  --nThreads=72 \
#  --walk-length=80 \
#  --walks-per-vertex=10 \
#  --p=1 \
#  --q=1 \
#  --ths-file=${ths_force} \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/twitter/trade-off/on-demand-load/2021-7-19/ths-force-0.01-${i}.log
#done

#for i in 1 2 3
#do
#  ${program} \
#  --file=/home/hongzheng/data/com-friendster.txt.unDir \
#  --nmblocks=2 \
#  --nThreads=72 \
#  --blocksize_kb=524288 \
#  --walk-length=80 \
#  --walks-per-vertex=10 \
#  --p=1 \
#  --q=1 \
#  --ths-file=${ths_friendster} \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/friendster/2021-7-22/on-demand-nopair-${i}.log
#done
#
#
#for i in 1 2 3
#do
#  ${program} \
#  --file=/home/hongzheng/data/com-friendster.txt.unDir \
#  --nmblocks=2 \
#  --nThreads=72 \
#  --blocksize_kb=524288 \
#  --walk-length=80 \
#  --walks-per-vertex=10 \
#  --p=1 \
#  --q=1 \
#  --ths-file=${ths_all_0} \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/friendster/2021-7-22/on-demand-th-0-${i}.log
#done
#
#for i in 1 2 3
#do
#  ${program} \
#  --file=/home/hongzheng/data/com-friendster.txt.unDir \
#  --nmblocks=2 \
#  --nThreads=72 \
#  --blocksize_kb=524288 \
#  --walk-length=20 \
#  --walks-per-vertex=5 \
#  --p=1 \
#  --q=1 \
#  --ths-file=${ths_friendster} \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/friendster/2021-7-22/light-task-20-5-on-demand-nopair-${i}.log
#done
#
#
#for i in 1 2 3
#do
#  ${program} \
#  --file=/home/hongzheng/data/com-friendster.txt.unDir \
#  --nmblocks=2 \
#  --nThreads=72 \
#  --blocksize_kb=524288 \
#  --walk-length=20 \
#  --walks-per-vertex=5 \
#  --p=1 \
#  --q=1 \
#  --ths-file=${ths_all_0} \
#  > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/friendster/2021-7-22/light-task-20-5-on-demand-th-0-${i}.log
#done