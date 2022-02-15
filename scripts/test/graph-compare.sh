#!/bin/bash

file[0]=/home/hongzheng/data/soc-LiveJournal1.txt.unDir
file[1]=/home/hongzheng/data/twitter_rv.net.unDir
file[2]=/home/hongzheng/data/uk_edge.txt.unDir
file[3]=/home/hongzheng/data/com-friendster.txt.unDir

graph_name[0]=LJ
graph_name[1]=twitter
graph_name[2]=uk
graph_name[3]=friendster

num_vertices[0]=4846609
num_vertices[1]=41652230
num_vertices[2]=105153952
num_vertices[3]=65608366

program[0]=/mnt/home/idmg/lhz/CLionProjects/IOE-SORW/build/plain/Plain_Node2Vec_TCI
program[1]=/mnt/home/idmg/lhz/CLionProjects/IOE-SORW/build/bucket_order/BucketOrder_Node2Vec_TCI
program[2]=/mnt/home/idmg/lhz/CLionProjects/IOE-SORW/build/dual_bucket/DualBucket_Node2Vec_TCI

program[3]=/mnt/home/idmg/lhz/CLionProjects/IOE-SORW/build/plain/Plain_Node2VecAR_TCI
program[4]=/mnt/home/idmg/lhz/CLionProjects/IOE-SORW/build/bucket_order/BucketOrder_Node2VecAR_TCI
program[5]=/mnt/home/idmg/lhz/CLionProjects/IOE-SORW/build/dual_bucket/DualBucket_Node2VecAR_TCI

program_name[0]=Plain_Node2Vec
program_name[1]=BucketOrder_Node2Vec
program_name[2]=DualBucket_Node2Vec
program_name[3]=Plain_Node2VecAR
program_name[4]=BucketOrder_Node2VecAR
program_name[5]=DualBucket_Node2VecAR

# Node2Vec
p[0]=1
p[1]=4
p[2]=0.25

q[0]=1
q[1]=0.25
q[2]=4

#for program_i in 1 2
#do
#  for fi in 1 2 3
#  do
#    for i in 0 1 2
#    do
#        ${program[program_i]} \
#        --file=${file[fi]} \
#        --num-vertices=${num_vertices[fi]} \
#        --nmblocks=2 \
#        --nThreads=16 \
#        --walk-length=80 \
#        --walks-per-vertex=10 \
#        --p=${p[i]} \
#        --q=${q[i]} \
#        > /mnt/home/idmg/lhz/CLionProjects/IOE-SORW/log/graph-compare/${program_name[program_i]}-${graph_name[fi]}-${p[i]}-${q[i]}.log
#    done
#  done
#done

#Autoregressive on Node2vec
alpha[0]=0.2
alpha[1]=0.8

#for fi in 0
#do
#  for program_i in 4
#  do
#    for i in 0
#    do
#        ${program[program_i]} \
#        --file=${file[fi]} \
#        --num-vertices=${num_vertices[fi]} \
#        --nmblocks=2 \
#        --nThreads=16 \
#        --walk-length=80 \
#        --walks-per-vertex=10 \
#        --alpha=${alpha[i]} \
#        --blocksize_kb=20000 \
#        --bucket-th=0 \
#        > /mnt/home/idmg/lhz/CLionProjects/IOE-SORW/log/graph-compare/${program_name[program_i]}-${graph_name[fi]}-${alpha[i]}-B=17.log
#    done
#  done
#done

blocksize=150000
while [ ${blocksize} -lt 1000000 ]
do
  for i in 0
  do
      /home/hongzheng/Codes/CLionProjects/IOE-SORW/cmake-build-idmg-monitor/IOE_SORW \
      --file=${file[1]} \
      --num-vertices=${num_vertices[1]} \
      --blocksize_kb=${blocksize} \
      --nmblocks=2 \
      --nThreads=64 \
      --walk-length=5 \
      --walks-per-vertex=2 \
      --p=${p[i]} \
      --q=${q[i]} \
      > /home/hongzheng/Codes/CLionProjects/IOE-SORW/log/load-tradeoff/NV-1-1-${graph_name[1]}.log
  done
  blocksize='expr ${blocksize} + 50000'
done
