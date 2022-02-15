#!/bin/bash
cd ../../
log_path=log/model-compare/
dual_bucket_path=build/dual_bucket/
bucket_order_path=build/bucket_order/
plain_path=build/plain/

dual_bucket_bins[0]=DualBucket_Autoregressive_TCI
dual_bucket_bins[1]=DualBucket_Node2VecAR_TCI

bucket_order_bins[0]=BucketOrder_Autoregressive_TCI
bucket_order_bins[1]=BucketOrder_Node2VecAR_TCI

plain_bins[0]=Plain_Autoregressive_TCI
plain_bins[1]=Plain_Node2VecAR_TCI

for bin in ${dual_bucket_bins[*]}
do
  ${dual_bucket_path}${bin} > ${log_path}${bin}.log
done

for bin in ${bucket_order_bins[*]}
do
  ${bucket_order_path}${bin} > ${log_path}${bin}.log
done

for bin in ${plain_bins[*]}
do
  ${plain_path}${bin} > ${log_path}${bin}.log
done