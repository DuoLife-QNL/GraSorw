#!/bin/bash
cd ../../
bin=cmake-build-idmg/IOE_SORW
log_path=log/bucket-order/

th=0
while [ $(echo "${th} < 15"|bc) -eq 1 ]; do
    ${bin} --bucket-th=${th} > ${log_path}th=${th}.log
    echo "Execute bucket test th = ${th} done."
    th=$(echo "${th} + 0.5" | bc)
done