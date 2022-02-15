# GraSorw

GraSorw is a disk-based graph processing system designed for scalable second-order random walk tasks.

## Quick Start

### Build the targets

CMake (version 3.5 or higher) is required to build the targets. The OpenMP should be installed to run in parallel.

It is recommended to build the targets in a "build" folder. In project root folder, do:
```
mkdir build
cd build
cmake ..
make
cd ..
```
Now the build targets are generated in folder `build/`, and you should now in the project root folder.

### Dataset preparation

To download the test dataset, in the project root folder, do:

```
mkdir datasets
cd datasets
wget https://snap.stanford.edu/data/soc-LiveJournal1.txt.gz
gzip -d soc-LiveJournal1.txt.gz
cd ..
```

Use the following command to generate an undirected version of this graph:

```
build/Convert2Undir datasets/soc-LiveJournal1.txt
```
A file called "soc-LiveJournal1.txt.unDir" should be generated in `datasets/`

### Configurations

To change the running configurations, modify file ```conf/GraSorw.cnf```

Meaning of arguments in `conf/GraSorw.cnf`:

* Graph information
  * file: the input graph file
  * num-vertices: the total number of vertices in the input graph. This value can be retrieved from the output log when converting the original graph to undirected.
  * blocksize_kb: the size of the current block and ancillary block in the engine in KiB. This is set for default sequential partitioning.
* Run-time configuration
  * nThreads: number of threads during execution. Default value is the maximum thread number of your machine.
* General random walk task settings
  * walk-length: the length of each walk. In PageRank tasks, this represents the maximum length of each walk.
  * walks-per-vertex: the number of walks starting from each activated vertex.
* Node2vec tasks
  * p: the hyper-parameter $p$ in Node2vec model
  * q: the hyper-parameter $q$ in Node2vec model
* PageRank tasks
  * start-vertices-file: the file path to indicate the query nodes in PageRank task, in which each line represents the ID of one query node.
  * num-start-vertices: the number of start vertices in the above start vertices file.
  * decay-factor: the decay factor in PageRank model.

Task switching:

Open the ` src/includes/engine/Settings.hpp` file,  and change the value of the micro definition of symbol `RWNV`and `PRNV`,  which represent the Random Walk on Nove2vec and the PageRank Query on Node2vec, respectively. The task whose corresponding symbol is set to one will be compiled and executed.

### Run

To process the task, in the project root folder, run:

```
build/GraSorw
```

## Configuration of learning-based block loading model

### Get the running logs

First build the executable under pure full load mode. Open the ` src/includes/engine/Settings.hpp` file and set `FULLY_LOAD` and `OUTPUT_FULLY_DATA` to 1, and set `ONDEMAND_LOAD` to 0. Then build the target. When running, add the following line in the configuration file `conf/GraSorw.cnf`:

```
load-test-output-file-dynamic=conf/on-demand-th-0.csv
```

This indicates the output path of the running logs under full load mode.

Then build the executable who accepts a given block loading method switching threshold. To do this, ensure that the `FULLY_LOAD` is set to 0 and `OUTPUT_ON_DEMAND_DATA` and `ONDEMAND_LOAD` is set to 1. Then build the target. Similarly, add the following line in the configuration file before running:

```
ths-file-dynamic=conf/ths/ths-all-1.txt
load-test-output-file-dynamic=conf/on-demand-th-1.csv
```

The first line indicates that the system would switch to on-demand loading when the ratio of $\mathbb{W} / N_v$ is less than 1, and the second line indicates where should the running logs be put.

After running the above two executables, the running logs should be output to folder `/conf`.

### Training

A python script to train the thresholds for each block is added in `scripts/Train.py` . You should change the `nBlocks` the script to the number of partitioned blocks for the graph you're processing. For the default block size set in our given configuration file for LiveJournal, the number of partioned blocks is 17. To generate the thresholds , run in project root path:

```  text
python scripts/Train.py
```

The thresholds should be generated as `conf/ths.txt`

### Run with the trained thresholds

Set the `ths-file-dynamic` in the configuration file to the generated threshold file. Following the above configuration, the file should be `conf/ths.txt`, so add the following line in `/conf`:

```
ths-file-dynamic=conf/ths.txt
```

To perform the learning-based loading model, run the target generated under`OUTPUT_ON_DEMAND_DATA`, `FULLY_LOAD` set to 0, and `ONDEMAND_LOAD` set to 1.

Note that the trained thresholds should be able to use in different tasks under this partition for the current graph. 

