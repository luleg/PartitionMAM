# Partitioning Directed Networks Based on their Motif Adjacency Matrices

This is a C++ software for partitioning directed networks based on their Motif Adjacency Matrices (MAM). Main steps of the software are illustrated below.

<img align="center" src="https://github.com/luleg/PartitionMAM/blob/main/Images/AlgoSimple.png" width="100%">

Namely,
1. The MAM of the directed network is built using the [buildMAM software](https://github.com/luleg/MotifAdjacencyMatrix).  The doc *GraphletIdentifiers.pdf* lists all the motifs upon which a MAM can be built.
2. This MAM then is partitioned via the [Louvain algorithm](https://arxiv.org/pdf/0803.0476.pdf).
3. The disconnected nodes from the MAM are postprocessed by a homemade adaptation of Louvain.

## A Word about the Postprocessing Step

The postprocessing is done by applying the Louvain algorithm on the directed network, symmetrised by forgetting edge directions, and in which nodes that belong to a same cluster have been merged into a unique meta-node. An additional constraint forbids the algorithm to put in a same cluster of the final partition two meta-nodes containing more than *k* nodes from the initial network.

## Requirements

A bash shell and a `g++` compiler are enough to compile and use the software.

Tested in an `Ubuntu 18.04` environment emulated via a `Windows Subsystem for Linux 1`, with `gcc version 7.5.0` as a compiler.

*The present release contains the required files from third-party software, and can be used as a standalone.*

### Third-Party Software

* This implementation uses the [SNAP Software v6.0](https://snap.stanford.edu/snap/download.html), e.g. for graph structure and time management.
* The partitioning and postprocessing steps have been derived fron the [implementation of generic Louvain](https://sourceforge.net/projects/louvain/files/GenericLouvain/).
* Building the MAM of directed networks is done using the [buildMAM software](https://github.com/luleg/MotifAdjacencyMatrix).

## Installation

On a bash command, at the root of the folder:

```bash
cd src/UtilsSNAP
make
cd ../Pipeline
make
cd ../..
```

## Usage

To keep it short, the software can be used for six different tasks, and each of the six following commands in the root folder explains how to run one of these tasks.

```bash
./src/Pipeline/pipeline -ma -h     # Task: Build a MAM.
./src/Pipeline/pipeline -pa -h     # Task: Partition a network/MAM.
./src/Pipeline/pipeline -po -h     # Task: Postprocess disconnected nodes.
./src/Pipeline/pipeline -mapa -h   # Task: Build a MAM and partition it.
./src/Pipeline/pipeline -papo -h   # Task: Partition a MAM and postprocess disconnected nodes.
./src/Pipeline/pipeline -h         # Task: Build a MAM, partition it, postprocess disconnected nodes.
```

### Detailed Usage

A number of arguments must/can be used for each task. The table below provides a summary of these arguments, with a brief description.
<img align="center" src="https://github.com/luleg/PartitionMAM/blob/main/Images/Params.png" width="100%">

More precisely:

* The flags ```-ma```, ```-pa```, ```-po```, are used to indicate the kind of atomic tasks one wants to perform. They can be merged to perform non atomic tasks.

<img align="right" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/toyGraph.png" width="15%" height="15%">

* Input argument ```-igraph PathToDirectedGraph``` provides the path to the directed graph, that must be an edgelist with integer nodes, with only two columns, as shown on the right.
* Input argument ```-isym PathToSymmetrisedGraph``` can be used instead of ```-igraph```, e.g. when the symmetrised graph has been computed, or when working solely on a MAM.
* Input argument ```-imam PathToMAM``` provides the path to an already computed MAM.
* Input argument ```-ipart PathToPartition``` provides the path to an already computed partition.

* Output argument ```-omam PathToMAM``` is the path to the file in which the computed MAM will be stored.
* Output argument ```-opart PathToPartition``` is the path to the file in which the computed partition will be stored.
* Output argument ```-oppart PathToPartition``` is the path to the file in which the partial partition (i.e. without postprocessing of disconnected nodes) will be stored, when both partial and full partitions are computed.

:warning: For the output arguments, if the file already exists, its content is destroyed.

* Parameter argument ```-m MotifIdentifier``` indicates the motif to use to build the MAM. See the doc *GraphletIdentifiers.pdf* for the list of all admissible motifs, along with their identifier.
* Parameter argument ```-nth NumberOfThreads``` is the number of threads to use for building the MAM.

:bulb: To use when the network is large and the motif is a quadrangle. Number of threads should be between 4 and 8.

* Parameter argument ```-l levelInLouvain``` is the Louvain hierarchical level used as the partition.

:bulb: Use *levelInLouvain=1* for the finest hierarchical level, which has the larger number of communities, and *levelInLouvain=-2* for the coarsest hierarchical level (i.e.  with the smallest number of communities).

* Parameter argument ```-c ResParamLouvain``` is the modularity resolution parameter to be used in Louvain.

:bulb: It is expected that the highest the resolution parameter, the smallest the communities. Classic values lie between 1 and 2.

* Parameter argument ```-cc ResParamPostproc``` is the modularity resolution parameter to be used in the postprocessing step.

:bulb: The default value is low (1e-3) to avoid the creation of new clusters in the final partition.

* Parameter argument ```-k SizeMergeableMetaNode``` is the maximum number of nodes that a meta-node can contain to be mergeable to another meta-node.

### Comparing Partitioning

The Python program *Analysis/CompParts.py* enables visual comparison of two partitionings. Namely, it relabels partitions in both partitionings using the [Munkres algorithm](https://software.clapper.org/munkres/index.html), and plots visual comparisons of the partitionings:

* In a staircase fashion (length of a segment of y-coordinate *i* being the number of nodes belonging to part *i*).
* With the confusion matrix
