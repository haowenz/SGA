# SGA
We provide a prototype implementation of two sequence to graph alignment algorithms [1] and [2]. *The implementation is straitforward without any code optimization.* A faster parallel version is still in progress.

## Installation

### Install dependencies
CMake, OpenMP, Protobuf 3.6.1 and GTest if you want to build the test.

### Download and install
First get the repo:
```
git clone git@github.com:haowenz/SGA.git
```
Then just run:
```
cd SGA && mkdir build && cd build
cmake -DCMAKE_SYSTEM_PREFIX_PATH="path_to_gtest;path_to_protobuf" -DSGAL_BUILD_TESTING=ON ..
make
```

## Usage
SGA is a header only C++ library. It supports loading graphs in the vg format and sequence files in fasta/q format. You can easily use its API to build your own applications. An example on how to use the library to align sequences to graphs is given. After the installation, you can simply test it by running:
```
./sga_example graph_file read_file
```

## References
[1] Jain, Chirag, Haowen Zhang, Yu Gao, and Srinivas Aluru. "On the complexity of sequence to graph alignment." In International Conference on Research in Computational Molecular Biology, pp. 85-100. Springer, Cham, 2019.

[2] Navarro, Gonzalo. "Improved approximate pattern matching on hypertext." Theoretical Computer Science 237, no. 1-2 (2000): 455-463.
