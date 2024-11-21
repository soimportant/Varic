Varic is a error correction tool for Third Generation Sequencing data (long read).
It utilized the variation graph for error pruning, and using De Bruijn graph for error correction.

## Requirements

It is recommanded to use virtual environment to install the dependencies.
The `requirements.txt` file contains all the dependencies needed to run the project.
  
```
conda create -n varic python=3.8
conda activate varic
pip install -r requirements.txt
```

- CMake >= 3.10
- compiler supporting C++20 (GCC >= 10.2, Clang >= 10.0)

## Build project

Varic is a C++ project, and it is built using CMake. To build the project, you need to create a build directory and run CMake to generate the build files.

```
mkdir build
cd build
cmake .. && make -j
```

The executable file `varic` will be generated in the `build` directory.

## Usage

```
> ./varic -h
  -h [ --help ]                      Show help message
  -r [ --raw_read ] arg              The path to TGS raw reads (.fasta or .fastq)
  -c [ --overlap_info ] arg          The path to the overlap information of raw
                                     reads (.paf)
  -o [ --output_path ] arg           The output path of corrected read
  -t [ --thread ] arg (=1)           The maximum number of threads
  --match arg (=5)                   The match score in alignment
  --mismatch arg (=-4)               The mismatch score in alignment
  --gap arg (=-8)                    The gap score in alignment
  --extend arg (=-6)                 The gap extension score in alignment
  --prune arg (=0.94999999999999996) The prune threshold for the variation
                                     graph, this value should set between 0 and
                                     1. The higher the value, the pruned graph
                                     will be more aggressive. If your data is
                                     noisy, which means has more haplotypes
                                     inside it, you should set this value
                                     higher
  -d [ --depth ] arg (=-1)           The maximum sequence to be used when
                                     building variation graph of a window, -1
                                     means take all sequences. This will be
                                     used to reduce the run-time of the
                                     program, but will effect the accuracy of
                                     the result.
  --seed arg                         The seed for random number generator
  -q [ --quiet ]                     disable all log
  --debug                            enable debug mode (print more verbose log to
                                     stderr)
```

### Example

For example, to correct the reads using the following command:

```
minimap2 -x ava-ont -t 8 reads.fasta reads.fasta > reads.paf
./varic -r reads.fasta -c reads.paf -o corrected.fasta -t 8
```