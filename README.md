# mantis
Mantis: A Fast, Small, and Exact Large-Scale Sequence-Search Index

Overview
--------

Mantis is a space-efficient data structure that can be used to index thousands of raw-
read experiments and facilitate large-scale sequence searches on those experiments. Mantis uses counting quotient
filters instead of Bloom filters, enabling rapid index builds and queries, small indexes, and exact results, i.e., no
false positives or negatives. Furthermore, Mantis is also a colored de Bruijn graph representation, so it supports fast
graph traversal and other topological analyses in addition to large-scale sequence-level searches.

Mantis was presented at RECOMB 2018, and a full journal paper is published in [Cell Systems](https://www.cell.com/cell-systems/abstract/S2405-4712(18)30239-4).  If you use Mantis, please cite the paper:

>Pandey, Prashant, Fatemeh Almodaresi, Michael A. Bender, Michael Ferdman, Rob Johnson, and Rob Patro. "Mantis: A Fast, Small, and Exact Large-Scale Sequence-Search Index." Cell Systems (2018).

A preprint of the paper is available [on bioRxiv](https://www.biorxiv.org/content/biorxiv/early/2017/11/10/217372.full.pdf).

API
--------
* `mantis build`: builds a mantis index from a collection of (squeakr) CQF files.
* `mantis query`: query k-mers in the mantis index.

Build
-------

Library dependencies (given version or higher):
 - [zlib](https://zlib.net/)
 - [sdsl-lite](https://github.com/simongog/sdsl-lite)
 
To build mantis, you will also need [CMake](https://cmake.org/) version 3.5 or higher.

The Counting Quotient Filter (CQF) code uses two new instructions to implement select on machine words
introduced in intel's Haswell line of CPUs. However, there is also an alternate
implementation of select on machine words to work on CPUs older than Haswell.
To build on an older hardware (older than Haswell) pass `-DNH=1` as a cmake argument.

```bash
 $ mkdir build
 $ cd build
 $ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../ ..
 $ make install
 $ cd ..
 $ ./bin/mantis build -s 20 -i raw/incqfs.lst -o raw/
```

If SDSL is not installed in a standard location, you can try and tell CMake where to look by adding
the following to the cmake command:

```
 -DSDSL_INSTALL_PATH=<path-to-sdsl-build-location>
```

The usage for this command are as follows:

Build
-------
`mantis build` creates a colored de Bruijn graph representation that can be used to query transcripts.

```bash
 $ ./bin/mantis build -s 28 -i txp_filtered_squeakrs.list -o txn_k23.idx/
```

```
SYNOPSIS
        mantis build [-e] -s <log-slots> -i <input_list> -o <build_output>

OPTIONS
        -e, --eqclass_dist
                    write the eqclass abundance distribution

        <log-slots> log of number of slots in the output CQF

        <input_list>
                    file containing list of input filters

        <build_output>
                    directory where results should be written
```

The following are the arguments to mantis build:
 - log-slots: The initial value for log of the number of slots in the CQF (i.e. the number of quotient bits).
 Mantis will automatically resize and go to the next value for log-slots very early on during the build process
 and it'll continue to resize until all the k-mers and their associated count/IDs fit.
 You would want to start from a reasonably small number so that you
 won't have a lot of resizing iterations and at the same time your final CQF
 if not larger than the size it could be. (for example starting from 30, while 
 even 28 was enough for log-slots)
 Suggested starting numbers based on the size of the squeakr you have is as following:
    - 28 for a small set of genomes like a bacterial genome
    - 30 for a large set of medium size read files.
    - 33 for a large set of big read files. 
 Notice that these are just suggestions. You can start with an arbitrarily small log-slot. 
 - input squeakr files: A list of input squeakr files (path to files) and cutoffs separated by tab. 
 A sample input squeakr file in provided in the raw dir. This is a list of squeakr output files that are generated after running squeakr on input experiments. We have provided two sample squeakr files in data dir.
 - build_output: The prefix filepath where all the output files will be written.

Build MST
-------
`mantis mst` encodes the Color information into a list of succinct 
int-vectors and bit-vectors. It creates a Color graph derived from de Bruijn graph of k-mers 
and encodes its Minimum Spanning Tree in a format to be able to retrieve the color classes.

```bash
 $ ./bin/mantis mst -p txp_k23.idx/ -t 8
```

The options and arguments are as follows:

```bash
SYNOPSIS
        mantis mst -p <index_prefix> [-t <num_threads>]

OPTIONS
        <index_prefix>
                    The directory where the index is stored.

        <num_threads>
                    number of threads
```
The command takes the following options :
 - `--input-prefix,-p`: the index directory.
 - `--threads,-t`: number of threads (default: `1`)

This step is will reduce the size of the color class 
significantly proportional to the original size of the color classes
(The larger the color class, the more the size reduction.)
So it is highly recommended that you run this step after `mantis build`
since this makes your query required memory much smaller and doesn't hurt
the query time. 

Query
-------

`mantis query` lets you query a mantis index with a set of sequences.

```bash
 $ ./bin/mantis query -p txp_k23.idx/ -o query.res input_txns.fa
```

The options and arguments are as follows:

```bash
SYNOPSIS
        mantis query [-b] [-1] [-j] [-k <kmer>] -p <query_prefix> [-o <output_file>] <query>

OPTIONS
        -b, --bulk  Process the whole input query file as a bulk.

        -1, --use-colorclasses
                    Use color classes as the color info representation instead of MST

        -j, --json  Write the output in JSON format
        <kmer>      size of k for kmer.

        <query_prefix>
                    Prefix of input files.

        <output_file>
                    Where to write query output.

        <query>     Prefix of input files.
```

 The only required option for the command is the following:
 - `--query-prefix,-p`: the directory where the output of coloreddbg command is present.
 
 additionally the command takes the following mandatory _positional_ argument :
 - query transcripts: input transcripts to be queried.

 There are also a couple of optional inputs:
 - `--use-colorclasses,-1`: This option runs a query over the list of color classes.
 - `--bulk,-b`: mantis supports two types of queries: bulk and serial.
 The default is `serial` unless this option is set.
 - `-k <kmer>`: mantis supports approximate queries for `k`
 larger than the `k` that the index and its de Bruijn graph was built with.
 `k` can only be larger than the `index k`. If not set, the default
 is providing exact query results for a `k` equal to the `index k`.
 
 **Note** that if you haven't run `mantis mst` and don't
 have the MST encoding of color information, the `--use-colorclasses,-1` option becomes
 mandatory, because the default behavior of query is to look for
 the MST encoding of the color information unless this option is set.
 
 Finally, rather than writing the results in the "simple" output format, they can be written in JSON if you
 provide the `--json,-j` flag to the `query` comamnd.
 
The output file in contains the list of experiments (i.e., hits) corresponding to each queried transcript.

Contributing
------------
Contributions via GitHub pull requests are welcome.

Authors
-------
- Prashant Pandey <ppandey@cs.stonybrook.edu>
- Fatemeh Almodaresi <falmodaresit@cs.stonybrook.edu>
- Michael Bender <bender@cs.stonybrook.edu>
- Mike Ferdman <mferdman@cs.stonybrook.edu>
- Rob Johnson <rob@cs.stonybrook.edu>
- Rob Patro <rob.patro@cs.stonybrook.edu>
