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
 - zlib1g-dev 1:1.2.8.dfsg-2ubuntu4
 - sdsl

The CQF code uses two new instructions to implement select on machine words
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

`mantis build` creates a colored de Bruijn graph representation that can be used to query transcripts.

The following are the arguments to mantis build:
 - log-slots: The log of the number of slots in the CQF (i.e. the number of quotient bits).
 - input squeakr files: A list of input squeakr files (path to files) and cutoffs separated by tab. A sample input squeakr file in provided in the raw dir. This is a list of squeakr output files that are generated after running squeakr on input experiments. We have provided two sample squeakr files in data dir.
 - build_output: The prefix filepath where all the output files will be written.

Query
-------

`mantis query` lets you query a mantis index with a set of sequences.

```bash
 $ ./bin/mantis query -p raw/ -o query.res raw/input_txns.fa
```

The options and arguments are as follows:

```bash
SYNOPSIS
        mantis query [-j] -p <query_prefix> [-o <output_file>] <query>

OPTIONS
        -j, --json  Write the output in JSON format

        <query_prefix>
                    Prefix of input files.

        <output_file>
                    Where to write query output.

        <query>     Prefix of input files.
```

 The command takes the following options :
 - `--query-prefix,-p`: the directory where the output of coloreddbg command is present.
 - `--output,-o`: the file where the query results should be written (default : `samples.output`).
 
 additionally the command takes the following mandatory _positional_ argument :
 - query transcripts: input transcripts to be queried.

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
