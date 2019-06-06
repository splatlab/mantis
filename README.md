# mantis
Mantis: A Fast, Small, and Exact Large-Scale Sequence-Search Index

Overview
--------

Mantis is a space-efficient data structure that can be used to index thousands of raw-
read experiments and facilitate large-scale sequence searches on those experiments. Mantis uses counting quotient
filters instead of Bloom filters, enabling rapid index builds and queries, small indexes, and exact results, i.e., no
false positives or negatives. Furthermore, Mantis is also a colored de Bruijn graph representation, so it supports fast
graph traversal and other topological analyses in addition to large-scale sequence-level searches.

Mantis was presented at RECOMB 2018, and a full journal paper is published in [Cell Systems](https://www.cell.com/cell-systems/abstract/S2405-4712(18)30239-4).
New version of Mantis that uses MST-based compression for equivalence class bit
vectors was presented at RECOMB 2019.
If you use Mantis, please cite these papers:
>Prashant Pandey, Fatemeh Almodaresi, Michael A. Bender, Michael Ferdman, Rob Johnson, and Rob Patro. "Mantis: A Fast, Small, and Exact Large-Scale Sequence-Search Index." Cell Systems (2018).
>Fatemeh Almodaresi, Prashant Pandey, Michael Ferdman, Rob Johnson, and Rob Patro. "An Efficient, Scalable and Exact Representation of High-Dimensional Color Information Enabled via de Bruijn Graph Search." RECOMB (2019).


A preprint of the paper is available [on bioRxiv](https://www.biorxiv.org/content/biorxiv/early/2017/11/10/217372.full.pdf).

Notes
--------

Mantis uses `mmap` to read input Squeakr files and write to output counting
quotient filter (CQF) file. During construction, input Squeakr files and output
CQF file are accessed sequentially. Each page is accessed only once and can be
removed from memory once it's used.  However, unless there is memory pressure
used pages are not cleared from the process memory which causes `/usr/bin/time`
tool to report very high max resident set size (RSS).

In the current release, we use `madvise` system call to explicitly let the
kernel know to remove used pages from process memory. This has negligible cost
to the overall runtime of the construction process and keeps max RSS in check.

Mantis should be used with the latest version of
[squeakr](https://github.com/splatlab/squeakr/tree/master), and we _highly
recommend_ running squeakr with the desired k-mer count threshold and the
`--no-counts` argument.  Early versions of mantis used _unfiltered_ squeakr
output to build the mantis data structure, which required considerable
intermediate disk-space, as those files represented the original k-mers and
their counts in each sample exactly.  When run with the `--no-counts` argument,
each squeakr file encodes the threshold only once in its metadata, and includes
only the k-mers that passed the abundance threshold; this can reduce the
intermediate storage requirements by over an order of magnitude.


API
--------
* `mantis build`: builds a mantis index from a collection of (squeakr) CQF files.
* `mantis mst`: builds a new encoding based on Minimum Spanning Trees for the color information.
* `mantis query`: query k-mers in the mantis index.

Build
-------

Library dependencies (given version or higher):
 - [zlib](https://zlib.net/)
 - [sdsl-lite](https://github.com/simongog/sdsl-lite)
 
To build mantis, you will also need [CMake](https://cmake.org/) version 3.9 or higher and C++17.

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
```

If SDSL is not installed in a standard location, you can try and tell CMake where to look by adding
the following to the cmake command:

```
 -DSDSL_INSTALL_PATH=<path-to-sdsl-build-location>
```

The usage for this command are as follows:

Build Mantis
-------
`mantis build` creates a colored de Bruijn graph representation that can be used to query transcripts.

``` bash
 $ ./bin/mantis build -s 20 -i raw/incqfs.lst -o raw/
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

'log-slots': The initial value for log of the number of slots in the CQF (i.e. the number of quotient bits).
 Mantis will automatically resize when the CQF reaches its capacity during the build process.
 Starting with a reasonable value is recommended so that the build process does not have to perform a bunch of resizes. Each resize operation will halt the build process and in-turn increase the overall build time.

Suggested starting values based on the size of input Squeakr files:
* 28 for a small set of genomes like a bacterial genomes.
* 30 for a large set of medium size read files.
* 33 for a large set of big read files.
Notice that these are just suggestions. You can start with a other smaller values as well.

Note: build process will open all input Squeakr files at the same time. So, please increase the limit on the number of open file handles to at least the number of input Squeakr files before running build.

Build MST
-------
`mantis mst` encodes the color information into a list of succinct 
int-vectors and bit-vectors. It creates a color graph derived from the de Bruijn graph of k-mers 
and encodes its minimum spanning tree (MST) in a format to be able to retrieve the color classes.

```bash
 $ ./bin/mantis mst -p raw/ -t 8 -k
```

The options and arguments are as follows:

```bash
SYNOPSIS
        mantis mst -p <index_prefix> [-t <num_threads>] (-k|-d)

OPTIONS
        <index_prefix>
                    The directory where the index is stored.

        <num_threads>
                    number of threads

        -k, --keep-RRR
                    Keep the previous color class RRR representation.

        -d, --delete-RRR
                    Remove the previous color class RRR representation.
```
This step is will further compress the color class representation.
It is highly recommended that you run this step after `mantis build`
since this makes your query required memory much smaller and doesn't hurt
the query time.

If you want to keep the RRR-compressed representation of color classes
after having the mst representation you require to use `-k` option
and if you want to delete this intermediate representation
you should use `-d`.

Query
-------

`mantis query` lets you query a mantis index with a set of sequences.

```bash
 $ ./bin/mantis query -p raw/ -o query.res raw/input_txns.fa
```

The options and arguments are as follows:

```bash
SYNOPSIS
        mantis query [-1] [-j] [-k <kmer>] -p <query_prefix> [-o <output_file>] <query>

OPTIONS
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
 
The output file contains the list of experiments (i.e., hits) corresponding to each queried transcript.

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
