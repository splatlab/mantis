# mantis
Mantis: A Fast, Small, and Exact Large-Scale Sequence-Search Index

Overview
--------

Mantis is a space-efficient data structure that can be used to index thousands of raw-
read experiments and facilitate large-scale sequence searches on those experiments. Mantis uses counting quotient
filters instead of Bloom filters, enabling rapid index builds and queries, small indexes, and exact results, i.e., no
false positives or negatives. Furthermore, Mantis is also a colored de Bruijn graph representation, so it supports fast
graph traversal and other topological analyses in addition to large-scale sequence-level searches.

A preprint of the paper describing mantis is available [here](https://www.biorxiv.org/content/biorxiv/early/2017/11/10/217372.full.pdf).

API
--------
* `mantis build`: builds a mantis index from a collection of (squeakr) CQF files.
* `mantis query`: query k-mers in the mantis index.

Build
-------

Library dependencies (given version or higher):
 - libboost-dev 1.58.0.1ubuntu1
 - zlib1g-dev 1:1.2.8.dfsg-2ubuntu4
 - sdsl

The CQF code uses two new instructions to implement select on machine words
introduced in intel's Haswell line of CPUs. However, there is also an alternate
implementation of select on machine words to work on CPUs older than Haswell.
To build on an older hardware (older than Haswell) use "NH=1" as a make argument.

`mantis build` creates a colored de Bruijn graph representation that can be used to query transcripts.

```bash
 $ make mantis
 $ ./mantis build -i raw/incqfs.lst -c raw/experiment_cutoffs.lst -o raw/
```
The usage for this command are as follows:

```
build the mantis index
Usage: ./src/mantis build [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -i,--input-list TEXT        file containing list of input filters
  -c,--cutoff-list TEXT       file containing list of experiment-specific cutoffs
  -o,--output TEXT            directory where results should be written
```

 Following are the arguments to coloreddbg:
 - input cqf files: a list of input cqf files. This is a list of squeakr output files that are generated after running squeakr on input experiments. We have provided two sample cqf files in data dir.
 - experiment cutoffs: The cutoff value for each input cqf file corresponding to the experiment. The cutoff value is the minimum count that a k-mer needs to be considered in the search.
 - prefix: prefix filepath where all the output files will be written.

Query
-------

`mantis query` lets you query a mantis index with a set of sequences.

```bash
 $ make mantis
 $ mantis query -i raw/ -o query.res raw/input_txns.fa
```

The options and arguments are as follows:

```bash
query the mantis index
Usage: ./src/mantis query [OPTIONS] query

Positionals:
  query TEXT                  Prefix of input files.

Options:
  -h,--help                   Print this help message and exit
  -j,--json                   Write the output in JSON format
  -p,--input-prefix TEXT      Prefix of input files.
  -o,--output TEXT            Where to write query output.
```

 The command takes the following options :
 - `--input-prefix,-p`: the directory where the output of coloreddbg command is present.
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
- Rob Johnson <rob@cs.stonybrook.edu>
- Rob Patro <rob.patro@cs.stonybrook.edu>
