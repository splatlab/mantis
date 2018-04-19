#!/bin/bash
. load_deps_spack.sh  # load mantis dependencies with spack and save environment variables

TOOLS=~/tools
mkdir -p $TOOLS/src
git clone https://github.com/simongog/sdsl-lite.git $TOOLS/src/sdsl-lite
$TOOLS/src/sdsl-lite/install.sh $TOOLS

git clone https://github.com/effect/mantis.git $TOOLS/mantis

LIBRARY_PATH=$TOOLS/lib:$LIBRARY_PATH
./patch_makefile.sh $TOOLS/mantis/Makefile
cd $TOOLS/mantis
MYINCL=$TOOLS/include make NH=1
