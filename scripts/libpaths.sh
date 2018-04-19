#!/bin/bash
IFS=":"
for line in $LIBRARY_PATH; do
  echo -L$line '\\'
done

