#!/bin/bash

FILES=path/to/conda/pkgs/*.tar.bz2
for f in $FILES
do
  anaconda upload $f
done