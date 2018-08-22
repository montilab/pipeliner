#!/bin/bash

# Must first log in to anaconda (anaconda login)

FILES=path/to/conda/pkgs/*.tar.bz2
for file in ${path}*.tar.bz2; do
    #echo ${path}${file##*/}
    anaconda upload ${path}${file##*/}
done