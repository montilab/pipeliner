#!/bin/bash

# Must first log in to anaconda (anaconda login)

path=path/to/conda/pkgs/
for file in ${path}*.tar.bz2; do
    #echo ${path}${file##*/}
    anaconda upload ${path}${file##*/}
done