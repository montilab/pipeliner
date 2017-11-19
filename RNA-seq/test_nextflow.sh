#!/bin/sh
for file in test_configs/*.config; do
    filename=${file##*/}
    testname=${filename%.*}
    pathtoconfig=${PWD}/${file}

    # Printing Commands
    printf "_%.0s" {1..40}
    printf "\nTest: ${testname}\n"
    echo 'Time:' `date +%Y-%m-%d\ %H:%M:%S`
    echo 'Status: Running...'
    
    # Running Nextflow
    mkdir -p test_logs
    rm -f -- test_logs/${testname}.txt
    ./nextflow main.nf -c $pathtoconfig > test_logs/${testname}.txt
done
