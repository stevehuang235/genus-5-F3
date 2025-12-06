#!/bin/bash

dir="./data_unfiltered"  # Replace with your directory
files=( "$dir"/* )        # Array of all files in the directory
max_jobs=5
num_files=${#files[@]}

for (( i=0; i<num_files; i+=max_jobs )); do
    # Launch a batch of up to max_jobs concurrently
    for (( j=i; j<i+max_jobs && j<num_files; j++ )); do
        base_file=$(basename "${files[j]}")
        magma -b Input:="${base_file}" curve_check.m &
    done
    # Wait for the current batch to finish before starting the next
    wait
done