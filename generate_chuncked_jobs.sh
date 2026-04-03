#!/bin/bash
# Generate chunked job input from data_input.txt
# Each directory is split into chunks of FILES_PER_JOB files
# Output: data_input_chunked.txt with columns: region base_dir job_index files_per_job

INPUT="${1:?Usage: $0 <input_file> [files_per_job]}"
FILES_PER_JOB=${2:-50}
OUTPUT="${INPUT%.txt}_chunked.txt"

> "$OUTPUT"
total_jobs=0

while IFS=' ' read -r region base_dir; do
    [[ -z "$region" || -z "$base_dir" ]] && continue
    nfiles=$(find "$base_dir" -name '*.root' 2>/dev/null | wc -l)
    njobs=$(( (nfiles + FILES_PER_JOB - 1) / FILES_PER_JOB ))
    if [[ $njobs -eq 0 ]]; then
        echo "Warning: no .root files in $base_dir, skipping"
        continue
    fi
    for (( j=0; j<njobs; j++ )); do
        echo "$region $base_dir $j $FILES_PER_JOB" >> "$OUTPUT"
    done
    echo "$base_dir: $nfiles files -> $njobs jobs"
    total_jobs=$((total_jobs + njobs))
done < "$INPUT"

echo "Written $total_jobs jobs to $OUTPUT"
