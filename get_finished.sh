#!/bin/bash

if [ "$#" -lt 1 ]; then
    "Usage: $0 <analysis_type>"
    exit 1;
fi

LOG_DIR="/storage/gpfs_small_files/kloe/gamrat/log/$1"
ROOT_DIR="/storage/gpfs_small_files/kloe/gamrat/root_files/${1^^}_SIGNAL_NoSmearing"

rm -f "finished/${1}_finished.txt"
touch "finished/${1}_finished.txt"

find "$LOG_DIR" -name "cut.analysis.log" -printf "%h\n" | sed 's|.*/||' | while read num; do
    ls -d "$ROOT_DIR"/*_"${num}.root" 2>/dev/null >> "finished/${1}_finished.txt"
done
