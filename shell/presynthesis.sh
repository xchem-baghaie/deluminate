#! /usr/bin/bash
for x in $@
do
    mkdir -p "/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/${x}/presynthesis/output/"
    echo "$x presynthesis enumeration started $(date +%m/%d/%Y--%H:%M:%S)"
    fn="/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/${x}/presynthesis/output/${x}_console_log.txt"
    echo > "$fn" #Create an empty log file
    python DELuminate/workflows/presynthesis.py run --library_name $x > "$fn"
    echo "$x presynthesis enumeration finished $(date +%m/%d/%Y--%H:%M:%S)"
done