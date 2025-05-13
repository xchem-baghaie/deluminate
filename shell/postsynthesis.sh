#! /usr/bin/bash
for x in $@
do
    mkdir -p "/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/${x}/postsynthesis/output/"
    echo "$x postsynthesis enumeration started $(date +%m/%d/%Y--%H:%M:%S)"
    fn="/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/${x}/postsynthesis/output/${x}_console_log.txt"
    echo > "$fn" #Create an empty log file
    python DELuminate/workflows/postsynthesis.py run --library_name $x > "$fn"
    echo "$x postsynthesis enumeration finished $(date +%m/%d/%Y--%H:%M:%S)"
done