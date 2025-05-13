#! /usr/bin/bash
for x in $@
do
    mkdir -p "/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/${x}/virtual/output/XBBs_only/"
    echo "$x virtual library enumeration started $(date +%m/%d/%Y--%H:%M:%S)"
    fn="/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/${x}/virtual/output/XBBs_only/${x}_console_log.txt"
    echo > "$fn" #Create an empty log file
    python DELuminate/workflows/virtual_presynthesis.py run --library_name $x --use_xbbs_only True --apply_mw_filter False --prohibit_exact_overlap_in_instance_enum False > "$fn"
    echo "$x virtual library enumeration finished $(date +%m/%d/%Y--%H:%M:%S)"
done