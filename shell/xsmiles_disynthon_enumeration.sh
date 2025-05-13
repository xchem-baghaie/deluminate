#! /usr/bin/bash
for x in $@
do
    mkdir -p "/mnt/d/XSMILES/Disynthon_Enumeration/${x}/"
    echo "$x disynthon enumeration for X-Viewer started $(date +%m/%d/%Y--%H:%M:%S)"
    fn="/mnt/d/XSMILES/Disynthon_Enumeration/${x}/${x}_console_log.txt"
    echo > "$fn" #Create an empty log file
    python DELuminate/workflows/xsmiles_disynthon_enumeration.py run --library_name $x --enum_destination 'X-Viewer' > "$fn"
    echo "$x disynthon enumeration for X-Viewer finished $(date +%m/%d/%Y--%H:%M:%S)"
    
    # mkdir -p "/mnt/d/Virtual_Libraries/DELs/logs/"
    # echo "$x disynthon enumeration for S3 started $(date +%m/%d/%Y--%H:%M:%S)"
    # fn="/mnt/d/Virtual_Libraries/DELs/logs/${x}_disynthon_enum_console_log.txt"
    # echo > "$fn" #Create an empty log file
    # python DELuminate/workflows/xsmiles_disynthon_enumeration.py run --library_name $x --enum_destination 'S3' > "$fn"
    # echo "$x disynthon enumeration for S3 finished $(date +%m/%d/%Y--%H:%M:%S)"
done