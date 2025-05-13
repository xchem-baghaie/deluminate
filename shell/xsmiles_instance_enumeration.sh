#! /usr/bin/bash
for x in $@
do
    mkdir -p "/mnt/d/Virtual_Libraries/DELs/logs/"
    echo "$x exhaustive instance enumeration started $(date +%m/%d/%Y--%H:%M:%S)"
    fn="/mnt/d/Virtual_Libraries/DELs/logs/${x}_instance_enum_console_log.txt"
    echo > "$fn" #Create an empty log file
    python DELuminate/workflows/xsmiles_instance_enumeration.py run --library_name $x --canonicalize_smiles False > "$fn"
    echo "$x exhaustive instance enumeration finished $(date +%m/%d/%Y--%H:%M:%S)"
done