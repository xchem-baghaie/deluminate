#! /usr/bin/bash
for x in $@
do
    echo "$x XSMILES export started $(date +%m/%d/%Y--%H:%M:%S)"
    python DELuminate/workflows/xsmiles_export_json.py run --library_name $x
    echo "$x XSMILES export finished $(date +%m/%d/%Y--%H:%M:%S)"
done