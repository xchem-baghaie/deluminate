#! /usr/bin/bash
for x in $@
do
    python DELuminate/workflows/enumerate_instances.py run --library_name $x
done