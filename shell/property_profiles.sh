#! /usr/bin/bash
for x in $@
do
    python DELuminate/workflows/property_profiles.py run --library_name $x
done