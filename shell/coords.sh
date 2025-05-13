#! /usr/bin/bash
for x in $@
do
    python DELuminate/visualization/calculate_coordinates.py run $x
done