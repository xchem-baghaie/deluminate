#! /usr/bin/bash
for x in $@
do
    python DELuminate/visualization/plot.py run $x
done