#! /usr/bin/bash
for x in $@
do
    python DELuminate/workflows/visualize_chemical_space.py run --filename $x --embedding ECFP2 --dimensionality_reduction UMAP --plot_config_string seaborn
    python DELuminate/workflows/visualize_chemical_space.py run --filename $x --embedding ECFP4 --dimensionality_reduction UMAP --plot_config_string seaborn
    python DELuminate/workflows/visualize_chemical_space.py run --filename $x --embedding ECFP6 --dimensionality_reduction UMAP --plot_config_string seaborn
done