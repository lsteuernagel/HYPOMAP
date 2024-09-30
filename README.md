# HYPOMAP

This repository contains the analysis scripts for the single-nucleus transcriptomics aspects of the human HYPOMAP dataset and paper.

See here for the [spatial analysis](https://github.com/georgiedowsett/HYPOMAP)

# Data

The Seurat objects used in this pipeline are available from : TBD

# Structure

- [paper_figures](paper_figures/) contains the scripts used to create the final publication plots as well as all original plots and supplementary and extended data

- [integration_pipeline](integration_pipeline/) contains the full integration & harmonization pipeline

- [merge_human_mouse_neurons](merge_human_mouse_neurons/) contains the pipeline to merge the human neuronal subsets with the neuronal cells of our previously published mouse [HypoMap](https://doi.org/10.1038/s42255-022-00657-y)

- [tadross_processing](tadross_processing/) contains the original scripts 01-06 that were used to prepare and QC the Tadross data sets

- [siletti_processing](siletti_processing/) contains the script used to load the Siletti et al data ( [Link to publication](https://doi.org/10.1126/science.add7046) )

- [scIntegration](scIntegration/) contains the json files with parameters used in the scIntegration pipeline, as well as the results of the model evaluation. See the [scIntegration](https://github.com/lsteuernagel/scIntegration) repository for the full pipeline

- [data](data/) contains additional data files required to run the pipeline or generate figures

# Usage

This docker image should has the required software and packages: [Image](https://hub.docker.com/r/lsteuernagel/r_scvi_docker_v3/tags)

The main scripts of the integration & processing pipeline are numbered from 0 to x and can be executed either interactively or using slurm (as commands via sbatch).
