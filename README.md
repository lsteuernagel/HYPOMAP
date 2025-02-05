# HYPOMAP

This repository contains the analysis scripts for the single-nucleus transcriptomics data of the human HYPOMAP publication:
[A comprehensive spatio-cellular map of the human hypothalamus](https://www.nature.com/articles/s41586-024-08504-8)
Tadross J.A., Steuernagel L., Dowsett G.K.C. et al. Nature, 2025 (doi.org/10.1038/s41586-024-08504-8)

See here for the [spatial analysis](https://github.com/georgiedowsett/HYPOMAP)

# Data

The Seurat objects (+anndata objects) used in this pipeline are deposited at University of Cambridge’s Apollo Repository [https://doi.org/10.17863/CAM.111988](https://doi.org/10.17863/CAM.111988). Newly generated raw human snRNA-seq read data ('Tadross') are deposited at the European Genome-Phenome Archive [https://ega-archive.org/](https://ega-archive.org/) under accession numbers EGAD50000000997.

The HYPOMAP snRNA-seq data is also available in an interactive [cellxgene viewer](https://cellxgene.cziscience.com/collections/d0941303-7ce3-4422-9249-cf31eb98c480).

# Structure

- [paper_figures](paper_figures/) contains the scripts used to create the final publication plots as well as all original plots and supplementary and extended data

- [integration_pipeline](integration_pipeline/) contains the full integration & harmonization pipeline

Parameters and highly variable features used can be found in integration_pipeline/parameters/

- [merge_human_mouse_neurons](merge_human_mouse_neurons/) contains the pipeline to merge the human neuronal subsets with the neuronal cells of our previously published mouse [HypoMap](https://doi.org/10.1038/s42255-022-00657-y)

This folder also includes the used parameters and cross-species highly variable features.

- [tadross_processing](tadross_processing/) contains the original scripts that were used to prepare and QC the Tadross data sets

- [siletti_processing](siletti_processing/) contains the script used to load the Siletti et al data ( [Link to publication](https://doi.org/10.1126/science.add7046) )

- [scIntegration](scIntegration/) contains the json files with parameters used in the scIntegration pipeline, as well as the results of the model evaluation. See the [scIntegration](https://github.com/lsteuernagel/scIntegration) repository for the full pipeline

- [data](data/) contains additional data files required to run the pipeline or generate figures

# Usage

This docker image has the required software and packages: [lsteuernagel/r_scvi_docker_v3:v7](https://hub.docker.com/r/lsteuernagel/r_scvi_docker_v3), used to generate the integration and data analysis in this repository. Please note that the original scvi model was trained on an older scvi-tools version which requires pandas <= 1.5.3 to re-run or to use when loading for projection.

The main scripts of the integration & processing pipeline are numbered from 0 to x and can be executed either interactively or using slurm (as commands via sbatch).

### Projection of new data 

We have included a notebook demonstrating how to project new data onto HYPOMAP here: [projection/hypoMap_scArches.ipynb](projection/hypoMap_scArches.ipynb).
This docker image has a compatible version of pandas (1.5.3), scvi (1.1.2) and juypter-notebook server to run the notebook [lsteuernagel/r_scvi_docker_rollback:v1](https://hub.docker.com/r/lsteuernagel/r_scvi_docker_rollback)


