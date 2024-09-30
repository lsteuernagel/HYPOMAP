
# sbatch run_slurm.sh ~/Documents/r_scvi_v3_42.simg merge_human_mouse_neurons/01_merge.R

# export params using 02a...

# sbatch -J crossspecies_scvi -o /beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_slurmlogs/crossspecies_scvi_slurm-%j.out -e /beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_slurmlogs/crossspecies_scvi_slurm-%j.err integration_pipeline/run_scripts/run_Python_slurm.sh ~/Documents/r_scvi_v3_42.simg merge_human_mouse_neurons/02_integrate_scVI_crossspecies.py /beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons_2/parameters_cross_species_neurons_v1.json

# sbatch -J umap_crossspecies_scvi -o /beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_slurmlogs/umap_crossspecies_scvi_slurm-%j.out -e /beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_slurmlogs/umap_crossspecies_scvi_slurm-%j.err integration_pipeline/run_scripts/run_Rscript_slurm.sh ~/Documents/r_scvi_v3_42.simg merge_human_mouse_neurons/03_calculate_umap.R 

# run 04 in rstudio

# run 05 in rstudio 