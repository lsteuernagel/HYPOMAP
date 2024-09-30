#!/bin/bash
#SBATCH --job-name='humanHypo'
#SBATCH --output='/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/other_logs/%j-script.out'
#SBATCH --error='/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/other_logs/%j-script.err'
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=56
#SBATCH --partition=blade-b
#SBATCH --nice

# need to get 3 variables from call
# singularity image
singularity_image=$1
# $method_script : script that should be executed (depend on method)
method_script=$2
# params
#param_file 

# Run script
echo "singularity exec "$singularity_image" Rscript "$method_script
srun singularity exec $singularity_image Rscript $method_script

# submit to slurm:
# singularity_path = "~/Documents/r_scvi_015.simg"
# script_path = "R/myscript.R"
# system(paste0("sbatch R/run_slurm.sh ",singularity_path," ",script_path),intern = TRUE)

# from terminal (when in project dir):
# sbatch run_slurm.sh ~/Documents/r_scvi_015.simg R/myscript.R
# sbatch run_slurm.sh ~/Documents/r_scvi_015.simg R/myscript.R
