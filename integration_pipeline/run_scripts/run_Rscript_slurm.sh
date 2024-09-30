#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=56
#SBATCH --partition=highmem
#SBATCH --nice

# need to get 3 variables from call
# singularity image
singularity_image=$1
# $method_script : script that should be executed (depend on method)
method_script=$2
# $param_file : specifies which object, assay etc are used and which parameters of the method
param_file=$3

# Run script
echo "singularity exec "$singularity_image" Rscript "$method_script" "$param_file
srun singularity exec $singularity_image Rscript $method_script $param_file
