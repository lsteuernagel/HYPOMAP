##########
### Load
##########

# source("R/harmonization_functions.R")
singularity_path = "~/Documents/r_scvi_411.simg"#"~/Documents/r_scvi_v3_42.simg"

# direct output and logs to some files on the local filesystem:
# where to store temporary json's with params for jobs:
param_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_params_2/"
# where to save log files --> use this path in the slurm.sh files!
log_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_slurmlogs_2/"

# load json file with all other information
params_harmonization = jsonlite::read_json("integration_pipeline/parameters/parameters_human_hypo_v2.json")
# if some fields are lists --> unlist
params_harmonization = lapply(params_harmonization,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

### try to creat dir if necessary:
system(paste0("mkdir -p ",paste0(param_path)))
system(paste0("mkdir -p ",paste0(log_path)))
system(paste0("mkdir -p ",paste0(params_harmonization$harmonization_folder_path)))


# define helper function
writeList_to_JSON = function (list_with_rows, filename){
  jsonfile = jsonlite::toJSON(list_with_rows, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  writeLines(jsonfile, con = paste0(filename))
}

##########
### [1] prepare
##########

# set params
param_set = params_harmonization
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"prepare_harmonization_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "integration_pipeline/run_scripts/01_prepare_harmonization.R"
# set sbatch params:
jobname = paste0("prepare_harmonization_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = ""
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_1 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [2] Integrate with scVI
##########

# set additional parameters for scvi
param_set = params_harmonization

# make unique id:
job_id=digest::digest(param_set)
param_set$job_id = job_id
# write to JSON as transfer file
param_file = paste0(param_path,"harmonization_scvi_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# execute job
script_path = "integration_pipeline/run_scripts/02_integrate_scVI.py"
# set sbatch params:
jobname = paste0("harmonization_scvi_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_1)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Python_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
#output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," integration_pipeline/run_scripts/run_Python_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_2 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [3] Basic harmonization
##########

# set params
param_set = params_harmonization
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"basic_harmonization_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "integration_pipeline/run_scripts/03_basic_harmonization.R"
# set sbatch params:
jobname = paste0("basic_harmonization_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_2)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
#output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_3 = stringr::str_remove(output_message,pattern = "Submitted batch job ")


##########
### [4] Initial clustering using python leiden
##########

# set additional parameters for clustering
param_set = params_harmonization
param_set$target_clusterN = param_set$target_clusterN_initial
param_set$start_res = param_set$start_res_initial
param_set$end_res = param_set$end_res_initial
param_set$step_size = param_set$step_size_initial
param_set$include_low_res = param_set$include_low_res_initial
param_set$clustering_key_name = param_set$clustering_key_name_initial 

# make unique id:
job_id=digest::digest(param_set)
param_set$job_id = job_id
# write to JSON as transfer file
param_file = paste0(param_path,"leiden_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# execute job
script_path = "integration_pipeline/run_scripts/04_basic_leiden_clustering.py"
# set sbatch params:
jobname = paste0("leiden_scanpy_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_3)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Python_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_4 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [5] Annotation based on inital harmonization
##########

# set params
param_set = params_harmonization
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"basic_annotation_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "integration_pipeline/run_scripts/05_annotate_full_data.R"
# set sbatch params:
jobname = paste0("annotation_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_4)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_5 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [6] Initial marker detection
##########

# set params
param_set = params_harmonization
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"basic_markers_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "integration_pipeline/run_scripts/06_basic_marker_detection.R"
# set sbatch params:
jobname = paste0("basic_markers_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_5)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_6 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [7] subset object
##########

# set params
param_set = params_harmonization
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"subset_celltypes_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "integration_pipeline/run_scripts/07_subset_object.R"
# set sbatch params:
jobname = paste0("subset_object_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = "-1" #c(slurm_id_5)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_7 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [MULTIPLE] RUN on different subsets
##########

# this is manually defined in the above script --> need to adjust that script if different subset are desired !! This will break if different names are used!
subset_names = c("AstroEpendymal","Oligodendrocytes","HypoNeuron","NonNeuron")





##########
### [MULTIPLE] Manually adjust on different subsets
##########

##########
### [9] Define cluster levels for mr_tree
##########


## Run selection script and save file in

# - load various files from split clusering job
# add to object
# make basic selection of 6 levels
# export selected levels
# also export mrged all

# parameter_list$clusters_for_mrtree_file

##########
### [MULTIPLE] Run all mrtree steps per subset
##########

for( subname in subset_names){
  
  message(subname)
  
  ##########
  ### [8] calculate KNN
  ##########
  
  slurm_id_8 = c()
  # set additional parameters for clustering
  param_set = params_harmonization
  
  # load subset specifc
  param_set$harmonization_folder_path = paste0(param_set$harmonization_folder_path,"/",subname,"/")
  param_set$new_name_suffix = paste0(param_set$new_name_suffix,"_",subname)
  param_set$clustering_key_name = paste0(subname,"_clusters")
  
  # make unique id:
  job_id=digest::digest(param_set)
  param_set$job_id = job_id
  # write to JSON as transfer file
  param_file = paste0(param_path,"knn_params_",subname,"_",job_id,".json")
  writeList_to_JSON(list_with_rows = param_set,filename = param_file)
  # execute job
  script_path = "integration_pipeline/run_scripts/08_calculate_knn.py"
  # set sbatch params:
  jobname = paste0("knn_scanpy_",job_id)
  outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
  errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
  dependency_ids = c(slurm_id_7)
  output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Python_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
  slurm_id_8 = stringr::str_remove(output_message,pattern = "Submitted batch job ")
  
  ##########
  ### [9] Clustering for hierarchical tree after curation
  ##########
  
  # which resolutions:
  if(subname=="HypoNeuron"){
    resolutions = params_harmonization$resolutions_neurons
  }else{
    resolutions = params_harmonization$resolutions_basic 
  }
  
  # set additional parameters for clustering
  param_set = params_harmonization
  
  # load subset specific
  param_set$harmonization_folder_path = paste0(param_set$harmonization_folder_path,"/",subname,"/")
  param_set$new_name_suffix = paste0(param_set$new_name_suffix,"_",subname)
  param_set$clustering_key_name = paste0(subname,"_clusters")
  param_set$n_cores =100
  
  # I split ino multiple jobs to be faster: Run a maximum of 10 resolutions per job.
  slurm_ids_9 = c()
  for( leiden_resolution in resolutions){
    ## full clustering extra
    param_set$leiden_resolution = leiden_resolution
    param_set$additional_clustering_suffix = paste0("res_",leiden_resolution,"_repeat_",param_set$leiden_repeat)
    
    # make unique id:
    job_id=digest::digest(param_set)
    param_set$job_id = job_id
    # write to JSON as transfer file
    param_file = paste0(param_path,"leiden_params_",job_id,".json")
    writeList_to_JSON(list_with_rows = param_set,filename = param_file)
    # execute job
    script_path = "integration_pipeline/run_scripts/09_multiple_leiden_clustering.py"
    # set sbatch params:
    jobname = paste0("leiden_scanpy_",param_set$leiden_resolution,"_",job_id)
    outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
    errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
    dependency_ids = c(slurm_id_8)
    output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Python_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
    slurm_ids_9[leiden_resolution] = stringr::str_remove(output_message,pattern = "Submitted batch job ")
  }
  
  ##########
  ### [10] Consensus clustering
  ##########
  #for( subname in subset_names){  
  # which resolutions:
  if(subname=="HypoNeuron"){
    resolutions = params_harmonization$resolutions_neurons
  }else{
    resolutions = params_harmonization$resolutions_basic 
  }
  
  # set additional parameters for clustering
  param_set = params_harmonization
  
  # load subset specific
  param_set$harmonization_folder_path = paste0(param_set$harmonization_folder_path,"/",subname,"/")
  param_set$new_name_suffix = paste0(param_set$new_name_suffix,"_",subname)
  param_set$clustering_key_name = paste0(subname,"_clusters")
  
  # I split ino multiple jobs to be faster: Run a maximum of 5 resolutions per job.
  n_res = 1
  res_groups = split(resolutions, ceiling(seq_along(resolutions)/n_res))
  counter = 10
  for( res_group in res_groups){
    message("Running ",paste0(res_group,collapse="|"))
    ## full clustering extra
    param_set$leiden_resolutions_to_run = res_group
    # make unique id:
    job_id=digest::digest(param_set)
    param_set$job_id = job_id
    # write to JSON as transfer file
    param_file = paste0(param_path,"consensus_params_",counter,"_",job_id,".json")
    writeList_to_JSON(list_with_rows = param_set,filename = param_file)
    # execute job
    script_path = "integration_pipeline/run_scripts/10_consensus_clustering.R"
    # set sbatch params:
    jobname = paste0("consensus_",subname,"_",counter,"_",job_id)
    outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
    errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
    dependency_ids = "-1"#c(slurm_ids_9)
    output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
    slurm_ids_9[counter] = stringr::str_remove(output_message,pattern = "Submitted batch job ")
    counter = counter+1
  }
  #}
  
  ##########
  ### [11] Hierarchical tree
  ##########
  
  for( subname in subset_names){  
    slurm_id_11 = c()
    # set additional parameters for subset
    param_set = params_harmonization
    
    # load subset specifc
    param_set$harmonization_folder_path = paste0(param_set$harmonization_folder_path,"/",subname,"/")
    param_set$new_name_suffix = paste0(param_set$new_name_suffix,"_",subname)
    param_set$clustering_key_name = paste0(subname,"_clusters")
    
    # set params for markers
    param_set$clustering_folder = paste0(param_set$harmonization_folder_path,"consensus/") 
    param_set$load_from_file =FALSE
    param_set$n_cores = 56
    
    # hardcoded which clusterings to use
    clusterings= list(
      "AstroEpendymal" = c("celltype_annotation","res_0.25_repeat_100","res_3_repeat_100"), # or res_0.25_repeat_100 instead of 0.75 ?
      "Oligodendrocytes" =c("res_0.25_repeat_100","res_0.75_repeat_100","res_1.5_repeat_100"),
      "HypoNeuron" = c("res_0.05_repeat_100","res_0.5_repeat_100","res_4_repeat_100","res_50_repeat_100"),
      "NonNeuron" =  c("res_0.01_repeat_100","res_0.1_repeat_100","res_2_repeat_100")
    )
    param_set$clusterings_to_use = clusterings[[subname]]
    
    # make unique id:
    job_id=digest::digest(param_set)
    # write to JSON as transfer file
    param_file = paste0(param_path,"mrtree_construction_params_",job_id,".json")
    writeList_to_JSON(list_with_rows = param_set,filename = param_file)
    # not a loop
    script_path = "integration_pipeline/run_scripts/11_build_mrtree.R"
    # set sbatch params:
    jobname = paste0("mrtree_construction_",subname,"_",job_id)
    outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
    errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
    dependency_ids = "-1"#c(slurm_ids_10) 
    output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
    slurm_id_11 = stringr::str_remove(output_message,pattern = "Submitted batch job ")
  
  ##########
  ### [12] Hierachical tree markers
  ##########

    # set additional parameters for subset
    param_set = params_harmonization
    
    # load subset specifc
    param_set$harmonization_folder_path = paste0(param_set$harmonization_folder_path,"/",subname,"/")
    param_set$new_name_suffix = paste0(param_set$new_name_suffix,"_",subname)
    param_set$clustering_key_name = paste0(subname,"_clusters")
    
    param_set$n_cores_markers = 4
    param_set$marker_suffix = "raw"
    param_set$start_node = "all" # "all" for everything
    
    # make unique id:
    job_id=digest::digest(param_set)
    # write to JSON as transfer file
    param_file = paste0(param_path,subname,"_mrtree_markers_params_",job_id,".json")
    writeList_to_JSON(list_with_rows = param_set,filename = param_file)
    # not a loop
    script_path = "integration_pipeline/run_scripts/12_markers_mrtree.R"
    # set sbatch params:
    jobname = paste0(subname,"_mrtree_markers_",job_id)
    outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
    errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
    dependency_ids = c(slurm_id_11) #"-1"#c(slurm_ids_11) 
    output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
    slurm_id_12 = stringr::str_remove(output_message,pattern = "Submitted batch job ")
  
  
  ##########
  ### [13] Hierachical tree pruning
  ##########
  prune_rounds = 5
 # for( subname in subset_names){  
  
  for( pr in 1:prune_rounds){
    # subname="AstroEpendymal"
    # set additional parameters for subset
    param_set = params_harmonization
    
    # load subset specifc
    param_set$harmonization_folder_path = paste0(param_set$harmonization_folder_path,"/",subname,"/")
    param_set$new_name_suffix = paste0(param_set$new_name_suffix,"_",subname)
    param_set$clustering_key_name = paste0(subname,"_clusters")
    if(pr == 1){
      param_set$mrtree_result_file = paste0(param_set$harmonization_folder_path,param_set$clustering_key_name,"_mrtree_results",".rds")
      param_set$marker_suffix = "raw"
      param_set$old_prefix = "K"
      dependency_ids = c(slurm_id_12)#"-1"#c(slurm_id_12) ## Might need to adjust this !
    }else{
      param_set$mrtree_result_file = paste0(param_set$harmonization_folder_path,param_set$clustering_key_name,"_pruned_",pr-1,"_mrtree_clustering_results",".rds")
      param_set$marker_suffix = paste0("pruned_",pr-1)
      param_set$old_prefix = "C"
      dependency_ids = slurm_id_14 # depend on last round markers
    }
    param_set$mrtree_result_file_output = paste0(param_set$harmonization_folder_path,param_set$clustering_key_name,"_pruned_",pr,"_mrtree_clustering_results",".rds")
    
    # move to parama script !
    param_set$pct_threshold = 90 # for donor
    param_set$min_cells = 50 # if below --> merge with neighbor
    param_set$min_specificity = 1.5 # min specificity for a sibling marker to count
    param_set$max_pvalue_prune = 0.001 # max pvalue for a sibling marker to count
    param_set$min_sibling_markers = 5 # how many sibling markers are required to not merge
    
    # set
    param_set$merge_marker_based =TRUE
    param_set$merge_sample_based=TRUE
    param_set$start_node = "all"
    
    # make unique id:
    job_id=digest::digest(param_set)
    # write to JSON as transfer file
    param_file = paste0(param_path,subname,"mrtree_pruning_",pr,"_params_",job_id,".json")
    writeList_to_JSON(list_with_rows = param_set,filename = param_file)
    # not a loop
    script_path = "integration_pipeline/run_scripts/13_prune_mrtree_flex.R"
    # set sbatch params:
    jobname = paste0("mrtree_pruning_",pr,"_",job_id)
    outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
    errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
    output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
    slurm_id_13 = stringr::str_remove(output_message,pattern = "Submitted batch job ")
  
    ##########
    ### [14] Hierachical tree markers
    ##########
    
    # set additional parameters for subset
    param_set = params_harmonization
    
    # load subset specifc
    param_set$harmonization_folder_path = paste0(param_set$harmonization_folder_path,"/",subname,"/")
    param_set$new_name_suffix = paste0(param_set$new_name_suffix,"_",subname)
    param_set$clustering_key_name = paste0(subname,"_clusters")
    
    # set up
    param_set$n_cores_markers = 4
    param_set$marker_suffix = paste0("pruned_",pr)
    param_set$start_node = "all" # "all" for everything
    param_set$mrtree_result_file = paste0(param_set$harmonization_folder_path,param_set$clustering_key_name,"_pruned_",pr,"_mrtree_clustering_results",".rds")
    
    # make unique id:
    job_id=digest::digest(param_set)
    # write to JSON as transfer file
    param_file = paste0(param_path,"mrtree_markers_params_",job_id,".json")
    writeList_to_JSON(list_with_rows = param_set,filename = param_file)
    # not a loop
    script_path = "integration_pipeline/run_scripts/12_markers_mrtree.R"
    # set sbatch params:
    jobname = paste0("mrtree_pruned_markers_",job_id)
    outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
    errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
    dependency_ids = c(slurm_id_13)
    # start
    output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
    slurm_id_14 = stringr::str_remove(output_message,pattern = "Submitted batch job ")
  }
  
}

##########
###  Put back together
##########


##########
### [15] Stitch together again
##########

# set params
param_set = params_harmonization
  
# set additional params
param_set$subset_names = subset_names
param_set$prune_round = prune_rounds

# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"merge_neuron_nonneuron_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "integration_pipeline/run_scripts/15_merge_subsets.R"
# set sbatch params:
jobname = paste0("merge_neuron_nonneuron_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = "-1"#c(slurm_id_13)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_15 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [16] C0 level marker detection
##########

# set params
param_set = params_harmonization

# set some relevant params:
param_set$seurat_object_markers = paste0(param_set$harmonization_folder_path, "human_hypo_combined" ,"/",  "human_hypo_combined" ,".rds")
param_set$cluster_column_markers = "C0"
param_set$harmonization_folder_path = paste0(param_set$harmonization_folder_path, "human_hypo_combined" ,"/")
param_set$basic_marker_filename = paste0("_markers_",param_set$cluster_column_markers)
  
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"c0_markers_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "integration_pipeline/run_scripts/06_basic_marker_detection.R"
# set sbatch params:
jobname = paste0("c0_markers_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = "-1"#c(slurm_id_16)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_16 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [17] Plot subsets and merged object
##########

# set params
param_set = params_harmonization

# set additional params
param_set$subset_names = subset_names
param_set$prune_round = prune_rounds

# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"plot_mrtree_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "integration_pipeline/run_scripts/17_plot_mrtree_results.R"
# set sbatch params:
jobname = paste0("plot_mrtree_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = "-1"#c(slurm_id_14)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_17 = stringr::str_remove(output_message,pattern = "Submitted batch job ")


##########
### [18] Annotate object
##########

# set params
param_set = params_harmonization
param_set$subname = "human_hypo_combined"
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"annotate_mrtree_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "integration_pipeline/run_scripts/18_annotate_mrtree.R"
# set sbatch params:
jobname = paste0("annotate_mrtree_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = "-1"#c(slurm_id_14)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_18 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [19] Build final object
##########

# set params
param_set = params_harmonization
# set
param_set$new_name_suffix = "human_hypo_combined"
param_set$harmonization_folder_path = paste0(param_set$harmonization_folder_path,param_set$new_name_suffix ,"/")
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"build_annotated_object_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "integration_pipeline/run_scripts/19_build_final_object.R"
# set sbatch params:
jobname = paste0("build_final_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = "-1"#c(slurm_id_13)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_19 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [20] Plot merged object
##########

# set params
param_set = params_harmonization

# set additional params
param_set$new_name_suffix = "human_hypo_combined"
param_set$harmonization_folder_path = paste0(param_set$harmonization_folder_path,param_set$new_name_suffix ,"/")
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"plot_anno_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "integration_pipeline/run_scripts/20_plot_annotated_results.R"
# set sbatch params:
jobname = paste0("plot_anno_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = "-1"#c(slurm_id_14)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes integration_pipeline/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_20 = stringr::str_remove(output_message,pattern = "Submitted batch job ")
