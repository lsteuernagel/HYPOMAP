library(Seurat)
library(ggplot2)

param_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons_2/parameters_cross_species_neurons_v2.json"
param_list = jsonlite::read_json(param_file)
param_list$batch_var
param_list$max_epochs
param_list$new_name_suffix
# make path for umaps
umap_pathname = "evaluated_umaps/"
umap_pathname = paste0(param_list$harmonization_folder_path,umap_pathname)

# load seurat object
hypothalamus_neurons_cross_species = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons_2/hypothalamus_neurons_cross_species.rds")

##########
### Add embeddings
##########

### Add scvi embedding
embed = c( "cross_species_neurons_scvi_1_scVI_reduction" )
embed_short = substr(embed, start = 1, stop = 45)

# make path for umaps
param_file = "merge_human_mouse_neurons/parameters_cross_species_neurons_scvi_v1.json"
param_list = jsonlite::read_json(param_file)
# find file
embed_file = c()
embed_file = list.files(paste0(param_list$harmonization_folder_path),pattern = embed,recursive = TRUE,full.names = TRUE)
embed_file = embed_file[!grepl("umap",embed_file)]
if(length(embed_file) != 1){
  message("Found ",length(embed_file)," files for embed. Please provide exactly one valid file." )
  next
}

# load
message("Load ",embed," ...")
embedding = read_embedding(filename_withpath = embed_file,seurat_object = hypothalamus_neurons_cross_species)

#add
dimred <- Seurat::CreateDimReducObject(
  embeddings = as.matrix(embedding),
  stdev = as.numeric(apply(embedding, 2, stats::sd)),
  assay = "RNA",
  key = embed
)
# add
hypothalamus_neurons_cross_species_updated@reductions[[embed]] = dimred

### Add umap of scvi embedding

emebeddings_to_load = c("cross_species_neurons_scvi_1_scVI_reduction")
#a1 =all_evaluation_results[all_evaluation_results$reduction %in% emebeddings_to_load,]
# iterate over embeddings
all_embed_short = c()
for ( embed in emebeddings_to_load){
  embed_short = substr(embed, start = 1, stop = 45)
  all_embed_short = c(all_embed_short,embed_short)
  hypothalamus_neurons_cross_species[[paste0("umap_",embed_short)]] = readRDS(paste0(umap_pathname,"umap_",embed,".rds"))
}

##########
### Add clustering results from human neurons using the final object
##########

human_hypo_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_3/human_hypo_combined/"
human_hypo_combined = readRDS(paste0(human_hypo_path,"human_hypo_combined.rds"))
combined_edgelist_mrtree = data.table::fread(paste0(human_hypo_path,"human_hypo_combined_edgelist_mrtree_annotated.txt"),data.table = F)

cluster_assignments = human_hypo_combined@meta.data %>% dplyr::select(Cell_ID,C0,C1,C2,C3,C4,C0_named,C1_named,C2_named,C3_named,C4_named)

# join
temp_meta = dplyr::left_join(hypothalamus_neurons_cross_species@meta.data[,c(1:11,14:17,55,58:95)],cluster_assignments,by="Cell_ID",suffix = c("mm","hs"))
nrow(temp_meta) == nrow(hypothalamus_neurons_cross_species@meta.data)
rownames(temp_meta) = temp_meta$Cell_ID

## add
hypothalamus_neurons_cross_species@meta.data = temp_meta

all_features = jsonlite::read_json("/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons_2/feature_set_cross_species.json")
all_features_df = data.frame(gene = unlist(all_features))

##########
### save updated object
##########



##########
### Make plots
##########

## make species plot
all_embed_short2 = all_embed_short[c(1,3,4)]
all_embed_short2 = all_embed_short[c(4,5)]
all_embed_short2 = all_embed_short
plot_list=list()
require(ggplot2)
for(embed_short in all_embed_short2){
  plot_list[[embed_short]] = Seurat::DimPlot(hypothalamus_neurons_cross_species,reduction = paste0("umap_",embed_short),
                                             cols = c("#E69F00","#009E73"),raster = TRUE,raster.dpi = c(1536,1536),shuffle = TRUE,group.by = "species")+
                                               NoLegend()+NoAxes()+ggtitle(embed_short)
}

cowplot::plot_grid(plotlist = plot_list,ncol=2)

featureplot_list=list()
for(embed_short in all_embed_short2){
  featureplot_list[[embed_short]] = Seurat::FeaturePlot(hypothalamus_neurons_cross_species,reduction = paste0("umap_",embed_short),
                                             raster = TRUE,raster.dpi = c(1536,1536),order = TRUE,features = "ANXA2")+
    NoLegend()+NoAxes()+ggtitle(embed_short)
}

cowplot::plot_grid(plotlist = featureplot_list,ncol=2)

####
names(hypothalamus_neurons_cross_species)
reduc="umap_cross_species_neurons_scvi_1_scVI_reduction"
p1=DimPlot(hypothalamus_neurons_cross_species,group.by = "species",reduction = reduc,shuffle = TRUE,raster.dpi = c(1024,1024),cols = c("#E69F00","#009E73" ))+ggtitle(NULL) #+NoLegend()
p2=FeaturePlot(hypothalamus_neurons_cross_species,features = "GHRH",reduction = reduc,raster.dpi = c(1024,1024),order = TRUE)
#p2
cowplot::plot_grid(p1,p2)

FeaturePlot(hypothalamus_neurons_cross_species,features = "HCRT",split.by = "species",order = TRUE,raster.dpi = c(1024,1024))

hvgs = unlist(jsonlite::read_json("/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons/feature_set_cross_species.json"))

FeaturePlot(harmonized_seurat_object,features = "SLC17A7",raster.dpi = c(1024,1024),order = TRUE)+NoAxes()

p1 = FeaturePlot(harmonized_seurat_object,features = "GHRH",raster.dpi = c(2048,2048),order = TRUE,pt.size = 1.3)+NoAxes()
# ggsave(filename = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/dump/human_ghrh.pdf"),
#        plot = p1, "pdf",dpi=300,width=320,height = 300,units="mm")


