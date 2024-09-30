
# join parents
#similarity_edgelist = dplyr::left_join(similarity_edgelist,combined_edgelist_mrtree %>% dplyr::filter(grepl("C5",from)) %>% dplyr::select(from_parent = from,to),by=c("from"="to"))
#similarity_edgelist = dplyr::left_join(similarity_edgelist,hypoMap_edgelist %>% dplyr::filter(grepl("C286",from_named)) %>% dplyr::select(to_parent = from_named,to_named),by=c("to"="to_named"))

similarity_edgelist = data.frame(from = "", to = "", from_parent = "", to_parent ="" , similarity = "")

# min_similarity_siblings and min_similarity are used in second round when adding back in.
# min_similarity_siblings checks that the matched sibling has at least this value but will add any sibling vwith value above min_similarity
prune_similarity_network = function(similarity_edgelist, adjust_similarity = FALSE,similarity_threshold = 0.6,  min_similarity_siblings = 0.75,min_similarity = 0.5) {
  
  # Input: edgelist with similarityrelatios/ auroc and parents joined
  # add edgenames
  similarity_edgelist = similarity_edgelist %>% dplyr::mutate(edgename = paste0(from,"_",to))
  
  if(adjust_similarity){
    ## get max values and subtract diff to max to weigh down multi-edges
    similarity_edgelist = similarity_edgelist %>%
      dplyr::group_by(from) %>%
      dplyr::mutate(from_max = max(similarity),from_max_diff = from_max - similarity, similarity_from = similarity - from_max_diff) %>%
      dplyr::group_by(to) %>%
      dplyr::mutate(to_max = max(similarity),to_max_diff = to_max - similarity, similarity_to = similarity - to_max_diff) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(similarity = pmin(similarity_from,similarity_to))
  }
  
  # 1 For each human cluster find the highest similarity
  from_max = similarity_edgelist %>%
    dplyr::group_by(from) %>%
    dplyr::slice_max(order_by = similarity,n = 1,with_ties = FALSE)
  
  # 2 For each mouse cluster find the highest similarity
  to_max = similarity_edgelist %>%
    dplyr::group_by(to) %>%
    dplyr::slice_max(order_by = similarity,n = 1,with_ties = FALSE)
  
  # 3 keep edges in both
  # edges_keep = intersect(from_max$edgename,to_max$edgename)
  all_edges = bind_rows(from_max,to_max) %>%
    dplyr::group_by(from,to) %>%
    dplyr::add_count(name="edge_occ") %>%
    dplyr::group_by(from) %>%
    dplyr::add_count(name="from_occ") %>%
    dplyr::group_by(to) %>%
    dplyr::add_count(name="to_occ") %>%
    dplyr::distinct(from,to,.keep_all = TRUE)
  
  # do I need to find a way to remove xY?
  all_edges_filt = all_edges %>%
    dplyr::filter(similarity >= similarity_threshold)
  message("Found ",nrow(all_edges_filt)," edges above main filtering threshold of: ",similarity_threshold)
  ## Key step: pply pruning and recount
  all_edges_filt_pruned = prune_edgelist(all_edges_filt) %>%
    dplyr::group_by(from,to) %>%
    dplyr::add_count(name="edge_occ") %>%
    dplyr::group_by(from) %>%
    dplyr::add_count(name="from_occ") %>%
    dplyr::group_by(to) %>%
    dplyr::add_count(name="to_occ") %>%
    dplyr::distinct(from,to,.keep_all = TRUE)
  
  ## this leaves us with a handul of many:many relations that are not resolved because the are between siblings and 
  ## hence kept because it is difficult to tell automatically how to sort them
  ## Next: remove  many:many relations by filtering out edges that have to_occ > 1 & from_occ > 1
  ## special case: some times both relations make it into the list (but after filtering the remaining one would actually be a valid 1:many or 1:1 relationship)
  ## resolve by always taking the minimum one
  edges_to_remove = all_edges_filt_pruned %>%
    dplyr::group_by(to) %>%
    dplyr::filter(to_occ > 1 & from_occ > 1) %>%
    dplyr::filter(similarity == min(similarity)) %>%
    dplyr::group_by(from) %>%
    dplyr::filter(to_occ > 1 & from_occ > 1) %>%
    dplyr::filter(similarity == min(similarity))
  
  ## final filtering and update occ
  all_edges_filt_pruned_rm = all_edges_filt_pruned %>%
    dplyr::filter(! edgename %in% edges_to_remove$edgename) %>%
    dplyr::group_by(from,to) %>%
    dplyr::add_count(name="edge_occ") %>%
    dplyr::group_by(from) %>%
    dplyr::add_count(name="from_occ") %>%
    dplyr::group_by(to) %>%
    dplyr::add_count(name="to_occ") %>%
    dplyr::distinct(from,to,.keep_all = TRUE)
  
  ## after all removal , cases exist where we have multiple subtypes in mouse and one in human (and vice versa) and only one ismapped and the otehr missed the cutoff (e.g. due to similarityrelation adjustement)
  ## we iterate over all nodes thathave zero edges and check if there is a good mapping of any of its direct siblings to excatly one other cluster
  ## if yes that edges is added back in
  
  nodes_from_without_edge = setdiff(similarity_edgelist$from, all_edges_filt_pruned_rm$from)
  nodes_to_without_edge = setdiff(similarity_edgelist$to, all_edges_filt_pruned_rm$to)
  
  
  edges_to_add = c()
  for(from_node in nodes_from_without_edge){
    from_node_parent = similarity_edgelist$from_parent[similarity_edgelist$from == from_node][1]
    siblings_in_pruned = all_edges_filt_pruned_rm[all_edges_filt_pruned_rm$from_parent %in% from_node_parent ,] %>% as.data.frame() # & all_edges_filt_pruned_rm$similarity >= min_similarity_siblings
    # if there are siblings with min_similarity
    if(nrow(siblings_in_pruned) > 0){
      # get sibling cluster targets
      sibling_target = unique(siblings_in_pruned$to)
      # if there is exactly one cluster the siblings are matched to AND # check that minc_similarity is fullfiled
      if(length(sibling_target) == 1 & siblings_in_pruned$similarity[1] > min_similarity_siblings){
        edge_to_add = similarity_edgelist$edgename[similarity_edgelist$from == from_node & similarity_edgelist$to == sibling_target & similarity_edgelist$similarity >= min_similarity ]
        edges_to_add = c(edges_to_add,edge_to_add)
      }
    }
  }
  # repeat for to nodes
  for(to_node in nodes_to_without_edge){
    to_node_parent = similarity_edgelist$to_parent[similarity_edgelist$to == to_node][1]
    siblings_in_pruned = all_edges_filt_pruned_rm[all_edges_filt_pruned_rm$to_parent %in% to_node_parent,] %>% as.data.frame()
    # if there are siblings with min_similarity
    if(nrow(siblings_in_pruned) > 0){
      # get sibling cluster targets
      sibling_target = unique(siblings_in_pruned$from)
      # if there is exactly one cluster the siblings are matched to
      if(length(sibling_target) == 1 & siblings_in_pruned$similarity[1] > min_similarity_siblings){
        edge_to_add = similarity_edgelist$edgename[similarity_edgelist$to == to_node & similarity_edgelist$from == sibling_target & similarity_edgelist$similarity >= min_similarity ]
        edges_to_add = c(edges_to_add,edge_to_add)
      }
    }
  }
  
  ## these edges back in
  edge_to_bind = similarity_edgelist[similarity_edgelist$edgename %in% edges_to_add & ! similarity_edgelist$edgename %in% edges_to_remove,]
  message("Adding ",nrow(edge_to_bind)," edges back in.")
  ## add and update numbers
  all_edges_filt_pruned_final = bind_rows(all_edges_filt_pruned_rm,edge_to_bind) %>%
    dplyr::group_by(from,to) %>%
    dplyr::add_count(name="edge_occ") %>%
    dplyr::group_by(from) %>%
    dplyr::add_count(name="from_occ") %>%
    dplyr::group_by(to) %>%
    dplyr::add_count(name="to_occ") %>%
    dplyr::distinct(from,to,.keep_all = TRUE)
  
  return(all_edges_filt_pruned_final)
  
}

##########
### Prune edgelist
##########

# helper for prune_similarity_network

# removes multi-edges :
# if multi-edges are clean 1:many connections they remain, because the neighbor will always be the max neighbor
# if the neighbor also has a connection to somewhere else, only the maximum connection will be kept (if that is the node itself the edge stays in)
# unless the maximum connection is to a direct sibling of the node. then it also stays in

prune_edgelist = function(edgelist_to_prune,weight_col = "similarity",edge_col = "edgename"){
  edgelist_to_prune = as.data.frame(edgelist_to_prune)
  edges_to_delete = c()
  all_from_nodes = unique(edgelist_to_prune$from)
  # iterate over all from nodes
  for(from_node in all_from_nodes){
    from_node_parent = edgelist_to_prune$from_parent[edgelist_to_prune$from == from_node][1]
    current_neighbors =  edgelist_to_prune$to[edgelist_to_prune$from == from_node]
    # if more than one neighbor
    if(length(current_neighbors) > 1){
      # check all neighbors
      for(current_neighbor in current_neighbors){
        # only keep neighbor if this is its highest weighted edge OR if the highest weight is a direct neighbor of current from node
        current_edges = edgelist_to_prune[edgelist_to_prune$to == current_neighbor,]
        max_weight = current_edges[,weight_col][current_edges[,weight_col] == max(current_edges[,weight_col])]
        max_node = current_edges$from[current_edges[,weight_col] == max(current_edges[,weight_col])]
        max_node_parent = current_edges$from_parent[current_edges[,weight_col] == max(current_edges[,weight_col])]
        # select edges to delete
        if(from_node != max_node ){#& max_node_parent != from_node_parent){
          edges_to_delete = c(edges_to_delete,edgelist_to_prune[edgelist_to_prune$to == current_neighbor & edgelist_to_prune$from == from_node,edge_col])
        }
      }
    }
  }
  
  # repeat for to
  all_to_nodes = unique(edgelist_to_prune$to)
  # iterate over all from nodes
  for(to_node in all_to_nodes){
    to_node_parent = edgelist_to_prune$to_parent[edgelist_to_prune$to == to_node][1]
    current_neighbors =  edgelist_to_prune$from[edgelist_to_prune$to == to_node]
    # if more than one neighbor
    if(length(current_neighbors) > 1){
      
      # check all neighbors
      for(current_neighbor in current_neighbors){
        # only keep neighbor if this is its highest weighted edge OR if the highest weight is a direct neighbor of current from node
        current_edges = edgelist_to_prune[edgelist_to_prune$from == current_neighbor,]
        max_weight = current_edges[,weight_col][current_edges[,weight_col] == max(current_edges[,weight_col])]
        max_node = current_edges$to[current_edges[,weight_col] == max(current_edges[,weight_col])]
        max_node_parent = current_edges$to_parent[current_edges[,weight_col] == max(current_edges[,weight_col])]
        # select edges to delete
        if(to_node != max_node ){#& max_node_parent != to_node_parent){
          edges_to_delete = c(edges_to_delete,edgelist_to_prune[edgelist_to_prune$from == current_neighbor & edgelist_to_prune$to == to_node,edge_col])
        }
      }
    }
  }
  edgelist_to_prune2 = edgelist_to_prune[!edgelist_to_prune[,edge_col] %in% edges_to_delete,]
  message("Removing ",length(edges_to_delete)," edges")
  return(edgelist_to_prune2)
}


#' Find reciprocal top hits, stratifying results by study.
#'
#' This function looks for reciprocal top hits in each target study
#' separately, allowing for as many reciprocal top hits as target studies.
#' This is the recommended function for extracting top hits.
#'
#' @param auroc matrix of celltype-to-celltype AUROC scores
#' (output from \code{\link{MetaNeighborUS}})
#' @param threshold AUROC threshold, must be between [0,1]. Default is 0.9.
#'  Only top hits above this threshold are included in the result table.
#' @param n_digits Number of digits for AUROC values in the result table. Set
#'  to "Inf" to skip rounding.
#' @param collapse_duplicates Collapse identical pairs of cell types (by
#'  default), effectively averaging AUROCs when reference and target roles are
#'  reversed. Setting this option to FALSE makes it easier to filter results
#'  by study or cell type.
#'  If collapse_duplicates is set to FALSE, "Celltype_1" is the
#'  reference cell type and "Celltype_2" is the target cell type (relevant
#'  if MetaNeighborUS was run with symmetric_output = FALSE). 
#'
#' @return Function returns a dataframe with cell types that are either reciprocal best 
#' matches, and/or those with AUROC values greater than or equal to threshold 
#' value
#'
#' @examples
#' data(mn_data)
#' var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
#' aurocs = MetaNeighborUS(var_genes = var_genes, 
#'                         dat = mn_data, 
#'                         study_id = mn_data$study_id,
#'                         cell_type = mn_data$cell_type)
#' top_hits = topHitsByStudy(aurocs)
#' top_hits
#'
#' @seealso \code{\link{topHits}}
#' @export
#'
topHitsByStudy = function(auroc, threshold = 0.9, n_digits = 2, collapse_duplicates = TRUE) {
  `%>%` <- dplyr::`%>%`
  #could be sped up by finding same study AUROC's when in matrix form and masking them to NA (as in the old topHits function)
  result <- tibble::as_tibble(auroc, rownames = "ref_cell_type") %>%
    tidyr::pivot_longer(cols = -ref_cell_type,
                        names_to = "target_cell_type",
                        values_to = "auroc") %>%
    dplyr::mutate(ref_study = getStudyId(ref_cell_type),
                  target_study = getStudyId(target_cell_type)) %>%
    dplyr::filter(ref_study != target_study) %>%
    dplyr::group_by(ref_cell_type, target_study) %>%
    dplyr::filter(auroc == max(auroc, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-ref_study, -target_study) %>%
    dplyr::mutate(is_reciprocal = is_reciprocal_top_hit(.)) %>%
    dplyr::filter(auroc >= threshold) 
  
  if (collapse_duplicates) {
    result <- result %>%
      dplyr::group_by(ref_cell_type, target_cell_type) %>%
      dplyr::mutate(pair_id = paste(sort(c(ref_cell_type, target_cell_type)),
                                    collapse = "")) %>%
      dplyr::group_by(pair_id) %>%
      dplyr::summarize(ref_cell_type = dplyr::first(ref_cell_type),
                       target_cell_type = dplyr::first(target_cell_type),
                       auroc = mean(auroc),
                       is_reciprocal = dplyr::first(is_reciprocal)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-pair_id)
  }
  
  # final formatting
  result <- result %>%
    dplyr::arrange(desc(auroc)) %>%
    dplyr::mutate(auroc = round(auroc, n_digits)) %>%    
    dplyr::mutate(Match_type = ifelse(is_reciprocal,
                                      "Reciprocal_top_hit",
                                      paste0("Above_", threshold))) %>%
    dplyr::select("Study_ID|Celltype_1" = ref_cell_type,
                  "Study_ID|Celltype_2" = target_cell_type,
                  "AUROC" = auroc,
                  Match_type)
  return(result)
}

# helper
is_reciprocal_top_hit = function(best_hits) {
  `%>%` <- dplyr::`%>%`
  best_hits <- best_hits %>%
    dplyr::select(-auroc)
  reverse_hits <- best_hits %>%
    dplyr::select(target_cell_type = ref_cell_type,
                  reciprocal_cell_type = target_cell_type)
  reciprocal_best_hits <- dplyr::inner_join(best_hits, reverse_hits,
                                            by = "target_cell_type") %>%
    dplyr::filter(ref_cell_type == reciprocal_cell_type) %>%
    dplyr::select(-reciprocal_cell_type) %>%
    tibble::add_column(is_reciprocal = TRUE)
  result <- dplyr::left_join(best_hits, reciprocal_best_hits,
                             by = c("ref_cell_type", "target_cell_type")) %>%
    tidyr::replace_na(replace = list(is_reciprocal = FALSE)) %>%
    dplyr::pull(is_reciprocal)
  return(result)
}