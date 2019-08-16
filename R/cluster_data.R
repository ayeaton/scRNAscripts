
meta_tbl <- 
pcs <- 

clusters <- calculate_clusters(pcs, 25)

clusters_metadata <- save_seurat_metadata(data = meta_tbl,
                                         metadata = clusters,
                                         out_path = out_path,
                                         proj_name = proj_name, 
                                         log_file = log_file,
                                         type = "unbiased.clust",
                                         write = TRUE)

calculate_clusters <- function(pcs, num_dim, num_neighbors = 30){
  
  message_str <- "========== Seurat::FindNeighbors() =========="
  write_message(mesage_str)
  
  if (num_dim < 5) stop("too few dims: ", num_dim)
  if (num_dim > 50) stop("too many dims: ", num_dim)
  
  snn_graph <- Seurat:::FindNeighbors.default(
    pcs[,1:num_dim],
    distance.matrix = FALSE,
    k.param = num_neighbors,
    compute.SNN = TRUE, 
    prune.SNN = 1/15,
    nn.eps = 0,
    force.recalc = TRUE)
  
  # I NEED TO CHECK THIS, but i think than snn is the one I want, not nn
  snn_graph <- as.Graph(snn_graph[["snn"]])
  
  message_str <- "\n\n ========== Seurat::FindClusters() ========== \n\n"
  write_message(mesage_str)
  
  res_range = seq(0.1, 2.5, 0.1)
  if (nrow(pcs) > 1000) res_range = c(res_range, 3, 4, 5, 6, 7, 8, 9)
  
  clusters <-  Seurat:::FindClusters.default(snn_graph,
                                algorithm = 3,
                                resolution = res_range,
                                verbose = FALSE)
  return(clusters)
}

# loop over resolutions
plot_scatter() 

# differential expression between defined groups
# input a list of groups to iterate over

data <- genes as rownames

list_groups <- list(c("hash5-sample-3762", "hash7-sample-0310"),
                    c("hash1-sample-0160","hash3-sample-3133"))

diff_exp <- function(data, metadata, metadata_column, list_groups, test.use = "wilcox"){
  
  diff_exp_stats <- tibble(
    p_val = numeric(),
    avg_logFC = numeric(),
    p_val_adj = numeric(), 
    group_1 = character(),
    group_2 = character()
  )
  
  for(current_group in list_groups){
    
   cell_group1 <-  metadata %>% 
      rownames_to_column("cell") %>% 
      filter(get(metadata_column) == current_group[1]) %>% 
      select("cell") 
   
   cell_group2 <-  metadata %>% 
     rownames_to_column("cell") %>% 
     filter(get(metadata_column) == current_group[2]) %>% 
     select("cell") 

   current_comparison <- Seurat:::FindMarkers.default(
     object = as.matrix(data),
     reduction = NULL,
     slot = "data",
     cells.1 = cell_group1$cell,
     cells.2 = cell_group2$cell,
     logfc.threshold = 0.001,
     test.use = test.use,
     min.pct =  0)
   
   current_comparison_filt <- current_comparison %>% 
     select(p_val, avg_logFC, p_val_adj) %>% 
     filter(avg_logFC > 1 | avg_logFC < -1) %>%   
     filter(p_val_adj < 0.1) %>% 
     mutate(group_1 = rep(current_group[1])) %>% 
     mutate(group_2 = rep(current_group[2]))
  
   diff_exp_stats <- rbind(diff_exp_stats, current_comparison_filt)
  }
  return(diff_exp_stats)
}



