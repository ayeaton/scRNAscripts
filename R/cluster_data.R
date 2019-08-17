
# read in metadata
metadata <- as.data.frame(data.table::fread("/home/ay1392/anna_beegfs/CITESeq/third_set/test_code/test.dim.log.metadata.csv"))
data <- as.data.frame(data.table::fread("/home/ay1392/anna_beegfs/CITESeq/third_set/test_code"))
#subset PC data
pcs <- metadata[,grep("^PC.log", colnames(metadata))]
rownames(pcs) <- metadata$cell

# calculate_clusters
res_range <- seq(0.1, 2.5, 0.1)
if (nrow(pcs) > 1000) res_range = c(res_range, 3, 4, 5, 6, 7, 8, 9)
clusters <- calculate_clusters(pcs, 25, log_file = "test",num_neighbors = 30, res = res_range)

# add to metadata
cluster_metadata <- save_seurat_metadata(data = metadata, 
                                         metadata = clusters,
                                         out_path = out_path,
                                         proj_name = proj_name,
                                         type = "cluster",
                                         log_file = log_file, 
                                         write = TRUE)
  
# plot scatter plots
loop_plot_scatter(metadata = cluster_metadata, 
                  out_path = out_path, 
                  proj_name = proj_name, 
                  log_file = log_file,
                  X = "UMAP.log1", 
                  Y = "UMAP.log2", 
                  colors_vect = paste0("res.", res_range))
 
# calculate cluster stats
cluster_stats_bar(metadata = cluster_metadata, 
                  group1 = "hash.ID",
                  group2 = "res.0.2",
                  write = TRUE)


# calculate cluster averages



# diff expression
list_groups <- list(c("hash5-sample-3762", "hash7-sample-0310"),
                    c("hash1-sample-0160","hash3-sample-3133"))


calculate_clusters <- function(pcs, num_dim, log_file, num_neighbors = 30, res = NULL){
  
  message_str <- "========== Seurat::FindNeighbors() =========="
  write_message(message_str, log_file)
  
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
  
  snn_graph <- as.Graph(snn_graph[["snn"]])
  
  message_str <- "\n\n ========== Seurat::FindClusters() ========== \n\n"
  write_message(message_str, log_file)
  
  if(is.null(res)){
    res_range <- seq(0.1, 2.5, 0.1)
    if (nrow(pcs) > 1000) res_range = c(res_range, 3, 4, 5, 6, 7, 8, 9)
  }else{
    res_range <- res
  }
  
  clusters <-  Seurat:::FindClusters.default(snn_graph,
                                algorithm = 3,
                                resolution = res_range,
                                verbose = FALSE)
  return(clusters)
}

loop_plot_scatter <- function(metadata, out_path, proj_name, log_file, X, Y, colors_vect){
  
  for(i in colors_vect){
    current_plot <- plot_scatter(metadata, 
                                 out_path, 
                                 proj_name,
                                 log_file,
                                 X, 
                                 Y, 
                                 i,
                                 write = TRUE)
  }

}

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

cluster_stats_bar <- function(metadata, group1, group2, write = FALSE){
  # TODO: pull out plots into new function
  # make barplots and output cluster stats
  summary_metadata <- metadata %>% 
    group_by(!!!syms(c(group1, group2))) %>% 
    summarize(n = n()) %>% 
    group_by(!!sym(group1)) %>% 
    mutate(pct_g2_in_g1 = n / sum(n)) %>% 
    group_by(!!sym(group2)) %>% 
    mutate(pct_g1_in_g2 = n / sum(n)) %>% 
    ungroup() 
  
  write_excel_csv(summary_metadata, path = glue("{out_path}/{proj_name}.summary.{group1}{group2}.csv"))
  
  if(write){
    
    # make both grouping variables factors
    summary_metadata %<>% mutate(!!group1 := as.factor(!!sym(group1)))
    summary_metadata %<>% mutate(!!group2 := as.factor(!!sym(group2)))
    
    mat_g1 = summary_metadata %>% 
      select(!!c(group1, group2, "pct_g1_in_g2")) %>% 
      spread(group2, 'pct_g1_in_g2', fill = 0) %>% 
      as.data.frame %>% 
      column_to_rownames(group1) %>% 
      as.matrix()
    
    hc_g1 = hclust(dist(mat_g1), method = 'ward.D2')  # clusters rows of mat
    levels_g1 = rownames(mat_g1)[order.dendrogram(as.dendrogram(hc_g1))]
    
    summary_metadata <- summary_metadata %>% 
      mutate(!!group1 := fct_relevel(!!sym(group1), levels_g1))
    
    mat_g2 = summary_metadata %>% 
      select(!!c(group1, group2, "pct_g2_in_g1")) %>% 
      spread(group1, 'pct_g2_in_g1', fill = 0) %>% 
      as.data.frame %>% 
      column_to_rownames(group2) %>% 
      as.matrix()
    
    hc_g2 = hclust(dist(mat_g2), method = 'ward.D2')  # clusters rows of mat
    levels_g2 = rownames(mat_g2)[order.dendrogram(as.dendrogram(hc_g2))]
    
    summary_metadata <- summary_metadata %>% 
      mutate(!!group2 := fct_relevel(!!sym(group2), levels_g2))

    # use levels to re-order factor
    
    group1_col <- create_color_vect(as.data.frame(summary_metadata[group1]))
    group2_col <- create_color_vect(as.data.frame(summary_metadata[group2]))
    
    summary_plots_g2 <- ggplot(summary_metadata) + 
      geom_col(aes_string(x = group2, y = "pct_g1_in_g2", fill = group1)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = group1_col,
                        name = group1) +
      ylab(glue("percent {group1} in {group2}"))
    
    summary_plots_g2_legend <- get_legend(summary_plots_g2)
    
    
    summary_plots_g1 <- ggplot(summary_metadata) + 
      geom_col(aes_string(x = group1, y = "pct_g2_in_g1", fill = group2)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = group2_col,
                        name = group2) +
      ylab(glue("percent {group2} in {group1}"))
    
    summary_plots_g1_legend <- get_legend(summary_plots_g1)
    
    
    summary_plots <- plot_grid(summary_plots_g2 + theme(legend.position = "none"),
                               summary_plots_g2_legend,
                               summary_plots_g1 + theme(legend.position = "none"),
                               summary_plots_g1_legend)
    
    ggsave(summary_plots, 
           file = glue("{out_path}/{proj_name}.{group1}{group2}.bar.png"),
           height = 10,
           width = 10)
  }
  return(summary_metadata)
}

calc_clust_averages <- function(metadata, data, group){
  # merge relevant metadata and data and row avg in group
  
  # get relevant metadata
  metadata <- metadata %>% 
    select("cell", group)
  
  # merge metadata and data on cell 
  if(nrow(metadata) != (ncol(data) - 1)) {stop("the number of cells in metadata is not the same as the number of cells in data")}
  # manitpulate data to merge with metadata
  
  data <- data %>% 
    column_to_rownames("gene") %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("cell")
    
  current_data <- full_join(metadata, data, by = "cell")
  t<- aggregate(current_data, 
            list(current_data$MULTI_classification),
            mean)
  
}

