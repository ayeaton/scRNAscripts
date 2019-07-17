
calculate_clusters = function(seurat_obj, num_dim, num_neighbors = 30) {
# perform graph-based clustering and tSNE
# specify neighbors for UMAP and FindNeighbors (default is 30 in Seurat 2 and 3 pre-release)
  # check if number of dimensions seems reasonable
  if (num_dim < 5) stop("too few dims: ", num_dim)
  if (num_dim > 50) stop("too many dims: ", num_dim)
  
  s_obj = seurat_obj
  
  message("\n\n ========== Seurat::FindNeighbors() ========== \n\n")
  
  message("assay: ", DefaultAssay(s_obj))
  message("num dims: ", num_dim)
  
  # construct a Shared Nearest Neighbor (SNN) Graph for a given dataset
  s_obj =
    FindNeighbors(
      s_obj, dims = 1:num_dim, k.param = num_neighbors,
      graph.name = "snn", compute.SNN = TRUE, force.recalc = TRUE
    )
  
  message("\n\n ========== Seurat::FindClusters() ========== \n\n")
  
  message("initial metadata fields: ", str_c(colnames(s_obj@meta.data), collapse = ", "))
  
  # resolutions for graph-based clustering
  # increased resolution values lead to more clusters (recommendation: 0.6-1.2 for 3K cells, 2-4 for 33K cells)
  res_range = seq(0.1, 2.5, 0.1)
  if (ncol(s_obj) > 1000) res_range = c(res_range, 3, 4, 5, 6, 7, 8, 9)
  
  # algorithm: 1 = original Louvain; 2 = Louvain with multilevel refinement; 3 = SLM
  # identify clusters of cells by SNN modularity optimization based clustering algorithm
  s_obj = FindClusters(s_obj, algorithm = 3, resolution = res_range, graph.name = "snn", verbose = FALSE)
  
  # remove "seurat_clusters" column that is added automatically (added in v3 late dev version)
  s_obj@meta.data = s_obj@meta.data %>% select(-seurat_clusters)
  
  message("new metadata fields: ", str_c(colnames(s_obj@meta.data), collapse = ", "))
  
  # create a separate sub-directory for cluster resolution plots
  clusters_dir = "clusters-resolutions"
  if (!dir.exists(clusters_dir)) dir.create(clusters_dir)
  
  # for calculated cluster resolutions: remove redundant (same number of clusters), rename, and plot
  res_cols = str_subset(colnames(s_obj@meta.data), "snn_res")
  res_cols = sort(res_cols)
  res_num_clusters_prev = 1
  for (res in res_cols) {
    
    # proceed if current resolution has more clusters than previous and less than the color scheme length
    res_vector = s_obj@meta.data[, res] %>% as.character()
    res_num_clusters_cur = res_vector %>% n_distinct()
    if (res_num_clusters_cur > res_num_clusters_prev && res_num_clusters_cur < length(colors_clusters)) {
      
      # check if the resolution still has original labels (characters starting with 0)
      if (min(res_vector) == "0") {
        
        # convert to character vector
        s_obj@meta.data[, res] = as.character(s_obj@meta.data[, res])
        # relabel identities so they start with 1 and not 0
        s_obj@meta.data[, res] = as.numeric(s_obj@meta.data[, res]) + 1
        # pad with 0s to avoid sorting issues
        s_obj@meta.data[, res] = str_pad(s_obj@meta.data[, res], width = 2, side = "left", pad = "0")
        # pad with "C" to avoid downstream numeric conversions
        s_obj@meta.data[, res] = str_c("C", s_obj@meta.data[, res])
        # encode as a factor
        s_obj@meta.data[, res] = factor(s_obj@meta.data[, res])
        
      }
      
      # resolution value based on resolution column name
      res_val = sub("snn_res\\.", "", res)
      
      # plot file name
      res_str = gsub("\\.", "", res)
      dr_filename = glue("{clusters_dir}/dr.{DefaultAssay(s_obj)}.{num_dim}.{res_str}.clust{res_num_clusters_cur}")
      
      s_obj = plot_clusters(seurat_obj = s_obj, resolution = res_val, filename_base = dr_filename)
      
      # add blank line to make output easier to read
      message(" ")
      
    } else {
      
      # remove resolution if the number of clusters is same as previous
      s_obj@meta.data = s_obj@meta.data %>% select(-one_of(res))
      
    }
    
    # update resolution cluster count for next iteration
    res_num_clusters_prev = res_num_clusters_cur
    
  }
  
  message("updated metadata fields: ", str_c(colnames(s_obj@meta.data), collapse = ", "))
  
  save_metadata(seurat_obj = s_obj)
  
  return(s_obj)
  
}

calculate_cluster_stats = function(seurat_obj, label) {
  
  message("\n\n ========== calculate cluster stats ========== \n\n")
  
  message("cluster names: ", str_c(levels(seurat_obj), collapse = ", "))
  
  # compile relevant cell metadata into a single table
  seurat_obj$cluster = Idents(seurat_obj)
  metadata_tbl = seurat_obj@meta.data %>% rownames_to_column("cell") %>% as_tibble() %>%
    select(cell, num_UMIs, num_genes, pct_mito, sample_name = orig.ident, cluster)
  tsne_tbl = seurat_obj[["tsne"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
  umap_tbl = seurat_obj[["umap"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
  cells_metadata = metadata_tbl %>% full_join(tsne_tbl, by = "cell") %>% full_join(umap_tbl, by = "cell")
  cells_metadata = cells_metadata %>% arrange(cell)
  write_excel_csv(cells_metadata, path = glue("metadata.{label}.csv"))
  
  # get number of cells split by cluster and by sample
  summary_cluster_sample =
    cells_metadata %>%
    select(cluster, sample_name) %>%
    mutate(num_cells_total = n()) %>%
    group_by(sample_name) %>%
    mutate(num_cells_sample = n()) %>%
    group_by(cluster) %>%
    mutate(num_cells_cluster = n()) %>%
    group_by(cluster, sample_name) %>%
    mutate(num_cells_cluster_sample = n()) %>%
    ungroup() %>%
    distinct() %>%
    mutate(
      pct_cells_cluster = num_cells_cluster / num_cells_total,
      pct_cells_cluster_sample = num_cells_cluster_sample / num_cells_sample
    ) %>%
    mutate(
      pct_cells_cluster = round(pct_cells_cluster * 100, 1),
      pct_cells_cluster_sample = round(pct_cells_cluster_sample * 100, 1)
    ) %>%
    arrange(cluster, sample_name)
  
  # get number of cells split by cluster (ignore samples)
  summary_cluster = summary_cluster_sample %>% select(-contains("sample")) %>% distinct()
  write_excel_csv(summary_cluster, path = glue("summary.{label}.csv"))
  
  # gene expression for an "average" cell in each identity class (averaging and output are in non-log space)
  cluster_avg_exp = AverageExpression(seurat_obj, assay = "RNA", verbose = FALSE)[["RNA"]]
  cluster_avg_exp = cluster_avg_exp %>% round(3) %>% rownames_to_column("gene") %>% arrange(gene)
  write_excel_csv(cluster_avg_exp, path = glue("expression.mean.{label}.csv"))
  
  Sys.sleep(1)
  
  # export results split by sample if multiple samples are present
  num_samples = cells_metadata %>% pull(sample_name) %>% n_distinct()
  if (num_samples > 1) {
    
    # number of cells split by cluster and by sample
    write_excel_csv(summary_cluster_sample, path = glue("summary.{label}.per-sample.csv"))
    
    # cluster averages split by sample
    sample_avg_exp = AverageExpression(seurat_obj, assay = "RNA", add.ident = "orig.ident", verbose = FALSE)[["RNA"]]
    sample_avg_exp = sample_avg_exp %>% round(3) %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
    write_excel_csv(sample_avg_exp, path = glue("expression.mean.{label}.per-sample.csv"))
    
  }
  
}
