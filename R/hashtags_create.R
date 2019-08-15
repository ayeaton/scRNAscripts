 
# create s obj with and without use of hashtags

# TODO: remove libraries 
suppressPackageStartupMessages({
  library(magrittr)
  library(glue)
  library(Seurat)
  library(future)
  library(Matrix)
  library(tidyverse)
  library(data.table)
  library(ggplot2); theme_set(theme_bw())
  library(cowplot)
  library(scales)
  library(pheatmap)
  library(RColorBrewer)
  library(ggsci)
  library(eulerr)
  library(UpSetR)
  library(GGally)
  library(compositions)
  library(reticulate)
})


assemble_seurat_obj_hto <- function(data_path, # path to 10x data /data_path/outs/
                                    sample_names, # name of samples, can be more than 1
                                    out_path, # path to deposit outputs
                                    num_dim, # number of dimensions to use for PCA, UMAP, and TSNE
                                    HTO_file, # path to HTO file,
                                    ADT_file, # path to ADT file
                                    sct = FALSE, # use ScTransform or not
                                    proj_name = NULL, # name of project
                                    log_file = NULL, # name of log file
                                    min_genes = NULL,
                                    max_genes = NULL,
                                    max_mt = NULL) {

  # Set project name to sample name if there is only one sample, and if projec
  # name is null
  if (is.null(proj_name) & length(sample_names) == 1) {
    proj_name = sample_names
  } else if (is.null(proj_name)) {
    proj_name = "proj"
    message("WARNING: Project name is proj. Is this what you want? 
            This is an easy wayfor files to be over-written")
  } else{
    proj_name = proj_name
  }

  # If log_file is null, set it to project name
  if (is.null(log_file)) {
    log_file = glue("{out_path}/{proj_name}_create.log")
  }

  # Set up log file
  write(glue("Starting analysis for {proj_name}"), 
        file = log_file, 
        append = TRUE)

  # Write message if proj_name is "proj"
  if (proj_name == "proj") {
    write(glue("WARNING: Project name is proj. Is this what you want?
               This is an easy way for files to be over-written"), 
          file = log_file, 
          append = TRUE)
  }

  # Log parameters
  write(glue("data_path = {data_path}
             sample_names = {sample_names}
             out_path = {out_path}
             num_dim = {num_dim}
             HTO_file = {HTO_file}
             sct = {sct}
             proj_name = {proj_name}
             log_file = {log_file}
             min_genes = {min_genes}
             max_genes = {max_genes}
             max_mt = {max_mt}"),
        file = log_file,
        append = TRUE)

  # unfiltered ----------------------
  # read in count data from 10x
  counts_mat <- load_sample_counts_matrix(sample_names = sample_names,
                                          data_path = data_path,
                                          log_file = log_file)

  # take the counts matrix and make a seurat object
  seurat_obj <- create_seurat_obj(counts_mat = counts_mat, 
                                    out_path = out_path, 
                                    proj_name = proj_name, 
                                    log_file = log_file)
  
   # save counts mat
   save_counts_matrix(seurat_obj = seurat_obj,
                      out_path = out_path,
                      proj_name = proj_name,
                      log_file = log_file,
                      type = "counts.unfiltered", 
                      assay = "RNA",
                      slot = "counts")
  
  # save unfiltered seurat metadata
  save_seurat_metadata(seurat_obj = seurat_obj, 
                       out_path = out_path,
                       proj_name = proj_name, 
                       log_file = log_file,
                       type = "unfiltered")
    
  # plot qc of the unfiltered seurat object
  plot_qc_seurat(seurat_obj = seurat_obj,
                 out_dir = out_dir,
                 proj_name = proj_name,
                 type = "unfiltered")

  # qc filter ----------------------------
  # filter seurat object for min genes, max genes and max mito pct
  seurat_obj <- filter_data(seurat_obj, 
                              out_dir = out_dir, 
                              proj_name = proj_name, 
                              log_file = log_file,
                              min_genes = min_genes, 
                              max_genes = max_genes, 
                              max_mt = max_mt)
  
  # plot qc plots for filtered seurat obj
  plot_qc_seurat(seurat_obj = seurat_obj,
                 out_dir = out_dir,
                 proj_name = proj_name,
                 type = "filtered")
  
  # save counts filtered data
  save_counts_matrix(seurat_obj = seurat_obj,
                     out_path = out_path,
                     proj_name = proj_name,
                     log_file = log_file,
                     type = "counts.filtered")
  
  #save filtered metadata
  save_seurat_metadata(seurat_obj = seurat_obj,
                       out_path = out_path,
                       log_file = log_file,
                       proj_name = proj_name, 
                       type = "filtered")
  
   
  # Hashtags -------------------------
  
  # clean hto data
  hto_data <- read_remove.unmapped(HTO_file)
  
  # here the proj name and the sample name have to match
  seurat_obj <- create_seurat_obj_hto(seurat_obj = seurat_obj, 
                                          HTO_counts = hto_data, 
                                      proj_name = proj_name, 
                                      out_dir = data_path,
                                      log_file = log_file, 
                                      hto_demux_quantile = 0.99,
                                      multi_demux_quantile = 0.7)

  # save normed HTO
  save_counts_matrix(seurat_obj,
                     out_path, 
                     proj_name,
                     type = "HTO.norm", 
                     log_file,
                     assay = "HTO",
                     slot = "data")
    
  # save metadata with HTO
  save_seurat_metadata(seurat_obj = seurat_obj,
                       out_path = out_path,
                       log_file = log_file,
                       proj_name = proj_name, 
                       type = "HTO")
  
  # plot qc plots for filtered seurat obj with HTO
  plot_qc_seurat(seurat_obj = seurat_obj,
                 out_dir = out_dir,
                 proj_name = proj_name,
                 type = "filtered.HTO",
                 group = "hash.ID")

  # plot HTO related plots
  seurat_plot_hto(seurat_obj = seurat_obj, 
                  out_path = out_path,
                  proj_name = proj_name)
  
  # Add ADT data -----------------------
  adt_data <- read_remove.unmapped(ADT_file)
  
  seurat_obj <- create_seurat_obj_adt(seurat_obj = seurat_obj, 
                                      out_dir = out_dir,
                                      ADT_counts = adt_data,
                                      proj_name = proj_name,
                                      log_file = log_file)
  
  # save normed ADTs
  save_counts_matrix(seurat_obj,
                     out_path, 
                     proj_name,
                     type = "ADT.norm", 
                     log_file,
                     assay = "ADT",
                     slot = "data")
  
  
  # Subset singlets ---------------------
  # manual_hto(seurat_obj, out_path, proj_name)
  
  # subset singlets
  seurat_obj <- subset_singlets(seurat_obj = seurat_obj,
                                    method = "HTO_demux", 
                                    log_file = log_file)
  
  # save normed ADTs singlet
  save_counts_matrix(seurat_obj,
                     out_path, 
                     proj_name,
                     type = "ADT.norm.singlet", 
                     log_file,
                     assay = "ADT",
                     slot = "data")
  
  # save normed HTOs singlet
  save_counts_matrix(seurat_obj,
                     out_path, 
                     proj_name,
                     type = "HTO.norm.singlet", 
                     log_file,
                     assay = "HTO",
                     slot = "data")
  
  # dimred based on ADTs ----------
  ADT_data <- seurat_obj@assays$ADT@data
  
  ADT_dimred <- run_dimensionality_reduction(ADT_data,num_pcs, num_dim, log_file, dim_red_suffix = ".ADT")
  
  seurat_obj_adt <-  add_dim_red_seurat(seurat_obj, ADT_dimred, dim_red_suffix = ".ADT")
  
  ADT_dim_metadata <- save_seurat_metadata(seurat_obj = seurat_obj_adt,
                       dim_red_list = ADT_dimred,
                       out_path = out_path,
                       proj_name = proj_name, 
                       log_file = log_file,
                       type = "dim.adt")
  
  plot_scatter(metadata = ADT_dim_metadata,
               out_path = out_path,
               proj_name = proj_name,
               log_file = log_file,
               X = "UMAP.ADT1",
               Y = "UMAP.ADT2",
               color = "orig.ident"))

  plot_scatter(metadata = ADT_dim_metadata,
             out_path = out_path,
             proj_name = proj_name,
             log_file = log_file,
             X = "tSNE.ADT1",
             Y = "tSNE.ADT2",
             color = "orig.ident"))

  plot_scatter(metadata = ADT_dim_metadata,
             out_path = out_path,
             proj_name = proj_name,
             log_file = log_file,
             X = "PC.ADT1",
             Y = "PC.ADT2",
             color = "orig.ident"))
  
  
  # log normalize data ----------
  seurat_obj_log <- log_normalize_data(seurat_obj = seurat_obj, 
                                       out_path = out_path,
                                       proj_name = proj_name,
                                       log_file = log_file)
  # save counts mat
  save_counts_matrix(seurat_obj = seurat_obj_log,
                     out_path = out_path,
                     proj_name = proj_name,
                     log_file = log_file,
                     type = "log.norm")
  
  # calculate variance and plot 
  seurat_obj_log <- calculate_variance(seurat_obj = seurat_obj_log,
                                   out_path = out_path,
                                   proj_name = proj_name,
                                   log_file = log_file)
  
  # run PCA, TSNE and UMAP
  seurat_log_dimred <- run_dimensionality_reduction(seurat_obj_log, num_pcs, num_dim, log_file, dim_red_suffix = ".log")
  
  # add dim red to seurat object 
  seurat_obj_log <-  add_dim_red_seurat(seurat_obj_log, seurat_log_dimred, dim_red_suffix = ".log")
  
  
  #save dim red in metadata 
  log_dim_metadata <- save_seurat_metadata(seurat_obj = seurat_obj_log,
                       dim_red_list = seurat_log_dimred,
                       out_path = out_path,
                       proj_name = proj_name, 
                       log_file = log_file,
                       type = "dim.log")
  
  log_dim_metadata <- save_seurat_metadata(seurat_obj = seurat_obj_log,
                                           dim_red_list = seurat_log_dimred,
                                           out_path = out_path,
                                           proj_name = proj_name, 
                                           log_file = log_file,
                                           type = "dim.adt")
  
  plot_scatter(metadata = log_dim_metadata,
                 out_path = out_path,
                 proj_name = proj_name,
                 log_file = log_file,
                 X = "UMAP.log1",
                 Y = "UMAP.log2",
                 color = "orig.ident"))
  
  plot_scatter(metadata = log_dim_metadata,
               out_path = out_path,
               proj_name = proj_name,
               log_file = log_file,
               X = "tSNE.log1",
               Y = "tSNE.log2",
               color = "orig.ident"))
  
  plot_scatter(metadata = log_dim_metadata,
               out_path = out_path,
               proj_name = proj_name,
               log_file = log_file,
               X = "PC.log1",
               Y = "PC.log2",
               color = "orig.ident"))

  
  saveRDS(seurat_obj_log,
          file = glue("{proj_name}.log.seurat_obj.rds"))
  
  # sct normalize data -----------------
  if(sct == TRUE){
    # sctransform data ( should save the sctransform in a new data slot)
    seurat_obj_sct <- sctransform_data(seurat_obj, out_path = out_path, log_file = log_file, proj_name = proj_name)
    
    # run PCA, TSNE and UMAP
    seurat_sct_dimred <- run_dimensionality_reduction(seurat_obj_sct, num_pcs, num_dim, log_file, dim_red_suffix = ".sct")
    
    # add dim red to seurat object 
    seurat_obj_sct <-  add_dim_red_seurat(seurat_obj_sct, seurat_sct_dimred, dim_red_suffix = ".sct")
    
    
    #save dim red in metadata 
    sct_dim_metadata <- save_seurat_metadata(seurat_obj = seurat_obj_sct,
                           dim_red_list = seurat_sct_dimred,
                           out_path = out_path,
                           proj_name = proj_name, 
                           log_file = log_file,
                           type = "dim.sct")
    
    
    plot_scatter(metadata = sct_dim_metadata,
                 out_path = out_path,
                 proj_name = proj_name,
                 log_file = log_file,
                 X = "UMAP.sct1",
                 Y = "UMAP.sct2",
                 color = "orig.ident"))

    plot_scatter(metadata = sct_dim_metadata,
                 out_path = out_path,
                 proj_name = proj_name,
                 log_file = log_file,
                 X = "tSNE.sct1",
                 Y = "tSNE.sct2",
                 color = "orig.ident"))
    
    plot_scatter(metadata = sct_dim_metadata,
                 out_path = out_path,
                 proj_name = proj_name,
                 log_file = log_file,
                 X = "PC.sct1",
                 Y = "PC.sct2",
                 color = "orig.ident"))
    
    
    saveRDS(seurat_obj_sct,
            file = glue("{proj_name}.sct.seurat_obj.rds"))
  }

}

write_message <- function(message_str, log_file) {
  # Small function to write to message and to log file
  message(message_str)
  write(message_str,
        file = log_file,
        append = TRUE)
}

load_sample_counts_matrix = function(sample_names, data_path, log_file) {
  # Reads in count data from 10x from one path or multiple paths
  #
  # Args:
  #   sample_names: Names to set each sample
  #   data_path: Paths to each sample /data_path/outs
  #   log_file: Name of log_file
  #
  # Returns:
  #   Counts matrix from 10x, merged from several samples or from one

  message_str <- "\n\n ========== import cell ranger counts matrix ========== \n\n"
  write_message(message_str, log_file)

  merged_counts_matrix = NULL

  for (i in 1:length(sample_names)) {

    sample_name = sample_names[i]

    message_str <- glue("loading counts matrix for sample: {sample_name}")
    write_message(message_str, log_file)
    

    # check if sample dir is valid
    if (!dir.exists(data_path)) stop(glue("dir {data_path} does not exist"))

    # determine counts matrix directory (HDF5 is not the preferred option)
    # "filtered_gene_bc_matrices" for single library
    # "filtered_gene_bc_matrices_mex" for aggregated
    # Cell Ranger 3.0: "genes" has been replaced by "features" to account for feature barcoding
    # Cell Ranger 3.0: the matrix and barcode files are now gzipped
    data_dir = glue("{data_path}/outs")
    if (!dir.exists(data_dir)) stop(glue("dir {data_path} does not contain outs directory"))
    data_dir = list.files(path = data_dir, pattern = "matrix.mtx", full.names = TRUE, recursive = TRUE)
    data_dir = str_subset(data_dir, "filtered_.*_bc_matri")[1]
    data_dir = dirname(data_dir)
    if (!dir.exists(data_dir)) stop(glue("dir {data_path} does not contain matrix.mtx"))

    message_str <- glue("loading counts matrix dir: {data_dir}")
    write_message(message_str, log_file)
    

    counts_matrix = Read10X(data_dir)

    message_str <- glue("library {sample_name} cells: {ncol(counts_matrix)}
                        library {sample_name} genes: {nrow(counts_matrix)}")
    write_message(message_str, log_file)
  
    # clean up counts matrix to make it more readable
    counts_matrix = counts_matrix[sort(rownames(counts_matrix)), ]
    colnames(counts_matrix) = str_c(sample_name, ":", colnames(counts_matrix))

    # combine current matrix with previous
    if (i == 1) {

      # skip if there is no previous matrix
      merged_counts_matrix = counts_matrix

    } else {

      # check if genes are the same for current and previous matrices
      if (!identical(rownames(merged_counts_matrix), rownames(counts_matrix))) {

        # generate a warning, since this is probably a mistake
        warning("counts matrix genes are not the same for different libraries")
        write("counts matrix genes are not the same for different libraries",
              file = log_file,
              append = TRUE)
        Sys.sleep(1)

        # get common genes
        common_genes = intersect(rownames(merged_counts_matrix), rownames(counts_matrix))
        common_genes = sort(common_genes)

        message_str <- glue("num genes for previous libraries: {length(rownames(merged_counts_matrix))}
                            num genes for current library: {length(rownames(counts_matrix))}
                            num genes in common: {length(common_genes)}")
        write_message(message_str, log_file)
        
        # exit if the number of overlapping genes is too few
        if (length(common_genes) < (length(rownames(counts_matrix)) * 0.9)) stop("libraries have too few genes in common")

        # subset current and previous matrix to overlapping genes
        merged_counts_matrix = merged_counts_matrix[common_genes, ]
        counts_matrix = counts_matrix[common_genes, ]

      }

      # combine current matrix with previous
      merged_counts_matrix = cbind(merged_counts_matrix, counts_matrix)
      Sys.sleep(1)

    }

  }

  return(merged_counts_matrix)

}

create_seurat_obj <- function(counts_matrix, out_path, proj_name, log_file, aggregated = NULL) {
  # convert a sparse matrix of counts to a Seurat object and generate some QC plots
  # 
  # Args:
  #   counts_matrix: Sparse matrix of gene counts 
  #   out_path: Output directory
  #   proj_name: Name of the project -- also names of output files
  #   aggregated: If using the aggregated counts from cellranger
  #
  # Returns:
  #   Seurat_object with light filtering min cells =5 and min features = 250
  

  message_str <- glue("\n\n ========== create seurat object ========== \n\n
                     input cells: {ncol(counts_matrix)}
                     input genes: {nrow(counts_matrix)}")
  write_message(message_str, log_file)
  

  if (is.null(aggregated)) {

    # if name is not set, then it's a manually merged counts matrix
    s_obj = CreateSeuratObject(counts = counts_matrix, min.cells = 5, min.features = 250, project = proj_name,
                               names.field = 1, names.delim = ":")
    rm(counts_matrix)

  } else if (aggregated == "aggregated") {

    # multiple libraries combined using Cell Ranger (cellranger aggr)
    # setup taking into consideration aggregated names delimiter
    s_obj = CreateSeuratObject(counts = counts_matrix, min.cells = 5, min.features = 250, project = proj_name,
                               names.field = 2, names.delim = "-")

    # import cellranger aggr sample sheet
    sample_sheet_csv = paste0(out_path, "/outs/aggregation_csv.csv")
    sample_sheet = read.csv(sample_sheet_csv, stringsAsFactors = FALSE)
    
    samples_names <- paste(sample_sheet[, 1], collapse=", ")
    message_str <- glue("samples: {samples_names}")
    

    # change s_obj@meta.data$orig.ident sample identities from numbers to names
    s_obj[["orig.ident"]][, 1] = factor(sample_sheet[s_obj[["orig.ident"]][, 1], 1])
    # set s_obj@ident to the new s_obj@meta.data$orig.ident
    s_obj = set_identity(seurat_obj = s_obj, group_var = "orig.ident")
    rm(counts_matrix)

  } else {

    stop("aggregated name set to unknown value")

  }

  message_str <- glue("imported cells: {ncol(s_obj)}
                      imported genes: {nrow(s_obj)}")
  write_message(message_str, log_file)
  
  # rename nCount_RNA and nFeature_RNA slots to make them more clear
  s_obj$num_UMIs = s_obj$nCount_RNA
  s_obj$num_genes = s_obj$nFeature_RNA

  # nFeature_RNA and nCount_RNA are automatically calculated for every object by Seurat
  # calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData
  mt_genes = grep("^MT-", rownames(GetAssayData(s_obj)), ignore.case = TRUE, value = TRUE)
  percent_mt = Matrix::colSums(GetAssayData(s_obj)[mt_genes, ]) / Matrix::colSums(GetAssayData(s_obj))
  percent_mt = round(percent_mt * 100, digits = 3)

  # add columns to object@meta.data, and is a great place to stash QC stats
  s_obj = AddMetaData(s_obj,
                      metadata = percent_mt,
                      col.name = "pct_mito")

  # check distribution of gene counts and mitochondrial percentage
  low_quantiles = c(0.05, 0.02, 0.01, 0.001)
  high_quantiles = c(0.95, 0.98, 0.99, 0.999)
  
  low_pct <- s_obj$num_genes %>%
    quantile(low_quantiles) %>%
    round(1) %>% 
    print()
  
  message_str <- glue("num genes low percentiles: {names(low_pct)}:{low_pct}       ")
  write_message(message_str, log_file)

  high_pct <- s_obj$num_genes %>%
    quantile(high_quantiles) %>%
    round(1) %>% 
    print()
  message_str <- glue("num genes high percentiles: {names(high_pct)}:{high_pct}       ")
  write_message(message_str, log_file)
  
  high_mito_pct <- s_obj$pct_mito %>%
    quantile(high_quantiles) %>%
    round(1) %>% 
    print()
  message_str <- glue("num genes high mito percentiles: {names(high_mito_pct)}:{high_mito_pct}       ")
  write_message(message_str, log_file)

  return(s_obj)
}

save_counts_matrix <- function(seurat_obj, out_path, proj_name, type, log_file, assay = "RNA", slot = "data") {
  # save counts matrix as a csv file (to be consistent with the rest of the tables)
  
  s_obj <- seurat_obj
  
  message_str <- "\n\n ========== saving counts ========== \n\n"
  write_message(message_str, log_file)
  
  # save counts matrix as a basic gzipped text file
  # object@data stores normalized and log-transformed single cell expression
  # used for visualizations, such as violin and feature plots, most diff exp tests, finding high-variance genes
  counts = GetAssayData(s_obj, assay = assay) %>%
    as.matrix() %>%
    round(3)
  
  counts = counts %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(gene)
  
  write_csv(counts,
            path = glue("{out_path}/{proj_name}.{type}.counts.csv.gz"))
  
  rm(counts)
}

save_seurat_metadata <- function(seurat_obj, dim_red_list = NULL, out_path, proj_name, type, log_file) {
  # save metadata from seurat object 
  
  message_str <- "\n\n ========== saving metadata ========== \n\n"
  write_message(message_str, log_file)
  
  s_obj <- seurat_obj
  
  if (!is.null(dim_red_list)) {
    if(is.null(dim_red_suffix)){
      message_str <- "Watch out: This dimensionality reduction will not get a unique name"
      write_message(message_str, log_file)
    }
    
    # compile all cell metadata into a single table
    metadata_tbl = s_obj@meta.data %>%
      rownames_to_column("cell") %>% 
      as_tibble() %>%
      mutate(sample_name = orig.ident)
    
    umap <- dim_red_list[[grep("^umap", names(dim_red_list))]] %>% 
      as.data.frame() %>% 
      rownames_to_column("cell") %>% 
      as_tibble() 
    
    tsne <- dim_red_list[[grep("^tsne", names(dim_red_list))]] %>% 
      as.data.frame() %>% 
      rownames_to_column("cell") %>% 
      as_tibble() 
    
    pca <- dim_red_list[[grep("^cell.embeddings", names(dim_red_list))]] %>% 
      as.data.frame() %>% 
      rownames_to_column("cell") %>% 
      as_tibble() 
    
    cells_metadata = metadata_tbl %>%
      full_join(umap ,by = "cell") %>% 
      full_join(tsne, by = "cell") %>% 
      full_join(pca, by = "cell")

    cells_metadata = cells_metadata %>%
      arrange(cell)
  } else {
    cells_metadata = s_obj@meta.data %>%
      rownames_to_column("cell") %>% 
      as_tibble() %>%
      mutate(sample_name = orig.ident) %>% 
      arrange(cell)
  }
  write_excel_csv(cells_metadata, path = glue("{out_path}/{proj_name}.{type}.metadata.csv"))
  return(cells_metadata)
}

create_color_vect <- function(seurat_obj, group = "orig.ident") {
  # create a vector of colors for the Idents of the s_obj
  sample_names <- switch(class(seurat_obj),
                 Seurat = s_obj[[group]] %>% unique() %>% arrange(get(group)),
                 data.frame = unique(seurat_obj))
  
  colors_samples = c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
  # create a named color scheme to ensure names and colors are in the proper order
  sample_names[] <- lapply(sample_names, as.character)
  colors_samples_named = colors_samples[1:nrow(sample_names)]
  names(colors_samples_named) = sample_names[,1]
  return(colors_samples_named)
}

plot_qc_seurat <- function(seurat_obj, out_dir, proj_name, type = "_", group = "orig.ident") {
  # plot qc plots from seurat obj like violin and scatter plots 
  
  s_obj <- seurat_obj
  
  colors_samples_named <-  create_color_vect(s_obj, group = group)
  
  # num genes violin
  num_genes_violin <- ggplot(s_obj@meta.data, aes(x = reorder(eval(as.name(group)), num_genes), y = num_genes, fill = eval(as.name(group)))) +
    geom_violin() +
    xlab(group) +
    ylab("Number of genes per cell") +
    scale_fill_manual(values = colors_samples_named,
                      name = group)
  
  num_umi_violin <- ggplot(s_obj@meta.data, aes(x = reorder(eval(as.name(group)), num_UMIs), y = num_UMIs, fill = eval(as.name(group)))) +
    geom_violin() +
    xlab(group) +
    ylab("Number of UMIs per cell") +
    scale_fill_manual(values = colors_samples_named,
                      name = group)
  
  pct_mito_violin <- ggplot(s_obj@meta.data, aes(x = reorder(eval(as.name(group)), pct_mito), y = pct_mito, fill = eval(as.name(group)))) +
    geom_violin() +
    xlab(group) +
    ylab("Percent Mitochondrial gene expression per cell") +
    scale_fill_manual(values = colors_samples_named,
                      name = group)

  qc_violin_plot = plot_grid(num_genes_violin,
                               num_umi_violin,
                               pct_mito_violin,
                                 ncol = 3)


  ggsave(file = glue("{out_path}/{proj_name}.{type}.qc.png"),
     plot = qc_violin_plot,
     width = 10,
     height = 6,
     units = "in")

  Sys.sleep(1)


  UMI_gene_scatter <- ggplot(s_obj@meta.data, aes(x = num_UMIs, y = num_genes, col = eval(as.name(group)))) +
    geom_point() +
    scale_color_manual(values = colors_samples_named, 
                       name = group) +
    coord_fixed(ratio = max(s_obj@meta.data$num_UMI)/max(s_obj@meta.data$num_genes))
    
  UMI_mito_scatter <- ggplot(s_obj@meta.data, aes(x = num_UMIs, y = pct_mito, col = eval(as.name(group)))) +
    geom_point() +
    scale_color_manual(values = colors_samples_named, 
                       name = group) +
    coord_fixed(ratio = max(s_obj@meta.data$num_UMI)/max(s_obj@meta.data$pct_mito))
  
  genes_mito_scatter <-  ggplot(s_obj@meta.data, aes(x = num_genes, y = pct_mito, col = eval(as.name(group)))) +
    geom_point() +
    scale_color_manual(values = colors_samples_named, 
                       name = group) +
    coord_fixed(ratio = max(s_obj@meta.data$num_genes)/max(s_obj@meta.data$pct_mito))
  
  qc_scatter_plots = plot_grid(UMI_gene_scatter,
                              UMI_mito_scatter,
                              genes_mito_scatter,
                              ncol = 3)


  ggsave(glue("{out_path}/{proj_name}.{type}.qc.correlations.png"),
         plot = qc_scatter_plots,
         width = 18,
         height = 5,
         units = "in")

  Sys.sleep(1)
}

get_dr_point_size <- function(seurat_obj) {
  # get point size for dim red plots
  
  pt_size = 1.8
  if (ncol(seurat_obj) > 1000) pt_size = 1.2
  if (ncol(seurat_obj) > 5000) pt_size = 1.0
  if (ncol(seurat_obj) > 10000) pt_size = 0.8
  if (ncol(seurat_obj) > 25000) pt_size = 0.6

  return(pt_size)

}

subset_singlets <- function(seurat_obj, method = "HTO_demux", log_file) {
  
  message_str <- glue("Keeping only HTO demux Singlets
                      
                      before: {table(seurat_obj$HTO_classification.global)}
                      ")
  write_message(message_str, log_file)
  
  if(method == "HTO_demux") {
    #subset singlets
    Idents(seurat_obj) <- "HTO_classification.global"
    seurat_obj <- subset(seurat_obj, cells = WhichCells(seurat_obj, idents = "Singlet"))
  }
  message_str <- glue("after: {table(seurat_obj$HTO_classification.global)}
                      ")
  write_message(message_str, log_file)
  
  return(seurat_obj)
}

filter_data <- function(seurat_obj, out_dir, proj_name, log_file, min_genes = NULL, max_genes = NULL, max_mt = 10) {
  # filter data by number of genes and mitochondrial percentage
  #
  # Args:
  #   seurat_obj: Seurat object
  #   out_dir: Output directory
  #   proj_name: Name or project and name of output files
  #   min_genes: Minimum number of genes
  #   max_genes: Maximum number of genes
  #   max_mt: Maximum mito pct
  #
  # Results:
  #   Filtered seurat object

  s_obj = seurat_obj

  message_str <- glue("\n\n ========== filter data matrix ========== \n\n
                      unfiltered min genes: {min(s_obj$num_genes)}
                      unfiltered max genes: {max(s_obj$num_genes)}
                      unfiltered mean num genes: {round(mean(s_obj$num_genes), 3)}
                      unfiltered median num genes: {median(s_obj$num_genes)}")
  write_message(message_str, log_file)

  # convert arguments to integers (command line arguments end up as characters)
  min_genes = as.numeric(min_genes)
  max_genes = as.numeric(max_genes)
  max_mt = as.numeric(max_mt)

  # default cutoffs (gene numbers rounded to nearest 10)
  # as.numeric() converts NULLs to 0 length numerics, so can't use is.null()
  if (!length(min_genes)) min_genes = s_obj$num_genes %>%
    quantile(0.02, names = FALSE) %>%
    round(-1)

  if (!length(max_genes)) max_genes = s_obj$num_genes %>%
    quantile(0.98, names = FALSE) %>%
    round(-1)

  if (!length(max_mt)) max_mt = 10

  message_str <- glue("min genes cutoff: {min_genes}
                      max genes cutoff: {max_genes}
                      max mitochondrial percentage cutoff: {max_mt}
                      imported cells: {ncol(s_obj)}
                      imported genes: {nrow(s_obj)}")
  write_message(message_str, log_file)

  # filter
  cells_subset =
    seurat_obj@meta.data %>%
    rownames_to_column("cell") %>%
    filter(nFeature_RNA > min_genes & nFeature_RNA < max_genes & pct_mito < max_mt) %>%
    pull(cell)
  s_obj = subset(s_obj, cells = cells_subset)

  message_str <- glue("filtered cells: {ncol(s_obj)}
                      filtered genes: {nrow(s_obj)}")
  return(s_obj)
}

log_normalize_data <- function(seurat_obj, out_path, proj_name, log_file) {
  # log normalize data
  
  message_str <- "\n\n ========== log normalize ========== \n\n"
  write_message(message_str, log_file)

  s_obj <- seurat_obj
  # after removing unwanted cells from the dataset, normalize the data
  # LogNormalize:
  # - normalizes the gene expression measurements for each cell by the total expression
  # - multiplies this by a scale factor (10,000 by default)
  # - log-transforms the result
  s_obj = NormalizeData(s_obj,
                        normalization.method = "LogNormalize",
                        assay = "RNA",
                        scale.factor = 10000,
                        verbose = FALSE)
  # log to file
  message_str <- glue("filtered cells: {ncol(s_obj)}
                      filtered genes: {nrow(s_obj)}
                      filtered mean num genes: {round(mean(s_obj$num_genes), 3)}
                      filtered median num genes: {median(s_obj$num_genes)}")
  write_message(message_str, log_file)
  
  return(s_obj)
}

calculate_variance <- function(seurat_obj, out_path, proj_name, log_file){
  # calculate variance of genes in a seurat object
  
  s_obj = seurat_obj

  message_str <- "\n\n ========== Seurat::FindVariableGenes() ========== \n\n"
  write_message(message_str, log_file)

  # identify features that are outliers on a 'mean variability plot'
  # Seurat v3 implements an improved method based on a variance stabilizing transformation ("vst")
  s_obj = FindVariableFeatures(s_obj,
                               selection.method = "vst",
                               nfeatures = 2000,
                               verbose = FALSE)

  # export highly variable feature information (mean, variance, variance standardized)
  hvf_tbl = HVFInfo(s_obj) %>%
    round(3) %>%
    rownames_to_column("gene") %>%
    arrange(-variance.standardized)

  write_excel_csv(hvf_tbl,
                  path = glue("{out_path}/{proj_name}.variance.csv"))

  # plot variance
  var_plot = VariableFeaturePlot(s_obj,
                                 pt.size = 0.5)

  var_plot = LabelPoints(var_plot,
                         points = head(hvf_tbl$gene, 30),
                         repel = TRUE,
                         xnudge = 0,
                         ynudge = 0)

  ggsave(glue("{out_path}/{proj_name}.variance.features.png"),
         plot = var_plot,
         width = 12,
         height = 5,
         units = "in")

  message_str <- "\n\n ========== Seurat::ScaleData() ========== \n\n"
  write_message(message_str, log_file)

  # regress out unwanted sources of variation
  # regressing uninteresting sources of variation can improve dimensionality reduction and clustering
  # could include technical noise, batch effects, biological sources of variation (cell cycle stage)
  # scaled z-scored residuals of these models are stored in scale.data slot
  # used for dimensionality reduction and clustering
  # RegressOut function has been deprecated, and replaced with the vars.to.regress argument in ScaleData
  # s_obj = ScaleData(s_obj, features = rownames(s_obj), vars.to.regress = c("num_UMIs", "pct_mito"), verbose = FALSE)
  s_obj = ScaleData(s_obj,
                    vars.to.regress = c("num_UMIs", "pct_mito"),
                    verbose = FALSE)

  return(s_obj)
}

sctransform_data <- function(seurat_obj, out_path, log_file, proj_name){
  # sc transform data
  
  s_obj <-seurat_obj

  message("\n\n ========== sc transform ========== \n\n")

  s_obj <- PercentageFeatureSet(s_obj,
                                pattern = "^MT-",
                                col.name = "percent.feature.set.mt")
  # run sctransform
  s_obj <- SCTransform(s_obj,
                       vars.to.regress = c("percent.feature.set.mt", "num_UMIs"),
                       verbose = FALSE)

  # save counts matrix as a basic gzipped text file
  # object@data stores normalized and log-transformed single cell expression
  # used for visualizations, such as violin and feature plots, most diff exp tests, finding high-variance genes
  counts_norm = GetAssayData(s_obj) %>%
    as.matrix() %>%
    round(3)
  counts_norm = counts_norm %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(gene)

  write_csv(counts_norm, path = glue("{out_path}/{proj_name}.counts.normalized_sc.csv.gz"))

  # log to file
  # log to file
  write(glue("filtered cells: {ncol(s_obj)}"), file = log_file, append = TRUE)
  write(glue("filtered genes: {nrow(s_obj)}"), file = log_file, append = TRUE)
  write(glue("filtered mean num genes: {round(mean(s_obj$num_genes), 3)}"), file = log_file, append = TRUE)
  write(glue("filtered median num genes: {median(s_obj$num_genes)}"), file = log_file, append = TRUE)

  return(s_obj)

}

run_dimensionality_reduction <- function(seurat_obj,  num_pcs, num_dim, log_file, dim_red_suffix = NULL){
  # taken from seurat 
  # I didnt 100% understand the Seurat code so I took the parts I needed and made it so that
  # I could calculate dim red on a matrix as well 
  
  message_str <- "\n\n ========== dimensionality reduction ========== \n\n"
  write_message(message_str, log_file)
  
  data <- switch(class(seurat_obj),
                 Seurat = Seurat:::PrepDR(seurat_obj, features = NULL),
                 matrix = seurat_obj)
  # if(!is.null(seurat_obj)){
  #   data <- Seurat:::PrepDR(object = seurat_obj,
  #          features = NULL,
  #          verbose = TRUE)
  # } else if(!is.null(matrix)){
  #   data <- matrix
  # } else {
  #   stop("ERROR: either input a seurat object or a matrix")
  # }
  
  # suffix is _ if not specified
  dim_red_suffix <- ifelse(is.null(dim_red_suffix), "_", dim_red_suffix)
  
  pca_out <- run_pca(data = data,
                     num_dim = num_pcs,
                     reduction.key = paste("PC",
                                           dim_red_suffix,
                                           sep = ""))
  
  feature.loadings <- pca_out[[1]]
  cell.embeddings <- pca_out[[2]]
  sdev = pca_out[[3]]
  
  
  tsne_out <- run_tsne(cell.embeddings[,1:num_dim],
                       reduction.key = paste("tSNE", 
                                             dim_red_suffix,
                                             sep = ""))
  
  umap_out <- run_umap(cell.embeddings[,1:num_dim],
                       reduction.key = paste("UMAP",
                                             dim_red_suffix,
                                             sep = ""))
  
  return(list(feature.loadings = feature.loadings,
              cell.embeddings = cell.embeddings, 
              sdev = sdev, 
              tsne_out = tsne_out,
              umap_out = umap_out))
}

run_pca <- function(data, num_dim, reduction.key = "PC_"){
  # copied from the Seurat package
  
  npcs <- min(num_dim, nrow(data))
  pca.results <- prcomp(t(data), rank. = npcs)
  feature.loadings <- pca.results$rotation
  sdev <- pca.results$sdev
  total.variance <- sum(sdev)
  cell.embeddings <- pca.results$x
  
  
  rownames(x = feature.loadings) <- rownames(data)
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:npcs)
  rownames(x = cell.embeddings) <- colnames(data)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  
  return(list(feature.loadings, cell.embeddings, sdev))
}

run_tsne <- function(data, seed.use = 1, tsne.method ='Rtsne',  dim.embed = 2, reduction.key = "tSNE_"){
  #copied from the Seurat package
    set.seed(seed = seed.use)
    tsne.data <- switch(
      EXPR = tsne.method,
      'Rtsne' = Rtsne(data, dims = dim.embed)$Y,
      'FIt-SNE' = fftRtsne(data, dims = dim.embed, rand_seed = seed.use),
      stop("Invalid tSNE method: please choose from 'Rtsne' or 'FIt-SNE'")
    )

    colnames(x = tsne.data) <- paste0(reduction.key, 1:ncol(x = tsne.data))
    rownames(x = tsne.data) <- rownames(x = data)

    return(tsne.data)
}

run_umap <- function(object, seed.use = 42L, n.neighbors = 30L, n.components = 2L, metric = "correlation",
                     n.epochs = NULL, learning.rate = 1, min.dist = 0.3, spread = 1, set.op.mix.ratio = 1,
                     local.connectivity = 1L, repulsion.strength = 1, negative.sample.rate = 5L, 
                     a = NULL, b = NULL, metric.kwds = NULL, verbose = TRUE,angular.rp.forest = FALSE, reduction.key = 'UMAP_'){
  #copied from the seurat package
  if (!py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
  }
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
    py_set_seed(seed = seed.use)
  }
  if (typeof(x = n.epochs) == "double") {
    n.epochs <- as.integer(x = n.epochs)
  }
  umap_import <- import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(x = n.neighbors),
    n_components = as.integer(x = n.components),
    metric = metric,
    n_epochs = n.epochs,
    learning_rate = learning.rate,
    min_dist = min.dist,
    spread = spread,
    set_op_mix_ratio = set.op.mix.ratio,
    local_connectivity = local.connectivity,
    repulsion_strength = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    a = a,
    b = b,
    metric_kwds = metric.kwds,
    angular_rp_forest = angular.rp.forest,
    verbose = verbose
  )
  umap_output <- umap$fit_transform(as.matrix(x = object))
  colnames(x = umap_output) <- paste0(reduction.key, 1:ncol(x = umap_output))
  rownames(x = umap_output) <- rownames(object)
  
  return(umap_output)
}

add_dim_red_seurat <- function(seurat_obj, dim_red_list, dim_red_suffix = NULL){
  # from seurat - removed bloat and warnings and changing my key for annoying reasons
  
  ## Turns out you cant add the PC names the way I wanted to anyways so 
  
  pca.dim.reduc <- new(Class = "DimReduc", 
                   cell.embeddings =  as.matrix(dim_red_list$cell.embeddings), 
                   feature.loadings = as.matrix(dim_red_list$feature.loadings), 
                   assay.used = "RNA", 
                   stdev = dim_red_list$sdev, 
                   key = paste0("PC", dim_red_suffix))
  
  tsne.dim.reduc <- new(Class = "DimReduc", 
                       cell.embeddings =  as.matrix(dim_red_list$tsne_out), 
                       assay.used = "RNA", 
                       key = paste0("tSNE", dim_red_suffix))
  
  umap.dim.reduc <- new(Class = "DimReduc", 
                        cell.embeddings =  as.matrix(dim_red_list$umap_out), 
                        assay.used = "RNA", 
                        key = paste0("UMAP", dim_red_suffix))

  
  dim_red_suffix <- ifelse(is.null(dim_red_suffix), "", dim_red_suffix)
  
  seurat_obj[[paste("pca", dim_red_suffix, sep = "")]] <- pca.dim.reduc
  seurat_obj[[paste("tsne", dim_red_suffix, sep = "")]] <- tsne.dim.reduc
  seurat_obj[[paste("umap", dim_red_suffix, sep = "")]] <- umap.dim.reduc
  
  return(seurat_obj)
}

plot_scatter<- function(metadata, out_path, proj_name, log_file, X, Y, color){
  
  colors_samples_named <- create_color_vect(as.data.frame(metadata[color]))
  
  current_plot <- ggplot(metadata, aes(x = eval(as.name(X)), y = eval(as.name(Y)), color = eval(as.name(color)))) +
    geom_point() +
    coord_fixed(ratio = max(metadata[X])/max(metadata[Y])) +
    xlab(X) + 
    ylab(Y) +
    scale_color_manual(values = colors_samples_named, 
                       name = color)

  
  ggsave(glue("{out_path}/{proj_name}.{X}.{Y}.{color}.png"),
         plot = current_plot,
         width = 8,
         height = 6,
         units = "in")
}


read_remove.unmapped <- function(HTO_file) {
  # remove unmapped row
  
  hto_matrix <- read.table(HTO_file)
  hto_matrix <- hto_matrix[-which(rownames(hto_matrix) == "unmapped"),]
  return(hto_matrix)
  
}

create_seurat_obj_hto <- function(seurat_obj, out_dir, HTO_counts, proj_name, log_file, 
                                  hto_demux_quantile = 0.99, multi_demux_quantile = 0.7) {
  # add HTO slot to an existing seurat object
  #
  # Args:
  #   seurat_obj: Seurat object
  #   out_dir: Output directory
  #   HTO_counts: HTO counts as a dataframe
  #   proj_name: Name of project and name of output files
  #
  # Results:
  #   A seurat object with HTO in HTO slot using HTO demux and MULTIseqDemux
  
  s_obj <- seurat_obj

  message_str <- "\n\n ========== creating HTO slot ========== \n\n"
  write_message(message_str, log_file)

  colnames(HTO_counts) <- str_c(proj_name, ":", colnames(HTO_counts))

  # check if the cells in the data are the same as the cells in the hashtag data
  cells_to_use <- intersect(colnames(seurat_obj), colnames(HTO_counts))

  if(length(s_obj) != length(cells_to_use)){
    message_str <- "some cells in scrna matrix not in hto matrix"
    write_message(message_str, log_file)
  }
  if(ncol(HTO_counts) != length(cells_to_use)){
    message_str <- "some cells in hto matrix not in scrna matrix"
    write_message(message_str, log_file)
  }

  # Subset RNA and HTO counts by joint cell barcodes
  HTO_counts <- as.matrix(HTO_counts[, cells_to_use])
  s_obj <- subset(s_obj, cells = cells_to_use)

  # add hashtag slot
  s_obj[["HTO"]] <- CreateAssayObject(counts = HTO_counts)

  # Normalize HTO data
  s_obj <- NormalizeData(s_obj, assay = "HTO", normalization.method = "CLR")
  # demultiplex htos
  s_obj <- HTODemux(s_obj, assay = "HTO", positive.quantile = hto_demux_quantile)
  s_obj <- MULTIseqDemux(s_obj, assay = "HTO", quantile = multi_demux_quantile)

  message_str <- paste0(capture.output(table(s_obj$HTO_classification.global)), collapse = "\n")
  write_message(message_str, log_file)
  
  return(s_obj)
}

seurat_plot_hto <- function(seurat_obj, out_path, proj_name) {
  # plot standard seurat plots for HTO analysis
  
  Idents(seurat_obj) <- "HTO_maxID"

  seurat_obj_ridge_plot <- RidgePlot(seurat_obj,
                                     assay = "HTO",
                                     features = rownames(seurat_obj[["HTO"]]),
                                     ncol = 2)
  ggsave(seurat_obj_ridge_plot,
         filename = glue("{out_path}/{proj_name}.htodemux_ridge_plot.png"),
         height = 7,
         width =7 ,
         units = "in")

  # heatmap
  seurat_obj_heatmap <- HTOHeatmap(seurat_obj,
                                   assay = "HTO",
                                   ncells = 800)
  ggsave(seurat_obj_heatmap, filename =  glue("{out_path}/{proj_name}.htodemux_heatmap.png"),
         height = 7,
         width =7 ,
         units = "in")

  # scatter plots
  HTO_seurat <- as.data.frame(t(as.data.frame(seurat_obj@assays$HTO@data)))
  hto_pairs <- ggpairs(HTO_seurat)
  ggsave(hto_pairs, filename =  glue("{out_path}/{proj_name}.htodemux_pairs.png"),
         height = 7,
         width =7 ,
         units = "in")

  # plot the doublet trends
  doublet_trend <- ggplot(as.data.frame(seurat_obj$HTO_classification), aes(seurat_obj$HTO_classification)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(doublet_trend, filename = glue("{out_path}/{proj_name}.htodemux_doubet.png"),
         height = 7,
         width =7 ,
         units = "in")

  doub_col <- c("red", "blue", "green")
  names(doub_col) <- c("Singlet", "Negative", "Doublet")
  doublet_bar <- data.frame(doub = seurat_obj$HTO_classification.global)
  doublet_bar <- doublet_bar %>%
    group_by(doub) %>%
    summarise (n = n()) %>%
    unique() %>%
    mutate(percentage = n /sum(n))
  doublet_bar <- doublet_bar %>%
    mutate(sample = rep(proj_name, nrow(doublet_bar)))
  plot_num_doublet <- ggplot(doublet_bar, aes(x = sample, y =percentage, fill = doub))+
    geom_col() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = doub_col) +
    theme_bw()
  ggsave(plot_num_doublet, filename =  glue("{out_path}/{proj_name}.htodemux_doublet_count.png"),
         height = 7,
         width =7 ,
         units = "in")


  # ridge plots
  Idents(seurat_obj) <- "MULTI_ID"
  seurat_obj_ridge_plot <- RidgePlot(seurat_obj,
                                     assay = "HTO",
                                     features = rownames(seurat_obj[["HTO"]]),
                                     ncol = 2)
  ggsave(seurat_obj_ridge_plot, filename = glue("{out_path}/{proj_name}.htomulti_ridge_plot.png"),
         height = 7,
         width =7 ,
         units = "in")

  # heatmap
  seurat_obj_heatmap <- HTOHeatmap(seurat_obj, assay = "HTO", ncells = 800)
  ggsave(seurat_obj_heatmap, filename =  glue("{out_path}/{proj_name}.htomulti_heatmap.png"),
                                              height = 7,
                                              width =7 ,
                                              units = "in")

  # scatter plots
  HTO_seurat <- as.data.frame(t(as.data.frame(seurat_obj@assays$HTO@data)))
  hto_pairs <- ggpairs(HTO_seurat)
  ggsave(hto_pairs, filename =  glue("{out_path}/{proj_name}.htomulti_pairs.png"),
         height = 7,
         width =7 ,
         units = "in")

  # plot the doublet trends
  doublet_trend <- ggplot(as.data.frame(seurat_obj$MULTI_classification), aes(seurat_obj$MULTI_classification)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(doublet_trend, filename = glue("{out_path}/{proj_name}.htomulti_doubet.png"),
         height = 7,
         width =7 ,
         units = "in")
}

create_seurat_obj_adt <- function(seurat_obj, out_dir, ADT_counts, proj_name, log_file) {
  # add ADT slot to an existing seurat object
  #
  # Args:
  #   seurat_obj: Seurat object
  #   out_dir: Output directory
  #   ADT_counts: ADT counts as a dataframe
  #   proj_name: Name of project and name of output files
  #
  # Results:
  #   A seurat object with ADT in ADT slot 
  
  s_obj <- seurat_obj
  
  message_str <- "\n\n ========== creating ADT slot ========== \n\n"
  write_message(message_str, log_file)
  
  colnames(ADT_counts) <- str_c(proj_name, ":", colnames(ADT_counts))
  
  # check if the cells in the data are the same as the cells in the hashtag data
  cells_to_use <- intersect(colnames(seurat_obj), colnames(ADT_counts))
  
  if(length(s_obj) != length(cells_to_use)){
    message_str <- "some cells in scrna matrix not in ADT matrix"
    write_message(message_str, log_file)
  }
  if(ncol(ADT_counts) != length(cells_to_use)){
    message_str <- "some cells in ADT matrix not in scrna matrix"
    write_message(message_str, log_file)
  }
  
  # Subset RNA and HTO counts by joint cell barcodes
  ADT_counts <- as.matrix(ADT_counts[, cells_to_use])
  s_obj <- subset(s_obj, cells = cells_to_use)
  
  # add ADT slot
  s_obj[["ADT"]] <- CreateAssayObject(counts = ADT_counts)
  
  # Normalize HTO data
  s_obj <- NormalizeData(s_obj, assay = "ADT", normalization.method = "CLR")
  s_obj <- ScaleData(s_obj, assay = "ADT")
  
  
  return(s_obj)
}

# try clr outside of Seurat bc I can't figure out what Seurat actually does
manual_hto <- function(HTO_counts, out_path, proj_name){

  hto_clr <- as.data.frame(t(as.data.frame(compositions::clr(HTO_counts))))

  hto_pairs <- ggpairs(hto_clr)
  ggsave(hto_pairs,
         filename =  glue("{out_dir}/{proj_name}.htomulti_manual_pairs.png"),
         height = 7,
         width =7 ,
         units = "in")
}



