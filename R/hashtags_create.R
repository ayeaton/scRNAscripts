
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
  library(cowplot)
  library(scales)
  library(pheatmap)
  library(RColorBrewer)
  library(ggsci)
  library(eulerr)
  library(UpSetR)
  library(GGally)
  library(compositions)
})


assemble_seurat_obj_hto <- function(data_path, # path to 10x data /data_path/outs/
                                    sample_names, # name of samples, can be more than 1
                                    out_path, # path to deposit outputs
                                    num_dim, # number of dimensions to use for PCA, UMAP, and TSNE
                                    HTO_file, # path to HTO file
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

  # read in count data from 10x
  counts_mat <- load_sample_counts_matrix(sample_names = sample_names,
                                          data_path = data_path,
                                          log_file = log_file)

  # take the counts matrix and make a seurat object
  seurat_obj_1 <- create_seurat_obj(counts_mat = counts_mat, 
                                    out_path = out_path, 
                                    proj_name = proj_name, 
                                    log_file = log_file)
  
  # save counts mat
  save_counts_matrix(seurat_obj = seurat_obj_1,
                     out_path = out_path,
                     proj_name = proj_name,
                     log_file = log_file,
                     type = "raw")
  
  # save unfiltered seurat metadata
  save_seurat_metadata(seurat_obj = seurat_obj_1, 
                       out_path = out_path,
                       proj_name = proj_name, 
                       log_file = log_file,
                       type = "unfiltered")
    
  # plot qc of the unfiltered seurat object
  plot_qc_seurat(seurat_obj = seurat_obj_1,
                 out_dir = out_dir,
                 proj_name = proj_name,
                 type = "unfiltered")

  # filter seurat object for min genes, max genes and max mito pct
  seurat_obj_2 <- filter_data(seurat_obj_1, 
                              out_dir = out_dir, 
                              proj_name = proj_name, 
                              log_file = log_file,
                              min_genes = min_genes, 
                              max_genes = max_genes, 
                             max_mt = max_mt)
  rm(seurat_obj_1)
  
  # plot qc plots for filtered seurat obj
  plot_qc_seurat(seurat_obj = seurat_obj_2,
                 out_dir = out_dir,
                 proj_name = proj_name,
                 type = "filtered")

  # add hto data
  hto_data <- clean_hto(HTO_file)

  # here the proj name and the sample name have to match
  seurat_obj_hto <- create_seurat_obj_hto(seurat_obj = seurat_obj_2, 
                                          HTO_counts = hto_data, 
                                      proj_name = proj_name, 
                                      out_dir = data_path,
                                      log_file = log_file, 
                                      hto_demux_quantile = 0.99,
                                      multi_demux_quantile = 0.7)
  rm(seurat_obj_2)

  # save metadata with HTO
  save_seurat_metadata(seurat_obj = seurat_obj_hto,
                       out_path = out_path,
                       log_file = log_file,
                       proj_name = proj_name, 
                       type = "HTO")
  
  # plot qc plots for filtered seurat obj with HTO
  plot_qc_seurat(seurat_obj = seurat_obj_hto,
                 out_dir = out_dir,
                 proj_name = proj_name,
                 type = "filtered.HTO",
                 group = "hash.ID")

  # plot HTO related plots
  seurat_plot_hto(seurat_obj = seurat_obj_hto, 
                  out_path = out_path,
                  proj_name = proj_name)
  
  # manual_hto(seurat_obj, out_path, proj_name)

  # log normalize data
  seurat_obj_log <- log_normalize_data(seurat_obj = seurat_obj_hto, 
                                       out_path = out_path,
                                       proj_name = proj_name,
                                       log_file = log_file)
  rm(seurat_obj_hto)
  
  # calculate variance and plot 
  seurat_obj_var <- calculate_variance(seurat_obj = seurat_obj_log,
                                   out_path = out_path,
                                   proj_name = proj_name,
                                   log_file = log_file)
  rm(seurat_obj_log)
  
  # run PCA, TSNE and UMAP
  seurat_obj_dim <- run_dimensionality_reduction(seurat_obj = seurat_obj_var, 
                                             assay = "RNA", 
                                             num_dim = num_dim,
                                             log_file = log_file)
  rm(seurat_obj_var)
  
  save_seurat_metadata(seurat_obj = seurat_obj_dim,
                       out_path = out_path,
                       proj_name = proj_name, 
                       log_file = log_file,
                       type = "dim")
  
  # plot PCA, UMAP, TSNE 
  plot_dimensionality_reduction(seurat_obj = seurat_obj_dim, 
                                out_path = out_path, 
                                proj_name = proj_name, 
                                log_file = log_file,
                                assay = "RNA",
                                num_pcs = num_dim)
  

  if(sct == T){
    # sctransform data ( should save the sctransform in a new data slot)
    seurat_obj <- sctransform_data(seurat_obj)
    seurat_obj <- run_dimensionality_reduction(seurat_obj, assay = "SCT", num_dim)
    saveRDS(seurat_obj,
            file = glue("{proj_name}.seurat_obj.rds"))
    plot_dimensionality_reduction(seurat_obj, out_path, proj_name, assay = "SCT", num_pcs = num_dim)
  }
  saveRDS(seurat_obj,
          file = glue("{proj_name}.seurat_obj.rds"))
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

save_counts_matrix <- function(seurat_obj, out_path, proj_name, type, log_file) {
  # save counts matrix as a csv file (to be consistent with the rest of the tables)
  
  s_obj <- seurat_obj
  
  message_str <- "\n\n ========== saving counts ========== \n\n"
  write_message(message_str, log_file)
  
  # save counts matrix as a basic gzipped text file
  # object@data stores normalized and log-transformed single cell expression
  # used for visualizations, such as violin and feature plots, most diff exp tests, finding high-variance genes
  counts = GetAssayData(s_obj) %>%
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

save_seurat_metadata <- function(seurat_obj, out_path, proj_name, type, log_file) {
  # save metadata from seurat object 
  
  message_str <- "\n\n ========== saving metadata ========== \n\n"
  write_message(message_str, log_file)
  
  s_obj <- seurat_obj
  
  if (length(which(names(s_obj@reductions) %in% "tsne")) > 0  & length(which(names(s_obj@reductions) %in% "umap")) > 0) {
    # compile all cell metadata into a single table
    metadata_tbl = s_obj@meta.data %>%
      rownames_to_column("cell") %>% 
      as_tibble() %>%
      mutate(sample_name = orig.ident)
    
    tsne_tbl = s_obj[["tsne"]]@cell.embeddings %>%
      round(3) %>%
      as.data.frame() %>%
      rownames_to_column("cell")
    
    umap_tbl = s_obj[["umap"]]@cell.embeddings %>%
      round(3) %>%
      as.data.frame() %>%
      rownames_to_column("cell")
    
    cells_metadata = metadata_tbl %>%
      full_join(tsne_tbl, by = "cell") %>%
      full_join(umap_tbl, by = "cell")
    
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
}

create_color_vect <- function(seurat_obj, group = "orig.ident") {
  # create a vector of colors for the Idents of the s_obj
  s_obj <- seurat_obj
  colors_samples = c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
  # create a named color scheme to ensure names and colors are in the proper order
  sample_names = s_obj[[group]] %>% unique() %>% arrange(get(group))
  sample_names[] <- lapply(sample_names, as.character)
  colors_samples_named = colors_samples[1:nrow(sample_names)]
  names(colors_samples_named) = sample_names[,1]
  return(colors_samples_named)
}

plot_qc_seurat <- function(seurat_obj, out_dir, proj_name, type = "_", group = "orig.ident") {
  # plot qc plots from seurat obj like violin and scatter plots 
  
  s_obj <- seurat_obj
  
  colors_samples_named <-  create_color_vect(s_obj, group = group)

  vln_theme =
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90,
                                 vjust = 0.5,
                                 hjust = 1),
      legend.position = "none"
    )

  suppressMessages({

    # num genes violin
    dist_unfilt_nft_plot =
      VlnPlot(
        s_obj,
        features = "num_genes",
        group.by = group,
        pt.size = 0.1,
        sort = TRUE,
        combine = TRUE,
        cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme

    #num umi violin
    dist_unfilt_nct_plot =
      VlnPlot(
        s_obj,
        features = "num_UMIs",
        group.by = group,
        pt.size = 0.1,
        sort = TRUE,
        combine = TRUE,
        cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme

    # percent mito violin
    dist_unfilt_pmt_plot =
      VlnPlot(
        s_obj,
        features = "pct_mito",
        group.by = group,
        pt.size = 0.1,
        sort = TRUE,
        combine = TRUE,
        cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme


    dist_unfilt_plot = plot_grid(dist_unfilt_nft_plot,
                                 dist_unfilt_nct_plot,
                                 dist_unfilt_pmt_plot,
                                 ncol = 3)


    ggsave(file = glue("{out_path}/{proj_name}.{type}.qc.png"),
           plot = dist_unfilt_plot,
           width = 10,
           height = 6,
           units = "in")
  })

  Sys.sleep(1)


  cor_ncr_nfr_plot =
    FeatureScatter(
      s_obj,
      feature1 = "num_UMIs",
      feature2 = "num_genes",
      group.by = group,
      cols = colors_samples_named
    ) +
    theme(aspect.ratio = 1)

  cor_ncr_pmt_plot =
    FeatureScatter(
      s_obj,
      feature1 = "num_UMIs",
      feature2 = "pct_mito",
      group.by = group,
      cols = colors_samples_named
    ) +
    theme(aspect.ratio = 1)


  cor_nfr_pmt_plot =
    FeatureScatter(
      s_obj,
      feature1 = "num_genes",
      feature2 = "pct_mito",
      group.by = group,
      cols = colors_samples_named
    ) +
    theme(aspect.ratio = 1)


  cor_unfilt_plot = plot_grid(cor_ncr_nfr_plot,
                              cor_ncr_pmt_plot,
                              cor_nfr_pmt_plot,
                              ncol = 3)


  ggsave(glue("{out_path}/{proj_name}.{type}.qc.correlations.png"),
         plot = cor_unfilt_plot,
         width = 18,
         height = 5,
         units = "in")

  Sys.sleep(1)
}

get_dr_point_size = function(seurat_obj) {
  # get point size for dim red plots
  
  pt_size = 1.8
  if (ncol(seurat_obj) > 1000) pt_size = 1.2
  if (ncol(seurat_obj) > 5000) pt_size = 1.0
  if (ncol(seurat_obj) > 10000) pt_size = 0.8
  if (ncol(seurat_obj) > 25000) pt_size = 0.6

  return(pt_size)

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

sctransform_data <- function(seurat_obj){
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
  write(glue("filtered cells: {ncol(s_obj)}"), file = "create.log", append = TRUE)
  write(glue("filtered genes: {nrow(s_obj)}"), file = "create.log", append = TRUE)
  write(glue("filtered mean num genes: {round(mean(s_obj$num_genes), 3)}"), file = "create.log", append = TRUE)
  write(glue("filtered median num genes: {median(s_obj$num_genes)}"), file = "create.log", append = TRUE)

  return(s_obj)

}

run_dimensionality_reduction <- function(seurat_obj, assay, num_dim, log_file) {
  # Runs PCA, UMAP, and TSNE - UMAP AND TSNE use all PCs
  
  s_obj <- seurat_obj
  
  message_str <- "\n\n ========== dimensionality reduction ========== \n\n"
  write_message(message_str, log_file)
  
  if (ncol(s_obj) < 100) num_dim = 20
  if (ncol(s_obj) < 25) num_dim = 5

  # PCA on the scaled data
  # PCA calculation stored in object[["pca"]]
  s_obj <- RunPCA(s_obj, 
                  assay = assay, 
                  features = VariableFeatures(s_obj),
                  npcs = num_dim, 
                  verbose = FALSE)

  s_obj = RunTSNE(s_obj, 
                  reduction = "pca",
                  assay = assay, 
                  dims.use = 1:num_dim)

  # runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
  s_obj = RunUMAP(s_obj,
                  reduction = "pca", 
                  assay = assay, 
                  dims = 1:num_dim, 
                  verbose = FALSE)

  return(s_obj)
}

plot_dimensionality_reduction <- function(seurat_obj, out_path, proj_name, assay, log_file, num_pcs = 30){

  s_obj <- seurat_obj

  Idents(s_obj) <- "orig.ident"

  colors_samples_named <- create_color_vect(seurat_obj)

  # reduce point size for larger datasets
  dr_pt_size = get_dr_point_size(s_obj)

  message_str <- "\n\n ========== plotting dimensionality reduction ========== \n\n"

  # plot the output of PCA analysis (shuffle cells so any one group does not appear overrepresented due to ordering)
  pca_plot =
    DimPlot(
      s_obj,
      cells = sample(colnames(s_obj)),
      group.by = "orig.ident",
      reduction = "pca",
      pt.size = 0.5,
      cols = colors_samples_named
    ) +
    theme(aspect.ratio = 1)
  ggsave(glue("{out_path}/{proj_name}.{assay}.variance.pca.png"),
         plot = pca_plot,
         width = 8,
         height = 6,
         units = "in")

  message_str <- "\n\n ========== Seurat::DimHeatmap() ========== \n\n"
  write_message(message_str, log_file)

  # PCHeatmap (former) allows for easy exploration of the primary sources of heterogeneity in a dataset
  png(glue("{out_path}/{proj_name}.{assay}.variance.pca.heatmap.png"),
      res = 300,
      width = 10,
      height = 16,
      units = "in")
  DimHeatmap(s_obj,
             reduction = "pca",
             dims = 1:15,
             nfeatures = 20,
             cells = 250,
             fast = TRUE)
  dev.off()

  message_str <- "\n\n ========== Seurat::PCElbowPlot() ========== \n\n"
  write_message(message_str, log_file)

  # a more ad hoc method for determining PCs to use, draw cutoff where there is a clear elbow in the graph
  elbow_plot = ElbowPlot(s_obj,
                         reduction = "pca",
                         ndims = num_pcs)

  ggsave(glue("{out_path}/{proj_name}.{assay}.variance.pca.elbow.png"),
         plot = elbow_plot,
         width = 8,
         height = 5,
         units = "in")

  # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
  plot_tsne =
    DimPlot(s_obj, reduction = "tsne",
            cells = sample(colnames(s_obj)),
            pt.size = dr_pt_size,
            cols = colors_samples_named) +
    theme(aspect.ratio = 1)

  ggsave(glue("{out_path}/{proj_name}.{assay}.tsne.{num_dim}.sample.png"),
         plot = plot_tsne,
         width = 10,
         height = 6,
         units = "in")
  Sys.sleep(1)

  ggsave(glue("{out_path}/{proj_name}.{assay}.tsne.{num_dim}.sample.pdf"),
         plot = plot_tsne,
         width = 10,
         height = 6,
         units = "in")
  Sys.sleep(1)
  
  features_plot <- names(s_obj[[]][which(names(s_obj[[]]) %in% c("num_UMIs", 
                                                                  "num_genes", 
                                                                  "pct_mito",
                                                                  "hash.ID", 
                                                                  "HTO_classification.global"))])
  
  feature_tsne =
    FeaturePlot(s_obj, reduction = "tsne",
                cells = sample(colnames(s_obj)),
                pt.size = dr_pt_size,
                features = features_plot) +
    theme(aspect.ratio = 1)
  
  ggsave(glue("{out_path}/{proj_name}.{assay}.tsne.{num_dim}.features.png"),
         plot = feature_tsne,
         width = 10,
         height = 6,
         units = "in")
  Sys.sleep(1)
  
  ggsave(glue("{out_path}/{proj_name}.{assay}.tsne.{num_dim}.features.pdf"),
         plot = feature_tsne,
         width = 10,
         height = 6,
         units = "in")
  Sys.sleep(1)
  
  # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
  plot_umap =
    DimPlot(s_obj, reduction = "umap",
            cells = sample(colnames(s_obj)),
            pt.size = dr_pt_size,
            cols = colors_samples_named) +
    theme(aspect.ratio = 1)

  ggsave(glue("{out_path}/{proj_name}.{assay}.umap.{num_dim}.sample.png"),
         plot = plot_umap,
         width = 10,
         height = 6,
         units = "in")
  Sys.sleep(1)

  ggsave(glue("{out_path}/{proj_name}.{assay}.umap.{num_dim}.sample.pdf"),
         plot = plot_umap,
         width = 10,
         height = 6,
         units = "in")
  Sys.sleep(1)

  return(s_obj)
}

clean_hto <- function(HTO_file) {
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
  s_obj <- SubsetData(s_obj, cells = cells_to_use)

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



