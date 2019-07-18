#combine_raw
source("/beegfs/ay1392/TET2/scRNAseq/Analysis/scRNAscripts/R/hashtags_create.R")


log_file <- "combine_raw.log"
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

out_path <- "/beegfs/ay1392/CITESeq/third_set/combine_raw_CDNA"
num_dim <- 30



proj_name <- "cdna3"
sample_names <- "cdna3"
data_path <- "../data/genome.med.nyu.edu/results/aifantislab/2019-runs/2019-07-01/cellranger/count-cDNA3"
HTO_file <- "../data/genome.med.nyu.edu/results/aifantislab/2019-runs/2019-07-01/cellranger/HTO3_Results/HTO_results_clean.tsv"
ADT_file <- "../data/genome.med.nyu.edu/results/aifantislab/2019-runs/2019-07-01/cellranger/ADT3_Results/ADT_results.tsv"
min_genes <- NULL
max_genes <- NULL
max_mt <- 20


cdna3 <- raw_seurat(data_path, sample_names, proj_name, HTO_file, ADT_file, log_file)
saveRDS(cdna3, file = "cdna3.seurat_obj_raw.rds")

combine_raw <- merge(cdna1, c(cdna2, cdna3))
proj_name <- "combined_raw"

plot_qc_seurat(seurat_obj = combine_raw,
               out_dir = out_dir,
               proj_name = proj_name,
               type = "unfiltered")

# plot HTO related plots
seurat_plot_hto(seurat_obj = combine_raw, 
                out_path = out_path,
                proj_name = proj_name)


# filter seurat object for min genes, max genes and max mito pct
seurat_obj_2 <- filter_data(combine_raw, 
                            out_dir = out_dir, 
                            proj_name = proj_name, 
                            log_file = log_file,
                            min_genes = min_genes, 
                            max_genes = max_genes, 
                            max_mt = max_mt)

# subset singlets
seurat_obj_hto_s <- subset_singlets(seurat_obj = seurat_obj_2,
                                    method = "HTO_demux", 
                                    log_file = log_file)

# plot qc plots for filtered seurat obj
plot_qc_seurat(seurat_obj = seurat_obj_hto_s,
               out_dir = out_dir,
               proj_name = proj_name,
               type = "filtered")

# log normalize data
seurat_obj_log <- log_normalize_data(seurat_obj = seurat_obj_hto_s, 
                                     out_path = out_path,
                                     proj_name = proj_name,
                                     log_file = log_file)


# calculate variance and plot 
seurat_obj_var <- calculate_variance(seurat_obj = seurat_obj_log,
                                     out_path = out_path,
                                     proj_name = proj_name,
                                     log_file = log_file)

seurat_obj_dim <- run_dimensionality_reduction(seurat_obj = seurat_obj_var, 
                                               assay = "RNA", 
                                               num_dim = num_dim,
                                               log_file = log_file)


# plot PCA, UMAP, TSNE 
plot_dimensionality_reduction(seurat_obj = seurat_obj_dim, 
                              out_path = out_path, 
                              proj_name = proj_name, 
                              log_file = log_file,
                              assay = "RNA",
                              num_pcs = num_dim)

# rescale HTO and ADT
seurat_obj_dim <- NormalizeData(seurat_obj_dim, assay = "HTO", normalization.method = "CLR")
seurat_obj_dim <- NormalizeData(seurat_obj_dim, assay = "ADT", normalization.method = "CLR")

saveRDS(seurat_obj_dim,
        file = glue("{proj_name}.seurat_obj.rds"))

save_seurat_metadata(seurat_obj = seurat_obj_dim,
                     out_path = out_path,
                     log_file = log_file,
                     proj_name = proj_name, 
                     type = "dim")
# save counts mat
save_counts_matrix(seurat_obj = seurat_obj_dim,
                   out_path = out_path,
                   proj_name = proj_name,
                   log_file = log_file,
                   type = "norm")

ADT <- c("CD3-CTCATTGTAACTCCT", 
         "CD4-GAGGTTAGTGATGGA", 
         "CD8-GCGCAACTTGATGAT", 
         "CD19-CTGGGCAATTACTCG", 
         "CD34-GCAGAAATCTCCCTT",
         "CD117-AGACTAATAGCTGAC",
         "CD64-AAGTATGCCCTACGA",
         "CD11b-GACAAGTGATCTGCA")

for( i in ADT) {
  current_plot <- FeaturePlot(s_obj, feature = i) + ggsci::scale_color_material("pink")
  ggsave(current_plot, file = paste(i, "_ADT.png", sep = ""))
}

match_ADT <- c("CD3E",
               "CD3G",
               "CD3D",
               "CD4",
               "CD8A",
               "CD8B",
               "CD19",
               "CD34",
               "KIT",
               "FCGR1A",
               "ITGAM")

for( i in match_ADT) {
  current_plot <- FeaturePlot(s_obj, feature = i) +ggsci::scale_color_material("red")
  ggsave(current_plot, file = paste(i, "_RNA.png", sep = ""))
}

Idents(seurat_obj) <- "snn_res.5" 

cluster_avg_exp <- read.csv("combine_raw_CDNA/expression.mean.{}.csv")
cluster_avg_exp <- cluster_avg_exp_t %>% 
  column_to_rownames("gene") %>% 
  as.matrix()

celltypes_csv = "/home/ay1392/anna_beegfs/TET2/scRNAseq/Analysis/scRNA_scripts/Igor_genes_wbm_cells_global.csv"
celltypes = data.table::fread(file = celltypes_csv, showProgress = FALSE, data.table = FALSE, nThread = 4)

mod_score <- module_score(module_tbl = celltypes,
             counts_norm = NULL, 
             counts_raw = cluster_avg_exp, 
             method = c("rsscore"))

raw_seurat <- function (data_path, sample_names, proj_name, HTO_file, ADT_file, log_file) {
  # Log parameters
  write(glue("data_path = {data_path}
             sample_names = {sample_names}
             out_path = {out_path}
             num_dim = {num_dim}
             HTO_file = {HTO_file}
             proj_name = {proj_name}
             log_file = {log_file})"),
        file = log_file,
        append = TRUE)
  
  counts_mat <- load_sample_counts_matrix(sample_names = sample_names,
                                          data_path = data_path,
                                          log_file = log_file)
  
  # take the counts matrix and make a seurat object
  seurat_obj_1 <- create_seurat_obj(counts_mat = counts_mat, 
                                    out_path = out_path, 
                                    proj_name = proj_name, 
                                    log_file = log_file)
  hto_data <- clean_hto(HTO_file)
  
  # here the proj name and the sample name have to match
  seurat_obj_hto <- create_seurat_obj_hto(seurat_obj = seurat_obj_1, 
                                          HTO_counts = hto_data, 
                                          proj_name = proj_name, 
                                          out_dir = data_path,
                                          log_file = log_file, 
                                          hto_demux_quantile = 0.99,
                                          multi_demux_quantile = 0.7)
  
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
  
  
  ADT_counts <- clean_hto(ADT_file)
  
  seurat_obj_adt <- create_seurat_obj_adt(seurat_obj_hto, out_dir, ADT_counts, proj_name, log_file)
 
}

create_seurat_obj_adt <- function(seurat_obj, out_dir, ADT_counts, proj_name, log_file, 
                                  hto_demux_quantile = 0.99, multi_demux_quantile = 0.7) {
  # add HTO slot to an existing seurat object
  #
  # Args:
  #   seurat_obj: Seurat object
  #   out_dir: Output directory
  #   ADT_counts: ADT counts as a dataframe
  #   proj_name: Name of project and name of output files
  #
  # Results:
  #   A seurat object with HTO in HTO slot using HTO demux and MULTIseqDemux
  
  s_obj <- seurat_obj
  
  message_str <- "\n\n ========== creating ADT slot ========== \n\n"
  write_message(message_str, log_file)
  
  colnames(ADT_counts) <- str_c(proj_name, ":", colnames(ADT_counts))
  
  # check if the cells in the data are the same as the cells in the hashtag data
  cells_to_use <- intersect(colnames(seurat_obj), colnames(ADT_counts))
  
  if(length(s_obj) != length(cells_to_use)){
    message_str <- "some cells in scrna matrix not in hto matrix"
    write_message(message_str, log_file)
  }
  if(ncol(ADT_counts) != length(cells_to_use)){
    message_str <- "some cells in hto matrix not in scrna matrix"
    write_message(message_str, log_file)
  }
  
  # Subset RNA and HTO counts by joint cell barcodes
  ADT_counts <- as.matrix(ADT_counts[, cells_to_use])
  s_obj <- SubsetData(s_obj, cells = cells_to_use)
  
  # add hashtag slot
  s_obj[["ADT"]] <- CreateAssayObject(counts = ADT_counts)
  
  # Normalize HTO data
  s_obj <- NormalizeData(s_obj, assay = "ADT", normalization.method = "CLR")
  s_obj <- ScaleData(s_obj, assay = "ADT")
  
  
  return(s_obj)
}


