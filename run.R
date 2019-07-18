source("/beegfs/ay1392/TET2/scRNAseq/Analysis/scRNAscripts/R/hashtags_create.R")

HTO_path <- "/beegfs/ay1392/CITESeq/third_set/data/genome.med.nyu.edu/results/aifantislab/2019-runs/2019-07-01/cellranger/HTO3_Results/HTO_results.tsv"
HTO = data.table::fread(file = HTO_path, showProgress = FALSE, data.table = FALSE, nThread = 4)
HTO$V1[which(HTO$V1 == "Hashtag1-GTCAACTCTTTAGCG")] <- "hash1_sample_0160"
HTO$V1[which(HTO$V1 == "Hashtag2-TGATGGCCTATTGGG")] <- "hash2_sample_0160"
HTO$V1[which(HTO$V1 == "Hashtag3-TTCCGCCTCTCTTTG")] <- "hash3_sample_3133"
HTO$V1[which(HTO$V1 == "Hashtag4-AGTAAGTTCAGCGTA")] <- "hash4_sample_3762"
HTO$V1[which(HTO$V1 == "Hashtag5-AAGTATCGTTTCGCA")] <- "hash5_sample_3762"
HTO$V1[which(HTO$V1 == "Hashtag9-CAGTAGTCACGGTCA")] <- "hash9_sample_0612"
HTO$V1[which(HTO$V1 == "Hashtag7-TGTCTTTCCTGCCAG")] <- "hash7_sample_0310"

HTO <- HTO %>%  
  column_to_rownames("V1")
write.table(HTO, 
file = "/beegfs/ay1392/CITESeq/third_set/data/genome.med.nyu.edu/results/aifantislab/2019-runs/2019-07-01/cellranger/HTO2_Results/HTO_results_clean.tsv",
row.names = T,
col.names = T)

data_path <- "/beegfs/ay1392/CITESeq/third_set/data/genome.med.nyu.edu/results/aifantislab/2019-runs/2019-07-01/cellranger/count-cDNA2/"
sample_names <- "cdna2"
out_path <- "/Users/anna/Desktop/Audrey_combine"
num_dim <- 30
HTO_file <- "/beegfs/ay1392/CITESeq/third_set/data/genome.med.nyu.edu/results/aifantislab/2019-runs/2019-07-01/cellranger/HTO2_Results/HTO_results_clean.tsv"
sct = FALSE
proj_name = NULL
log_file = NULL
min_genes = NULL
max_genes = NULL
max_mt = NULL



seurat_obj <- assemble_seurat_obj_hto(data_path = data_path,
                                              sample_names = sample_names, # name of samples, can be more than 1
                                              out_path = out_path, # path to deposit outputs
                                              num_dim = num_dim, # number of dimensions to use for PCA, UMAP, and TSNE
                                              HTO_file = HTO_file, # path to HTO file
                                              sct = FALSE, # use ScTransform or not
                                              proj_name = NULL, # name of project
                                              log_file = NULL, # name of log file
                                              min_genes = NULL,
                                              max_genes = NULL,
                                              max_mt = 20) 

# add adt

seurat_obj <- readRDS("")

ADT_file<- "/beegfs/ay1392/CITESeq/third_set/data/genome.med.nyu.edu/results/aifantislab/2019-runs/2019-07-01/cellranger/ADT3_Results/ADT_results.tsv"
ADT_counts <- clean_hto(ADT_file)


proj_name = "cdna3"
log_file = NULL
out_path <- "/home/ay1392/anna_beegfs/CITESeq/third_set/CDNA3"


if (is.null(log_file)) {
  log_file = glue("{out_path}/{proj_name}_create.log")
}

# Set up log file
write(glue("Starting analysis for {proj_name}"), 
      file = log_file, 
      append = TRUE)


# call module scores
module_score_list <- module_score(module_tbl, counts_norm = NULL, counts_raw = NULL, 
                         method = c("iscore", "sscore","zscore","rsscore"))

# plot UMAPs, density maps, barplots
