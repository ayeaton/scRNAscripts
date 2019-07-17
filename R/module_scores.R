
# module scores 

# scoring methods are from Igor and Teo

# scoring qc - silhouette

# compare scoring 

# scoring plots


module_score <- function(module_tbl, counts_norm = NULL, counts_raw = NULL, 
                         method = c("iscore", "sscore","zscore","rsscore")) {
  # returns module scores in a uniform way for the method specified
  # 
  # Args:
  #   module_tbl: a tibble in long format with celltype and gene
  #   counts_norm: normalized counts table for use with Teo's scores
  #   count_raw: raw counts table for use with Igor's scores
  #   method: character vector with the score 
  #
  # Results:
  #   list of dataframes of module scores
  
  scoring_methods <- c(iscore = "run_iscore",
                       sscore = "run_sscore",
                       zscore = "run_zscore",
                       rsscore = "run_rsscore")
  
  current_methods <- scoring_methods[method]
  
  module_scores <- lapply(current_methods, function(x){
    scores <- eval(as.name(x))(module_tbl= module_tbl, 
                                             counts_norm = counts_norm,
                                             counts_raw = counts_raw)
    
    scores_max <- max_scores(scores = scores, method = names(x))
  })
}

max_scores <- function(scores, method, threshold = 0) {
  # gets max pos or zero 
  # idk how to do Igors yet
  scores$unknown <- rep(threshold, nrow(scores))
  scores <- scores %>% 
    column_to_rownames("cell")
  scores$module <- colnames(scores)[apply(scores,1,which.max)]
  scores <- scores %>% 
    separate(module, c("module", NA))
  var <- sym(paste(eval(method), "_module", sep = ""))
  scores <- scores %>% 
    rename(!!var := module) %>% 
    select(-contains("unknown"))
  return(scores)
}

run_iscore <- function(module_tbl, counts_norm, counts_raw = NULL) {
  # module tbl should be in long format with celltype and gene
  
  module_tbl <- module_tbl %>% 
    filter(.$gene %in% rownames(counts_norm)) 
  
  counts_sub <- counts_norm[which(rownames(counts_norm) %in% module_tbl$gene),]
  counts_sub <- scale_rows(counts_sub)

  iscores <- module_tbl %>% 
    group_by(celltype) %>% 
    do(
      s <- simple_score(counts_sub, .$gene)
    ) %>% 
    spread(celltype, scores) %>% 
    rename_at(vars(-contains("cell")), list(~paste0(., "_Iscore")))
  
  return(iscores)
}

run_sscore <- function(module_tbl, counts_norm, counts_raw = NULL) {
  # module tbl should be in long format with celltype and gene

  module_list <- module_tbl %>%
    filter(.$gene %in% rownames(counts_norm)) %>% 
    with(split(.$gene, celltype))
  
  celltypes <- unique(module_tbl$celltype)
  
  rmodule_binned <- lapply(module_list, get_binned_module, 
                           binned_genes = bin_genes(counts_norm, nbins = 25L))
  
  score = vapply(celltypes, calc_boot_score, numeric(ncol(counts_norm)), score = boot_seurat, rmodule = rmodule_binned)
  sscore <- score %>% 
    as.data.frame() %>% 
    rownames_to_column("cell") %>% 
    rename_at(vars(-contains("cell")), list(~paste0(., "_Sscore")))
  return(sscore)
}

run_zscore <- function(module_tbl, counts_norm, counts_raw = NULL) {
  # module tbl should be in long format with celltype and gene
  
  module_list <- module_tbl %>%
    filter(.$gene %in% rownames(counts_norm)) %>% 
    with(split(.$gene, celltype))
  
  celltypes <- unique(module_tbl$celltype)
  
  rmodule_smooth = lapply(module_list, get_smooth_module, ave = rowMeans(counts_norm))
  
  score = vapply(celltypes, calc_boot_score, numeric(ncol(counts_norm)), score = boot_zval, rmodule = rmodule_smooth)
  zscore <- score %>% 
    as.data.frame() %>% 
    rownames_to_column("cell") %>% 
    rename_at(vars(-contains("cell")), list(~paste0(., "_Zscore")))
  return(zscore)
}

run_rsscore = function(module_tbl, counts_raw, min_cpm = 0, limit_pct = 1, counts_norm = NULL) {
  # perform the cell type enrichment calculation based on rescaled values
  
  module_list <- module_tbl %>%
    filter(.$gene %in% rownames(counts_norm)) %>% 
    with(split(.$gene, celltype))
  
  if (class(counts_raw) != "matrix") { stop("expression matrix is not a matrix") }
  if (max(counts_raw) < 100) { stop("expression values appear to be log-scaled") }
  
  # filter matrix for expressed genes only
  counts_raw = filter_mat_by_cpm(counts_raw = counts_raw, min_cpm = min_cpm)
  
  # rescale matrix for expressed genes only
  counts_raw_subs = normalize_mat_by_gene(counts_raw = counts_raw[unlist(module_list), ], limit_pct = limit_pct)
  

  # check if enough genes pass filter
  if (min(lengths(module_list)) < 3) { stop("too few genes per celltype") }
  
  # calculate average z-score per celltype
  celltype_scores_tbl = tibble()
  for (ct in names(module_list)) {
    celltype_scores_tbl =
      bind_rows(
        celltype_scores_tbl,
        tibble(
          cell = colnames(counts_raw_subs),
          celltype = ct,
          score = colMeans(counts_raw_subs[module_list[[ct]], ])
        )
      )
    ct_scores = colnames(counts_raw_subs)
  }

  celltype_scores_tbl <- celltype_scores_tbl %>% 
    spread(celltype, score) %>% 
    rename_at(vars(-contains("cell")), list(~paste0(., "_RSscore")))
  
  return(celltype_scores_tbl)
}

filter_mat_by_cpm = function(counts_raw, min_cpm = 0) {
  # filter matrix by a specified CPM value (higher in at least one sample/column for each gene/row)
  
  if (class(counts_raw) != "matrix") { stop("expression matrix is not a matrix") }
  if (max(counts_raw) < 100) { stop("expression values appear to be log-scaled") }
  if (nrow(counts_raw) < 10000) { stop("expression matrix has too few genes") }
  
  # expression level equivalent to 1 CPM (1 for 1m, 0.01 for 10k)
  exp_cpm1 = (1 / 1000000) * median(colSums(counts_raw))
  # expression level equivalent to the requested CPM
  min_exp = exp_cpm1 * min_cpm
  # filtered expression matrix
  counts_raw = counts_raw[matrixStats::rowMaxs(counts_raw) > min_exp, ]
  
  return(counts_raw)
  
}

rescale_vector = function(x, limit_pct = 1) {
  x / quantile(x, limit_pct)
}

normalize_mat_by_gene = function(counts_raw, limit_pct = 1) {
  
  if (limit_pct > 1) { stop("percentile should be expressed as a fraction") }
  if (class(counts_raw) != "matrix") { stop("expression matrix is not a matrix") }
  if (max(counts_raw) < 100) { stop("expression values appear to be log-scaled") }
  
  counts_raw = apply(counts_raw, MARGIN = 1, FUN = rescale_vector, limit_pct = limit_pct)
  counts_raw = t(counts_raw)
  counts_raw[counts_raw > 1] = 1
  
  return(counts_raw)
  
}
calc_boot_score = function(celltype, score, rmodule) 
  score(counts_norm, module_list[[celltype]], rmodule[[celltype]])

bin_genes = function(mat, nbins = 25L) {
  # given a matrix it splits the genes into `bins`
  # and also returns a vector (`gene2bin`) mapping genes to bins
  ave = rowMeans(mat)
  ave_bins = as.integer(Hmisc::cut2(ave, m = round(nrow(mat) / nbins)))
  names(ave_bins) = rownames(mat)
  list(bins = split(rownames(mat), ave_bins), gene2bin = ave_bins)
}

get_binned_module = function(module, binned_genes) {
  # given a list of binned genes and a gene-module, it returns a function
  # that generates `n` "similar" random modules
  pools = with(binned_genes, bins[gene2bin[module]])
  pools = lapply(pools, setdiff, module)
  function(n) 
    vapply(pools, sample, rep("", n), size = n, replace = TRUE, USE.NAMES = FALSE)
}
get_smooth_module = function(module, ave, p = 2, scale. = 1) {
  # given a gene `module` and a named vector of gene averages (`ave`)
  # it returns a function that generates `n` "similar" random modules
  # `p` is used to scale the negative exponential distances
  # p = 2 -> gaussian kernel around each gene, p = 1 -> laplacian kernel etc
  ix = names(ave) %in% module
  module_mean = ave[ix]
  ave = ave[!ix]
  genes = names(ave)
  function(n)
    vapply(module_mean, function(x) {
      d = abs(x - ave) / scale.
      prob = exp(-d^p)
      sample(genes, n, prob = prob, replace = TRUE)
    }, rep("", n), USE.NAMES = FALSE)
}
get_boot_scores = function(mat, rmodule, boot_size = 100L) 
  apply(rmodule(boot_size), 1L, function(x) colMeans(mat[x, , drop = FALSE]))

scale_rows = function(mat, epsilon = 1e-5) (mat - rowMeans(mat)) / (rowSds(mat) + epsilon)
scale_rows_01 = function(mat) mat / rowMaxs(mat)
scale_rows_rank = function(mat) t(apply(mat[module_genes,], 1L, rank))

simple_score = function(mat, module) {
  # assuming mat is scaled by row
  ix = rownames(mat) %in% module
  out <- as.data.frame(colMeans(mat[ix, , drop = FALSE]) - colMeans(mat[!ix, , drop = FALSE]))
  out$cell <- rownames(out)
  colnames(out) <- c("scores", "cell")
  return(out)
}
boot_seurat = function(mat, module, rmodule, boot_size = 100L, seed = 1234L) {
  set.seed(seed)
  boot = colMeans(mat[which(rownames(mat) %in% module), ]) - get_boot_scores(mat, rmodule, boot_size)
  rowMeans(boot)
}
boot_zval = function(mat, module, rmodule, boot_size = 100L, seed = 1234L) {
  set.seed(seed)
  boot = cbind(colMeans(mat[module, ]), get_boot_scores(mat, rmodule, boot_size))
  boot = scale_rows(boot)
  boot[, 1L]
}
boot_pval = function(mat, module, rmodule, boot_size = 100L, seed = 1234L) {
  set.seed(seed)
  boot = colMeans(mat[module, ]) > get_boot_scores(mat, rmodule, boot_size)
  rowMeans(boot)
}
softmax = function(x, mask = NULL) {
  y = exp(x)
  if (!is.null(mask))
    y[x < mask] = 0
  mass = rowSums(y)
  mass[mass == 0] = 1  # avoid division by 0
  y / mass
}

