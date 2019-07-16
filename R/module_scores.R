
# module scores 

# scoring methods are from Igor and Teo

# scoring qc - silhouette

# compare scoring 

# scoring plots


module_score <- function(module_tbl, counts, method = c("Iscore", "Yscore","Sscore","Zscore", "Pscore", "RSscore")) {
  # returns module scores in a uniform way for the method specified
  # 
  # module tbl should be in long format celltype and gene

  scoring_methods <- c(Iscore = "run_Iscore")
  
  # run over scoring methods and get outputs and merge them to one output
  eval(as.name(scoring_methods))(module_tbl, counts)
  

    Yscore = lapply(module_list, function(x, mat) colMeans(mat[x, ]), mat = scale_rows_01(counts[module_genes,])),  # Igor's score
    Sscore = lapply(celltypes, calc_boot_score, score = boot_seurat, rmodule = rmodule_binned),
    Zscore = lapply(celltypes, calc_boot_score, score = boot_zval, rmodule = rmodule_smooth),
    Pscore = lapply(celltypes, calc_boot_score, score = boot_pval, rmodule = rmodule_binned),
    RSscore = calculate_celltypes_rescaled(counts, celltypes_tbl)
  
}

run_Iscore <- function(module_tbl, counts) {
  # module tbl should be in long format with celltype and gene
  
  # subset counts to only have genes to be used
  counts_sub <- semi_join(counts, module_tbl, by = "gene")
  counts_sub = counts_sub %>% as.data.frame() %>% column_to_rownames("gene") %>% as.matrix()

  iscores <- module_tbl %>% 
    group_by(celltype) %>% 
    do(
      s <- simple_score(counts_sub, .$gene)
    ) %>% 
    spread(celltype, score) %>% 
    rename_at(vars(-contains("cell")), list(~paste0(., "_Iscore")))
  
  return(iscores)
}

run_Yscore <- function(module_tbl, counts) {
  # module tbl should be in long format with celltype and gene
  
  # subset counts to only have genes to be used
  counts_sub <- semi_join(counts, ., by = "gene"))
  counts_sub = counts_sub %>% as.data.frame() %>% column_to_rownames("gene") %>% as.matrix()

  yscores <- module_tbl %>% 
    group_by(celltype)
  
}


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
  out = as.data.frame(colMeans(mat[ix, , drop = FALSE]) - colMeans(mat[!ix, , drop = FALSE]))
  out$cell = colnames(mat)
  colnames(out) <- c("score", "cell")
  return(out)
}

boot_seurat = function(mat, module, rmodule, boot_size = 100L, seed = 1234L) {
  set.seed(seed)
  boot = colMeans(mat[module, ]) - get_boot_scores(mat, rmodule, boot_size)
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

calc_boot_score = function(celltype, score, rmodule) 
  score(counts, module_list[[celltype]], rmodule[[celltype]])

softmax = function(x, mask = NULL) {
  y = exp(x)
  if (!is.null(mask))
    y[x < mask] = 0
  mass = rowSums(y)
  mass[mass == 0] = 1  # avoid division by 0
  y / mass
}


calculate_celltypes_rescaled = function(exp_mat, celltypes_tbl, min_cpm = 0, limit_pct = 1) {
  # perform the cell type enrichment calculation based on rescaled values
  
  if (class(exp_clust_mat) != "matrix") { stop("expression matrix is not a matrix") }
  if (max(exp_mat) < 100) { stop("expression values appear to be log-scaled") }
  
  # filter matrix for expressed genes only
  exp_mat = filter_mat_by_cpm(exp_mat = exp_mat, min_cpm = min_cpm)
  
  celltype_genes = celltypes_tbl %>% pull(gene)
  celltype_genes = intersect(celltype_genes, rownames(exp_mat))
  
  # rescale matrix for expressed genes only
  exp_mat = normalize_mat_by_gene(exp_mat = exp_mat[celltype_genes, ], limit_pct = limit_pct)
  
  # convert celltypes to list
  celltypes_tbl = celltypes_tbl %>% filter(gene %in% celltype_genes) %>% arrange(celltype)
  celltypes_list = celltypes_tbl %>% split(x = .$gene, f = .$celltype)
  
  # check if enough genes pass filter
  if (min(lengths(celltypes_list)) < 3) { stop("too few genes per celltype") }
  
  # calculate average z-score per celltype
  celltype_scores_tbl = tibble()
  for (ct in names(celltypes_list)) {
    celltype_scores_tbl =
      bind_rows(
        celltype_scores_tbl,
        tibble(
          cluster = colnames(exp_mat),
          celltype = ct,
          num_celltype_genes = length(celltypes_list[[ct]]),
          score = colMeans(exp_mat[celltypes_list[[ct]], ])
        )
      )
    ct_scores = colnames(exp_mat)
  }
  
  # keep only the top scoring celltype for each cluster
  celltype_scores_top_tbl =
    celltype_scores_tbl %>%
    arrange(desc(score)) %>%
    group_by(cluster) %>%
    top_n(1, score) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(score = round(score, 3))
  
  return(celltype_scores_top_tbl)
  
}

# filter matrix by a specified CPM value (higher in at least one sample/column for each gene/row)
filter_mat_by_cpm = function(exp_mat, min_cpm = 0) {
  
  if (class(exp_clust_mat) != "matrix") { stop("expression matrix is not a matrix") }
  if (max(exp_mat) < 100) { stop("expression values appear to be log-scaled") }
  if (nrow(exp_mat) < 10000) { stop("expression matrix has too few genes") }
  
  # expression level equivalent to 1 CPM (1 for 1m, 0.01 for 10k)
  exp_cpm1 = (1 / 1000000) * median(colSums(exp_mat))
  # expression level equivalent to the requested CPM
  min_exp = exp_cpm1 * min_cpm
  # filtered expression matrix
  exp_mat = exp_mat[matrixStats::rowMaxs(exp_mat) > min_exp, ]
  
  return(exp_mat)
  
}

rescale_vector = function(x, limit_pct = 1) {
  x / quantile(x, limit_pct)
}

normalize_mat_by_gene = function(exp_mat, limit_pct = 1) {
  
  if (limit_pct > 1) { stop("percentile should be expressed as a fraction") }
  if (class(exp_clust_mat) != "matrix") { stop("expression matrix is not a matrix") }
  if (max(exp_mat) < 100) { stop("expression values appear to be log-scaled") }
  
  exp_mat = apply(exp_mat, MARGIN = 1, FUN = rescale_vector, limit_pct = limit_pct)
  exp_mat = t(exp_mat)
  exp_mat[exp_mat > 1] = 1
  
  return(exp_mat)
  
}

