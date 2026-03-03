# Statistical analysis functions (DE, Fisher, KW + DSCF)
# Place in R/stats.R

run_deseq_contrasts <- function(all_counts, coldata, ref_level = "Basal", outdir) {
  coldata$subtype <- factor(coldata$subtype, levels = unique(coldata$subtype))
  levels(coldata$subtype) <- gsub("[^A-Za-z0-9_.]", "_", levels(coldata$subtype))
  dds_de <- DESeqDataSetFromMatrix(round(all_counts), coldata, design = ~ subtype)
  dds_de <- DESeq(dds_de)
  present_lvls <- levels(coldata$subtype)
  targets <- setdiff(present_lvls, ref_level)
  de_results <- list()
  for (sub in targets) {
    res <- results(dds_de, contrast = c("subtype", sub, ref_level), alpha = 0.05)
    de_results[[sub]] <- as.data.frame(res)
  }
  dir.create(file.path(outdir, "de_results"), showWarnings = FALSE, recursive = TRUE)
  lapply(names(de_results), function(s) write.csv(de_results[[s]], file.path(outdir, "de_results", sprintf("de_%s_vs_%s.csv", s, ref_level))))
  de_results
}

run_fisher_tests <- function(cluster_list, dataset_lab, tcga_se = NULL, V_tcga = NULL, tcga_sub_col = NULL, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  fisher_results <- list()
  # generic cluster vs dataset tests
  for (nm in names(cluster_list)) {
    ct <- table(cluster_list[[nm]], dataset_lab)
    if (all(dim(ct) >= 2) && sum(ct) >= 5) {
      ft <- fisher.test(ct, simulate.p.value = nrow(ct) > 2)
      fisher_results[[nm]] <- data.frame(test = sprintf("Fisher_%s_vs_dataset", nm), p_value = ft$p.value)
    }
  }
  # Cluster vs PAM50 on TCGA if available
  if (!is.null(tcga_se) && !is.null(tcga_sub_col) && !is.null(V_tcga)) {
    tcga_samples <- colnames(V_tcga)
    lab_tcga <- as.character(colData(tcga_se)[tcga_samples, tcga_sub_col])
    for (nm in names(cluster_list)) {
      ct <- table(cluster_list[[nm]][tcga_samples], lab_tcga)
      if (all(dim(ct) >= 2) && sum(ct) >= 5) {
        ft <- fisher.test(ct, simulate.p.value = nrow(ct) > 2)
        fisher_results[[paste0(nm, "_pam50")]] <- data.frame(test = sprintf("Fisher_%s_vs_PAM50", nm), p_value = ft$p.value)
      }
    }
  }
  if (length(fisher_results)) {
    fisher_df <- do.call(rbind, fisher_results)
    write.csv(fisher_df, file.path(outdir, "fisher_tests.csv"), row.names = FALSE)
    invisible(fisher_df)
  } else {
    invisible(NULL)
  }
}

run_kruskal_and_dscf <- function(mat, cluster_factor, genes, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  kw_results <- list(); dscf_results <- list()
  for (gene in genes) {
    if (!(gene %in% rownames(mat))) next
    expr <- mat[gene, ]
    keep <- !is.na(expr)
    expr_ok <- expr[keep]
    grp_ok <- droplevels(factor(cluster_factor[keep]))
    if (length(unique(as.character(grp_ok))) < 2) next
    kw <- kruskal.test(expr_ok ~ grp_ok)
    kw_results[[gene]] <- data.frame(gene = gene, p_value = kw$p.value, statistic = kw$statistic)
    if (kw$p.value < 0.05) {
      dscf <- try(PMCMRplus::dscfAllPairsTest(expr_ok, grp_ok, p.adjust.method = "BH"), silent = TRUE)
      if (!inherits(dscf, "try-error") && !is.null(dscf$p.value)) {
        p_mat <- as.matrix(dscf$p.value)
        lev_present <- levels(grp_ok)
        if (is.null(rownames(p_mat)) || all(rownames(p_mat) == "")) rownames(p_mat) <- lev_present
        if (is.null(colnames(p_mat)) || all(colnames(p_mat) == "")) colnames(p_mat) <- lev_present
        pairs_idx <- which(upper.tri(p_mat, diag = FALSE), arr.ind = TRUE)
        if (nrow(pairs_idx) > 0) {
          subtype1 <- rownames(p_mat)[pairs_idx[,1]]
          subtype2 <- colnames(p_mat)[pairs_idx[,2]]
          pvals <- p_mat[cbind(pairs_idx[,1], pairs_idx[,2])]
          keep2 <- is.finite(pvals)
          dscf_results[[gene]] <- data.frame(gene = gene, subtype1 = subtype1[keep2], subtype2 = subtype2[keep2], p_value = pvals[keep2], stringsAsFactors = FALSE)
        }
      }
    }
  }
  if (length(kw_results)) write.csv(do.call(rbind, kw_results), file.path(outdir, "kruskal_wallis_tests.csv"), row.names = FALSE)
  if (length(dscf_results)) write.csv(do.call(rbind, dscf_results), file.path(outdir, "dscf_tests.csv"), row.names = FALSE)
  list(kw = kw_results, dscf = dscf_results)
}