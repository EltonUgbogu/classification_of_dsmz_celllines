#!/usr/bin/env Rscript
# NTP-style subtype assignment using TCGA subtype templates, then applied to DSMZ
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(tidyverse); library(SummarizedExperiment); library(limma)
})

set.seed(7)

# -------- config (edit) ----------
TCGA_SE_RDS     <- "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment.rds"
TCGA_SUBTYPE_CSV<- "/home/chu25/data/tcga/pancan_subtypes.csv"  # columns: sample, cohort, subtype
DSMZ_RDS        <- "/home/chu25/data/dsmz/DSMZ_count_gene.rds"
OUTDIR          <- "/home/chu25/dsmz/results/subtypes_ntp"
MIN_GENES_PER_TMPL <- 20
ACC_THRESHOLD      <- 0.80
dir.create(OUTDIR, FALSE, TRUE)

# -------- load & harmonise ------
tcga_se <- readRDS(TCGA_SE_RDS)
X       <- assay(tcga_se); rownames(X) <- sub("\\..*$","",rownames(X))
if (any(duplicated(rownames(X)))) X <- rowsum(X, rownames(X), reorder=TRUE)
stb     <- read.csv(TCGA_SUBTYPE_CSV)
stb$sample <- as.character(stb$sample); stb$cohort <- sub("^TCGA-","", as.character(stb$cohort))
stb <- stb %>% filter(sample %in% colnames(X))

dsmz_raw <- readRDS(DSMZ_RDS)
stopifnot(all(c("Ensembl_ID","gene_name") %in% colnames(dsmz_raw)))
scols <- setdiff(colnames(dsmz_raw), c("Ensembl_ID","gene_name","Ensembl_ID_with_version"))
dsmz_raw[scols] <- lapply(dsmz_raw[scols], \(x) as.numeric(as.character(x)))
Xd <- as.matrix(dsmz_raw[, scols, drop=FALSE]); rownames(Xd) <- dsmz_raw$Ensembl_ID
if (any(duplicated(rownames(Xd)))) Xd <- rowsum(Xd, rownames(Xd), reorder=TRUE)

# shared genes
g <- intersect(rownames(X), rownames(Xd)); X <- X[g, , drop=FALSE]; Xd <- Xd[g, , drop=FALSE]

# -------- helpers ----------------
mk_templates <- function(Xcohort, y) {
  # voom+limma DE: each subtype vs rest; return list of signatures (up genes only), as centroids
  dge <- DGEList(counts = round(Xcohort + 0.1))  # small offset; counts-like
  v   <- voom(dge, plot=FALSE)
  lev <- sort(unique(y))
  sigs <- list()
  for (s in lev) {
    design <- model.matrix(~ 0 + factor(ifelse(y==s, s, "OTHER")))
    colnames(design) <- c("OTHER", s)
    fit <- lmFit(v, design)
    cont <- makeContrasts(ctr = get(s) - OTHER, levels = design)
    fit2 <- eBayes(contrasts.fit(fit, cont))
    tt  <- topTable(fit2, number = Inf, sort.by="P")
    sig <- tt %>% filter(adj.P.Val < 0.01 & logFC > 1)
    gsig <- rownames(sig)
    if (length(gsig) >= MIN_GENES_PER_TMPL) {
      # centroid on z-scored expression within cohort
      Z <- t(scale(t(Xcohort[gsig, , drop=FALSE])))
      tmpl <- rowMeans(Z[, y==s, drop=FALSE])
      sigs[[s]] <- list(genes=gsig, centroid=tmpl)
    }
  }
  sigs
}

predict_ntp <- function(Xquery, sigs) {
  # score = Pearson correlation between query z-score and subtype centroid
  if (length(sigs)==0) return(tibble(sample=colnames(Xquery), subtype=NA, score=NA))
  # common support
  G <- sort(unique(unlist(lapply(sigs, `[[`, "genes"))))
  Zq <- t(scale(t(Xquery[G, , drop=FALSE]))); Zq[is.na(Zq)] <- 0
  scores <- sapply(names(sigs), function(s) {
    cg <- sigs[[s]]$genes
    ct <- sigs[[s]]$centroid
    common <- intersect(cg, rownames(Zq))
    if (length(common) < MIN_GENES_PER_TMPL) return(rep(NA_real_, ncol(Zq)))
    cor(Zq[common, , drop=FALSE], ct[common], method="pearson")
  })
  if (is.null(dim(scores))) scores <- cbind(scores) # one subtype edge case
  pred <- apply(scores, 1, function(v) names(which.max(v)))
  scr  <- apply(scores, 1, max, na.rm=TRUE)
  tibble(sample = colnames(Xquery), subtype = pred, score = scr)
}

cv80_acc <- function(Xcohort, y, sigs) {
  set.seed(1)
  idx <- sample(seq_along(y), floor(0.8*length(y)))
  ytr <- y[idx]; yte <- y[-idx]
  Xtr <- Xcohort[, idx, drop=FALSE]; Xte <- Xcohort[, -idx, drop=FALSE]
  sigs_tr <- mk_templates(Xtr, ytr); if (!length(sigs_tr)) return(NA_real_)
  pd <- predict_ntp(Xte, sigs_tr)
  mean(pd$subtype == yte, na.rm=TRUE)
}

# -------- run per cohort --------------
out_all <- list(); k <- 1
for (ct in sort(unique(stb$cohort))) {
  samp <- stb %>% filter(cohort==ct)
  if (n_distinct(samp$subtype) < 2) next
  keep <- samp$sample
  Xc   <- X[, keep, drop=FALSE]
  y    <- as.character(samp$subtype); names(y) <- keep

  sigs <- mk_templates(Xc, y)
  if (!length(sigs)) next
  acc  <- tryCatch(cv80_acc(Xc, y, sigs), error=function(e) NA_real_)
  if (!is.finite(acc) || acc < ACC_THRESHOLD) next

  # apply to DSMZ
  pred_d <- predict_ntp(Xd, sigs) %>% mutate(cohort = ct)
  out_all[[k]] <- pred_d; k <- k + 1
  write.csv(pred_d, file.path(OUTDIR, sprintf("dsmz_ntp_%s.csv", ct)), row.names=FALSE)
}
pred_all <- bind_rows(out_all)
write.csv(pred_all, file.path(OUTDIR, "dsmz_ntp_all_cohorts.csv"), row.names = FALSE)

cat("[OK] Subtype classification complete → ", OUTDIR, "\n")
