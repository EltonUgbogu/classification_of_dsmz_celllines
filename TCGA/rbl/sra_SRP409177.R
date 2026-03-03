#!/usr/bin/env Rscript

# === Config ===
BASE      <- Sys.getenv("BASE",      "/home/chu25/data/rbl")
SRA_DIR   <- file.path(BASE, "sra")
FASTQ_DIR <- file.path(BASE, "fastq")
TMPDIR    <- Sys.getenv("TMPDIR", file.path(BASE, "tmp"))
THREADS   <- as.integer(Sys.getenv("THREADS", "8"))

# Optional: read SRRs from a file passed on the command line (one per line).
# If no file given, use the 47 SRR IDs from SRP409177 (your neuroblastoma set).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1 && file.exists(args[1])) {
  SRRs <- scan(args[1], what = character(), quiet = TRUE)
} else {
  SRRs <- c(
    "SRR22373268","SRR22373269","SRR22373270","SRR22373271","SRR22373272",
    "SRR22373273","SRR22373274","SRR22373275","SRR22373276","SRR22373277",
    "SRR22373278","SRR22373279","SRR22373280","SRR22373281","SRR22373282",
    "SRR22373283","SRR22373284","SRR22373285","SRR22373286","SRR22373287",
    "SRR22373288","SRR22373289","SRR22373290","SRR22373291","SRR22373292",
    "SRR22373293","SRR22373294","SRR22373295","SRR22373296","SRR22373297",
    "SRR22373298","SRR22373299","SRR22373300","SRR22373301","SRR22373302",
    "SRR22373303","SRR22373304","SRR22373305","SRR22373306","SRR22373307",
    "SRR22373308","SRR22373309","SRR22373310","SRR22373311","SRR22373312",
    "SRR22373313","SRR22373314"
  )
}

# === Pre-flight checks ===
dir.create(SRA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FASTQ_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TMPDIR, recursive = TRUE, showWarnings = FALSE)

prefetch_bin     <- Sys.which("prefetch")
fasterq_dump_bin <- Sys.which("fasterq-dump")
pigz_bin         <- Sys.which("pigz")
gzip_bin         <- if (nzchar(pigz_bin)) pigz_bin else Sys.which("gzip")

if (!nzchar(prefetch_bin) || !nzchar(fasterq_dump_bin)) {
  stop("sra-tools not found on PATH. Activate your env (e.g., 'conda activate sra3').")
}

# TLS CA bundle (helps avoid certificate errors on some clusters)
if (file.exists("/etc/ssl/certs/ca-certificates.crt")) {
  Sys.setenv(SSL_CERT_FILE = "/etc/ssl/certs/ca-certificates.crt")
} else if (file.exists("/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem")) {
  Sys.setenv(SSL_CERT_FILE = "/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem")
}

run_cmd <- function(cmd, args, workdir = NULL) {
  cat(sprintf(">> %s %s\n", cmd, paste(shQuote(args), collapse = " ")))
  status <- system2(cmd, args = args, stdout = "", stderr = "")
  if (status != 0) stop(sprintf("Command failed: %s", paste(c(cmd, args), collapse = " ")))
}

locate_src <- function(srr) {
  d <- file.path(SRA_DIR, srr)
  f <- file.path(SRA_DIR, paste0(srr, ".sra"))
  if (dir.exists(d)) return(d)
  if (file.exists(f)) return(f)
  return(srr) # let fasterq-dump fetch on the fly if needed
}

# === Download & convert ===
cat(sprintf("SRRs: %d\n", length(SRRs)))
cat(sprintf("SRA_DIR=%s\nFASTQ_DIR=%s\nTMPDIR=%s\nTHREADS=%d\n\n", SRA_DIR, FASTQ_DIR, TMPDIR, THREADS))

for (srr in SRRs) {
  cat(sprintf("==== %s ====\n", srr))

  # Prefetch (skip if already present)
  has_dir <- dir.exists(file.path(SRA_DIR, srr))
  has_sra <- file.exists(file.path(SRA_DIR, paste0(srr, ".sra")))
  if (!has_dir && !has_sra) {
    run_cmd(prefetch_bin, c("--max-size", "0", "-O", SRA_DIR, srr))
  } else {
    cat(".. SRA object already present; skipping prefetch\n")
  }

  # Convert to FASTQ (skip if fastqs already gzipped)
  fq_pattern <- file.path(FASTQ_DIR, paste0(srr, "_*.fastq.gz"))
  already_done <- length(Sys.glob(fq_pattern)) > 0
  if (!already_done) {
    src <- locate_src(srr)
    run_cmd(fasterq_dump_bin, c(
      src,
      "--split-files",
      "--skip-technical",
      "--threads", THREADS,
      "--temp", TMPDIR,
      "--outdir", FASTQ_DIR
    ))

    # compress any produced .fastq
    fq_plain <- Sys.glob(file.path(FASTQ_DIR, paste0(srr, "_*.fastq")))
    if (length(fq_plain) > 0) {
      if (nzchar(pigz_bin)) {
        run_cmd(pigz_bin, c("-p", THREADS, fq_plain))
      } else if (nzchar(gzip_bin)) {
        # gzip sequentially
        for (f in fq_plain) run_cmd(gzip_bin, f)
      } else {
        cat("!! No pigz/gzip found; leaving .fastq uncompressed\n")
      }
    } else {
      cat("!! No .fastq produced; check logs above.\n")
    }
  } else {
    cat(".. FASTQ(.gz) already exists; skipping conversion\n")
  }

  cat(sprintf("---- done %s ----\n\n", srr))
}

# === Summary ===
gz_files <- unlist(lapply(SRRs, function(srr) Sys.glob(file.path(FASTQ_DIR, paste0(srr, "_*.fastq.gz")))))
cat("Completed. Gzipped FASTQs found:\n")
if (length(gz_files)) {
  cat(paste0(" - ", basename(gz_files)), sep = "\n")
} else {
  cat(" - none found\n")
}
