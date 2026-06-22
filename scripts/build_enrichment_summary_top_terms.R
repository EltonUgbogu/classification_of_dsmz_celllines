#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table)})
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- sub("^--file=", "", script_arg)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
pipeline_root <- Sys.getenv("PIPELINE_ROOT", unset = repo_root)
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[[i+1]] }
infile <- get_arg('--input')
outfile <- get_arg('--output', 'results/unsupervised/enrichment_summary_top_terms.tsv')
manifest <- get_arg('--query-manifest', NA_character_)
if (is.na(infile) || !file.exists(infile)) stop('--input is required and must exist')
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
dt <- fread(infile)
rename_first <- function(x, candidates, target) {
  if (target %in% names(x)) return(x)
  hit <- intersect(candidates, names(x))[1]
  if (!is.na(hit) && hit != target) setnames(x, hit, target)
  x
}
dt <- rename_first(dt, c('query','query_name','query_id'), 'query_id')
dt <- rename_first(dt, c('name','term_name','term'), 'term_name')
dt <- rename_first(dt, c('p_value','pval','p.value'), 'p_value')
dt <- rename_first(dt, c('adjusted_p_value','p_value_adjusted','padj','p_adj'), 'adjusted_p_value')
dt <- rename_first(dt, c('source','domain'), 'source')
required <- c('query_id','term_name','source','p_value')
missing <- setdiff(required, names(dt)); if (length(missing)) stop('Missing required columns after normalisation: ', paste(missing, collapse=', '))
if (!'adjusted_p_value' %in% names(dt)) dt[, adjusted_p_value := as.numeric(p_value)]
parts <- tstrsplit(dt$query_id, '__', fixed = TRUE, fill = NA_character_)
if (length(parts) < 2L) parts[[2L]] <- rep(NA_character_, nrow(dt))
if (!'cohort' %in% names(dt)) dt[, cohort := fifelse(!is.na(parts[[1]]), parts[[1]], NA_character_)]
dt[, marker_direction := fifelse(grepl('__down$|_down$', query_id), 'down', fifelse(grepl('__up$|_up$', query_id), 'up', 'all'))]
if (!'contrast_group' %in% names(dt)) dt[, contrast_group := fifelse(!is.na(parts[[2]]), parts[[2]], NA_character_)]
if (!is.na(manifest) && file.exists(manifest)) {
  man <- fread(manifest)
  if (!'query_id' %in% names(man) && 'query_name' %in% names(man)) setnames(man, 'query_name', 'query_id')
  if ('query_id' %in% names(man)) {
    if (!'category' %in% names(man) && 'query_family' %in% names(man)) man[, category := query_family]
    if (!'cohort' %in% names(man) && 'owner_profile' %in% names(man)) man[, cohort := owner_profile]
    if (!'contrast' %in% names(man) && 'source_contrast' %in% names(man)) man[, contrast := source_contrast]
    if (!'source_marker_table' %in% names(man) && 'source_marker_table_path' %in% names(man)) man[, source_marker_table := source_marker_table_path]
    if (!'source_marker_table' %in% names(man) && 'source_table' %in% names(man)) man[, source_marker_table := source_table]
    if (!'gene_list_path' %in% names(man) && 'genes_path' %in% names(man)) man[, gene_list_path := genes_path]
    if (!'gene_list_path' %in% names(man) && 'genes_tsv' %in% names(man)) man[, gene_list_path := genes_tsv]
    if (!'background_path' %in% names(man) && 'background_tsv' %in% names(man)) man[, background_path := background_tsv]
    keep <- intersect(c('query_id','category','cohort','contrast','direction',
                        'source_marker_table','gene_list_path','background_path',
                        'profile','disease','group_id','query_family',
                        'gene_count','background_count'), names(man))
    meta <- unique(man[, ..keep], by = 'query_id')
    merge_cols <- setdiff(keep, 'query_id')
    existing <- intersect(merge_cols, names(dt))
    for (col in existing) setnames(dt, col, paste0(col, '.result'))
    dt <- merge(dt, meta, by='query_id', all.x=TRUE, sort=FALSE)
    for (col in existing) {
      result_col <- paste0(col, '.result')
      dt[, (col) := fifelse(!is.na(get(result_col)) & nzchar(as.character(get(result_col))),
                            as.character(get(result_col)), as.character(get(col)))]
      dt[, (result_col) := NULL]
    }
  }
}
for (col in c('category','cohort','contrast','direction','source_marker_table')) {
  if (!col %in% names(dt)) dt[, (col) := NA_character_]
}
setorder(dt, query_id, adjusted_p_value, p_value, source, term_name)
fwrite(dt, outfile, sep='\t')
prov <- data.table(figure_name='Fig_enrichment_top_terms_heatmap.pdf', script='scripts/build_enrichment_summary_top_terms.R', command=paste(c('Rscript','scripts/build_enrichment_summary_top_terms.R','--input',infile,'--output',outfile,'--query-manifest',manifest), collapse=' '), git_commit=Sys.getenv('GIT_COMMIT', unset='unavailable_not_git_worktree'), timestamp=format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z'), input_files=paste(na.omit(c(normalizePath(infile, mustWork=FALSE), if (!is.na(manifest)) normalizePath(manifest, mustWork=FALSE))), collapse=';'), output_files=normalizePath(outfile, mustWork=FALSE), upstream_tables=normalizePath(infile, mustWork=FALSE), key_parameters='deterministic_column_normalisation;query_metadata_preserved', software_versions=paste0('R=', getRversion(), ';data.table=', packageVersion('data.table')), figure_type='enrichment', source_pipeline_root=pipeline_root, copied_to_figure_export_path='', legacy_source_path='', notes='Builds enrichment_summary_top_terms.tsv from gprofiler top_terms.tsv')
fwrite(prov, sub('\\.tsv$', '_provenance.tsv', outfile), sep='\t')
cat('[OK] wrote ', outfile, '\n', sep='')
