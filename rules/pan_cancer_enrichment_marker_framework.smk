# Explicit current pan-cancer marker-framework enrichment targets.
# These targets prepare local query/background files and optional heatmap inputs.
# They do not run g:Profiler, DESeq2, or pan-cancer feature construction.

MARKER_FRAMEWORK_QUERY_OUTDIR_REL = os.path.join(
    "results", "unsupervised", "pan_cancer", "enrichment", "query_sets", "marker_framework"
)
MARKER_FRAMEWORK_QUERY_MANIFEST = os.path.join(MARKER_FRAMEWORK_QUERY_OUTDIR_REL, "query_manifest.tsv")
MARKER_FRAMEWORK_QUERY_SCRIPT = os.path.join(
    SCRIPTS_DIR, "build_pan_cancer_enrichment_marker_framework_query_sets.py"
)
MARKER_FRAMEWORK_GPROFILER_EXPORT = os.path.join(
    "results", "unsupervised", "pan_cancer", "enrichment", "gprofiler_exports",
    "marker_framework_gprofiler_export.csv"
)
MARKER_FRAMEWORK_GPROFILER_EXPORT_PARSER = os.path.join(
    SCRIPTS_DIR, "parse_gprofiler_export_to_enrichment_summary.py"
)
MARKER_FRAMEWORK_HEATMAP_OUTDIR_REL = os.path.join(
    "results", "unsupervised", "pan_cancer", "enrichment", "figures", "marker_framework"
)
MARKER_FRAMEWORK_HEATMAP_SUMMARY = os.path.join(
    MARKER_FRAMEWORK_HEATMAP_OUTDIR_REL, "marker_framework_enrichment_summary_top_terms.tsv"
)

rule build_pan_cancer_enrichment_marker_framework_query_sets:
    """
    Build current pan-cancer marker-framework g:Profiler query/background files.
    This prepares local inputs only; it does not run g:Profiler.
    """
    input:
        features = GPROFILER_FEATURE_TABLE,
        clean_features = os.path.join(
            "results", "unsupervised", "pan_cancer", "feature_space", "pan_cancer_features_clean.txt"
        ),
        retained_source_contrasts = (
            [GPROFILER_RETAINED_SOURCE_CONTRASTS]
            if GPROFILER_RETAINED_SOURCE_CONTRASTS else []
        ),
        cohort_manifest = os.path.join(GPROFILER_QUERY_OUTDIR_REL, "query_manifest.tsv"),
        script = MARKER_FRAMEWORK_QUERY_SCRIPT
    output:
        manifest = MARKER_FRAMEWORK_QUERY_MANIFEST,
        query_counts = os.path.join(MARKER_FRAMEWORK_QUERY_OUTDIR_REL, "query_set_counts.tsv"),
        background_validation = os.path.join(MARKER_FRAMEWORK_QUERY_OUTDIR_REL, "background_validation.tsv"),
        category_availability = os.path.join(MARKER_FRAMEWORK_QUERY_OUTDIR_REL, "category_availability.tsv"),
        validation_summary = os.path.join(MARKER_FRAMEWORK_QUERY_OUTDIR_REL, "validation_summary.tsv")
    params:
        outdir = MARKER_FRAMEWORK_QUERY_OUTDIR_REL,
        cohort_query_dir = GPROFILER_QUERY_OUTDIR_REL,
        profiles = "brca,nbl,rbl"
    log:
        os.path.join(LOGROOT, "build_pan_cancer_enrichment_marker_framework_query_sets.log")
    conda:
        CONDA_ENV_PY
    shell:
        r"""
        set -euo pipefail
        python "{input.script}" \
          --features "{input.features}" \
          --clean-features "{input.clean_features}" \
          --retained-source-contrasts "{input.retained_source_contrasts}" \
          --cohort-query-dir "{params.cohort_query_dir}" \
          --outdir "{params.outdir}" \
          --profiles "{params.profiles}" \
          > "{log}" 2>&1
        test -s "{output.manifest}"
        test -s "{output.query_counts}"
        test -s "{output.background_validation}"
        test -s "{output.category_availability}"
        test -s "{output.validation_summary}"
        """

rule build_pan_cancer_enrichment_selected_term_heatmap:
    """
    Build a current selected-term heatmap from an explicit manual g:Profiler export.
    The export must be downloaded separately; this rule does not run g:Profiler.
    """
    input:
        export = MARKER_FRAMEWORK_GPROFILER_EXPORT,
        query_manifest = MARKER_FRAMEWORK_QUERY_MANIFEST,
        feature_table = GPROFILER_FEATURE_TABLE,
        parser = MARKER_FRAMEWORK_GPROFILER_EXPORT_PARSER,
        plot_script = os.path.join(SCRIPTS_DIR, "plot_enrichment_top_terms_heatmap.R")
    output:
        summary = MARKER_FRAMEWORK_HEATMAP_SUMMARY,
        summary_provenance = MARKER_FRAMEWORK_HEATMAP_SUMMARY + ".provenance.tsv",
        pdf = os.path.join(MARKER_FRAMEWORK_HEATMAP_OUTDIR_REL, "Fig_enrichment_top_terms_heatmap.pdf"),
        png = os.path.join(MARKER_FRAMEWORK_HEATMAP_OUTDIR_REL, "Fig_enrichment_top_terms_heatmap.png"),
        matrix = os.path.join(MARKER_FRAMEWORK_HEATMAP_OUTDIR_REL, "Fig_enrichment_top_terms_heatmap_matrix.tsv"),
        selected = os.path.join(MARKER_FRAMEWORK_HEATMAP_OUTDIR_REL, "Fig_enrichment_top_terms_heatmap_selected_terms.tsv"),
        excluded = os.path.join(MARKER_FRAMEWORK_HEATMAP_OUTDIR_REL, "Fig_enrichment_top_terms_heatmap_excluded_terms.tsv"),
        provenance = os.path.join(MARKER_FRAMEWORK_HEATMAP_OUTDIR_REL, "Fig_enrichment_top_terms_heatmap_provenance.tsv")
    params:
        outdir = MARKER_FRAMEWORK_HEATMAP_OUTDIR_REL
    log:
        os.path.join(LOGROOT, "build_pan_cancer_enrichment_selected_term_heatmap.log")
    conda:
        CONDA_ENV_R
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        python "{input.parser}" \
          --input "{input.export}" \
          --query-manifest "{input.query_manifest}" \
          --output "{output.summary}" \
          > "{log}" 2>&1
        Rscript "{input.plot_script}" \
          --input "{output.summary}" \
          --outdir "{params.outdir}" \
          --feature-table "{input.feature_table}" \
          >> "{log}" 2>&1
        test -s "{output.summary}"
        test -s "{output.summary_provenance}"
        test -s "{output.pdf}"
        test -s "{output.provenance}"
        """
