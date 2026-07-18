# Canonical graph-derived functional-enrichment workflow.

if FUNCTIONAL_ENRICHMENT_ENABLED and MARKER_POST_ENABLED:
    FE_QUERY_CFG = FUNCTIONAL_ENRICHMENT_CFG.get("query_construction", {})
    FE_GPROFILER_CFG = FUNCTIONAL_ENRICHMENT_CFG.get("gprofiler", {})
    FE_ORDERED_SENSITIVITY_CFG = FUNCTIONAL_ENRICHMENT_CFG.get("ordered_query_sensitivity", {})
    FE_HEATMAP_CFG = FUNCTIONAL_ENRICHMENT_CFG.get("heatmap", {})

    FE_QUERY_FAMILIES_CFG = FE_QUERY_CFG.get("query_families", {})
    FE_ENABLED_QUERY_FAMILIES = [
        family for family, enabled in FE_QUERY_FAMILIES_CFG.items()
        if bool(enabled)
    ]
    if not FE_ENABLED_QUERY_FAMILIES:
        raise ValueError("functional_enrichment.query_construction.query_families enables no query families")

    FE_GPROFILER_SOURCES = FE_GPROFILER_CFG.get("sources", [])
    if not FE_GPROFILER_SOURCES:
        raise ValueError("functional_enrichment.gprofiler.sources must be non-empty")

    FE_DISPLAY_MAPPING = FE_HEATMAP_CFG.get("display_group_mapping", "")
    FE_EXCLUDED_TERM_MAPPING = FE_HEATMAP_CFG.get("excluded_term_mapping", "")
    if not FE_DISPLAY_MAPPING:
        raise ValueError("functional_enrichment.heatmap.display_group_mapping is required")
    if not FE_EXCLUDED_TERM_MAPPING:
        raise ValueError("functional_enrichment.heatmap.excluded_term_mapping is required")

    def _functional_enrichment_contrast_manifest_args():
        args = []
        for prof in PAN_PROFILES:
            args.append(
                "--contrast-manifest "
                + shlex.quote(f"{prof}={profile_contrast_marker_manifest_abs(prof)}")
            )
        return " ".join(args)

    rule construct_graph_derived_enrichment_queries:
        input:
            features = PAN_FEATURES_TSV,
            selected_feature_genes = PAN_FEATURES_CLEAN,
            eligible_background = os.path.join(
                PAN_FEATURES_OUTDIR,
                "pan_cancer_feature_selection_eligible_gene_background.tsv",
            ),
            marker_readiness = [
                profile_marker_readiness_rel(p)
                for p in PAN_PROFILES
            ],
            script = os.path.join(SCRIPTS_DIR, "construct_graph_derived_enrichment_queries.py")
        output:
            manifest = FUNCTIONAL_ENRICHMENT_QUERY_MANIFEST,
            query_summary = os.path.join(FUNCTIONAL_ENRICHMENT_QUERY_DIR_REL, "query_summary.tsv"),
            skipped_queries = os.path.join(FUNCTIONAL_ENRICHMENT_QUERY_DIR_REL, "skipped_queries.tsv"),
            integrity = FUNCTIONAL_ENRICHMENT_QUERY_INTEGRITY,
            recurrence = os.path.join(FUNCTIONAL_ENRICHMENT_QUERY_DIR_REL, "enrichment_contrast_marker_recurrence.tsv"),
            done = touch(os.path.join(FUNCTIONAL_ENRICHMENT_QUERY_DIR_REL, ".done"))
        params:
            pipeline_root = PIPE_ROOT,
            outdir = FUNCTIONAL_ENRICHMENT_QUERY_DIR_ABS,
            minimum_query_gene_count = FE_QUERY_CFG.get("minimum_query_gene_count", 5),
            enabled_query_families = ",".join(FE_ENABLED_QUERY_FAMILIES),
            recurrence_minimum = 2,
            contrast_manifest_args = _functional_enrichment_contrast_manifest_args()
        log:
            os.path.join(LOGROOT, "construct_graph_derived_enrichment_queries.log")
        conda:
            CONDA_ENV_PY
        shell:
            r'''
            set -euo pipefail
            python "{input.script}" \
              --pipeline-root "{params.pipeline_root}" \
              --features "{input.features}" \
              --selected-feature-genes "{input.selected_feature_genes}" \
              --eligible-gene-background "{input.eligible_background}" \
              {params.contrast_manifest_args} \
              --outdir "{params.outdir}" \
              --minimum-query-gene-count {params.minimum_query_gene_count} \
              --enabled-query-families "{params.enabled_query_families}" \
              --recurrence-minimum-retained-contrast-count {params.recurrence_minimum} \
              > "{log}" 2>&1
            test -s "{output.manifest}" || (echo "ERROR: missing {output.manifest}" >&2; exit 1)
            test -s "{output.integrity}" || (echo "ERROR: missing {output.integrity}" >&2; exit 1)
            '''

    rule run_gprofiler_query_batch:
        input:
            manifest = FUNCTIONAL_ENRICHMENT_QUERY_MANIFEST,
            query_done = os.path.join(FUNCTIONAL_ENRICHMENT_QUERY_DIR_REL, ".done"),
            script = os.path.join(SCRIPTS_DIR, "run_gprofiler_from_manifest.R")
        output:
            success = touch(FUNCTIONAL_ENRICHMENT_GPROFILER_BATCH_SUCCESS)
        params:
            pipeline_root = PIPE_ROOT,
            outdir = FUNCTIONAL_ENRICHMENT_GPROFILER_DIR_ABS,
            organism = FE_GPROFILER_CFG.get("organism", "hsapiens"),
            sources = ",".join(FE_GPROFILER_SOURCES),
            correction_method = FE_GPROFILER_CFG.get("correction_method", "g_SCS"),
            alpha = FE_GPROFILER_CFG.get("alpha", 0.05),
            endpoint = FE_GPROFILER_CFG.get("endpoint", "https://biit.cs.ut.ee/gprofiler/api"),
            primary_exclude_iea = "TRUE" if FE_GPROFILER_CFG.get("primary_exclude_iea", False) else "FALSE",
            run_iea_sensitivity = "TRUE" if FE_GPROFILER_CFG.get("run_iea_sensitivity", True) else "FALSE",
            reuse_matching_query_results = "TRUE" if FE_GPROFILER_CFG.get("reuse_matching_query_results", False) else "FALSE",
            fail_on_query_error = "TRUE" if FE_GPROFILER_CFG.get("fail_on_query_error", True) else "FALSE",
            batch_start = 1,
            batch_end = 999999
        log:
            os.path.join(LOGROOT, "run_gprofiler_query_batch.log")
        conda:
            CONDA_ENV_R
        shell:
            r'''
            set -euo pipefail
            Rscript "{input.script}" \
              --pipeline-root "{params.pipeline_root}" \
              --query-manifest "{input.manifest}" \
              --outdir "{params.outdir}" \
              --organism "{params.organism}" \
              --sources "{params.sources}" \
              --alpha {params.alpha} \
              --correction-method "{params.correction_method}" \
              --endpoint "{params.endpoint}" \
              --primary-exclude-iea {params.primary_exclude_iea} \
              --run-iea-sensitivity {params.run_iea_sensitivity} \
              --reuse-matching-query-results {params.reuse_matching_query_results} \
              --fail-on-query-error {params.fail_on_query_error} \
              --batch-start {params.batch_start} \
              --batch-end {params.batch_end} \
              > "{log}" 2>&1
            test -s "{output.success}" || (echo "ERROR: missing {output.success}" >&2; exit 1)
            '''

    rule aggregate_gprofiler_results:
        input:
            manifest = FUNCTIONAL_ENRICHMENT_QUERY_MANIFEST,
            query_success = FUNCTIONAL_ENRICHMENT_GPROFILER_BATCH_SUCCESS,
            script = os.path.join(SCRIPTS_DIR, "aggregate_gprofiler_results.R")
        output:
            all_terms = os.path.join(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, "gprofiler_all_term_results.tsv"),
            significant_terms = os.path.join(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, "gprofiler_significant_term_results.tsv"),
            run_manifest = os.path.join(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, "gprofiler_query_run_manifest.tsv"),
            iea_sensitivity = os.path.join(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, "gprofiler_iea_sensitivity.tsv"),
            run_summary = os.path.join(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, "gprofiler_run_summary.tsv"),
            success = touch(FUNCTIONAL_ENRICHMENT_AGGREGATE_SUCCESS)
        params:
            pipeline_root = PIPE_ROOT,
            gprofiler_dir = FUNCTIONAL_ENRICHMENT_GPROFILER_DIR_ABS,
            outdir = FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_ABS,
            alpha = FE_GPROFILER_CFG.get("alpha", 0.05),
            run_iea_sensitivity = "TRUE" if FE_GPROFILER_CFG.get("run_iea_sensitivity", True) else "FALSE"
        log:
            os.path.join(LOGROOT, "aggregate_gprofiler_results.log")
        conda:
            CONDA_ENV_R
        shell:
            r'''
            set -euo pipefail
            Rscript "{input.script}" \
              --pipeline-root "{params.pipeline_root}" \
              --query-manifest "{input.manifest}" \
              --gprofiler-dir "{params.gprofiler_dir}" \
              --outdir "{params.outdir}" \
              --alpha {params.alpha} \
              --run-iea-sensitivity {params.run_iea_sensitivity} \
              > "{log}" 2>&1
            test -s "{output.success}" || (echo "ERROR: missing {output.success}" >&2; exit 1)
            '''

    rule construct_enrichment_reporting_terms:
        input:
            significant_terms = os.path.join(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, "gprofiler_significant_term_results.tsv"),
            manifest = FUNCTIONAL_ENRICHMENT_QUERY_MANIFEST,
            script = os.path.join(SCRIPTS_DIR, "build_enrichment_summary_top_terms.R")
        output:
            reporting_terms = os.path.join(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, "gprofiler_top_terms_reporting.tsv"),
            provenance = os.path.join(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, "gprofiler_top_terms_reporting_provenance.tsv")
        params:
            top_terms_per_query = FE_HEATMAP_CFG.get("top_terms_per_query", 10)
        log:
            os.path.join(LOGROOT, "construct_enrichment_reporting_terms.log")
        conda:
            CONDA_ENV_R
        shell:
            r'''
            set -euo pipefail
            Rscript "{input.script}" \
              --input "{input.significant_terms}" \
              --query-manifest "{input.manifest}" \
              --output "{output.reporting_terms}" \
              --top-terms-per-query {params.top_terms_per_query} \
              > "{log}" 2>&1
            test -s "{output.reporting_terms}" || (echo "ERROR: missing {output.reporting_terms}" >&2; exit 1)
            test -s "{output.provenance}" || (echo "ERROR: missing {output.provenance}" >&2; exit 1)
            '''

    rule summarize_contrast_term_support:
        input:
            significant_terms = os.path.join(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, "gprofiler_significant_term_results.tsv"),
            manifest = FUNCTIONAL_ENRICHMENT_QUERY_MANIFEST,
            script = os.path.join(SCRIPTS_DIR, "summarize_enrichment_contrast_term_support.R")
        output:
            support = os.path.join(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, "enrichment_contrast_term_support.tsv"),
            provenance = os.path.join(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, "enrichment_contrast_term_support_provenance.tsv")
        log:
            os.path.join(LOGROOT, "summarize_contrast_term_support.log")
        conda:
            CONDA_ENV_R
        shell:
            r'''
            set -euo pipefail
            Rscript "{input.script}" \
              --significant-terms "{input.significant_terms}" \
              --query-manifest "{input.manifest}" \
              --output "{output.support}" \
              > "{log}" 2>&1
            test -s "{output.support}" || (echo "ERROR: missing {output.support}" >&2; exit 1)
            test -s "{output.provenance}" || (echo "ERROR: missing {output.provenance}" >&2; exit 1)
            '''

    rule plot_primary_enrichment_heatmap:
        input:
            significant_terms = os.path.join(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, "gprofiler_significant_term_results.tsv"),
            manifest = FUNCTIONAL_ENRICHMENT_QUERY_MANIFEST,
            display_mapping = FE_DISPLAY_MAPPING,
            excluded_term_mapping = FE_EXCLUDED_TERM_MAPPING,
            script = os.path.join(SCRIPTS_DIR, "plot_enrichment_top_terms_heatmap.R")
        output:
            pdf = FUNCTIONAL_ENRICHMENT_PRIMARY_HEATMAP_PDF,
            png = os.path.join(FUNCTIONAL_ENRICHMENT_FIGURE_DIR_REL, "Fig_enrichment_top_terms_heatmap.png"),
            matrix = os.path.join(FUNCTIONAL_ENRICHMENT_FIGURE_DIR_REL, "Fig_enrichment_top_terms_heatmap_matrix.tsv"),
            selected = os.path.join(FUNCTIONAL_ENRICHMENT_FIGURE_DIR_REL, "Fig_enrichment_top_terms_heatmap_selected_terms.tsv"),
            excluded = os.path.join(FUNCTIONAL_ENRICHMENT_FIGURE_DIR_REL, "Fig_enrichment_top_terms_heatmap_excluded_terms.tsv"),
            provenance = os.path.join(FUNCTIONAL_ENRICHMENT_FIGURE_DIR_REL, "Fig_enrichment_top_terms_heatmap_provenance.tsv")
        params:
            outdir = FUNCTIONAL_ENRICHMENT_FIGURE_DIR_ABS,
            top_terms_per_query = FE_HEATMAP_CFG.get("top_terms_per_query", 10),
            maximum_terms = FE_HEATMAP_CFG.get("maximum_terms", 30),
            score_cap = FE_HEATMAP_CFG.get("score_cap", 30),
            term_wrap_width = FE_HEATMAP_CFG.get("term_wrap_width", 50),
            figure_width = FE_HEATMAP_CFG.get("figure_width", 16),
            minimum_figure_height = FE_HEATMAP_CFG.get("minimum_figure_height", 7),
            row_height_increment = FE_HEATMAP_CFG.get("row_height_increment", 0.34),
            correction_method = FE_GPROFILER_CFG.get("correction_method", "g_SCS")
        log:
            os.path.join(LOGROOT, "plot_primary_enrichment_heatmap.log")
        conda:
            CONDA_ENV_R
        shell:
            r'''
            set -euo pipefail
            Rscript "{input.script}" \
              --input "{input.significant_terms}" \
              --query-manifest "{input.manifest}" \
              --display-mapping "{input.display_mapping}" \
              --excluded-term-mapping "{input.excluded_term_mapping}" \
              --outdir "{params.outdir}" \
              --top-terms-per-query {params.top_terms_per_query} \
              --maximum-terms {params.maximum_terms} \
              --score-cap {params.score_cap} \
              --term-wrap-width {params.term_wrap_width} \
              --figure-width {params.figure_width} \
              --minimum-figure-height {params.minimum_figure_height} \
              --row-height-increment {params.row_height_increment} \
              --correction-method "{params.correction_method}" \
              > "{log}" 2>&1
            test -s "{output.pdf}" || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
            test -s "{output.provenance}" || (echo "ERROR: missing {output.provenance}" >&2; exit 1)
            '''

    rule plot_contrast_term_support_heatmap:
        input:
            support = os.path.join(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, "enrichment_contrast_term_support.tsv"),
            script = os.path.join(SCRIPTS_DIR, "plot_enrichment_contrast_term_support_heatmap.R")
        output:
            pdf = FUNCTIONAL_ENRICHMENT_CONTRAST_SUPPORT_HEATMAP_PDF,
            png = os.path.join(FUNCTIONAL_ENRICHMENT_FIGURE_DIR_REL, "Fig_enrichment_contrast_term_support_heatmap.png"),
            matrix = os.path.join(FUNCTIONAL_ENRICHMENT_FIGURE_DIR_REL, "Fig_enrichment_contrast_term_support_heatmap_matrix.tsv"),
            selected = os.path.join(FUNCTIONAL_ENRICHMENT_FIGURE_DIR_REL, "Fig_enrichment_contrast_term_support_heatmap_selected_terms.tsv"),
            provenance = os.path.join(FUNCTIONAL_ENRICHMENT_FIGURE_DIR_REL, "Fig_enrichment_contrast_term_support_heatmap_provenance.tsv")
        params:
            outdir = FUNCTIONAL_ENRICHMENT_FIGURE_DIR_ABS,
            maximum_terms = FE_HEATMAP_CFG.get("maximum_terms", 30),
            term_wrap_width = FE_HEATMAP_CFG.get("term_wrap_width", 50),
            figure_width = FE_HEATMAP_CFG.get("figure_width", 16),
            minimum_figure_height = FE_HEATMAP_CFG.get("minimum_figure_height", 7),
            row_height_increment = FE_HEATMAP_CFG.get("row_height_increment", 0.34)
        log:
            os.path.join(LOGROOT, "plot_contrast_term_support_heatmap.log")
        conda:
            CONDA_ENV_R
        shell:
            r'''
            set -euo pipefail
            Rscript "{input.script}" \
              --input "{input.support}" \
              --outdir "{params.outdir}" \
              --maximum-terms {params.maximum_terms} \
              --term-wrap-width {params.term_wrap_width} \
              --figure-width {params.figure_width} \
              --minimum-figure-height {params.minimum_figure_height} \
              --row-height-increment {params.row_height_increment} \
              > "{log}" 2>&1
            test -s "{output.pdf}" || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
            test -s "{output.provenance}" || (echo "ERROR: missing {output.provenance}" >&2; exit 1)
            '''
