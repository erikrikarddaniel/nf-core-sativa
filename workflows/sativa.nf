/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_sativa_pipeline'
include { CHECKNAMECONSISTENCY   } from '../modules/local/checknameconsistency/main'
include { EMBOSS_SEQRET          } from '../modules/nf-core/emboss/seqret/main'
include { ENSURE_ALIGNED         } from '../subworkflows/local/ensure_aligned'
include { GAPFILTER              } from '../modules/local/gapfilter/main'
include { RAXTAX_PREFILTER       } from '../subworkflows/local/raxtax_prefilter'
include { SATIVA as SWF_SATIVA   } from '../subworkflows/local/sativa'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SATIVA {

    take:
    ch_taxonomy    // channel: taxonomy file
    ch_alignment   // channel: alignment file
    skip_raxtax    // value:   skip the raxtax prefilter?
    skip_gapfilter // value:   skip the gap filter?
    hmm            // value:   path to an HMM profile database, or null/empty if not needed
    hmm_name       // value:   name of a specific profile within hmm, or null/empty
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'sativa_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: Validate that taxonomy and alignment name the same sequences, and
    // rewrite characters that are difficult for downstream tools (e.g. parens) in
    // both. Runs first, as a process (not inline Nextflow code) so a large input
    // doesn't inflate the head job's memory/CPU footprint.
    //
    CHECKNAMECONSISTENCY(
        ch_taxonomy.combine(ch_alignment).map { tax, aln -> [ [ id: 'user-alignment' ], tax, aln ] }
    )
    def ch_taxonomy_checked  = CHECKNAMECONSISTENCY.out.checked.map { _meta, tax, _aln -> tax }
    def ch_alignment_checked = CHECKNAMECONSISTENCY.out.checked.map { _meta, _tax, aln -> aln }

    //
    // MODULE: Normalise the alignment to FASTA once, here, rather than separately
    // inside RAXTAX_PREFILTER and SWF_SATIVA (previously duplicated). Also gives
    // ENSURE_ALIGNED a single canonical format to inspect for the unaligned-input
    // support below.
    //
    EMBOSS_SEQRET(ch_alignment_checked.map { [ [ id: 'user-alignment' ], it ] }, 'fasta')
    def ch_alignment_fasta = EMBOSS_SEQRET.out.outseq.map { _meta, aln -> aln }

    //
    // SUBWORKFLOW: ENSURE_ALIGNED
    //
    // Transparently accepts unaligned input too, with no separate mode-switch param:
    // already-aligned content passes straight through; unaligned content is aligned
    // via hmmalign against the hmm/hmm_name profile before continuing.
    //
    ENSURE_ALIGNED(ch_alignment_fasta, hmm, hmm_name)
    def ch_alignment_aligned = ENSURE_ALIGNED.out.alignment

    //
    // MODULE: GAPFILTER (optional, skip_gapfilter to disable)
    //
    // Drops sequences too short/gappy to place reliably (below params.min_nongap
    // non-gap proportion in the -- possibly hmmalign-realigned -- alignment),
    // reporting them separately rather than silently discarding them.
    //
    def ch_taxonomy_for_raxtax
    def ch_alignment_for_raxtax
    // Coerce explicitly: a CLI-supplied `--skip_gapfilter false` arrives as the
    // *string* "false" -- see the analogous skip_raxtax coercion below for why
    // .toString().toBoolean() is needed even with nf-schema's cli_typecast enabled.
    def run_gapfilter = !skip_gapfilter.toString().toBoolean()
    if (run_gapfilter) {
        GAPFILTER(
            ch_taxonomy_checked.combine(ch_alignment_aligned).map { tax, aln -> [ [ id: 'user-alignment' ], tax, aln ] }
        )
        ch_taxonomy_for_raxtax  = GAPFILTER.out.taxonomy.map { _meta, tax -> tax }
        ch_alignment_for_raxtax = GAPFILTER.out.alignment.map { _meta, aln -> aln }
    } else {
        ch_taxonomy_for_raxtax  = ch_taxonomy_checked
        ch_alignment_for_raxtax = ch_alignment_aligned
    }

    //
    // SUBWORKFLOW: RAXTAX_PREFILTER (optional, skip_raxtax to disable)
    //
    // Fast raxtax self-classification prefilter ahead of the expensive EPA-ng-based
    // placement below. Sequences it flags never reach SWF_SATIVA -- they're reported
    // directly via ch_raxtax_mislabels instead.
    //
    def ch_taxonomy_for_sativa
    def ch_alignment_for_sativa
    def ch_raxtax_mislabels
    // Coerce explicitly: a CLI-supplied `--skip_raxtax false` arrives as the *string*
    // "false", and Groovy's `!"false"` is false (any non-empty string is truthy) --
    // .toBoolean() parses both real Booleans and "true"/"false" strings correctly.
    // Confirmed empirically that nf-schema's cli_typecast (enabled just above, in
    // PIPELINE_INITIALISATION) validates the string against the boolean schema type but
    // does not itself replace params.skip_raxtax with a real Boolean, so this is still
    // needed even with cli_typecast on.
    def run_raxtax = !skip_raxtax.toString().toBoolean()
    if (run_raxtax) {
        RAXTAX_PREFILTER(ch_taxonomy_for_raxtax, ch_alignment_for_raxtax)
        ch_taxonomy_for_sativa  = RAXTAX_PREFILTER.out.taxonomy
        ch_alignment_for_sativa = RAXTAX_PREFILTER.out.alignment
        ch_raxtax_mislabels     = RAXTAX_PREFILTER.out.mislabels
    } else {
        ch_taxonomy_for_sativa  = ch_taxonomy_for_raxtax
        ch_alignment_for_sativa = ch_alignment_for_raxtax
        ch_raxtax_mislabels     = channel.empty()
    }

    //
    // SUBWORKFLOW: SATIVA
    //
    // This implements all the logic in the workflow.
    //
    // The later two params are meant to pass a reference tree and a model file respectively. Not implemented yet.
    //
    SWF_SATIVA(ch_taxonomy_for_sativa, ch_alignment_for_sativa, [], [])

    //
    // Merge raxtax-flagged mislabels (skipped placement entirely) with SATIVASCORE's own
    // into one final report. Both share the same TSV schema (method column distinguishes
    // detection source), so collectFile with keepHeader can concatenate them directly --
    // no bridging process needed just to reshape/combine two files.
    //
    ch_raxtax_mislabels
        .mix(SWF_SATIVA.out.mislabels)
        .map { _meta, tsv -> tsv }
        .collectFile(name: 'user-alignment.mislabels.tsv', storeDir: "${outdir}/mislabels", keepHeader: true, skip: 1)

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'sativa'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
