/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SATIVA - Semi-Automatic Taxonomy Improvement and Validation Algorithm
    Reverse-engineered from https://github.com/amkozlov/sativa

    Workflow:
      1. Build a taxonomy-constrained ML reference tree   (epa_trainer.py)
      2. Classify every sequence via leave-one-out EPA    (epa_classifier.py)
      3. Score placements, report mismatches              (mislabels_handler.py)

main.nf
  └── PIPELINE_INITIALISATION   (subworkflows/local/utils_nfcore_sativa_pipeline/main.nf)
        validates params, parses samplesheet → ch_samplesheet channel
  └── NFCORE_SATIVA
        └── SATIVA               (workflows/sativa.nf)  ← main logic lives here
              ├── FASTQC          (modules/nf-core/fastqc/)
              └── MULTIQC         (modules/nf-core/multiqc/)
  └── PIPELINE_COMPLETION        (subworkflows/local/utils_nfcore_sativa_pipeline/main.nf)
        sends email / completion summary

    Required nf-core modules (install before use):
      nf-core modules install raxmlng/search
      nf-core modules install raxmlng/evaluate
      nf-core modules install epang/hmmbuild
      nf-core modules install epang/place
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TAXONOMYTREE } from '../../../modules/local/taxonomytree/main'
//include { RAXMLNG_SEARCH   } from '../../../modules/nf-core/raxmlng/search/main'
//include { RAXMLNG_EVALUATE } from '../../../modules/nf-core/raxmlng/evaluate/main'
//include { EPANG_HMMBUILD   } from '../../../modules/nf-core/epang/hmmbuild/main'
include { EPANG_PLACE         } from '../../../modules/nf-core/epang/place/main'

// Decompose the labeled MSA into N (query, reference) pairs, one per sequence.
// For each sequence i: query_i = seq_i alone; ref_i = all other N-1 sequences.
// File names use the sequence ID as the basename so downstream processes can
// correlate queries with references after .transpose() scatter.
// Corresponds to the per-sequence loop inside sativa.py::LeaveOneTest.run().
process SATIVA_LOO_SPLIT {
    label 'process_low'

    conda "conda-forge::python=3.11 bioconda::biopython=1.84"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("queries/*.fasta"),    emit: queries     // N files, one per seq
    tuple val(meta), path("references/*.fasta"), emit: references  // N files, complement sets
    path "versions.yml",                         emit: versions

    script:
    // TODO: implement bin/sativa_loo_split.py
    // For each record in ${alignment}:
    //   write record         → queries/<seq_id>.fasta
    //   write all others     → references/<seq_id>.fasta
    // Gaps-only columns must be stripped from each reference alignment.
    """
    mkdir -p queries references
    sativa_loo_split.py \\
        --alignment ${alignment} \\
        --queries   queries/ \\
        --refs      references/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p queries references
    echo ">stub" > queries/stub.fasta
    echo "A"    >> queries/stub.fasta
    echo ">stub" > references/stub.fasta
    echo "A"    >> references/stub.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}

// Parse the per-sequence jplace placements, score them against original taxonomy
// labels using likelihood-weighted voting over placement edges, and emit a TSV
// of putative mislabels ranked by confidence.
// Corresponds to mislabels_handler.py in the original Sativa.
//
// Scoring logic:
//   For each jplace file → identify the set of pendant/neighbouring leaves
//   → majority-vote the taxonomy at each rank → compare to original label
//   → flag as mislabel if they disagree and LW confidence > cutoff (-C in original)
process SATIVA_SCORE {
    label 'process_low'

    conda "conda-forge::python=3.11 bioconda::biopython=1.84"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'biocontainers/biopython:1.84' }"

    input:
    // All N jplace files for one job are staged under placements/
    tuple val(meta), path(jplace_files, stageAs: "placements/*"), path(taxonomy)

    output:
    tuple val(meta), path("*.mislabels.tsv"), emit: mislabels
    tuple val(meta), path("*.summary.txt"),   emit: summary
    path "versions.yml",                      emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ""  // e.g. "--confidence 0.5 --min-lwr 0.0"
    // TODO: implement bin/sativa_score.py
    // Output TSV columns:
    //   seq_name, original_label, predicted_label, confidence, lwr, mislabel_rank
    """
    sativa_score.py \\
        --placements placements/ \\
        --taxonomy   ${taxonomy} \\
        --output     ${prefix}.mislabels.tsv \\
        --summary    ${prefix}.summary.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mislabels.tsv ${prefix}.summary.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}

// ─── Subworkflow ──────────────────────────────────────────────────────────────

workflow SATIVA {

    take:
    ch_alignment  // channel: [ val(meta), path(alignment.fasta) ]
                  //   Aligned, labeled sequences (FASTA or PHYLIP).
                  //   Sequence IDs must match the first column of ch_taxonomy.

    ch_taxonomy   // channel: [ val(meta), path(taxonomy.tsv) ]
                  //   Tab-separated: seq_name <TAB> Kingdom;Phylum;Class;...
                  //   The taxonomic code (BAC/BOT/ZOO/VIR) is the first token.

    ch_ref_tree   // channel: [ val(meta), path(tree.nwk) ]
                  //   Pre-built reference tree. Pass Channel.empty() to build one.

    ch_ref_model  // channel: [ val(meta), path(model.txt) ]
                  //   RAxML-NG model file matching ch_ref_tree. Channel.empty() if none.

    main:
    def ch_versions = channel.empty()

    // ── Phase 1: Reference tree construction (epa_trainer) ────────────────────
    //
    // Build a multifurcating guide tree from taxonomy strings, then run RAxML-NG
    // with that tree as a topology constraint.  Multiple independent searches are
    // controlled via ext.args (e.g. "--searches 10").  The best-scoring tree is
    // then model-optimised.  The resulting tree + model are reusable across runs
    // (pass via ch_ref_tree / ch_ref_model to skip this phase).

    TAXONOMYTREE(ch_taxonomy)

//    // Combine alignment with guide tree; pass guide as topology constraint.
//    // TODO: RAXMLNG_SEARCH needs ext.args = "--tree-constraint <guide.nwk>"
//    //       until the nf-core module exposes a dedicated input for constraint files,
//    //       stage guide_tree alongside the alignment and reference it in ext.args.
//    def ch_search_input = ch_alignment
//        .join(TAXONOMYTREE.out.guide_tree)
//        .map { meta, fasta, guide -> [ meta, fasta, guide, [] ] }
//
//    RAXMLNG_SEARCH(ch_search_input)
//    ch_versions = ch_versions.mix(RAXMLNG_SEARCH.out.versions)
//
//    // Optimise model parameters on the winning tree (raxml-ng --evaluate)
//    def ch_eval_input = ch_alignment
//        .join(RAXMLNG_SEARCH.out.bestTree)
//        .map { meta, fasta, tree -> [ meta, fasta, tree, [] ] }
//
//    RAXMLNG_EVALUATE(ch_eval_input)
//    ch_versions = ch_versions.mix(RAXMLNG_EVALUATE.out.versions)
//
//    // Any externally supplied tree/model takes precedence over what we just built
//    def ch_tree  = ch_ref_tree .mix(RAXMLNG_EVALUATE.out.tree)
//    def ch_model = ch_ref_model.mix(RAXMLNG_EVALUATE.out.model)
//
//    // ── Phase 2: HMM profile ──────────────────────────────────────────────────
//    //
//    // EPA-ng uses an HMM profile built from the reference MSA to re-align each
//    // LOO query sequence before placement.  Corresponds to the hmmbuild call in
//    // epa_trainer.py.  Not needed if query sequences are already aligned.
//
//    EPANG_HMMBUILD(ch_alignment)
//    ch_versions = ch_versions.mix(EPANG_HMMBUILD.out.versions)
//
//    // ── Phase 3: Leave-one-out scatter (epa_classifier) ───────────────────────
//    //
//    // Split the full MSA into N independent (query, reference) pairs and classify
//    // each sequence in parallel.  This is the computationally dominant phase; the
//    // fan-out means Nextflow schedules up to N EPA-ng jobs simultaneously.
//
//    SATIVA_LOO_SPLIT(ch_alignment)
//    ch_versions = ch_versions.mix(SATIVA_LOO_SPLIT.out.versions)
//
//    // Scatter: flatten the list outputs into one channel item per sequence.
//    // Both channels use the sequence ID as the file basename so they remain
//    // paired after the transpose; we embed seq_id in meta to keep them aligned
//    // through the join and to label the per-sequence jplace output files.
//    def ch_queries = SATIVA_LOO_SPLIT.out.queries
//        .transpose()
//        .map { meta, q -> [ meta + [seq_id: q.baseName], q ] }
//
//    def ch_refs = SATIVA_LOO_SPLIT.out.references
//        .transpose()
//        .map { meta, r -> [ meta + [seq_id: r.baseName], r ] }
//
//    // Each EPA-ng call gets: query.fasta, ref.fasta, ref.nwk, ref.model, ref.hmm
//    // The tree/model/hmm inputs are keyed on the base meta (without seq_id),
//    // so we broadcast them via a cross-join on meta.id.
//    // TODO: verify EPANG_PLACE module input signature and adjust tuple structure
//    def ch_place_input = ch_queries
//        .join(ch_refs)
//        // Attach tree: strip seq_id for the join, then re-attach
//        .map { meta, q, r -> [ meta.subMap(meta.keySet() - ['seq_id']), meta.seq_id, q, r ] }
//        .join(ch_tree)
//        .join(ch_model)
//        .join(EPANG_HMMBUILD.out.hmm)
//        .map { meta, seq_id, q, r, tree, model, hmm ->
//            [ meta + [seq_id: seq_id], q, r, tree, model, hmm ]
//        }
//
//    EPANG_PLACE(ch_place_input)
//    ch_versions = ch_versions.mix(EPANG_PLACE.out.versions)
//
//    // ── Phase 4: Gather and score (mislabels_handler) ─────────────────────────
//    //
//    // Collect all N per-sequence jplace files back into one item per input dataset,
//    // then compare each EPA classification to the original taxonomy label.
//
//    def ch_score_input = EPANG_PLACE.out.jplace
//        // Strip seq_id before grouping so all N results land in one tuple
//        .map { meta, jplace -> [ meta.subMap(meta.keySet() - ['seq_id']), jplace ] }
//        .groupTuple()
//        .join(ch_taxonomy)
//
//    SATIVA_SCORE(ch_score_input)
//    ch_versions = ch_versions.mix(SATIVA_SCORE.out.versions)

    emit:
//    mislabels = SATIVA_SCORE.out.mislabels  // [ meta, tsv ]  putative mislabels, ranked
//    summary   = SATIVA_SCORE.out.summary    // [ meta, txt ]  run statistics
//    tree      = ch_tree                     // [ meta, nwk ]  reference tree (cache for reuse)
//    model     = ch_model                    // [ meta, txt ]  RAxML-NG model  (cache for reuse)
}
