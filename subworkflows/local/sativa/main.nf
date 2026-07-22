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
      nf-core modules install emboss/seqret
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TAXONOMYTREE   } from '../../../modules/local/taxonomytree/main'
include { EMBOSS_SEQRET  } from '../../../modules/nf-core/emboss/seqret/main'
include { IQTREE         } from '../../../modules/nf-core/iqtree/main'
include { SATIVALOOSPLIT } from '../../../modules/local/sativaloosplit/main'
include { SATIVASCORE    } from '../../../modules/local/sativascore/main'
include { EPANG_PLACE    } from '../../../modules/nf-core/epang/place/main'

/**
process CHECKNAMECONSISTENCY {
    label 'process_low'

    conda "conda-forge::python=3.11 bioconda::biopython=1.84"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'biocontainers/biopython:1.84' }"

    input:
    path taxonomy
    path alignment

    output:

    script:

    """
    """

    stub:

    """
    """
}
**/

// ─── Subworkflow ──────────────────────────────────────────────────────────────

workflow SATIVA {

    take:
    ch_taxonomy   // channel: [ val(meta), path(taxonomy.tsv) ]
                  //   Tab-separated: seq_name <TAB> Kingdom;Phylum;Class;...
                  //   The taxonomic code (BAC/BOT/ZOO/VIR) is the first token.

    ch_alignment  // channel: [ val(meta), path(alignment) ]
                  //   Aligned, labeled sequences. FASTA, Clustal or PHYLIP; format
                  //   is auto-detected and normalised to PHYLIP by EMBOSS_SEQRET
                  //   below, since SATIVALOOSPLIT (our own Python code) only trusts
                  //   PHYLIP rather than trying to parse every format itself.
                  //   Sequence IDs must match the first column of ch_taxonomy.

    ch_ref_tree   // channel: [ val(meta), path(tree.nwk) ]
                  //   Pre-built reference tree. Pass Channel.empty() to build one.

    ch_ref_model  // channel: [ val(meta), path(model.txt) ]
                  //   RAxML-NG model file matching ch_ref_tree. Channel.empty() if none.

    main:
//    def ch_versions = channel.empty()

    //
    // Check that ch_taxonomy and ch_alignment have the same set of unique names
    //
    //CHECKNAMECONSISTENCY(ch_taxonomy, ch_alignment)

    // Normalise the input alignment to PHYLIP regardless of whether it arrived as
    // FASTA, Clustal or PHYLIP. EMBOSS auto-detects the input format from content,
    // so no format-sniffing of our own is needed here. IQTREE (a compiled tool) can
    // likely handle any of the three directly, but SATIVALOOSPLIT's hand-rolled
    // Python parser only supports PHYLIP, so every consumer downstream is kept on
    // one guaranteed format instead.
    EMBOSS_SEQRET(ch_alignment.map { [ [ id: 'user-alignment' ], it ] }, 'phylip')
    def ch_alignment_phylip = EMBOSS_SEQRET.out.outseq

    // ── Phase 1: Reference tree construction (epa_trainer) ────────────────────
    //
    // Build a multifurcating guide tree from taxonomy strings, then run RAxML-NG
    // with that tree as a topology constraint.  Multiple independent searches are
    // controlled via ext.args (e.g. "--searches 10").  The best-scoring tree is
    // then model-optimised.  The resulting tree + model are reusable across runs
    // (pass via ch_ref_tree / ch_ref_model to skip this phase).

    TAXONOMYTREE(ch_taxonomy.map { it -> [ [ id: 'guide-tree' ], it ] })

    IQTREE(
        ch_alignment_phylip.map { meta, aln -> [ meta, aln, [] ] },         // Alignment
        [],                                                                 // tree_te
        [],                                                                 // lmclust
        [],                                                                 // mdef
        [],                                                                 // partitions_equal
        [],                                                                 // partitions_proportional
        [],                                                                 // partitions_unlinked
        [],                                                                 // guide_tree
        [],                                                                 // sitefreq_in
        TAXONOMYTREE.out.guide_tree.map { _meta, tree -> tree },            // constraint_tree
        [],                                                                 // trees_z
        [],                                                                 // suptree
        []                                                                  // trees_rf
    )

    // TODO: give externally supplied ch_ref_tree / ch_ref_model precedence over what
    // we just built, once main.nf actually exposes a way to pass them in (currently
    // always called with `[]`, so mixing them in here would inject a spurious
    // empty-list item into the channel).
    def ch_tree = IQTREE.out.phylogeny

    // IQTREE doesn't emit a separate model file; the chosen substitution model is
    // only reported in its run log (e.g. "Best-fit model: GTR+F+I chosen according
    // to BIC"), so parse it out of there instead. Absent under -stub-run, where the
    // log is just an empty touched file.
    def ch_model = IQTREE.out.log.map { meta, log ->
        def matcher = log.text =~ /Best-fit model: (.*) chosen according to/
        [ meta, matcher.find() ? matcher.group(1) : null ]
    }

//    // ── Phase 2: HMM profile ──────────────────────────────────────────────────
//    //
//    // EPA-ng uses an HMM profile built from the reference MSA to re-align each
//    // LOO query sequence before placement.  Corresponds to the hmmbuild call in
//    // epa_trainer.py.  Skipped for now: pipeline input is already an aligned MSA,
//    // so EPA-ng can be run directly (--query is pre-aligned) without a profile.
//    // Revisit if/when unaligned input becomes a supported entry point.
//
//    EPANG_HMMBUILD(ch_alignment)
//    ch_versions = ch_versions.mix(EPANG_HMMBUILD.out.versions)
//
    // ── Phase 3: Leave-one-out scatter (epa_classifier) ────────────────────────
    //
    // Split the full MSA + reference tree into N independent (query, reference,
    // tree) triples, one per held-out sequence.  This is the computationally
    // dominant phase; the fan-out means Nextflow schedules up to N EPA-ng jobs
    // simultaneously once wired to EPANG_PLACE below.

    SATIVALOOSPLIT(ch_alignment_phylip.join(ch_tree))

    // SATIVALOOSPLIT emits one tuple per input dataset, with the query/reference/tree
    // outputs as same-length file lists (one entry per held-out sequence).  .transpose()
    // unpacks that into one channel item per sequence, but every unpacked item still
    // carries the same meta as the parent call — so a running counter (embedded by the
    // module in each file's basename) is folded into meta.id here to keep the N items
    // distinct downstream (e.g. for EPANG_PLACE and later grouping/joins).
    def ch_loo = SATIVALOOSPLIT.out.loo
        .transpose()
        .map { meta, queryaln, referencealn, referencetree ->
            def counter = queryaln.baseName.tokenize('_')[0]
            [ meta + [ id: "${meta.id}_${counter}" ], queryaln, referencealn, referencetree ]
        }

    // epa-ng refuses to run without an explicit --model (see ext.args in
    // conf/modules.config); fold the IQTREE-derived model string into each split's
    // meta so the config closure can read it. ch_model holds a single item per
    // input dataset, so .combine() broadcasts it across all N ch_loo items.
    // NB: combine directly on the [meta, model] tuples rather than unwrapping model
    // into its own .map() first — under -stub-run (or any run where the regex finds
    // no match) model is null, and a bare null returned from .map() is silently
    // dropped by Nextflow, which would empty out this channel entirely.
    def ch_loo_with_model = ch_loo
        .combine(ch_model)
        .map { meta, queryaln, referencealn, referencetree, _model_meta, model ->
            [ meta + [ model: model ], queryaln, referencealn, referencetree ]
        }

    // Place each held-out sequence back into its pruned reference tree.  No HMM
    // profile needed: query/reference alignments are both subsets of the same
    // input MSA, so they already share the same column coordinate space.
    EPANG_PLACE(ch_loo_with_model, [], [])

    // ── Phase 4: Gather and score (mislabels_handler) ──────────────────────────
    //
    // Collect all N per-sequence jplace files back into one item per input dataset,
    // then compare each EPA classification to the original taxonomy label.

    def ch_score_input = EPANG_PLACE.out.jplace
        // meta carries both the per-sequence counter suffix added after SATIVALOOSPLIT's
        // transpose (e.g. "user-alignment_0001") and the 'model' key folded in above for
        // epa-ng's ext.args; rebuild a bare [id:...] meta (strip the counter and drop
        // 'model') rather than merging, so it matches ch_taxonomy's meta below exactly
        // and .join() doesn't silently match nothing.
        .map { meta, jplace -> [ [ id: meta.id.tokenize('_')[0..-2].join('_') ], jplace ] }
        .groupTuple()
        // ch_taxonomy is a bare file channel (see take: above); give it the same
        // 'user-alignment' meta used elsewhere so it lines up with ch_score_input.
        .join(ch_taxonomy.map { [ [ id: 'user-alignment' ], it ] })

    SATIVASCORE(ch_score_input)

    emit:
    mislabels = SATIVASCORE.out.mislabels   // [ meta, tsv ]  putative mislabels, ranked
    summary   = SATIVASCORE.out.summary     // [ meta, txt ]  run statistics
    tree      = ch_tree                     // [ meta, nwk ]  reference tree (cache for reuse)
    model     = ch_model                    // [ meta, txt ]  IQTREE model  (cache for reuse)
    loo       = ch_loo                      // [ meta, queryaln, referencealn, referencetree ]  per-sequence LOO triples
    jplace    = EPANG_PLACE.out.jplace      // [ meta, jplace.gz ]  per-sequence EPA-ng placement result
}
