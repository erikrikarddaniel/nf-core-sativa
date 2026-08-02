/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENSURE_ALIGNED - transparently align unaligned input via hmmalign

    Sequences that are already a multiple sequence alignment (every record the same
    length) pass through unchanged. Sequences that aren't are aligned against an HMM
    profile (hmm, optionally a specific named profile within it via hmm_name) using
    hmmalign, then converted back to plain FASTA -- no separate mode-switch param
    needed, detection is by content.

    Workflow:
      1. Detect aligned vs unaligned (equal sequence lengths?)     (CHECKALIGNED)
      2. unaligned only: extract the named profile, if given       (HMMER_HMMFETCH)
      3. unaligned only: align against the profile                 (HMMER_HMMALIGN)
      4. unaligned only: mask down to the profile's own match-state
         columns, discarding insert-state columns                  (HMMER_ESLALIMASK)
      5. unaligned only: convert Stockholm back to aligned FASTA   (HMMER_ESLREFORMAT)
      6. unaligned only: decompress (eslreformat always gzips)     (GUNZIP)

    Step 4 matters beyond tidiness: hmmalign's raw output width is match-state
    columns *plus* the union of every sequence's own insert-state columns, so a
    single divergent sequence's insertion pads every other sequence with extra
    gap columns there. Left unmasked, that badly deflates a downstream non-gap-
    proportion filter for genuinely full-length sequences too (confirmed
    empirically: masking took the median non-gap proportion across a real
    121-sequence archaeal 16S set from 0.61 to 0.99). Masking first makes the
    alignment width exactly the HMM's match-state count.

    The two branches are emitted SEPARATELY (alignment_passthrough,
    alignment_from_hmm), not merged into one, because the caller applies a
    different gap/coverage filter -- with its own threshold default -- to each:
    a directly-provided alignment (e.g. from MAFFT) is naturally wide with real
    cross-lineage indels even for full-length sequences, whereas this subworkflow's
    masked hmmalign output has a fixed-length, much tighter distribution (see
    GAPFILTER vs PROFILECOVER in workflows/sativa.nf).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHECKALIGNED      } from '../../../modules/local/checkaligned/main'
include { HMMER_HMMFETCH    } from '../../../modules/nf-core/hmmer/hmmfetch/main'
include { HMMER_HMMALIGN    } from '../../../modules/nf-core/hmmer/hmmalign/main'
include { HMMER_ESLALIMASK  } from '../../../modules/nf-core/hmmer/eslalimask/main'
include { HMMER_ESLREFORMAT } from '../../../modules/nf-core/hmmer/eslreformat/main'
include { GUNZIP            } from '../../../modules/nf-core/gunzip/main'

workflow ENSURE_ALIGNED {

    take:
    ch_sequences // channel: sequences file, aligned or not, already FASTA (format-normalised by the caller)
    hmm          // value:   path to an HMM profile database, or null/empty if not needed
    hmm_name     // value:   name of a specific profile within hmm, or null/empty

    main:
    CHECKALIGNED(ch_sequences.map { [ [ id: 'user-alignment' ], it ] })

    // Fail clearly, but only if there's actually something unaligned to align --
    // aligned-input users should never be forced to supply --hmm. Deferred into this
    // .map() (rather than checked eagerly above) so it only fires if an unaligned item
    // genuinely flows through; error() halts the pipeline with the message.
    def ch_unaligned = CHECKALIGNED.out.unaligned
        .map { meta, fasta ->
            if (!hmm) {
                error("Unaligned input detected, but --hmm was not provided. Supply " +
                    "--hmm (path to an HMM profile database) to align it -- and " +
                    "--hmm_name too, if that database holds more than one profile.")
            }
            [ meta, fasta ]
        }

    // hmm/hmm_name are plain values, known before execution starts, so it's safe to
    // branch on them here at compose time rather than inside a channel operator.
    def ch_hmm_profile
    if (hmm) {
        def ch_hmm_db = channel.fromPath(hmm, checkIfExists: true)
            .map { [ [ id: 'hmm' ], it ] }
        if (hmm_name) {
            // hmmfetch works directly against a raw (unindexed) multi-profile
            // database when given an explicit key -- no separate --index step needed.
            HMMER_HMMFETCH(ch_hmm_db, hmm_name, [], [])
            ch_hmm_profile = HMMER_HMMFETCH.out.hmm
        } else {
            ch_hmm_profile = ch_hmm_db
        }
    } else {
        ch_hmm_profile = channel.empty()
    }

    HMMER_HMMALIGN(ch_unaligned, ch_hmm_profile.map { _meta, hmm_file -> hmm_file })

    // --rf-is-mask derives the mask from the alignment's own #=GC RF (reference/
    // match-state) annotation -- no separate maskfile needed, unlike the tool's
    // other usage form. The six false values are the module's optional
    // mask-report-file toggles (fmask/gmask/pmask), none of which are needed here.
    HMMER_ESLALIMASK(
        HMMER_HMMALIGN.out.sto.map { meta, sto -> [ meta, sto, false, false, false, false, false, false ] },
        []
    )

    HMMER_ESLREFORMAT(HMMER_ESLALIMASK.out.maskedaln, '')

    // eslreformat's output is unconditionally gzipped (hardcoded in the vendored
    // module's script); decompress so it matches the plain-text passthrough branch
    // below before merging the two back into one channel.
    GUNZIP(HMMER_ESLREFORMAT.out.seqreformated)

    emit:
    alignment_passthrough = CHECKALIGNED.out.aligned.map { _meta, fasta -> fasta } // channel: alignment file, already aligned on input
    alignment_from_hmm    = GUNZIP.out.gunzip.map { _meta, fasta -> fasta }        // channel: alignment file, aligned via hmmalign then masked to match-state columns
}
