/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENSURE_ALIGNED - transparently align unaligned input via hmmalign

    Sequences that are already a multiple sequence alignment (every record the same
    length) pass through unchanged. Sequences that aren't are aligned against an HMM
    profile (params.hmm, optionally a specific named profile within it via
    params.hmm_name) using hmmalign, then converted back to plain FASTA -- no separate
    mode-switch param needed, detection is by content.

    Workflow:
      1. Detect aligned vs unaligned (equal sequence lengths?)     (CHECKALIGNED)
      2. unaligned only: extract the named profile, if given       (HMMER_HMMFETCH)
      3. unaligned only: align against the profile                 (HMMER_HMMALIGN)
      4. unaligned only: convert Stockholm back to aligned FASTA   (HMMER_ESLREFORMAT)
      5. unaligned only: decompress (eslreformat always gzips)     (GUNZIP)
      6. Merge the two branches back into one `alignment` emit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHECKALIGNED      } from '../../../modules/local/checkaligned/main'
include { HMMER_HMMFETCH    } from '../../../modules/nf-core/hmmer/hmmfetch/main'
include { HMMER_HMMALIGN    } from '../../../modules/nf-core/hmmer/hmmalign/main'
include { HMMER_ESLREFORMAT } from '../../../modules/nf-core/hmmer/eslreformat/main'
include { GUNZIP            } from '../../../modules/nf-core/gunzip/main'

workflow ENSURE_ALIGNED {

    take:
    ch_alignment  // channel: alignment file, already FASTA (format-normalised by the caller)

    main:
    CHECKALIGNED(ch_alignment.map { [ [ id: 'user-alignment' ], it ] })

    // Fail clearly, but only if there's actually something unaligned to align --
    // aligned-input users should never be forced to supply --hmm. Deferred into this
    // .map() (rather than checked eagerly above) so it only fires if an unaligned item
    // genuinely flows through; error() halts the pipeline with the message.
    def ch_unaligned = CHECKALIGNED.out.unaligned
        .map { meta, fasta ->
            if (!params.hmm) {
                error("Unaligned input detected, but --hmm was not provided. Supply " +
                    "--hmm (path to an HMM profile database) to align it -- and " +
                    "--hmm_name too, if that database holds more than one profile.")
            }
            [ meta, fasta ]
        }

    // params.hmm/params.hmm_name are plain pipeline params, known before execution
    // starts, so it's safe to branch on them here at compose time rather than inside
    // a channel operator.
    def ch_hmm_profile
    if (params.hmm) {
        def ch_hmm_db = channel.fromPath(params.hmm, checkIfExists: true)
            .map { [ [ id: 'hmm' ], it ] }
        if (params.hmm_name) {
            // hmmfetch works directly against a raw (unindexed) multi-profile
            // database when given an explicit key -- no separate --index step needed.
            HMMER_HMMFETCH(ch_hmm_db, params.hmm_name, [], [])
            ch_hmm_profile = HMMER_HMMFETCH.out.hmm
        } else {
            ch_hmm_profile = ch_hmm_db
        }
    } else {
        ch_hmm_profile = channel.empty()
    }

    HMMER_HMMALIGN(ch_unaligned, ch_hmm_profile.map { _meta, hmm -> hmm })

    HMMER_ESLREFORMAT(HMMER_HMMALIGN.out.sto, '')

    // eslreformat's output is unconditionally gzipped (hardcoded in the vendored
    // module's script); decompress so it matches the plain-text passthrough branch
    // below before merging the two back into one channel.
    GUNZIP(HMMER_ESLREFORMAT.out.seqreformated)

    emit:
    alignment = CHECKALIGNED.out.aligned
        .mix(GUNZIP.out.gunzip)
        .map { _meta, fasta -> fasta }  // channel: alignment file, aligned either way
}
