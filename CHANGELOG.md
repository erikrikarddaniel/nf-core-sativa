# nf-core/taxmarker: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/taxmarker, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Optional `raxtax`-based prefilter ahead of the EPA-ng placement stage: quickly self-classifies the reference set and reports severely mislabeled sequences directly, skipping the more expensive placement step for them ([#NN](https://github.com/nf-core/taxmarker/pull/NN))
- `test_gtdb` profile and pipeline-level tests using a curated, real archaeal 16S dataset from GTDB, exercising the pipeline on full-length real-world sequences rather than the small structural fixtures used elsewhere ([#NN](https://github.com/nf-core/taxmarker/pull/NN))
- Unaligned `--sequences` input is now supported: detected automatically (no separate mode-switch parameter) and aligned via `hmmalign` against an HMM profile (`--hmm`, optionally `--hmm_name` to pick one profile out of a multi-profile database) before continuing through the rest of the pipeline as normal ([#NN](https://github.com/nf-core/taxmarker/pull/NN))
- Sequences too short/incomplete to place reliably are now filtered out before placement, reported separately rather than silently dropped. Already-aligned input is filtered by non-gap proportion (disable with `--skip_gapfilter`; tune with `--min_nongap`, default `0.3` -- kept low since a taxonomically broad reference alignment is naturally wide, with many columns real only for a handful of divergent taxa). hmmalign-derived input is filtered separately, by HMM profile coverage (disable with `--skip_profile_cover`; tune with `--min_profile_cover`, default `0.8`) -- its own threshold since a masked hmmalign alignment's width is fixed by the profile itself, so real full-length sequences cluster tightly near full coverage rather than the broader spread seen in a directly-provided alignment ([#NN](https://github.com/nf-core/taxmarker/pull/NN))
- `--taxonomy` is now optional: if omitted, taxonomy is derived from `--sequences` record headers instead (GTDB's own single-file convention, `>id taxonomy;string`); if both are given, the file wins, with a warning logged rather than the header text being silently ignored ([#NN](https://github.com/nf-core/taxmarker/pull/NN))
- `--skip_sativa` skips the phylogenetic placement subworkflow (reference tree, leave-one-out placement, scoring) entirely, turning the pipeline into a general-purpose taxonomy-resolution/alignment/prefilter QC tool: taxonomy resolution, name-consistency checking, alignment, and gap/profile-coverage/raxtax filtering all still run as configured, without the much more expensive EPA-ng-based placement step ([#NN](https://github.com/nf-core/taxmarker/pull/NN))

### `Changed`

- `--alignment` renamed to `--sequences`, reflecting that it may be aligned or unaligned (auto-detected; see the unaligned-input entry above) ([#NN](https://github.com/nf-core/taxmarker/pull/NN))

### `Fixed`

- Conda environments now pull `biopython` from `conda-forge` instead of `bioconda`, whose build stops at 1.70: `-profile conda` failed resolving `bioconda::biopython=1.84` since bioconda never published a matching version ([#NN](https://github.com/nf-core/taxmarker/pull/NN))
- `CHECKNAMECONSISTENCY` now rewrites any character outside a safe set (was a small, growing blocklist), preventing real-world sequence identifiers (e.g. GTDB's `ACCESSION~CONTIG` names) from desyncing between the alignment/taxonomy and the tree IQTREE builds, which silently mangles the same characters in leaf names ([#NN](https://github.com/nf-core/taxmarker/pull/NN))
- `IQTREE`'s model search is now restricted to the GTR family (`-mset GTR`): ModelFinder could otherwise pick a model name (e.g. `K2P`) that EPA-ng's `--model` doesn't recognise, aborting placement ([#NN](https://github.com/nf-core/taxmarker/pull/NN))
- `SATIVALOOSPLIT` now consumes FASTA instead of PHYLIP: EMBOSS's phylip writer truncates sequence names to 10 characters, silently colliding for longer real-world identifiers ([#NN](https://github.com/nf-core/taxmarker/pull/NN))
- Unaligned input realigned via `hmmalign` no longer loses nearly every sequence to gap/coverage filtering on real, diverse datasets: `--trim` (previously assumed to restrict `hmmalign`'s output to the HMM's match-state columns) actually only trims each sequence's own unaligned terminal tails, so insert-state columns triggered by individual divergent sequences padded every other sequence with extra gap columns, badly deflating the filter's non-gap-proportion metric even for genuinely full-length sequences. `HMMER_ESLALIMASK` (`--rf-is-mask`) now masks the alignment down to true match-state columns before filtering (see the new `PROFILECOVER` entry above) ([#NN](https://github.com/nf-core/taxmarker/pull/NN))
- `IQTREE` could fail with "Taxon ... in constraint tree does not appear in full tree" on real-world runs: when `--skip_gapfilter`/`--skip_profile_cover` disabled the filtering step for the branch (already-aligned vs. hmmalign-derived) that wasn't even active for a given run, its `else` fallback paired an always-present taxonomy channel with an alignment channel that could be empty, silently desyncing their item counts and letting a stray, unfiltered taxonomy set leak through to build the reference tree against a smaller, filtered alignment ([#NN](https://github.com/nf-core/taxmarker/pull/NN))

### `Dependencies`

| Tool  | Previous version | New version |
| ----- | ---------------- | ----------- |
| HMMER |                  | 3.4         |

### `Deprecated`
