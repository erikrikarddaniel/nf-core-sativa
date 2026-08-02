# nf-core/sativa: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/sativa, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Optional `raxtax`-based prefilter ahead of the EPA-ng placement stage: quickly self-classifies the reference set and reports severely mislabeled sequences directly, skipping the more expensive placement step for them ([#NN](https://github.com/nf-core/sativa/pull/NN))
- `test_gtdb` profile and pipeline-level tests using a curated, real archaeal 16S dataset from GTDB, exercising the pipeline on full-length real-world sequences rather than the small structural fixtures used elsewhere ([#NN](https://github.com/nf-core/sativa/pull/NN))
- Unaligned `--alignment` input is now supported: detected automatically (no separate mode-switch parameter) and aligned via `hmmalign` against an HMM profile (`--hmm`, optionally `--hmm_name` to pick one profile out of a multi-profile database) before continuing through the rest of the pipeline as normal ([#NN](https://github.com/nf-core/sativa/pull/NN))
- Sequences with too high a proportion of alignment gaps to place reliably are now filtered out before placement, reported separately rather than silently dropped (disable with `--skip_gapfilter`; tune the threshold with `--min_nongap`, default `0.3` -- kept low since a taxonomically broad reference alignment is naturally wide, with many columns real only for a handful of divergent taxa) ([#NN](https://github.com/nf-core/sativa/pull/NN))

### `Fixed`

- `CHECKNAMECONSISTENCY` now rewrites any character outside a safe set (was a small, growing blocklist), preventing real-world sequence identifiers (e.g. GTDB's `ACCESSION~CONTIG` names) from desyncing between the alignment/taxonomy and the tree IQTREE builds, which silently mangles the same characters in leaf names ([#NN](https://github.com/nf-core/sativa/pull/NN))
- `IQTREE`'s model search is now restricted to the GTR family (`-mset GTR`): ModelFinder could otherwise pick a model name (e.g. `K2P`) that EPA-ng's `--model` doesn't recognise, aborting placement ([#NN](https://github.com/nf-core/sativa/pull/NN))
- `SATIVALOOSPLIT` now consumes FASTA instead of PHYLIP: EMBOSS's phylip writer truncates sequence names to 10 characters, silently colliding for longer real-world identifiers ([#NN](https://github.com/nf-core/sativa/pull/NN))

### `Dependencies`

| Tool  | Previous version | New version |
| ----- | ---------------- | ----------- |
| HMMER |                  | 3.4         |

### `Deprecated`
