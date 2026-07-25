# AGENTS.md

This file provides guidance to coding agents (e.g. Claude Code) when working with code in this repository.

## Project Overview

**nf-core/sativa** is a Nextflow bioinformatics pipeline that checks sequences for their phylogenetic signal against their taxonomy. It is built from the nf-core template (v4.0.2) and uses Nextflow DSL2. The pipeline is currently in early development (v1.0.0dev) — many `TODO nf-core:` comments mark where domain-specific logic still needs to be added.

Requires Nextflow ≥ 25.10.4.

## Commands

**Run the pipeline (with Docker):**

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results
```

**Run minimal test suite:**

```bash
nextflow run main.nf -profile test,docker --outdir results
```

**Run nf-test (unit/integration tests):**

```bash
nf-test test tests/default.nf.test
# Run all tests (respects nf-test.config ignore rules):
nf-test test
```

**Lint with nf-core tools:**

```bash
nf-core pipelines lint
```

**Format code (Prettier + Nextflow lint via pre-commit):**

```bash
pre-commit run --all-files
```

**Update nf-core modules:**

```bash
nf-core modules update <module-name>
```

## Architecture

### Execution flow

```
main.nf
  └── PIPELINE_INITIALISATION   (subworkflows/local/utils_nfcore_sativa_pipeline/main.nf)
        validates params, parses samplesheet → ch_samplesheet channel
  └── NFCORE_SATIVA
        └── SATIVA               (workflows/sativa.nf)  ← main logic lives here
              ├── FASTQC          (modules/nf-core/fastqc/)
              └── MULTIQC         (modules/nf-core/multiqc/)
  └── PIPELINE_COMPLETION        (subworkflows/local/utils_nfcore_sativa_pipeline/main.nf)
        sends email / completion summary
```

### Key conventions

- **Module arguments**: Pass extra CLI flags to tools via `ext.args` in `conf/modules.config`, not in the module itself.
- **Output paths**: Default publish rule in `conf/modules.config` derives directory from the process name (e.g., `FASTQC` → `outdir/fastqc/`). Override per-process with a `publishDir` block.
- **Samplesheet input**: Validated against `assets/schema_input.json`. Required columns: `sample`, `fastq_1`; optional: `fastq_2`. Single-end vs. paired-end is inferred from the presence of `fastq_2`.
- **Parameter schema**: `nextflow_schema.json` defines all pipeline parameters and is used for CLI validation (via nf-schema plugin) and help text generation.
- **Software versions**: Collected via a `channel.topic("versions")` stream and written to `pipeline_info/nf_core_sativa_software_mqc_versions.yml` for MultiQC.
- **nf-core modules**: Modules under `modules/nf-core/` and subworkflows under `subworkflows/nf-core/` are managed by nf-core tools — do not edit them directly. Custom/local code goes in `subworkflows/local/`.

### Container registries

All container profiles (`docker`, `singularity`, `apptainer`, etc.) default to `quay.io` as the registry. The `wave` profile enables on-demand container building via Seqera Wave, required for ARM64.

### Test infrastructure

- `nf-test.config` defines test directories and triggers (files that force a full test run when changed).
- Tests run with `-profile test` by default; the test profile uses nf-core's public test datasets hosted on GitHub.
- Snapshot files (`*.snap`) track expected outputs — update them with `nf-test test --update-snapshot` after intentional output changes.
