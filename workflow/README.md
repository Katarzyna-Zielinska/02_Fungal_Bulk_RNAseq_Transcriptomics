# Workflow

## Purpose

This directory is reserved for workflow-level orchestration and execution definitions for the *Candida auris* bulk RNA-seq project.

The current analysis is implemented as a collection of R and Python scripts in `scripts/`, supported by command-line tools and two Conda environments (`p2_ngs` and `p2_rmats`). The project audit did not identify an implemented workflow-engine entry point in this directory.

## Current status

**No executable workflow definition is currently provided in this directory.**

In particular, the audited project copy did not contain:
- `main.nf`;
- `nextflow.config`;
- a Snakemake `Snakefile`;
- workflow-specific rule or process files.

Nextflow is installed in the analysis environment, but its installation alone does not mean that this project already has a runnable Nextflow pipeline. This README documents the current status and intended role of the directory; it does not claim that workflow orchestration has been implemented or tested.

## Current analysis organization

The analysis is currently organized into script-based stages:

1. Read QC, preprocessing, alignment, and transcript quantification.
2. Baseline transcript import and sample QC.
3. Baseline differential expression with edgeR and limma-voom.
4. Functional annotation and baseline pathway analysis.
5. rMATS alternative-splicing analysis and integration with expression.
6. Duplicate-aware transcript quantification and grouped differential expression.
7. Grouped fgsea/GSVA analysis and leading-edge integration.
8. GO enrichment, core-gene selection, and annotation refinement.

The detailed script inventory, dependencies, and stage relationships are documented in [`../docs/pipeline_overview.md`](../docs/pipeline_overview.md).

## Intended future contents

If workflow orchestration is added, this directory may contain:

| File | Intended role |
|---|---|
| `main.nf` or `Snakefile` | Main workflow entry point, depending on the selected workflow engine |
| `nextflow.config` or workflow configuration | Parameters, profiles, resource settings, and execution options |
| `modules/` or `rules/` | Reusable process or rule definitions |
| `README.md` | Workflow-specific usage and execution instructions |

Only add files for the workflow engine actually selected and implemented. Do not keep placeholder definitions that imply the pipeline is executable.

## Execution notes

Until a workflow engine is implemented, run the existing scripts in their documented dependency order from the project root. Check that each required input exists before running a stage and verify the expected outputs afterwards.

Some scripts in the current working copy contain machine-specific absolute paths. These need to be made portable before the project can be described as reproducible from a fresh clone.

See:
- [`../docs/pipeline_overview.md`](../docs/pipeline_overview.md) for script order and dependencies;
- [`../docs/reproducibility.md`](../docs/reproducibility.md) for environment versions and reproducibility limitations;
- [`../README.md`](../README.md) for the project overview.
