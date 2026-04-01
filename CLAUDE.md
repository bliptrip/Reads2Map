# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Reads2Map is a collection of WDL (Workflow Description Language) workflows for constructing linkage maps from sequencing reads, primarily GBS (Genotyping-by-Sequencing) data. It benchmarks combinations of 4 SNP callers (GATK, Freebayes, TASSEL, STACKs), 3 genotyping tools (updog, polyRAD, SuperMASSA), and 3 map builders (OneMap, GUSMap, MAPpoly) across both diploid and polyploid organisms.

**Citation:** Taniguti et al. (2023) GigaScience 12, giad092

## Pipeline Architecture

The system is hierarchical — top-level pipelines compose subworkflows, which call tasks:

```
EmpiricalFullPipeline (entry point)
  ├── PreprocessingReads (demux + adapter trim via STACKs/cutadapt)
  ├── GenomeIndex (samtools faidx + picard dict + BWA index)
  └── EmpiricalReads2Map
        ├── EmpiricalSNPCalling (4 callers in parallel → Array[VCF])
        └── EmpiricalMaps (per-population scatter)
              ├── ApplyRandomFiltersArray (VCF filtering + INFO/DP computation)
              ├── SplitMarkers (biallelic/multiallelic)
              ├── RemoveNonInformative (optional)
              └── scatter per chromosome × genotyper:
                    ├── ReGenotyping (updog/polyrad/supermassa)
                    ├── JointMarkers → SetProbs
                    ├── FiltersReportEmp → MapsReportEmp (OneMap)
                    ├── GUSMap maps
                    └── MAPpoly maps (polyploid only)
```

A parallel `SimulatedReads2Map` pipeline uses RADinitio/pirs/SimuSCoP to simulate reads for validation.

## Key Directories

- `pipelines/` — 7 top-level workflow WDL files (entry points for Cromwell)
- `subworkflows/` — 18 reusable intermediate workflows
- `tasks/` — 21 task modules (the actual tool invocations)
- `structs/` — WDL struct definitions (`ReferenceFasta`, `PopulationSpec`, `Specifications`, etc.)
- `.configurations/` — Cromwell backend configs (local, PostgreSQL, MySQL, SGE, Slurm)
- `.dockerfiles/` — Docker image build contexts
- `tests/` — CI test inputs and pytest-based subworkflow tests

## Common Commands

### Validate WDL syntax (all pipelines)
```bash
for i in pipelines/*/*.wdl; do java -jar womtool.jar validate $i; done
```

### Run tests (requires miniwdl + pytest)
```bash
pip3 install -r tests/requirements.txt
pytest --git-aware tests/subworkflows/gatk_genotyping/
pytest --git-aware tests/subworkflows/stacks_genotyping/
pytest --git-aware tests/subworkflows/freebayes_genotyping/
```

### Submit a pipeline via Cromwell (using pumbaa wrapper)
```bash
pumbaa wf submit -w pipelines/EmpiricalFullPipeline/EmpiricalFullPipeline.wdl \
  -i <inputs.json> -o <options.json> -d <Reads2Map.zip>
```

### Query running workflows
```bash
pumbaa wf q
```

### Cromwell API (local server at localhost:8000)
```bash
# Status
curl -s http://localhost:8000/api/workflows/v1/<UUID>/status
# Full metadata
curl -s "http://localhost:8000/api/workflows/v1/<UUID>/metadata?expandSubWorkflows=true"
```

## WDL Conventions

- **Dual container specs**: Every task declares both `docker:` and `singularity:` for portability across cloud and HPC.
- **Dual runtime attributes**: Tasks specify both cloud-style (`memory: "X MiB"`, `disks:`) and Slurm-style (`mem:`, `time:`) attributes.
- **Dynamic resource sizing**: Memory and disk are computed from input file sizes, e.g. `ceil(size(vcf, "MiB") * 12 + 4000)`.
- **Scatter parallelism**: Heavy use of `scatter` over BAM chunks, VCFs, chromosomes, families, and random seeds.
- **Conditional execution**: SNP callers are gated by booleans (`if(run_gatk)`, `if(run_freebayes)`, etc.).
- **Bash blocks should use `set -euo pipefail`** to prevent silent failures that produce empty outputs and cascade downstream.

## Critical Task Files

- `tasks/utils.wdl` — VCF manipulation: `ApplyRandomFiltersArray` (filtering + INFO/DP), `SubsetVcfToPopulation`, `SplitMarkers`, `JointMarkers`, `RemoveNonInformative`. Uses `lifebitai/bcftools:1.13`.
- `tasks/utilsR.wdl` — R-based tasks: `ReGenotyping`, `SetProbs`, `FiltersReportEmp`, `MapsReportEmp`, `CheckDepths`. Uses `cristaniguti/reads2map:0.0.8`.
- `tasks/gatk.wdl` — GATK HaplotypeCaller → GenomicsDBImport → GenotypeGVCFs pipeline.

## Data Structs

The `structs/` directory defines the typed data model. Key structs:
- `ReferenceFasta` — reference genome + all index files (`.fai`, `.dict`, BWA indices)
- `PopulationSpec` — dataset name, parents, cross type, chromosomes, pedigree aliases
- `Specifications` — barcode key files, FASTQ files, enzyme/adapter info for preprocessing

## Known Pitfalls

- **`bcftools +fill-tags` does not support `INFO/DP`** as a tag name. Use a `bcftools query + awk + bcftools annotate` pipeline instead.
- **`onemap::est_rf_out` is an internal (non-exported) function** — access it with `:::` not `::`.
- **R matrix dimension-drop**: When subsetting a matrix to 1 column/row, R silently drops to a vector. OneMap's `est_rf_out` is vulnerable to this with small linkage groups.
- **PSOCK cluster memory**: Each R PSOCK worker serializes a full copy of the data. With 4 workers, peak memory can be 5-6x the input object size. Keep workers at 2 for large objects.
- **Cromwell call caching**: Cromwell hashes command strings for cache keys. If a task produces wrong output (e.g., empty VCF) and is cached, fixing the WDL command text will invalidate the cache automatically. If only inputs change, you may need `write_to_cache: false` in options.
- **Stacks VCFs lack `##contig` headers**, which causes `bcftools view --samples` to abort on bcftools >=1.13. Inject contig headers via `bcftools reheader` before subsetting.
- **Cromwell execution root**: `/mnt/share/maule2/Zalapa/Data/Cranberries/cromwell-executions` — task logs are deeply nested under `call-<TaskName>/execution/stderr`.

## Agent Infrastructure

Five specialized Claude agents are defined in `.claude/agents/` for supervised pipeline execution:
1. **empirical-pipeline-supervisor** — orchestrates launch → monitor → debug cycle
2. **empirical-pipeline-launcher** — runs make targets, extracts workflow UUID
3. **cromwell-workflow-monitor** — polls Cromwell API for status
4. **cromwell-workflow-debugger-lite** — rapid failure triage
5. **reads2map-pipeline-debugger** — deep root cause analysis with suggested fixes

The supervisor is invoked via `make sup` from the project Makefile.

## CI/CD

CircleCI pipeline (`.circleci/config.yml`):
1. `wdl-validate` — womtool syntax check on all pipelines
2. `run-tasks-tests` — pytest on subworkflow tests (GATK, STACKs, Freebayes)
3. `release-to-github` — automated release bundling (dev on `develop` branch, prod on tags)
