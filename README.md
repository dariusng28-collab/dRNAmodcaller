# dRNA Methylation Pipeline

A Snakemake pipeline for RNA base modification calling from Oxford Nanopore direct RNA sequencing (dRNA-seq) data. Based on the [DOGME](https://github.com/mortazavilab/dogme) workflow and the approach of Laurens Lambrechts, adapted for SGE cluster execution.

---
## Repository structure

dRNA_methylation_pipeline/
├── README.md               — this file
├── Snakefile               — main pipeline
├── cluster.yaml            — SGE resource requests per rule
├── submit.sh               — cluster submission script
├── cluster_qsub.sh         — SGE jobscript template
├── config/
│   ├── config.yaml         — paths, containers, samplesheet, modifications, thresholds
│   └── samples.csv         — sample metadata and pod5 directories
└── env/
    └── minimap.yaml        — conda environment (samtools, minimap2, bedtools)


## Workflow

```
pod5 files (per sample)
    └── 0. dorado_models_download  — download models once (skipped if present)
    └── 1. dorado_basecall         — basecall each pod5 → unaligned BAM
    └── 2. merge_unaligned_bam     — merge per-pod5 BAMs into one per sample
    └── 3. minimap2_align          — genome and/or transcriptome alignment
            ├── 4. sample_probs    — QC: modification probability distributions
            ├── 5. modkit_pileup   — call modifications → raw BED
            └── 6. modkit_summary  — QC: global modification summary
                    └── 7. filterbed   — filter by coverage and mod percentage
                            └── 8. splitbed    — split BED by modification type
                                    └── 9. qc_report   — final QC report
```

### Supported modifications

| Modification | Description | Modkit `--modified-bases` | BED code |
|---|---|---|---|
| m5C | 5-Methylcytosine | `C:m` | `m` |
| m6A | N6-Methyladenosine | `A:a` | `a` |
| inosine | Inosine (A-to-I editing) | `A:17596` | `17596` |
| pseU | Pseudouridine | `T:17802` | `17802` |
| 2OmeC | 2'-O-methylcytidine | `C:19228` | `19228` |
| 2OmeA | 2'-O-methyladenosine | `A:69426` | `69426` |
| 2OmeG | 2'-O-methylguanosine | `G:19229` | `19229` |
| 2OmeU | 2'-O-methyluridine | `T:19227` | `19227` |

---

## Requirements

### Software

- [Snakemake](https://snakemake.readthedocs.io) ≥ 7.0
- [Singularity](https://sylabs.io/singularity/) ≥ 3.8 or Apptainer
- [Conda](https://docs.conda.io) or [Mamba](https://github.com/mamba-org/mamba)
- SGE cluster with GPU nodes

### Singularity containers

```bash
# Dorado basecaller
singularity pull ontresearch-dorado-1.4.0.sif docker://ontresearch/dorado:1.4.0

# modkit modification calling
singularity pull ontresearch-modkit-0.6.3.sif docker://ontresearch/modkit:0.6.3
```

### Conda environments

```bash
conda env create -f env/minimap.yaml   # samtools, minimap2, bedtools
```

### Reference files

Download from Gencode:

```bash
# Reference genome
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.p14.genome.fa.gz

# GTF annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz

# Transcriptome FASTA, if using transcriptome alignment
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz
```

Generate the junction BED file once from the GTF (required for splice-aware alignment):

```bash
paftools.js gff2bed gencode.v49.primary_assembly.annotation.gtf > gencode.v49.basic.annotation.junc.bed
```

`paftools.js` is included with minimap2.

---

## Installation

```bash
git clone https://github.com/your-username/dRNA_methylation_pipeline.git
cd dRNA_methylation_pipeline
```

---

## Configuration

### 1. Edit `config/config.yaml`

Set all file paths and thresholds:

```yaml
# Output root for all workflow results and logs
outdir: "results"

# Container runtime for Dorado/modkit rules: "singularity" or "apptainer"
container_engine: "singularity"

# Singularity/Apptainer bind mount — adjust to your cluster mount point.
# Leave empty ("") if no bind mount is needed.
singularity_bind: "/your/data/mount:/your/data/mount"

# Container paths (dorado v.1.4.0; modkit v.0.6.3 recommended)
dorado_sif: "/path/to/ontresearch-dorado-1.4.0.sif"
modkit_sif: "/path/to/ontresearch-modkit-0.6.3.sif"

# Dorado models — full names required with --models-directory.
# These are downloaded automatically and passed to dorado basecaller.
dorado_models:
  - "rna004_130bps_sup@v5.3.0"
  - "rna004_130bps_sup@v5.3.0_inosine_m6A_2OmeA@v1"
  - "rna004_130bps_sup@v5.3.0_m5C_2OmeC@v1"

# Local directory for pre-downloaded dorado models
dorado_models_dir: "/path/to/dorado_models"

# Alignment branches.
# Genome alignment gives genomic coordinates and uses splice-aware minimap2.
# Transcriptome alignment is useful for transcript/isoform-level dRNA analyses.
alignments:
  genome:
    enabled: true
  transcriptome:
    enabled: false

# Modkit pileup settings
modkit:
  pileup:
    filter_thresholds:
      - "0.7"
    mod_thresholds: []     # optional override; otherwise use modifications.*.mod_threshold
    extra_args: []

# Enabled modifications.
# primary_base + code builds Modkit v0.6+ --modified-bases entries.
# code also identifies the modification in bedMethyl column 4.
modifications:
  m5C:
    enabled: true
    primary_base: "C"
    code: "m"
    output_suffix: "m5C"
    mod_threshold: 0.99
  m6A:
    enabled: true
    primary_base: "A"
    code: "a"
    output_suffix: "m6A"
    mod_threshold: 0.99

# Reference files
genome:  "/path/to/GRCh38.p14.genome.fa"
transcriptome: "/path/to/gencode.v49.transcripts.fa"
gtf:     "/path/to/gencode.v49.primary_assembly.annotation.gtf"
gtf_bed: "/path/to/gencode.v49.basic.annotation.junc.bed"

# Filterbed thresholds (can be changed and re-applied without re-running pileup)
min_coverage: 5   # minimum reads covering a site (col 5 of modkit BED)
mod_pct:      5   # minimum % modified reads per site (col 11 of modkit BED)

# Samplesheet with sample metadata and pod5 directories.
# The CSV filename can be anything.
samplesheet: "config/samples.csv"
```

Pod5 files must be organised in **separate directories per sample**.

The samplesheet must be a CSV and can be named anything. It must include a recognised sample-name column and a recognised pod5-directory column:

- Sample-name headers: `sample`, `sample_id`, `sampleid`, `id`, `name`
- Pod5-directory headers: `pod5_dir`, `pod5_directory`, `pod5_path`, `pod5`, `raw_dir`, `raw_path`, `path`

Extra metadata columns such as `group`, `condition`, and `replicate` are preserved in provenance and can be used by future downstream analysis steps.

```csv
sample,group,condition,replicate,pod5_dir
WT_01,WT,control,1,/path/to/WT_01/pod5/
KO_01,KO,knockout,1,/path/to/KO_01/pod5/
```

### 2. Choose modifications and thresholds in `config/config.yaml`

Disable a modification by setting `enabled: false`. Final split BED outputs and QC counts are generated only for enabled modifications.

Each enabled modification must define:

- `primary_base`: reference base used by Modkit `--modified-bases` (`A`, `C`, `G`, or `T`). For direct RNA, uridine is represented as `T` in FASTA references.
- `code`: Modkit/SAM modification code, also used to split bedMethyl column 4.
- `output_suffix`: suffix for `{sample}.{suffix}.filtered.bed`.
- `mod_threshold`: optional per-modification threshold used to build repeated `--mod-threshold` arguments.

For Modkit v0.6.0+, the Snakefile automatically passes `--reference` and `--modified-bases` using the enabled modification definitions, e.g. `A:a C:m A:17596`.

**`modkit.pileup.filter_thresholds`** — values emitted as repeated `--filter-threshold` arguments. Use entries like `"0.7"` or `"C:0.8"`.

**`modifications.*.mod_threshold`** — per-modification thresholds emitted as repeated `--mod-threshold` arguments for `modkit pileup`. Use `modkit.pileup.mod_thresholds` only if you want to override the generated list directly.

### 3. Choose alignment branches

Use `alignments.genome.enabled` and `alignments.transcriptome.enabled` to control which reference spaces are processed.

Genome alignment uses `minimap2 -ax splice -uf` with the configured `gtf_bed` junction BED, producing genomic coordinates. Transcriptome alignment uses `minimap2 -ax map-ont -k14` against the transcript FASTA, and Modkit pileup receives `--preload-references`, which Modkit recommends for transcriptome-aligned direct RNA.

Both branches share the same downstream Modkit/QC rules and write into separate output folders, so the results do not collide.



---

## Running the pipeline


## Output structure

```
results/  # or the configured outdir
├── basecalled/
│   └── {sample}/
│       └── {sample}.merged.unaligned.bam     ← merged unaligned BAM
├── bams/
│   └── {alignment}/
│       └── {sample}/
│       ├── {sample}.bam                       ← aligned BAM
│       └── {sample}.bam.bai
├── modkit/
│   └── {alignment}/
│       ├── {sample}.raw.bed                   ← all sites, unfiltered
│       └── {sample}.summary.tsv               ← modkit summary QC
├── bedMethyl/
│   └── {alignment}/
│       ├── {sample}.filtered.bed              ← all sites, filtered
│       └── {sample}.{enabled_modification_suffix}.filtered.bed
├── sample_probs/
│   └── {alignment}/
│       └── {sample}/
│           ├── probabilities.tsv
│           ├── thresholds.tsv
│           ├── counts.html                    ← open in browser
│           └── proportion.html                ← open in browser
├── qc/
│   └── {alignment}/
│       └── {sample}.qc_summary.txt            ← final QC report
└── provenance/
    ├── config.resolved.yaml                  ← config used for the run
    └── run_metadata.txt                      ← models, tool versions, reference checksums
```

If a modification is absent in the data, the corresponding BED file will be empty — the pipeline does not fail.

---
