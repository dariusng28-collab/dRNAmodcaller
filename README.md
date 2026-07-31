# dRNA Methylation Pipeline

A Snakemake workflow for Oxford Nanopore direct RNA sequencing (dRNA-seq) basecalling, alignment, RNA modification pileup, and QC.

The workflow is structured after the official Snakemake workflow-template style: workflow code lives under `workflow/`, user configuration lives under `config/`, and runtime behavior is controlled through config files, samplesheets, conda environments, containers, and profiles.

## Workflow

```text
pod5 files per sample
  -> Dorado model download
  -> Dorado basecalling (whole sample or in chunks)
  -> merge unaligned BAMs per sample
  -> build minimap2 index (per alignment)
  -> genome and/or transcriptome alignment
  -> modkit sample-probs
  -> modkit pileup
  -> modkit summary
  -> bedMethyl filtering
  -> split bedMethyl by modification code
  -> per-sample QC, NanoPlot, aggregate MultiQC report
  -> provenance
```

Genome and transcriptome outputs are separated by `{alignment}` so coordinate spaces do not collide.

## Repository Layout

```text
dRNA_methylation/
|-- README.md
|-- LICENSE
|-- CITATION.cff
|-- submit.sh                       # SGE convenience wrapper
|-- config/
|   |-- config.yaml                 # analysis configuration
|   `-- samples.csv                 # sample metadata and pod5 paths
|-- workflow/
|   |-- Snakefile                   # workflow entrypoint (auto-discovered)
|   |-- rules/
|   |   |-- common.smk
|   |   |-- reference.smk
|   |   |-- basecalling.smk
|   |   |-- alignment.smk
|   |   |-- modkit.smk
|   |   `-- qc.smk
|   |-- schemas/
|   |   |-- config.schema.yaml
|   |   `-- samples.schema.yaml
|   |-- scripts/
|   |   |-- split_bed_by_mod.py
|   |   `-- write_provenance.py
|   `-- profiles/
|       `-- sge/
|           `-- profile.yaml
`-- .tests/
    `-- unit/
        `-- test_split_bed_by_mod.sh
```

There is no top-level `Snakefile`: Snakemake automatically discovers
`workflow/Snakefile`. Always invoke Snakemake from the repository root so that
`config/`, `workflow/`, and the samplesheet paths resolve correctly.

## Requirements

- Snakemake with container support.
- Singularity or Apptainer. **Every rule runs in a container — no conda
  required.** Dorado and Modkit use local `.sif` images (or ONT image URIs); the
  generic tools (minimap2/samtools, NanoPlot, MultiQC, Python) use pinned public
  images that Snakemake pulls automatically on first run.
- Optional, for SGE cluster execution: Snakemake 8+ and
  `snakemake-executor-plugin-cluster-generic`. The SGE profile uses the
  cluster-generic executor with a `qsub` wrapper and a `qstat`-based status
  check (`workflow/scripts/sge_submit.py` / `sge_status.sh`), so it works on
  clusters where `qacct` accounting is unavailable.

Dorado and Modkit are the only images you supply. Pull them once, e.g.:

```bash
singularity pull ontresearch-dorado-1.4.0.sif docker://ontresearch/dorado:1.4.0
singularity pull ontresearch-modkit.sif docker://ontresearch/modkit:latest
```

The other images are defined (and overridable) in `config/config.yaml` under
`minimap_container`, `nanoplot_container`, `multiqc_container`, and
`python_container`, with pinned defaults so nothing else needs configuring.

### Container image caching on HPC

Converting a docker image to a Singularity `.sif` unpacks tens of thousands of
small files, which is slow on network storage and fails if the site's default
`SINGULARITY_TMPDIR` points somewhere that does not exist on the submit host.
`submit.sh` therefore builds images on **node-local disk** and stores the built
images in a **persistent, shared prefix** so they are pulled once and reused
across runs (compute nodes read them over the shared filesystem). Override the
defaults in `config/config.yaml` if needed:

- `singularity_tmpdir`: fast local scratch for the SIF build (default `/tmp`).
- `singularity_prefix`: shared directory of built images (default
  `<outdir>/singularity`). Point this at a stable location to reuse images
  across working directories.

## Configuration

Edit `config/config.yaml`.

Important keys:

- `outdir`: root directory for all results and workflow logs.
- `samplesheet`: CSV file with sample metadata and pod5 directories.
- `dorado_container`, `modkit_container`: local `.sif` paths or container URIs.
- `dorado_models`: Dorado canonical model plus modification models.
- `basecall_chunk_size`: pod5 files per GPU job. `0` (default) basecalls each
  sample's whole pod5 directory in one dorado invocation (models load onto the
  GPU once per sample); a positive integer splits large samples across GPUs.
- `alignments`: enable `genome`, `transcriptome`, or both.
- `modifications`: enabled modification definitions used for Modkit and split BED outputs.
- `modkit.pileup`: filter thresholds, per-mod thresholds, and advanced extra args.

The older `dorado_sif` and `modkit_sif` keys are still accepted for compatibility, but new configs should use `dorado_container` and `modkit_container`.

## Samplesheet

The samplesheet can be named anything as long as it is a CSV. It must include one recognized sample column and one recognized pod5 directory column.

Accepted sample headers:

```text
sample, sample_id, sampleid, id, name
```

Accepted pod5 directory headers:

```text
pod5_dir, pod5_directory, pod5_path, pod5, raw_dir, raw_path, path
```

Extra columns such as `group`, `condition`, and `replicate` are preserved in provenance.

Example:

```csv
sample,group,condition,replicate,pod5_dir
WT_01,WT,control,1,/path/to/WT_01/pod5/
KO_01,KO,knockout,1,/path/to/KO_01/pod5/
```

## Biological Interpretation

Genome alignment uses `minimap2 -ax splice -uf` with a junction BED and produces genomic coordinates. This is the better default when the downstream question is locus-level modification evidence, alignment QC, or integration with genome annotations.

Transcriptome alignment uses `minimap2 -ax map-ont -k14` against transcript FASTA and passes `--preload-references` to Modkit pileup. This is useful for transcript or isoform-level direct RNA analyses.

For Modkit v0.6.0 and later, the workflow passes `--reference` and config-driven `--modified-bases`. Each enabled modification defines:

- `primary_base`: reference base for Modkit, using `T` for uridine/pseudouridine because FASTA uses DNA alphabet.
- `code`: Modkit/SAM modification code and bedMethyl column-4 split key.
- `output_suffix`: output BED suffix.
- `mod_threshold`: optional per-modification call threshold.

Default supported RNA modification entries:

| Modification | Modkit modified base | BED code |
| --- | --- | --- |
| m5C | `C:m` | `m` |
| m6A | `A:a` | `a` |
| inosine | `A:17596` | `17596` |
| pseU | `T:17802` | `17802` |
| 2OmeC | `C:19228` | `19228` |
| 2OmeA | `A:69426` | `69426` |
| 2OmeG | `G:19229` | `19229` |
| 2OmeU | `T:19227` | `19227` |

## Running

Dry-run locally:

```bash
snakemake -n --configfile config/config.yaml
```

Run locally:

```bash
snakemake --cores 8 --use-singularity --configfile config/config.yaml
```

For Apptainer installations, replace `--use-singularity` with `--use-apptainer`.

Override output root:

```bash
snakemake -n --configfile config/config.yaml --config outdir=test_results
```

Run on SGE with the workflow profile:

```bash
bash submit.sh
```

The SGE wrapper stores scheduler logs under `${outdir}/logs/cluster/<timestamp>` and forwards any extra arguments to Snakemake:

```bash
bash submit.sh -n
```

## Output Structure

```text
results/
|-- references/
|   `-- index/{alignment}.mmi        # prebuilt minimap2 index (per alignment)
|-- basecalled/
|   `-- {sample}/
|       |-- chunks/{chunk}.unaligned.bam   # per-chunk dorado output (temp)
|       `-- {sample}.merged.unaligned.bam
|-- bams/
|   `-- {alignment}/{sample}/{sample}.bam
|-- modkit/
|   `-- {alignment}/{sample}.raw.bed
|-- bedMethyl/
|   |-- {alignment}/{sample}.filtered.bed
|   `-- {alignment}/{sample}.{mod}.filtered.bed
|-- sample_probs/
|   `-- {alignment}/{sample}/
|-- qc/
|   |-- {alignment}/{sample}.qc_summary.txt   # per-sample text summary
|   |-- samtools/{alignment}.{sample}.*.txt    # stats / flagstat / idxstats
|   |-- nanoplot/{sample}/                      # read-length & quality plots
|   `-- multiqc/multiqc_report.html            # aggregate report (all samples)
|-- provenance/
|   |-- config.resolved.yaml
|   `-- run_metadata.txt
`-- logs/
```

## Quality Control

The workflow produces QC at three levels:

- **Per-sample text summary** (`qc/{alignment}/{sample}.qc_summary.txt`): read
  counts, `samtools flagstat`, top reference sequences, estimated poly(A) tail
  lengths, modification-probability histogram, and per-modification site counts.
- **NanoPlot** (`qc/nanoplot/{sample}/`): alignment-independent read-length and
  quality distributions (N50, mean/median read length and quality) computed from
  the raw basecalled reads.
- **MultiQC** (`qc/multiqc/multiqc_report.html`): a single interactive report
  aggregating `samtools stats`/`flagstat`/`idxstats` and NanoPlot across every
  sample and alignment. This is the recommended starting point for run-level QC.

## Reproducibility

Every run writes `provenance/config.resolved.yaml` (the fully resolved
configuration) and `provenance/run_metadata.txt`. Every tool runs in a
version-pinned container, so the image references recorded in the provenance
metadata are the authoritative software-version record — a run can be reproduced
from the provenance file plus those images.

## Developer Checks

Run these before publishing or sharing workflow changes:

```bash
snakemake --lint
snakefmt --check workflow/Snakefile workflow/rules/*.smk
bash .tests/unit/test_split_bed_by_mod.sh
```

Useful dry-run scenarios:

```bash
snakemake -n --configfile config/config.yaml
snakemake -n --configfile config/config.yaml --config outdir=test_results
```

To dry-run with both alignments, enable `transcriptome` under `alignments` in
`config/config.yaml` (nested `--config` overrides are not reliably parsed by the
CLI, so edit the config file for structured changes).

## Troubleshooting

- If config validation fails, fix the first reported missing path or invalid sample entry before rerunning.
- If Modkit produces empty split BEDs, inspect `results/modkit/{alignment}/{sample}.raw.bed` and confirm column 4 contains the configured modification code.
- If transcriptome alignment is enabled, `transcriptome` must point to the transcript FASTA used for downstream interpretation.
- If SGE submission fails, first verify the SGE executor plugin is installed and that `qsub`, `qstat`, and `qacct` are available on the login node.

## Citations

If you use this workflow, please cite Snakemake and the underlying tools:

- Mölder F, *et al.* Sustainable data analysis with Snakemake. *F1000Research* 2021;10:33. doi:10.12688/f1000research.29032.2
- Oxford Nanopore Technologies. Dorado basecaller. https://github.com/nanoporetech/dorado
- Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 2018;34(18):3094–3100. doi:10.1093/bioinformatics/bty191
- Danecek P, *et al.* Twelve years of SAMtools and BCFtools. *GigaScience* 2021;10(2):giab008. doi:10.1093/gigascience/giab008
- Oxford Nanopore Technologies. Modkit: a bioinformatics tool for working with modified bases. https://github.com/nanoporetech/modkit
- De Coster W, Rademakers R. NanoPack2: population-scale evaluation of long-read sequencing data. *Bioinformatics* 2023;39(5):btad311. doi:10.1093/bioinformatics/btad311
- Ewels P, *et al.* MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics* 2016;32(19):3047–3048. doi:10.1093/bioinformatics/btw354

## License

Released under the MIT License. See [LICENSE](LICENSE).
