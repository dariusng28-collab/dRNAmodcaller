import csv
import os
import re
from pathlib import Path

# ================================================================
#
#   dRNA METHYLATION PIPELINE — Snakemake / SGE cluster
#   Based on DOGME (nanoporeModule.nf v1.2.3)
#   + Laurens Lambrechts approach (dorado v1.4.0 + junc.bed)
#
#   WORKFLOW:
#   0.  dorado_models_download — download dorado models once
#   1.  dorado_basecall        — basecall each pod5 → unaligned BAM
#   2.  merge_unaligned_bam    — merge BAMs per sample
#   3.  minimap2_align         — genome/transcriptome alignment
#   4.  sample_probs           — modkit sample-probs QC
#   5.  modkit_pileup          — call modifications
#   6.  modkit_summary         — modkit summary QC
#   7.  filterbed              — filter by coverage and mod percentage
#   8.  splitbed               — split BED by modification type
#   9.  qc_report              — final QC report
#
# ================================================================
#
#   USER CONFIGURATION
#
#   Paths, models, modifications, thresholds, and samples are set in
#   config/config.yaml. Ordinary runs should not require Snakefile edits.
#
# ================================================================

configfile: "config/config.yaml"

DORADO_SIF = config.get("dorado_sif", "")
MODKIT_SIF = config.get("modkit_sif", "")
CONTAINER_ENGINE = config.get("container_engine", "singularity")
BIND       = config.get("singularity_bind", "")
SINGULARITY_BIND_ARG = f"-B {BIND}" if BIND else ""

OUTDIR = config.get("outdir", "results").rstrip("/")
LOGDIR = f"{OUTDIR}/logs"
WORKFLOW_DIR = Path(workflow.basedir)

def out_path(path):
    return f"{OUTDIR}/{path}"

def log_path(path):
    return f"{LOGDIR}/{path}"

def resolve_path(path):
    p = Path(str(path)).expanduser()
    if not p.is_absolute():
        p = WORKFLOW_DIR / p
    return p

def config_list(value):
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    return list(value or [])

def config_args(value):
    if isinstance(value, str):
        return value.strip()
    return " ".join(str(item) for item in list(value or []))

def validate_token(value, field):
    if not re.match(r"^[A-Za-z0-9_.+-]+$", str(value)):
        raise ValueError(
            f"Invalid {field!r} value {value!r}. Use letters, numbers, '.', '_', '+', or '-'."
        )

def parse_modifications():
    mods = config.get("modifications")
    if not isinstance(mods, dict) or not mods:
        raise ValueError("Configure at least one entry under 'modifications' in config/config.yaml.")

    parsed = []
    seen_suffixes = set()
    for name, settings in mods.items():
        settings = settings or {}
        if not settings.get("enabled", True):
            continue

        primary_base = str(
            settings.get("primary_base", settings.get("canonical_base", ""))
        ).strip().upper()
        code = str(settings.get("code", "")).strip()
        output_suffix = str(settings.get("output_suffix", name)).strip()
        if primary_base not in {"A", "C", "G", "T"}:
            raise ValueError(
                f"Modification {name!r} must set primary_base to one of A, C, G, or T."
            )
        if not code:
            raise ValueError(f"Modification {name!r} is enabled but has no 'code'.")
        if not output_suffix:
            raise ValueError(f"Modification {name!r} is enabled but has no 'output_suffix'.")

        validate_token(name, "modification name")
        validate_token(code, f"modification {name} code")
        validate_token(output_suffix, f"modification {name} output_suffix")
        if output_suffix in seen_suffixes:
            raise ValueError(f"Duplicate modification output_suffix {output_suffix!r}.")
        seen_suffixes.add(output_suffix)

        parsed.append({
            "name": name,
            "primary_base": primary_base,
            "code": code,
            "output_suffix": output_suffix,
            "mod_threshold": settings.get("mod_threshold", settings.get("threshold"))
        })

    if not parsed:
        raise ValueError("At least one modification must be enabled in config/config.yaml.")
    return parsed

def check_existing_path(errors, key, path, kind="file"):
    p = resolve_path(path)
    if kind == "dir" and not p.is_dir():
        errors.append(f"{key}: directory does not exist: {path}")
    elif kind == "file" and not p.is_file():
        errors.append(f"{key}: file does not exist: {path}")

SAMPLE_COLUMN_ALIASES = ("sample", "sample_id", "sampleid", "id", "name")
POD5_DIR_COLUMN_ALIASES = (
    "pod5_dir", "pod5_directory", "pod5_path", "pod5", "raw_dir", "raw_path", "path"
)

def find_column(fieldnames, aliases, label, samplesheet):
    lower_to_original = {field.lower().strip(): field for field in fieldnames or []}
    for alias in aliases:
        if alias in lower_to_original:
            return lower_to_original[alias]
    raise ValueError(
        f"samplesheet {samplesheet!r} is missing a {label} column. "
        f"Accepted headers: {', '.join(aliases)}"
    )

def parse_samples():
    samplesheet = str(config.get("samplesheet", "") or "").strip()
    if not samplesheet:
        raise ValueError("Configure a CSV samplesheet with the 'samplesheet' key in config/config.yaml.")

    sheet_path = resolve_path(samplesheet)
    if not sheet_path.is_file():
        raise ValueError(f"Configured samplesheet does not exist: {samplesheet}")
    if sheet_path.suffix.lower() != ".csv":
        raise ValueError(f"Configured samplesheet must be a .csv file: {samplesheet}")

    with open(sheet_path, newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise ValueError(f"samplesheet {samplesheet!r} has no header row.")
        sample_col = find_column(reader.fieldnames, SAMPLE_COLUMN_ALIASES, "sample name", samplesheet)
        pod5_col = find_column(reader.fieldnames, POD5_DIR_COLUMN_ALIASES, "pod5 directory", samplesheet)

        rows = []
        for line_number, row in enumerate(reader, start=2):
            sample = (row.get(sample_col) or "").strip()
            pod5_dir = (row.get(pod5_col) or "").strip()
            if not sample and not pod5_dir:
                continue
            if not sample or not pod5_dir:
                raise ValueError(
                    f"samplesheet {samplesheet!r} line {line_number} must include sample and pod5 directory values."
                )
            row = {k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
            row["sample"] = sample
            row["pod5_dir"] = pod5_dir
            rows.append(row)
        return rows

DORADO_MODELS = config_list(config.get("dorado_models", config.get("dorado_model")))
if not DORADO_MODELS:
    raise ValueError("Configure at least one Dorado model with 'dorado_models' in config/config.yaml.")

DORADO_MODELS_DIR = config.get("dorado_models_dir", "")
DORADO_MODEL_ARG = ",".join(DORADO_MODELS)
DORADO_MODEL_DIR_PATHS = [
    f"{DORADO_MODELS_DIR}/{model}"
    for model in DORADO_MODELS
]
DORADO_MODEL_DIRS = [
    directory(path)
    for path in DORADO_MODEL_DIR_PATHS
]

GTF_BED_CONFIG = str(config.get("gtf_bed", "") or "").strip()
GENERATE_GTF_BED = not bool(GTF_BED_CONFIG)
JUNC_BED = GTF_BED_CONFIG if GTF_BED_CONFIG else out_path("references/annotation.junc.bed")

def parse_alignments():
    configured = config.get("alignments")
    if configured is None:
        configured = {"genome": {"enabled": True}}
    if isinstance(configured, list):
        configured = {name: {"enabled": True} for name in configured}
    if not isinstance(configured, dict) or not configured:
        raise ValueError("Configure at least one alignment under 'alignments'.")

    parsed = {}
    for name, settings in configured.items():
        settings = settings or {}
        if not settings.get("enabled", True):
            continue
        validate_token(name, "alignment name")
        if name == "genome":
            parsed[name] = {
                "reference_key": "genome",
                "reference": config.get("genome", ""),
                "juncbed": JUNC_BED,
                "mode": "genome",
                "preload_references": False,
            }
        elif name == "transcriptome":
            parsed[name] = {
                "reference_key": "transcriptome",
                "reference": config.get("transcriptome", ""),
                "juncbed": "",
                "mode": "transcriptome",
                "preload_references": True,
            }
        else:
            reference_key = settings.get("reference_key", name)
            mode = settings.get("mode", "transcriptome")
            if mode not in {"genome", "transcriptome"}:
                raise ValueError(
                    f"Alignment {name!r} has unsupported mode {mode!r}. Use 'genome' or 'transcriptome'."
                )
            parsed[name] = {
                "reference_key": reference_key,
                "reference": settings.get("reference", config.get(reference_key, "")),
                "juncbed": settings.get("juncbed", ""),
                "mode": mode,
                "preload_references": bool(settings.get("preload_references", False)),
            }

    if not parsed:
        raise ValueError("At least one alignment must be enabled in config/config.yaml.")
    return parsed

ALIGNMENTS = parse_alignments()
ALIGNMENT_NAMES = list(ALIGNMENTS.keys())

MODIFICATIONS = parse_modifications()
MOD_OUTPUT_SUFFIXES = [m["output_suffix"] for m in MODIFICATIONS]
MOD_SPLIT_OUTPUTS = expand(out_path("bedMethyl/{{alignment}}/{{sample}}.{mod}.filtered.bed"), mod=MOD_OUTPUT_SUFFIXES)
MOD_TABLE = " ".join(
    f"{m['name']}|{m['code']}|{m['output_suffix']}"
    for m in MODIFICATIONS
)
MODKIT_MODIFIED_BASES = [
    f"{m['primary_base']}:{m['code']}"
    for m in MODIFICATIONS
]
MODKIT_MODIFIED_BASES_ARGS = "--modified-bases " + " ".join(MODKIT_MODIFIED_BASES)

MODKIT_CONFIG = config.get("modkit", {})
MODKIT_PILEUP_CONFIG = MODKIT_CONFIG.get("pileup", {})
MODKIT_FILTER_THRESHOLDS = config_list(MODKIT_PILEUP_CONFIG.get("filter_thresholds", ["0.7"]))
MODKIT_MOD_THRESHOLDS = config_list(MODKIT_PILEUP_CONFIG.get("mod_thresholds", []))
if not MODKIT_MOD_THRESHOLDS:
    MODKIT_MOD_THRESHOLDS = [
        f"{m['code']}:{m['mod_threshold']}"
        for m in MODIFICATIONS
        if m["mod_threshold"] is not None
    ]
MODKIT_FILTER_ARGS = " ".join(
    f"--filter-threshold {threshold}"
    for threshold in MODKIT_FILTER_THRESHOLDS
)
MODKIT_MOD_THRESHOLD_ARGS = " ".join(
    f"--mod-threshold {threshold}"
    for threshold in MODKIT_MOD_THRESHOLDS
)
MODKIT_EXTRA_ARGS = config_args(MODKIT_PILEUP_CONFIG.get("extra_args", []))
MODKIT_PILEUP_ARGS = " ".join(
    arg for arg in [MODKIT_FILTER_ARGS, MODKIT_MOD_THRESHOLD_ARGS, MODKIT_EXTRA_ARGS] if arg
)

SAMPLE_ROWS = parse_samples()
SAMPLE_METADATA = {row["sample"]: row for row in SAMPLE_ROWS}
SAMPLE_POD5_DIRS = {row["sample"]: row["pod5_dir"] for row in SAMPLE_ROWS}

def validate_config():
    errors = []
    required_keys = [
        "dorado_sif", "modkit_sif", "dorado_models_dir",
    ]
    for key in required_keys:
        if key not in config or config[key] in (None, ""):
            errors.append(f"Missing required config key: {key}")

    if not SAMPLE_ROWS:
        errors.append("Configure at least one sample in samplesheet.")

    for key in ["dorado_sif", "modkit_sif"]:
        if config.get(key):
            check_existing_path(errors, key, config[key], "file")

    for alignment, settings in ALIGNMENTS.items():
        reference = settings["reference"]
        reference_key = settings["reference_key"]
        if not reference:
            errors.append(f"Alignment {alignment!r}: missing reference config key: {reference_key}")
        else:
            check_existing_path(errors, reference_key, reference, "file")

    if "genome" in ALIGNMENTS:
        if GENERATE_GTF_BED:
            if not config.get("gtf"):
                errors.append("Genome alignment requires either 'gtf_bed' or 'gtf' in config/config.yaml.")
            else:
                check_existing_path(errors, "gtf", config["gtf"], "file")
        else:
            check_existing_path(errors, "gtf_bed", JUNC_BED, "file")

    seen_samples = set()
    for row in SAMPLE_ROWS:
        sample = row["sample"]
        pod5_dir = row["pod5_dir"]
        validate_token(sample, "sample name")
        if sample in seen_samples:
            errors.append(f"Duplicate sample name: {sample}")
            continue
        seen_samples.add(sample)
        p = resolve_path(pod5_dir)
        if not p.is_dir():
            errors.append(f"Sample {sample!r}: pod5 directory does not exist: {pod5_dir}")
            continue
        if not list(p.glob("*.pod5")):
            errors.append(f"Sample {sample!r}: no .pod5 files found in {pod5_dir}")

    if errors:
        raise ValueError("Config validation failed:\n  - " + "\n  - ".join(errors))

validate_config()


# =========================================================
# FILE DISCOVERY
# =========================================================
def collect_jobs():
    jobs = []
    for sample, path in SAMPLE_POD5_DIRS.items():
        p = resolve_path(path).resolve()
        for f in sorted(p.glob("*.pod5")):
            jobs.append({
                "sample"  : sample,
                "pod5"    : str(f),
                "basename": f.stem
            })
    return jobs

JOBS    = collect_jobs()
SAMPLES = [row["sample"] for row in SAMPLE_ROWS]

JOB_POD5 = {
    (j["sample"], j["basename"]): j["pod5"]
    for j in JOBS
}


# =========================================================
# FINAL OUTPUT
# =========================================================
rule all:
    input:
        expand(out_path("bedMethyl/{alignment}/{sample}.{mod}.filtered.bed"), alignment=ALIGNMENT_NAMES, sample=SAMPLES, mod=MOD_OUTPUT_SUFFIXES),
        expand(out_path("sample_probs/{alignment}/{sample}/probabilities.tsv"), alignment=ALIGNMENT_NAMES, sample=SAMPLES),
        expand(out_path("qc/{alignment}/{sample}.qc_summary.txt"), alignment=ALIGNMENT_NAMES, sample=SAMPLES),
        out_path("provenance/config.resolved.yaml"),
        out_path("provenance/run_metadata.txt")


# =========================================================
# OPTIONAL REFERENCE PREPARATION
# =========================================================
if "genome" in ALIGNMENTS and GENERATE_GTF_BED:
    rule gtf_to_juncbed:
        input:
            gtf=config["gtf"]

        output:
            bed=JUNC_BED

        log:
            log_path("references/gtf_to_juncbed.log")

        conda:
            "env/minimap.yaml"

        shell:
            """
            mkdir -p $(dirname {output.bed}) $(dirname {log})
            paftools.js gff2bed {input.gtf} > {output.bed} 2> {log}
            """


# =========================================================
# RUN PROVENANCE
# =========================================================
rule run_provenance:
    input:
        references=lambda wc: [settings["reference"] for settings in ALIGNMENTS.values()],
        juncbed=JUNC_BED if "genome" in ALIGNMENTS else []

    output:
        config_copy=out_path("provenance/config.resolved.yaml"),
        metadata=out_path("provenance/run_metadata.txt")

    log:
        logfile=log_path("provenance/run_metadata.log")

    run:
        import datetime
        import hashlib
        import subprocess
        import yaml

        Path(str(output.metadata)).parent.mkdir(parents=True, exist_ok=True)
        Path(str(log.logfile)).parent.mkdir(parents=True, exist_ok=True)

        with open(output.config_copy, "w") as fh:
            yaml.safe_dump(dict(config), fh, sort_keys=False)

        def sha256_file(path):
            h = hashlib.sha256()
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                    h.update(chunk)
            return h.hexdigest()

        def command_output(args):
            try:
                completed = subprocess.run(
                    args,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
                text = completed.stdout.strip()
                if completed.returncode != 0:
                    return f"exit {completed.returncode}: {text}"
                return text
            except Exception as exc:
                return f"unavailable: {exc}"

        bind_args = ["-B", BIND] if BIND else []
        dorado_version = command_output(
            [CONTAINER_ENGINE, "exec"] + bind_args + [DORADO_SIF, "dorado", "--version"]
        )
        modkit_version = command_output(
            [CONTAINER_ENGINE, "exec"] + bind_args + [MODKIT_SIF, "modkit", "--version"]
        )
        minimap2_version = command_output(["minimap2", "--version"])
        samtools_version = command_output(["samtools", "--version"])

        lines = [
            f"generated_at: {datetime.datetime.now().isoformat(timespec='seconds')}",
            f"outdir: {OUTDIR}",
            f"container_engine: {CONTAINER_ENGINE}",
            f"dorado_sif: {DORADO_SIF}",
            f"modkit_sif: {MODKIT_SIF}",
            f"dorado_version: {dorado_version}",
            f"modkit_version: {modkit_version}",
            f"minimap2_version: {minimap2_version.splitlines()[0] if minimap2_version else ''}",
            f"samtools_version: {samtools_version.splitlines()[0] if samtools_version else ''}",
            "dorado_models:",
        ]
        lines.extend(f"  - {model}" for model in DORADO_MODELS)
        lines.append("enabled_modifications:")
        lines.extend(
            f"  - name: {m['name']}, primary_base: {m['primary_base']}, code: {m['code']}, output_suffix: {m['output_suffix']}, mod_threshold: {m['mod_threshold']}"
            for m in MODIFICATIONS
        )
        lines.append("alignments:")
        for alignment, settings in ALIGNMENTS.items():
            lines.append(
                f"  - name: {alignment}, mode: {settings['mode']}, reference: {settings['reference']}, preload_references: {settings['preload_references']}"
            )
        lines.append("samples:")
        for sample in SAMPLES:
            sample_row = SAMPLE_METADATA[sample]
            metadata = ", ".join(
                f"{key}: {value}" for key, value in sample_row.items()
                if key not in {"sample", "pod5_dir"} and value not in (None, "")
            )
            suffix = f", {metadata}" if metadata else ""
            lines.append(f"  - sample: {sample}, pod5_dir: {sample_row['pod5_dir']}{suffix}")
        for reference in input.references:
            lines.append(f"reference: {reference}")
            lines.append(f"reference_sha256: {sha256_file(reference)}")
        if input.juncbed:
            lines.append(f"junction_bed: {input.juncbed}")
            lines.append(f"junction_bed_sha256: {sha256_file(input.juncbed)}")
        lines.append(f"resolved_config: {output.config_copy}")

        with open(output.metadata, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        with open(log.logfile, "w") as fh:
            fh.write("Wrote run provenance metadata.\n")


# =========================================================
# 0. DORADO MODELS DOWNLOAD
#
# Downloads all required models to dorado_models_dir.
# Skipped automatically by Snakemake if model directories
# already exist on disk.
# =========================================================
rule dorado_models_download:
    output:
        DORADO_MODEL_DIRS

    params:
        model_dir=DORADO_MODELS_DIR,
        models=" ".join(DORADO_MODELS)

    log:
        log_path("dorado_download/download.log")

    shell:
        """
        mkdir -p {params.model_dir} $(dirname {log})

        for model in {params.models}; do

            echo "Downloading $model..." >> {log}
            {CONTAINER_ENGINE} exec {SINGULARITY_BIND_ARG} {DORADO_SIF} \
                dorado download \
                    --model $model \
                    --directory {params.model_dir} 2>> {log}
        done

        echo "Done:" >> {log}
        ls {params.model_dir} >> {log}
        """


# =========================================================
# 1. DORADO BASECALL (1 POD5 → 1 unaligned BAM)
#
# --emit-moves:        required to preserve modification tags
# --estimate-poly-a:   polyA tail estimation
# --batchsize 64:      conservative, safe across GPU types
# --device cuda:0:     use GPU assigned by SGE
# --models-directory:  use pre-downloaded local models
# =========================================================
rule dorado_basecall:
    input:
        pod5=lambda wc: JOB_POD5[(wc.sample, wc.basename)],
        # ensures models are downloaded before basecalling starts
        models=DORADO_MODEL_DIR_PATHS

    output:
        bam=temp(out_path("basecalled/{sample}/{basename}.unaligned.bam"))

    log:
        log_path("dorado/{sample}_{basename}.log")

    params:
        model     =DORADO_MODEL_ARG,
        models_dir=DORADO_MODELS_DIR

    threads: 2

    resources:
        gpu=1

    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log})

        {CONTAINER_ENGINE} exec {SINGULARITY_BIND_ARG} --nv {DORADO_SIF} \
        dorado basecaller {params.model} {input.pod5} \
            --emit-moves \
            --estimate-poly-a \
            --batchsize 64 \
            --device cuda:0 \
            --models-directory {params.models_dir} \
        > {output.bam} 2>> {log}
        """


# =========================================================
# 2. MERGE UNALIGNED BAM PER SAMPLE
# =========================================================
rule merge_unaligned_bam:
    input:
        lambda wc: [
            out_path(f"basecalled/{j['sample']}/{j['basename']}.unaligned.bam")
            for j in JOBS if j["sample"] == wc.sample
        ]

    output:
        bam=out_path("basecalled/{sample}/{sample}.merged.unaligned.bam")

    log:
        log_path("merge_unaligned/{sample}.log")

    conda:
        "env/minimap.yaml"

    threads: 8

    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        samtools merge -f --threads {threads} {output.bam} {input} 2> {log}
        """


# =========================================================
# 3. MINIMAP2 ALIGNMENT
#
# Genome mode uses splice-aware alignment and a junction BED.
# Transcriptome mode uses map-ont against transcript FASTA and passes
# --preload-references to modkit pileup downstream.
#
# -y:             transfers MM/ML/pt tags from unaligned BAM
#                 WITHOUT -y modification tags are lost
# -L:             long CIGAR format (required by modkit)
# --secondary=no: primary alignments only
# --MD:           MD tag for mismatches (QC)
# =========================================================
rule minimap2_align:
    input:
        bam    =out_path("basecalled/{sample}/{sample}.merged.unaligned.bam"),
        reference=lambda wc: ALIGNMENTS[wc.alignment]["reference"],
        juncbed=lambda wc: ALIGNMENTS[wc.alignment]["juncbed"] if ALIGNMENTS[wc.alignment]["juncbed"] else []

    output:
        bam=out_path("bams/{alignment}/{sample}/{sample}.bam"),
        bai=out_path("bams/{alignment}/{sample}/{sample}.bam.bai")

    log:
        log_path("minimap2/{alignment}/{sample}.log")

    conda:
        "env/minimap.yaml"

    threads: 4

    params:
        mode=lambda wc: ALIGNMENTS[wc.alignment]["mode"]

    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log})

        if [ "{params.mode}" = "genome" ]; then
            samtools bam2fq --threads {threads} -T MM,ML,pt \
                {input.bam} 2>> {log} | \
            minimap2 \
                -ax splice \
                -uf \
                -G 500000 \
                -L \
                --secondary=no \
                --MD \
                -y \
                --junc-bed {input.juncbed} \
                -t {threads} \
                {input.reference} - 2>> {log} | \
            samtools sort --threads {threads} -o {output.bam} 2>> {log}
        elif [ "{params.mode}" = "transcriptome" ]; then
            samtools bam2fq --threads {threads} -T MM,ML,pt \
                {input.bam} 2>> {log} | \
            minimap2 \
                -ax map-ont \
                -k14 \
                -L \
                --secondary=no \
                --MD \
                -y \
                -t {threads} \
                {input.reference} - 2>> {log} | \
            samtools sort --threads {threads} -o {output.bam} 2>> {log}
        else
            echo "Unsupported alignment mode: {params.mode}" > {log}
            exit 1
        fi

        samtools index -@ {threads} {output.bam} 2>> {log}
        """


# =========================================================
# 4. MODKIT SAMPLE-PROBS (QC)
#
# Checks MM/ML tags are present and probability distributions
# look sensible before running the full pileup.
# Inspect counts.html and proportion.html in a browser.
# =========================================================
rule sample_probs:
    input:
        bam=out_path("bams/{alignment}/{sample}/{sample}.bam"),
        bai=out_path("bams/{alignment}/{sample}/{sample}.bam.bai")

    output:
        probs      =out_path("sample_probs/{alignment}/{sample}/probabilities.tsv"),
        thresholds =out_path("sample_probs/{alignment}/{sample}/thresholds.tsv"),
        counts_html=out_path("sample_probs/{alignment}/{sample}/counts.html"),
        prop_html  =out_path("sample_probs/{alignment}/{sample}/proportion.html")

    log:
        log_path("sample_probs/{alignment}/{sample}.log")

    params:
        outdir=out_path("sample_probs/{alignment}/{sample}")

    threads: 4

    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})

        {CONTAINER_ENGINE} exec {SINGULARITY_BIND_ARG} {MODKIT_SIF} \
        modkit sample-probs \
            --hist \
            --threads {threads} \
            --out-dir {params.outdir} \
            {input.bam} \
        2> {log}
        """


# =========================================================
# 5. MODKIT PILEUP
#
# Runs on the merged aligned BAM (no strand separation).
# Strand info is preserved in col 6 of the output BED.
#
# For Modkit v0.6.0+, pileup is run with --modified-bases for
# clearer output and better performance. The argument is generated
# from modifications.*.primary_base and modifications.*.code.
#
# Pileup thresholds and extra arguments are configured under:
#   modkit:
#     pileup:
#       filter_thresholds
#       mod_thresholds
#       extra_args
# =========================================================
rule modkit_pileup:
    input:
        bam=out_path("bams/{alignment}/{sample}/{sample}.bam"),
        bai=out_path("bams/{alignment}/{sample}/{sample}.bam.bai")

    output:
        bed=out_path("modkit/{alignment}/{sample}.raw.bed")

    log:
        log_path("modkit/{alignment}/{sample}.log")

    threads: 8

    params:
        reference=lambda wc: ALIGNMENTS[wc.alignment]["reference"],
        modified_bases_args=MODKIT_MODIFIED_BASES_ARGS,
        preload_references=lambda wc: "--preload-references" if ALIGNMENTS[wc.alignment]["preload_references"] else "",
        pileup_args=MODKIT_PILEUP_ARGS

    shell:
        """
        mkdir -p $(dirname {output.bed}) $(dirname {log})

        {CONTAINER_ENGINE} exec {SINGULARITY_BIND_ARG} {MODKIT_SIF} \
        modkit pileup \
            -t {threads} \
            --reference {params.reference} \
            {params.preload_references} \
            {params.modified_bases_args} \
            {params.pileup_args} \
            {input.bam} {output.bed} \
            --log-filepath {log}
        """


# =========================================================
# 6. MODKIT SUMMARY (QC)
# =========================================================
rule modkit_summary:
    input:
        bam=out_path("bams/{alignment}/{sample}/{sample}.bam"),
        bai=out_path("bams/{alignment}/{sample}/{sample}.bam.bai")

    output:
        tsv=out_path("modkit/{alignment}/{sample}.summary.tsv")

    log:
        log_path("modkit/{alignment}/{sample}.summary.log")

    threads: 4

    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})

        {CONTAINER_ENGINE} exec {SINGULARITY_BIND_ARG} {MODKIT_SIF} \
        modkit summary \
            --threads {threads} \
            --tsv \
            {input.bam} > {output.tsv} 2> {log}
        """


# =========================================================
# 7. FILTERBED  <-- thresholds in config/config.yaml
#
# modkit pileup BED columns:
#   col 5  ($5)  = valid_coverage
#   col 11 ($11) = fraction_modified (0-100)
#
# Raw BED is always preserved. To re-filter with different
# thresholds without re-running pileup:
#   snakemake --forcerun filterbed
# =========================================================
rule filterbed:
    input:
        bed=out_path("modkit/{alignment}/{sample}.raw.bed")

    output:
        bed=out_path("bedMethyl/{alignment}/{sample}.filtered.bed")

    log:
        log_path("filterbed/{alignment}/{sample}.log")

    wildcard_constraints:
        alignment="|".join(re.escape(a) for a in ALIGNMENT_NAMES),
        sample="|".join(re.escape(s) for s in SAMPLES)

    params:
        min_coverage=config["min_coverage"],
        mod_pct     =config["mod_pct"]

    shell:
        """
        mkdir -p $(dirname {output.bed}) $(dirname {log})

        awk 'NR==1 || $1~/^#/ || ($5 >= {params.min_coverage} && $11 >= {params.mod_pct})' \
            {input.bed} > {output.bed} 2> {log}

        echo "Raw:      $(grep -vc '^#' {input.bed}  || echo 0)" >> {log}
        echo "Filtered: $(grep -vc '^#' {output.bed} || echo 0)" >> {log}
        """


# =========================================================
# 8. SPLITBED
#
# Splits filtered BED by enabled modification code using column 4
# of the modkit pileup BED output. If a modification is absent,
# the corresponding output BED is empty.
# =========================================================
rule splitbed:
    input:
        bed=out_path("bedMethyl/{alignment}/{sample}.filtered.bed")

    output:
        beds=MOD_SPLIT_OUTPUTS

    log:
        log_path("splitbed/{alignment}/{sample}.log")

    wildcard_constraints:
        alignment="|".join(re.escape(a) for a in ALIGNMENT_NAMES),
        sample="|".join(re.escape(s) for s in SAMPLES)

    params:
        mod_table=MOD_TABLE

    shell:
        """
        mkdir -p {OUTDIR}/bedMethyl/{wildcards.alignment} $(dirname {log})
        : > {log}

        for item in {params.mod_table}; do
            name="${{item%%|*}}"
            rest="${{item#*|}}"
            code="${{rest%%|*}}"
            suffix="${{rest##*|}}"
            outfile="{OUTDIR}/bedMethyl/{wildcards.alignment}/{wildcards.sample}.${{suffix}}.filtered.bed"

            awk -v code="${{code}}" 'BEGIN {{ FS=OFS="\t" }} $1 !~ /^#/ && $4 == code {{ print }}' \
                {input.bed} > "${{outfile}}" 2>> {log}

            echo "${{name}}: $(wc -l < "${{outfile}}")" >> {log}
        done
        """


# =========================================================
# 9. FINAL QC REPORT
# =========================================================
rule qc_report:
    input:
        bam_unaligned=out_path("basecalled/{sample}/{sample}.merged.unaligned.bam"),
        bam          =out_path("bams/{alignment}/{sample}/{sample}.bam"),
        bai          =out_path("bams/{alignment}/{sample}/{sample}.bam.bai"),
        summary      =out_path("modkit/{alignment}/{sample}.summary.tsv"),
        sample_probs =out_path("sample_probs/{alignment}/{sample}/probabilities.tsv"),
        raw          =out_path("modkit/{alignment}/{sample}.raw.bed"),
        filtered     =out_path("bedMethyl/{alignment}/{sample}.filtered.bed"),
        mod_beds     =MOD_SPLIT_OUTPUTS

    output:
        report=out_path("qc/{alignment}/{sample}.qc_summary.txt")

    log:
        log_path("qc/{alignment}/{sample}.log")

    conda:
        "env/minimap.yaml"

    threads: 2

    params:
        mod_table=MOD_TABLE

    shell:
        """
        mkdir -p $(dirname {output.report}) $(dirname {log})

        {{
        echo "===== QC SUMMARY — {wildcards.sample} ({wildcards.alignment}) ====="
        echo "Date: $(date)"
        echo ""

        echo "--- Total reads (unaligned BAM) ---"
        echo "Reads: $(samtools view -c -@ {threads} {input.bam_unaligned})"
        echo ""

        echo "--- Alignment {wildcards.alignment} (flagstat) ---"
        samtools flagstat --threads {threads} {input.bam}
        echo ""

        echo "--- Reads per chromosome (idxstats, top 25) ---"
        samtools idxstats {input.bam} | sort -k3 -rn | head -25 || true
        echo ""

        echo "--- Modification probabilities (sample-probs) ---"
        cat {input.sample_probs}
        echo ""

        echo "--- Modified sites ---"
        echo "Raw:      $(grep -vc '^#' {input.raw}      || true)"
        echo "Filtered: $(grep -vc '^#' {input.filtered} || true)"
        echo ""

        echo "--- Sites per modification type ---"
        for item in {params.mod_table}; do
            name="${{item%%|*}}"
            rest="${{item#*|}}"
            suffix="${{rest##*|}}"
            bed="{OUTDIR}/bedMethyl/{wildcards.alignment}/{wildcards.sample}.${{suffix}}.filtered.bed"
            echo "${{name}}: $(wc -l < "${{bed}}")"
        done
        echo ""

        echo "--- modkit summary ---"
        cat {input.summary}
        }} > {output.report} 2> {log}
        """
