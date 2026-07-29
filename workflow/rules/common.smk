import csv
import os
import re
from pathlib import Path

from snakemake.utils import validate


REPO_DIR = Path.cwd()
OUTDIR = config.get("outdir", "results").rstrip("/")
LOGDIR = f"{OUTDIR}/logs"

CONTAINER_ENGINE = config.get("container_engine", "singularity")
BIND = config.get("singularity_bind", "")
DORADO_CONTAINER = config.get("dorado_container", config.get("dorado_sif", ""))
MODKIT_CONTAINER = config.get("modkit_container", config.get("modkit_sif", ""))

# Public container images for the generic tools. Overridable via config; the
# defaults are pinned, widely-mirrored images so every rule is fully
# containerised (no conda required).
MINIMAP_CONTAINER = config.get(
    "minimap_container",
    "docker://quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:365b17b986c1a60c1b82c6066a9345f38317b763-0",
)
NANOPLOT_CONTAINER = config.get(
    "nanoplot_container",
    "docker://quay.io/biocontainers/nanoplot:1.44.1--pyhdfd78af_0",
)
MULTIQC_CONTAINER = config.get(
    "multiqc_container",
    "docker://multiqc/multiqc:v1.25",
)
PYTHON_CONTAINER = config.get(
    "python_container",
    "docker://python:3.12-slim",
)


def out_path(path):
    return f"{OUTDIR}/{path}"


def log_path(path):
    return f"{LOGDIR}/{path}"


def workflow_path(path):
    return str(REPO_DIR / path)


def resolve_path(path):
    p = Path(str(path)).expanduser()
    if not p.is_absolute():
        p = REPO_DIR / p
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


def is_local_container(path):
    value = str(path or "")
    return value and "://" not in value


def check_existing_path(errors, key, path, kind="file"):
    p = resolve_path(path)
    if kind == "dir" and not p.is_dir():
        errors.append(f"{key}: directory does not exist: {path}")
    elif kind == "file" and not p.is_file():
        errors.append(f"{key}: file does not exist: {path}")


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
            "mod_threshold": settings.get("mod_threshold", settings.get("threshold")),
        })

    if not parsed:
        raise ValueError("At least one modification must be enabled in config/config.yaml.")
    return parsed


SAMPLE_COLUMN_ALIASES = ("sample", "sample_id", "sampleid", "id", "name")
POD5_DIR_COLUMN_ALIASES = (
    "pod5_dir",
    "pod5_directory",
    "pod5_path",
    "pod5",
    "raw_dir",
    "raw_path",
    "path",
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
DORADO_MODELS_DIR = config.get("dorado_models_dir", "")
# Number of pod5 files per basecalling job. 0 (default) basecalls the whole
# sample in a single dorado invocation, which minimises GPU model-reload
# overhead. Set to a positive integer to split large samples into chunks that
# can run on separate GPUs.
BASECALL_CHUNK_SIZE = int(config.get("basecall_chunk_size", 0) or 0)
WHOLE_SAMPLE_BASECALL = BASECALL_CHUNK_SIZE <= 0
DORADO_MODEL_ARG = ",".join(DORADO_MODELS)
DORADO_MODEL_DIR_PATHS = [f"{DORADO_MODELS_DIR}/{model}" for model in DORADO_MODELS]
DORADO_MODEL_DIRS = [directory(path) for path in DORADO_MODEL_DIR_PATHS]

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
                "index_preset": "-x splice",
                "preload_references": False,
            }
        elif name == "transcriptome":
            parsed[name] = {
                "reference_key": "transcriptome",
                "reference": config.get("transcriptome", ""),
                "juncbed": "",
                "mode": "transcriptome",
                "index_preset": "-x map-ont -k14",
                "preload_references": True,
            }
        else:
            reference_key = settings.get("reference_key", name)
            mode = settings.get("mode", "transcriptome")
            if mode not in {"genome", "transcriptome"}:
                raise ValueError(
                    f"Alignment {name!r} has unsupported mode {mode!r}. Use 'genome' or 'transcriptome'."
                )
            index_preset = "-x splice" if mode == "genome" else "-x map-ont -k14"
            parsed[name] = {
                "reference_key": reference_key,
                "reference": settings.get("reference", config.get(reference_key, "")),
                "juncbed": settings.get("juncbed", ""),
                "mode": mode,
                "index_preset": settings.get("minimap2_index_preset", index_preset),
                "preload_references": bool(settings.get("preload_references", False)),
            }

    if not parsed:
        raise ValueError("At least one alignment must be enabled in config/config.yaml.")
    return parsed


validate(config, workflow_path("workflow/schemas/config.schema.yaml"))

ALIGNMENTS = parse_alignments()
ALIGNMENT_NAMES = list(ALIGNMENTS.keys())

MODIFICATIONS = parse_modifications()
MOD_OUTPUT_SUFFIXES = [m["output_suffix"] for m in MODIFICATIONS]
MOD_SPLIT_OUTPUTS = expand(
    out_path("bedMethyl/{{alignment}}/{{sample}}.{mod}.filtered.bed"),
    mod=MOD_OUTPUT_SUFFIXES,
)
MOD_TABLE = " ".join(
    f"{m['name']}|{m['code']}|{m['output_suffix']}"
    for m in MODIFICATIONS
)
MODKIT_MODIFIED_BASES = [f"{m['primary_base']}:{m['code']}" for m in MODIFICATIONS]
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
    f"--filter-threshold {threshold}" for threshold in MODKIT_FILTER_THRESHOLDS
)
MODKIT_MOD_THRESHOLD_ARGS = " ".join(
    f"--mod-threshold {threshold}" for threshold in MODKIT_MOD_THRESHOLDS
)
MODKIT_EXTRA_ARGS = config_args(MODKIT_PILEUP_CONFIG.get("extra_args", []))
MODKIT_PILEUP_ARGS = " ".join(
    arg for arg in [MODKIT_FILTER_ARGS, MODKIT_MOD_THRESHOLD_ARGS, MODKIT_EXTRA_ARGS] if arg
)

SAMPLE_ROWS = parse_samples()
validate({"samples": SAMPLE_ROWS}, workflow_path("workflow/schemas/samples.schema.yaml"))
SAMPLE_METADATA = {row["sample"]: row for row in SAMPLE_ROWS}
SAMPLE_POD5_DIRS = {row["sample"]: row["pod5_dir"] for row in SAMPLE_ROWS}


def validate_config_semantics():
    errors = []
    required_keys = ["dorado_models_dir"]
    for key in required_keys:
        if key not in config or config[key] in (None, ""):
            errors.append(f"Missing required config key: {key}")

    if not DORADO_MODELS:
        errors.append("Configure at least one Dorado model with 'dorado_models'.")
    if not DORADO_CONTAINER:
        errors.append("Configure 'dorado_container' or legacy 'dorado_sif'.")
    if not MODKIT_CONTAINER:
        errors.append("Configure 'modkit_container' or legacy 'modkit_sif'.")
    if is_local_container(DORADO_CONTAINER):
        check_existing_path(errors, "dorado_container", DORADO_CONTAINER, "file")
    if is_local_container(MODKIT_CONTAINER):
        check_existing_path(errors, "modkit_container", MODKIT_CONTAINER, "file")
    if not SAMPLE_ROWS:
        errors.append("Configure at least one sample in samplesheet.")

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
                errors.append("Genome alignment requires either 'gtf_bed' or 'gtf'.")
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


validate_config_semantics()


def collect_chunks():
    """Group each sample's pod5 files into one or more basecalling chunks.

    Returns:
      chunk_pod5:      {(sample, chunk): [pod5 file paths]}
      sample_chunks:   {sample: [chunk names]}
      sample_pod5_dir: {sample: resolved pod5 directory}
    """
    chunk_pod5 = {}
    sample_chunks = {}
    sample_pod5_dir = {}
    for sample, path in SAMPLE_POD5_DIRS.items():
        p = resolve_path(path).resolve()
        sample_pod5_dir[sample] = str(p)
        files = [str(f) for f in sorted(p.glob("*.pod5"))]
        if not files:
            sample_chunks[sample] = []
            continue
        if WHOLE_SAMPLE_BASECALL:
            groups = [files]
        else:
            groups = [
                files[i:i + BASECALL_CHUNK_SIZE]
                for i in range(0, len(files), BASECALL_CHUNK_SIZE)
            ]
        width = max(3, len(str(len(groups) - 1)))
        names = []
        for idx, group in enumerate(groups):
            name = f"chunk{idx:0{width}d}"
            chunk_pod5[(sample, name)] = group
            names.append(name)
        sample_chunks[sample] = names
    return chunk_pod5, sample_chunks, sample_pod5_dir


CHUNK_POD5, SAMPLE_CHUNKS, SAMPLE_POD5_RESOLVED = collect_chunks()
SAMPLES = [row["sample"] for row in SAMPLE_ROWS]

wildcard_constraints:
    alignment="|".join(re.escape(a) for a in ALIGNMENT_NAMES),
    sample="|".join(re.escape(s) for s in SAMPLES),
    chunk=r"chunk\d+"
