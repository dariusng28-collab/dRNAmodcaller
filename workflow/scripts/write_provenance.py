from pathlib import Path
import datetime
import hashlib
import subprocess

import yaml


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
            check=False,
        )
        text = completed.stdout.strip()
        if completed.returncode != 0:
            return f"exit {completed.returncode}: {text}"
        return text
    except Exception as exc:
        return f"unavailable: {exc}"


def container_command(container):
    bind = snakemake.params.bind
    engine = snakemake.params.container_engine
    bind_args = ["-B", bind] if bind else []
    return [engine, "exec"] + bind_args + [container]


Path(str(snakemake.output.metadata)).parent.mkdir(parents=True, exist_ok=True)
Path(str(snakemake.log.logfile)).parent.mkdir(parents=True, exist_ok=True)

with open(snakemake.output.config_copy, "w") as fh:
    yaml.safe_dump(dict(snakemake.params.config_data), fh, sort_keys=False)

dorado_version = command_output(
    container_command(snakemake.params.dorado_container) + ["dorado", "--version"]
)
modkit_version = command_output(
    container_command(snakemake.params.modkit_container) + ["modkit", "--version"]
)
minimap2_version = command_output(["minimap2", "--version"])
samtools_version = command_output(["samtools", "--version"])

lines = [
    f"generated_at: {datetime.datetime.now().isoformat(timespec='seconds')}",
    f"outdir: {snakemake.params.outdir}",
    f"container_engine: {snakemake.params.container_engine}",
    f"dorado_container: {snakemake.params.dorado_container}",
    f"modkit_container: {snakemake.params.modkit_container}",
    f"dorado_version: {dorado_version}",
    f"modkit_version: {modkit_version}",
    f"minimap2_version: {minimap2_version.splitlines()[0] if minimap2_version else ''}",
    f"samtools_version: {samtools_version.splitlines()[0] if samtools_version else ''}",
    "dorado_models:",
]
lines.extend(f"  - {model}" for model in snakemake.params.dorado_models)

lines.append("enabled_modifications:")
lines.extend(
    "  - name: {name}, primary_base: {primary_base}, code: {code}, "
    "output_suffix: {output_suffix}, mod_threshold: {mod_threshold}".format(**mod)
    for mod in snakemake.params.modifications
)

lines.append("alignments:")
for alignment, settings in snakemake.params.alignments.items():
    lines.append(
        f"  - name: {alignment}, mode: {settings['mode']}, "
        f"reference: {settings['reference']}, preload_references: {settings['preload_references']}"
    )

lines.append("samples:")
for sample_row in snakemake.params.samples:
    metadata = ", ".join(
        f"{key}: {value}"
        for key, value in sample_row.items()
        if key not in {"sample", "pod5_dir"} and value not in (None, "")
    )
    suffix = f", {metadata}" if metadata else ""
    lines.append(f"  - sample: {sample_row['sample']}, pod5_dir: {sample_row['pod5_dir']}{suffix}")

for reference in snakemake.input.references:
    lines.append(f"reference: {reference}")
    lines.append(f"reference_sha256: {sha256_file(reference)}")

if snakemake.input.juncbed:
    lines.append(f"junction_bed: {snakemake.input.juncbed}")
    lines.append(f"junction_bed_sha256: {sha256_file(snakemake.input.juncbed)}")

lines.append(f"resolved_config: {snakemake.output.config_copy}")

with open(snakemake.output.metadata, "w") as fh:
    fh.write("\n".join(lines) + "\n")

with open(snakemake.log.logfile, "w") as fh:
    fh.write("Wrote run provenance metadata.\n")
