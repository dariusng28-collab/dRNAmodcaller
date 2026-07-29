from pathlib import Path
import datetime
import hashlib
import json


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


Path(str(snakemake.output.metadata)).parent.mkdir(parents=True, exist_ok=True)
Path(str(snakemake.log.logfile)).parent.mkdir(parents=True, exist_ok=True)

# Written as JSON (a valid YAML subset) so the provenance rule needs no PyYAML,
# letting it run in a minimal python container.
with open(snakemake.output.config_copy, "w") as fh:
    json.dump(dict(snakemake.params.config_data), fh, indent=2, default=str)
    fh.write("\n")

# Every tool runs in a version-pinned container, so the image references are the
# authoritative software-version record for the run.
lines = [
    f"generated_at: {datetime.datetime.now().isoformat(timespec='seconds')}",
    f"outdir: {snakemake.params.outdir}",
    f"container_engine: {snakemake.params.container_engine}",
    "containers:",
    f"  dorado:            {snakemake.params.dorado_container}",
    f"  modkit:            {snakemake.params.modkit_container}",
    f"  minimap2_samtools: {snakemake.params.minimap_container}",
    f"  nanoplot:          {snakemake.params.nanoplot_container}",
    f"  multiqc:           {snakemake.params.multiqc_container}",
    f"  python:            {snakemake.params.python_container}",
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
