#!/usr/bin/env bash
set -euo pipefail


yaml_value() {
    awk -F: -v key="$1" '$1 ~ "^[[:space:]]*" key "[[:space:]]*$" { value=$2; sub(/^[[:space:]]*/, "", value); sub(/[[:space:]]*$/, "", value); gsub(/^"|"$/, "", value); print value; exit }' config/config.yaml
}

OUTDIR=$(yaml_value outdir)
OUTDIR=${OUTDIR:-results}
FOLDER="${OUTDIR%/}/logs/cluster/$(date +"%Y%m%d%H%M")"
BIND=$(yaml_value singularity_bind)
CONTAINER_ENGINE=$(yaml_value container_engine)
CONTAINER_ENGINE=${CONTAINER_ENGINE:-singularity}

mkdir -p "${FOLDER}"

CONTAINER_ARGS=(--use-singularity)
if [ "${CONTAINER_ENGINE}" = "apptainer" ]; then
    CONTAINER_ARGS=(--use-apptainer)
fi

SINGULARITY_ARGS=()
if [ -n "${BIND}" ]; then
    SINGULARITY_ARGS=(--singularity-args "-B ${BIND}")
fi

snakemake \
    --profile workflow/profiles/sge \
    --configfile config/config.yaml \
    --sge-logdir "${FOLDER}" \
    "${CONTAINER_ARGS[@]}" \
    "${SINGULARITY_ARGS[@]}" \
    "$@"
