#!/usr/bin/env bash
set -euo pipefail


yaml_value() {
    # Return the value of a top-level key, preserving colons in the value
    # (e.g. bind mounts like "/a:/a,/b:/b" must not be truncated).
    awk -v key="$1" '
        $0 ~ "^[[:space:]]*" key "[[:space:]]*:" {
            sub("^[[:space:]]*" key "[[:space:]]*:[[:space:]]*", "")
            sub(/[[:space:]]*$/, "")
            gsub(/^"|"$/, "")
            print
            exit
        }' config/config.yaml
}

OUTDIR=$(yaml_value outdir)
OUTDIR=${OUTDIR:-results}
FOLDER="${OUTDIR%/}/logs/cluster/$(date +"%Y%m%d%H%M")"
BIND=$(yaml_value singularity_bind)
CONTAINER_ENGINE=$(yaml_value container_engine)
CONTAINER_ENGINE=${CONTAINER_ENGINE:-singularity}

# ---------------------------------------------------------------------------
# Singularity/Apptainer image hygiene (matters a lot on shared-FS clusters).
#   * Build SIFs on fast node-local disk. Converting conda-based images to SIF
#     touches tens of thousands of tiny files; doing that on network storage
#     (or a site default like /scratch0/tmp that is absent on the login node)
#     is glacially slow or fails outright.
#   * Store the built images in a persistent, shared prefix so they are pulled
#     and converted ONCE and reused across runs. Compute nodes read them over
#     the shared filesystem.
# Both are overridable in config.yaml via `singularity_tmpdir` (default /tmp,
# node-local) and `singularity_prefix` (default <outdir>/singularity).
# ---------------------------------------------------------------------------
SING_TMPDIR=$(yaml_value singularity_tmpdir)
SING_TMPDIR=${SING_TMPDIR:-/tmp}
SING_PREFIX=$(yaml_value singularity_prefix)
SING_PREFIX=${SING_PREFIX:-${OUTDIR%/}/singularity}

mkdir -p "${FOLDER}" "${SING_PREFIX}" "${SING_PREFIX}/cache"

export SINGULARITY_TMPDIR="${SING_TMPDIR}" APPTAINER_TMPDIR="${SING_TMPDIR}"
export SINGULARITY_CACHEDIR="${SING_PREFIX}/cache" APPTAINER_CACHEDIR="${SING_PREFIX}/cache"

# Where the cluster-generic qsub wrapper (workflow/scripts/sge_submit.py) writes
# scheduler stdout/stderr for each job.
export SGE_LOGDIR="${FOLDER}"

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
    --apptainer-prefix "${SING_PREFIX}" \
    "${CONTAINER_ARGS[@]}" \
    "${SINGULARITY_ARGS[@]}" \
    "$@"
