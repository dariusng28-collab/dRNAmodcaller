#!/usr/bin/env bash
# Status check for the cluster-generic executor on SGE, WITHOUT qacct
# (this cluster's accounting file is not readable, which breaks the dedicated
# SGE plugin). Completion is inferred from qstat: gone from the queue = done,
# and Snakemake then verifies the job's output files actually exist.
set -u
jobid="$1"

# State of this job in the queue (column 5), if present.
state=$(qstat 2>/dev/null | awk -v id="$jobid" '$1 == id { print $5; exit }')

if [ -z "$state" ]; then
    echo "success"          # no longer in the queue -> finished
elif [[ "$state" == *E* ]]; then
    echo "failed"           # Eqw / error state
else
    echo "running"          # qw, r, t, hqw, etc.
fi
