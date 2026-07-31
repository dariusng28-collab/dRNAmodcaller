#!/usr/bin/env python3
"""Submit a Snakemake job to SGE for the cluster-generic executor.

Invoked by ``cluster-generic-submit-cmd`` as::

    sge_submit.py <rule> <threads> <mem_mb> <runtime_min> <sge_pe> <gpu> <jobscript>

Maps Snakemake resources onto this cluster's qsub conventions (UCL CS style:
per-slot ``tmem``/``h_vmem`` memory, ``h_rt`` walltime, ``-pe smp|gpu``) and
prints the submitted job id (via ``qsub -terse``) for the status checker.
"""
import math
import os
import subprocess
import sys

rule = sys.argv[1]
threads = max(1, int(float(sys.argv[2])))
mem_mb = max(1, int(float(sys.argv[3])))
runtime = max(1, int(float(sys.argv[4])))  # minutes
sge_pe = sys.argv[5] or "smp"
gpu = int(float(sys.argv[6]))
jobscript = sys.argv[-1]

# tmem/h_vmem are enforced PER SLOT on this cluster, so divide the rule's total
# memory across the number of slots actually requested (threads for an smp job,
# the GPU count for a gpu job).
slots = gpu if gpu > 0 else threads
per_slot_mb = max(1, math.ceil(mem_mb / max(1, slots)))
hours, minutes = divmod(runtime, 60)
h_rt = f"{hours:d}:{minutes:02d}:00"

logdir = os.environ.get("SGE_LOGDIR", "results/logs/cluster")
os.makedirs(logdir, exist_ok=True)

cmd = [
    "qsub", "-cwd", "-V", "-S", "/bin/bash", "-j", "y", "-terse",
    "-N", f"smk.{rule}",
    "-o", logdir,
    "-l", f"tmem={per_slot_mb}M,h_vmem={per_slot_mb}M,h_rt={h_rt}",
]
if gpu > 0:
    cmd += ["-pe", "gpu", str(gpu), "-l", "gpu=true", "-R", "y"]
else:
    cmd += ["-pe", sge_pe, str(threads), "-R", "y"]
cmd.append(jobscript)

proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
if proc.returncode != 0:
    sys.stderr.write(proc.stdout)
    sys.stderr.write(proc.stderr)
    sys.exit(1)

# -terse prints just the job id (single job; cluster-generic does not use arrays).
print(proc.stdout.strip().splitlines()[-1].strip())
