#!/usr/bin/env python3
import argparse
import json
from pathlib import Path


def split_bed_by_mod(input_bed, outputs_by_suffix, modifications, log_path):
    counts = {mod["name"]: 0 for mod in modifications}
    handles = {}
    try:
        for mod in modifications:
            suffix = mod["output_suffix"]
            output_path = Path(outputs_by_suffix[suffix])
            output_path.parent.mkdir(parents=True, exist_ok=True)
            handles[suffix] = output_path.open("w")

        with open(input_bed) as bed:
            for line in bed:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                for mod in modifications:
                    if fields[3] == str(mod["code"]):
                        handles[mod["output_suffix"]].write(line)
                        counts[mod["name"]] += 1
                        break
    finally:
        for handle in handles.values():
            handle.close()

    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as log:
        for mod in modifications:
            log.write(f"{mod['name']}: {counts[mod['name']]}\n")


def run_from_snakemake():
    modifications = list(snakemake.params.modifications)
    outputs_by_suffix = {
        mod["output_suffix"]: output
        for mod, output in zip(modifications, snakemake.output.beds)
    }
    split_bed_by_mod(
        snakemake.input.bed,
        outputs_by_suffix,
        modifications,
        snakemake.log[0],
    )


def run_from_cli():
    parser = argparse.ArgumentParser(description="Split a modkit bedMethyl file by column-4 modification code.")
    parser.add_argument("--input-bed", required=True)
    parser.add_argument("--modifications-json", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    modifications = json.loads(args.modifications_json)
    outputs_by_suffix = {
        mod["output_suffix"]: str(Path(args.output_dir) / f"{args.sample}.{mod['output_suffix']}.filtered.bed")
        for mod in modifications
    }
    split_bed_by_mod(args.input_bed, outputs_by_suffix, modifications, args.log)


if "snakemake" in globals():
    run_from_snakemake()
else:
    run_from_cli()
