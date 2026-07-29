rule dorado_models_download:
    output:
        DORADO_MODEL_DIRS
    params:
        model_dir=lambda w, output: os.path.dirname(output[0]),
        models=" ".join(DORADO_MODELS)
    log:
        log_path("dorado_download/download.log")
    container:
        DORADO_CONTAINER
    resources:
        mem_mb=4000,
        runtime=60
    shell:
        """
        mkdir -p {params.model_dir} $(dirname {log})

        for model in {params.models}; do
            echo "Downloading $model..." >> {log}
            dorado download --model $model --directory {params.model_dir} 2>> {log}
        done

        echo "Done:" >> {log}
        ls {params.model_dir} >> {log}
        """


rule dorado_basecall:
    input:
        pod5=lambda wc: CHUNK_POD5[(wc.sample, wc.chunk)],
        models=DORADO_MODEL_DIR_PATHS
    output:
        bam=temp(out_path("basecalled/{sample}/chunks/{chunk}.unaligned.bam"))
    log:
        log_path("dorado/{sample}_{chunk}.log")
    params:
        model=DORADO_MODEL_ARG,
        models_dir=lambda w, input: os.path.dirname(input.models[0]),
        whole_sample="1" if WHOLE_SAMPLE_BASECALL else "0",
        pod5_dir=lambda w, input: os.path.dirname(input.pod5[0])
    threads: 4
    # A failed GPU basecall is usually deterministic (OOM, bad model/pod5), so
    # cap retries below the global restart-times to avoid burning GPU hours.
    retries: 1
    resources:
        mem_mb=32000,
        runtime=2400,
        gpu=1,
        sge_pe="gpu"
    container:
        DORADO_CONTAINER
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log})

        # dorado basecaller takes a single reads path (file or directory). For a
        # whole-sample job pass the pod5 directory directly; for a chunk, stage
        # this chunk's pod5 files as symlinks in a temp directory and pass that.
        STAGE=""
        cleanup() {{ if [ -n "$STAGE" ]; then rm -rf "$STAGE"; fi }}
        trap cleanup EXIT

        if [ "{params.whole_sample}" = "1" ]; then
            POD5_INPUT="{params.pod5_dir}"
        else
            STAGE=$(mktemp -d)
            for f in {input.pod5}; do
                ln -s "$f" "$STAGE/"
            done
            POD5_INPUT="$STAGE"
        fi

        dorado basecaller {params.model} "$POD5_INPUT" \
            --estimate-poly-a \
            --device cuda:all \
            --models-directory {params.models_dir} \
        > {output.bam} 2>> {log}
        """


rule merge_unaligned_bam:
    input:
        lambda wc: [
            out_path(f"basecalled/{wc.sample}/chunks/{chunk}.unaligned.bam")
            for chunk in SAMPLE_CHUNKS[wc.sample]
        ]
    output:
        bam=out_path("basecalled/{sample}/{sample}.merged.unaligned.bam")
    log:
        log_path("merge_unaligned/{sample}.log")
    conda:
        workflow_path("workflow/envs/minimap.yaml")
    threads: 2
    resources:
        mem_mb=8000,
        runtime=120,
        sge_pe="smp"
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        # Chunks share an identical header, so concatenation (cat) is faster than
        # a header-merging samtools merge.
        samtools cat -o {output.bam} {input} 2> {log}
        """
