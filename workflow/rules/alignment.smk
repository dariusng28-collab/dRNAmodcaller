rule minimap2_index:
    input:
        reference=lambda wc: ALIGNMENTS[wc.alignment]["reference"]
    output:
        mmi=out_path("references/index/{alignment}.mmi")
    log:
        log_path("minimap2_index/{alignment}.log")
    container:
        MINIMAP_CONTAINER
    threads: 4
    resources:
        mem_mb=32000,
        runtime=120,
        sge_pe="smp"
    params:
        preset=lambda wc: ALIGNMENTS[wc.alignment]["index_preset"]
    shell:
        """
        mkdir -p $(dirname {output.mmi}) $(dirname {log})
        minimap2 -t {threads} {params.preset} -d {output.mmi} {input.reference} 2> {log}
        """


rule minimap2_align:
    input:
        bam=out_path("basecalled/{sample}/{sample}.merged.unaligned.bam"),
        reference=out_path("references/index/{alignment}.mmi"),
        juncbed=lambda wc: ALIGNMENTS[wc.alignment]["juncbed"] if ALIGNMENTS[wc.alignment]["juncbed"] else []
    output:
        bam=out_path("bams/{alignment}/{sample}/{sample}.bam"),
        bai=out_path("bams/{alignment}/{sample}/{sample}.bam.bai")
    log:
        log_path("minimap2/{alignment}/{sample}.log")
    container:
        MINIMAP_CONTAINER
    threads: 4
    resources:
        mem_mb=48000,
        runtime=720,
        sge_pe="smp"
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
                -y \
                --junc-bed {input.juncbed} \
                -t {threads} \
                {input.reference} - 2>> {log} | \
            samtools sort --threads {threads} -m 2G -o {output.bam} 2>> {log}
        elif [ "{params.mode}" = "transcriptome" ]; then
            samtools bam2fq --threads {threads} -T MM,ML,pt \
                {input.bam} 2>> {log} | \
            minimap2 \
                -ax map-ont \
                -k14 \
                -L \
                --secondary=no \
                -y \
                -t {threads} \
                {input.reference} - 2>> {log} | \
            samtools sort --threads {threads} -m 2G -o {output.bam} 2>> {log}
        else
            echo "Unsupported alignment mode: {params.mode}" > {log}
            exit 1
        fi

        samtools index -@ {threads} {output.bam} 2>> {log}
        """
