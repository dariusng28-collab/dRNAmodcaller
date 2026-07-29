rule sample_probs:
    input:
        bam=out_path("bams/{alignment}/{sample}/{sample}.bam"),
        bai=out_path("bams/{alignment}/{sample}/{sample}.bam.bai")
    output:
        probs=out_path("sample_probs/{alignment}/{sample}/probabilities.tsv"),
        thresholds=out_path("sample_probs/{alignment}/{sample}/thresholds.tsv"),
        counts_html=out_path("sample_probs/{alignment}/{sample}/counts.html"),
        prop_html=out_path("sample_probs/{alignment}/{sample}/proportion.html")
    log:
        log_path("sample_probs/{alignment}/{sample}.log")
    params:
        outdir=lambda w, output: os.path.dirname(output.probs)
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
        sge_pe="smp"
    container:
        MODKIT_CONTAINER
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})

        modkit sample-probs \
            --hist \
            --threads {threads} \
            --out-dir {params.outdir} \
            {input.bam} \
        2> {log}
        """


rule modkit_pileup:
    input:
        bam=out_path("bams/{alignment}/{sample}/{sample}.bam"),
        bai=out_path("bams/{alignment}/{sample}/{sample}.bam.bai")
    output:
        bed=out_path("modkit/{alignment}/{sample}.raw.bed")
    log:
        log_path("modkit/{alignment}/{sample}.log")
    threads: 8
    resources:
        mem_mb=32000,
        runtime=1440,
        sge_pe="smp"
    params:
        reference=lambda wc: ALIGNMENTS[wc.alignment]["reference"],
        modified_bases_args=MODKIT_MODIFIED_BASES_ARGS,
        preload_references=lambda wc: "--preload-references" if ALIGNMENTS[wc.alignment]["preload_references"] else "",
        pileup_args=MODKIT_PILEUP_ARGS
    container:
        MODKIT_CONTAINER
    shell:
        """
        mkdir -p $(dirname {output.bed}) $(dirname {log})

        modkit pileup \
            -t {threads} \
            --reference {params.reference} \
            {params.preload_references} \
            {params.modified_bases_args} \
            {params.pileup_args} \
            {input.bam} {output.bed} \
            --log-filepath {log}
        """


rule modkit_summary:
    input:
        bam=out_path("bams/{alignment}/{sample}/{sample}.bam"),
        bai=out_path("bams/{alignment}/{sample}/{sample}.bam.bai")
    output:
        tsv=out_path("modkit/{alignment}/{sample}.summary.tsv")
    log:
        log_path("modkit/{alignment}/{sample}.summary.log")
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120,
        sge_pe="smp"
    container:
        MODKIT_CONTAINER
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})

        modkit summary \
            --threads {threads} \
            --tsv \
            {input.bam} > {output.tsv} 2> {log}
        """


rule filterbed:
    input:
        bed=out_path("modkit/{alignment}/{sample}.raw.bed")
    output:
        bed=out_path("bedMethyl/{alignment}/{sample}.filtered.bed")
    log:
        log_path("filterbed/{alignment}/{sample}.log")
    container:
        PYTHON_CONTAINER
    params:
        min_coverage=config["min_coverage"],
        mod_pct=config["mod_pct"]
    resources:
        mem_mb=4000,
        runtime=60
    shell:
        """
        mkdir -p $(dirname {output.bed}) $(dirname {log})

        # modkit pileup bedMethyl has no header row. Column 10 is Nvalid_cov
        # (valid read coverage) and column 11 is percent modified (0-100).
        awk '$1 ~ /^#/ || ($10 >= {params.min_coverage} && $11 >= {params.mod_pct})' \
            {input.bed} > {output.bed} 2> {log}

        echo "Raw:      $(grep -vc '^#' {input.bed}  || echo 0)" >> {log}
        echo "Filtered: $(grep -vc '^#' {output.bed} || echo 0)" >> {log}
        """


rule splitbed:
    input:
        bed=out_path("bedMethyl/{alignment}/{sample}.filtered.bed")
    output:
        beds=MOD_SPLIT_OUTPUTS
    log:
        log_path("splitbed/{alignment}/{sample}.log")
    container:
        PYTHON_CONTAINER
    params:
        modifications=MODIFICATIONS
    resources:
        mem_mb=4000,
        runtime=60
    script:
        workflow_path("workflow/scripts/split_bed_by_mod.py")
