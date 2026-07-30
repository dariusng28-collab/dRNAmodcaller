rule run_provenance:
    input:
        references=lambda wc: [settings["reference"] for settings in ALIGNMENTS.values()],
        juncbed=JUNC_BED if "genome" in ALIGNMENTS else []
    output:
        config_copy=out_path("provenance/config.resolved.yaml"),
        metadata=out_path("provenance/run_metadata.txt")
    log:
        logfile=log_path("provenance/run_metadata.log")
    container:
        PYTHON_CONTAINER
    params:
        config_data=dict(config),
        alignments=ALIGNMENTS,
        modifications=MODIFICATIONS,
        samples=SAMPLE_ROWS,
        dorado_models=DORADO_MODELS,
        outdir=lambda w, output: os.path.dirname(os.path.dirname(output.metadata)),
        container_engine=CONTAINER_ENGINE,
        bind=BIND,
        dorado_container=DORADO_CONTAINER,
        modkit_container=MODKIT_CONTAINER,
        minimap_container=MINIMAP_CONTAINER,
        nanoplot_container=NANOPLOT_CONTAINER,
        multiqc_container=MULTIQC_CONTAINER,
        python_container=PYTHON_CONTAINER
    resources:
        mem_mb=1000,
        runtime=30
    script:
        workflow_path("workflow/scripts/write_provenance.py")


rule qc_report:
    input:
        bam_unaligned=out_path("basecalled/{sample}/{sample}.merged.unaligned.bam"),
        bam=out_path("bams/{alignment}/{sample}/{sample}.bam"),
        bai=out_path("bams/{alignment}/{sample}/{sample}.bam.bai"),
        summary=out_path("modkit/{alignment}/{sample}.summary.tsv"),
        sample_probs=out_path("sample_probs/{alignment}/{sample}/probabilities.tsv"),
        raw=out_path("modkit/{alignment}/{sample}.raw.bed"),
        filtered=out_path("bedMethyl/{alignment}/{sample}.filtered.bed"),
        mod_beds=MOD_SPLIT_OUTPUTS
    output:
        report=out_path("qc/{alignment}/{sample}.qc_summary.txt")
    log:
        log_path("qc/{alignment}/{sample}.log")
    container:
        MINIMAP_CONTAINER
    threads: 2
    resources:
        mem_mb=8000,
        runtime=120,
        sge_pe="smp"
    params:
        mod_table=MOD_TABLE,
        bedmethyl_dir=lambda w, input: os.path.dirname(input.filtered),
        idxstats_label=lambda wc: "Reads per reference sequence" if wc.alignment == "transcriptome" else "Reads per chromosome"
    shell:
        """
        mkdir -p $(dirname {output.report}) $(dirname {log})

        {{
        echo "===== QC SUMMARY - {wildcards.sample} ({wildcards.alignment}) ====="
        echo "Date: $(date)"
        echo ""

        echo "--- Total reads (unaligned BAM) ---"
        echo "Reads: $(samtools view -c -@ {threads} {input.bam_unaligned})"
        echo ""

        echo "--- Poly(A) tail length (pt tag, estimated reads only) ---"
        samtools view -@ {threads} {input.bam_unaligned} \
            | grep -oP 'pt:i:\\K[0-9]+' \
            | awk 'BEGIN{{n=0;s=0;mn=1e18;mx=0}}
                   $1>0{{n++;s+=$1;if($1<mn)mn=$1;if($1>mx)mx=$1}}
                   END{{if(n>0) printf "Reads with estimate: %d\\nMean: %.1f\\nMin: %d\\nMax: %d\\n",n,s/n,mn,mx;
                        else print "Reads with estimate: 0"}}' || true
        echo ""

        echo "--- Alignment {wildcards.alignment} (flagstat) ---"
        samtools flagstat --threads {threads} {input.bam}
        echo ""

        echo "--- {params.idxstats_label} (idxstats, top 25) ---"
        samtools idxstats {input.bam} | sort -k3 -rn | head -25 || true
        echo ""

        echo "--- Modification probabilities (sample-probs) ---"
        cat {input.sample_probs}
        echo ""

        echo "--- Modified sites ---"
        echo "Raw:      $(grep -vc '^#' {input.raw}      || true)"
        echo "Filtered: $(grep -vc '^#' {input.filtered} || true)"
        echo ""

        echo "--- Sites per modification type ---"
        for item in {params.mod_table}; do
            name="${{item%%|*}}"
            rest="${{item#*|}}"
            suffix="${{rest##*|}}"
            bed="{params.bedmethyl_dir}/{wildcards.sample}.${{suffix}}.filtered.bed"
            echo "${{name}}: $(wc -l < "${{bed}}")"
        done
        echo ""

        echo "--- modkit summary ---"
        cat {input.summary}
        }} > {output.report} 2> {log}
        """


rule samtools_stats:
    input:
        bam=out_path("bams/{alignment}/{sample}/{sample}.bam"),
        bai=out_path("bams/{alignment}/{sample}/{sample}.bam.bai")
    output:
        stats=out_path("qc/samtools/{alignment}.{sample}.stats.txt"),
        flagstat=out_path("qc/samtools/{alignment}.{sample}.flagstat.txt"),
        idxstats=out_path("qc/samtools/{alignment}.{sample}.idxstats.txt")
    log:
        log_path("samtools_stats/{alignment}/{sample}.log")
    container:
        MINIMAP_CONTAINER
    threads: 2
    resources:
        mem_mb=4000,
        runtime=60,
        sge_pe="smp"
    shell:
        """
        mkdir -p $(dirname {output.stats}) $(dirname {log})
        samtools stats    --threads {threads} {input.bam} > {output.stats}    2> {log}
        samtools flagstat --threads {threads} {input.bam} > {output.flagstat} 2>> {log}
        samtools idxstats {input.bam}                     > {output.idxstats} 2>> {log}
        """


rule nanoplot:
    input:
        bam=out_path("basecalled/{sample}/{sample}.merged.unaligned.bam")
    output:
        stats=out_path("qc/nanoplot/{sample}/{sample}.NanoStats.txt")
    log:
        log_path("nanoplot/{sample}.log")
    container:
        NANOPLOT_CONTAINER
    threads: 4
    resources:
        mem_mb=16000,
        runtime=720,
        sge_pe="smp"
    params:
        outdir=lambda w, output: os.path.dirname(output.stats),
        prefix="{sample}.",
        downsample=config.get("nanoplot_downsample", 500000)
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        # Alignment-independent read QC (length N50, quality) from the raw
        # (unaligned) basecalled reads. dRNA BAMs can be tens of millions of
        # reads; QC distributions are unchanged by random-sampling a subset,
        # which keeps this fast instead of iterating the whole (multi-GB) file.
        NanoPlot \
            --ubam {input.bam} \
            --threads {threads} \
            --downsample {params.downsample} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            2> {log}
        """


rule multiqc:
    input:
        samtools=expand(
            out_path("qc/samtools/{alignment}.{sample}.stats.txt"),
            alignment=ALIGNMENT_NAMES,
            sample=SAMPLES,
        ),
        flagstat=expand(
            out_path("qc/samtools/{alignment}.{sample}.flagstat.txt"),
            alignment=ALIGNMENT_NAMES,
            sample=SAMPLES,
        ),
        idxstats=expand(
            out_path("qc/samtools/{alignment}.{sample}.idxstats.txt"),
            alignment=ALIGNMENT_NAMES,
            sample=SAMPLES,
        ),
        nanoplot=expand(
            out_path("qc/nanoplot/{sample}/{sample}.NanoStats.txt"),
            sample=SAMPLES,
        )
    output:
        report=out_path("qc/multiqc/multiqc_report.html")
    log:
        log_path("multiqc/multiqc.log")
    container:
        MULTIQC_CONTAINER
    resources:
        mem_mb=8000,
        runtime=60
    params:
        scan_dir=lambda w, output: os.path.dirname(os.path.dirname(output.report)),
        outdir=lambda w, output: os.path.dirname(output.report)
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        multiqc \
            --force \
            --outdir {params.outdir} \
            --filename multiqc_report.html \
            --ignore "*/multiqc/*" \
            {params.scan_dir} \
        2> {log}
        """
