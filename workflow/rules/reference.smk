if "genome" in ALIGNMENTS and GENERATE_GTF_BED:
    rule gtf_to_juncbed:
        input:
            gtf=config["gtf"]
        output:
            bed=JUNC_BED
        log:
            log_path("references/gtf_to_juncbed.log")
        conda:
            workflow_path("workflow/envs/minimap.yaml")
        resources:
            mem_mb=4000,
            runtime=60
        shell:
            """
            mkdir -p $(dirname {output.bed}) $(dirname {log})
            paftools.js gff2bed {input.gtf} > {output.bed} 2> {log}
            """
