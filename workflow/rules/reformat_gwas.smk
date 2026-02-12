# Combine GWAS subset with lifted positions
rule reformat_sumstat:
    input:
        gwas  = rules.subset_gwas.output.sumstat,
        pos38 = rules.liftover_vcf.output.pos38,
        freq  = rules.compute_maf.output,
    output:
        report = temp(ws_path("report/{locuseq}.report"))
    params:
        gwas  = lambda wildcards: get_ofilename(wildcards.locuseq),
        seqid = lambda wildcards: get_column(wildcards.locuseq, "SeqID"),
        locus = lambda wildcards: get_column(wildcards.locuseq, "locus"),
        snp   = lambda wildcards: get_column(wildcards.locuseq, "SNPID"),
        gene  = lambda wildcards: get_column(wildcards.locuseq, "HARMONIZED_GENE_NAME"),
        uniprot = lambda wildcards: get_column(wildcards.locuseq, "UniProt_ID"),
        target  = lambda wildcards: get_column(wildcards.locuseq, "Target_Name"),
        targetf = lambda wildcards: get_column(wildcards.locuseq, "Target_Full_Name"),
        fixed_n = config.get("constants").get("sample_size"),
    conda:
        "../envs/liftover.yml"
    resources:
        runtime=lambda wc, attempt: 60 + attempt * 30
    script:
        "../scripts/reformat_sumstat.R"
