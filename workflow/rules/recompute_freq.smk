
rule compute_maf:
    input:
        rules.subset_gwas.output.snplist
    output:
        temp(ws_path("freq/{locuseq}.afreq"))
    params:
        locus = lambda wildcards: get_column(wildcards.locuseq, "locus"),
        ofile = lambda wildcards, output: output[0].replace(".afreq", ""),
        pgen  = config.get("path_geno"),
    resources:
        runtime=lambda wc, attempt: 120 + attempt * 60,
    shell:
        """
    source /exchange/healthds/singularity_functions
        
    # take chromosomal number from locus string
    chrom=$(echo {params.locus} | cut -d'_' -f1)

    plink2 \
      --pfile {params.pgen}"$chrom" \
      --extract {input[0]} \
      --freq \
      --out {params.ofile} \
      --memory 4000
        """
