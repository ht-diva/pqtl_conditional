
rule subset_gwas:
    input:
        gwas = lambda wildcards: get_gwas(wildcards.locuseq),
    output:
        sumstat = ws_path("gwas/{locuseq}.tsv"),
        snplist = temp(ws_path("gwas/{locuseq}.snps")),
    params:
        locus = lambda wildcards: get_column(wildcards.locuseq, "locus"),
    resources:
        runtime=lambda wc, attempt: 120 + attempt * 60,
    shell:
        """
        source /exchange/healthds/singularity_functions

        echo "Genomic region: {params.locus}"

        # take region bounaries from locus string
        chr=$(echo {params.locus} | cut -d'_' -f1)
        beg=$(echo {params.locus} | cut -d'_' -f2)
        end=$(echo {params.locus} | cut -d'_' -f3)
        
        # reformat locus to be readable by tabix
        region=${{chr}}:${{beg}}-${{end}}
        
        echo "region is: $region"

        tabix  {input.gwas} $region -h > {output.sumstat}
        tail -n+2 {output.sumstat} | cut -f3 > {output.snplist}
        """
