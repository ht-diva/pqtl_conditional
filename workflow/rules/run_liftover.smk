
rule convert_to_vcf:
    input:
        gwas = rules.subset_gwas.output.sumstat
    output:
        vcf = temp(ws_path("VCF/{locuseq}.vcf"))
    # log:
    #     ws_path("logs/VCF/{locuseq}.log")
    conda:
        "../envs/liftover.yml"
    resources:
        runtime=lambda wc, attempt: 120 + attempt * 60,
    script:
        "../scripts/convert_to_vcf.R"


# Lift over VCF files from GRCh37 to GRCh38
rule liftover_vcf:
    input:
        ws_path("VCF/{locuseq}.vcf"),
    output:
        pos38 = temp(ws_path("VCF_lifted/{locuseq}.txt"))
    # log:
    #     ws_path("logs/liftover/{locuseq}.log"),
    params:
        chain_37to38 = config.get("chain_37to38"),
        chain_38to37 = config.get("chain_38to37"),
        fasta_hg37   = config.get("fasta_hg37"),
        fasta_hg38   = config.get("fasta_hg38"),
        start_build  = config.get("start_build"),
        target_build = config.get("target_build"),
    script:
        "../scripts/liftover_vcf.sh"
