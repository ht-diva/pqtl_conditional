#!/bin/bash


########################################
# USER SETTINGS: REFERENCES & CHAIN FILES
########################################

vcf_input=${snakemake_input}
txt_out=${snakemake_output}
chain_37to38="${snakemake_params[chain_37to38]}"
chain_38to37="${snakemake_params[chain_38to37]}"
fasta_hg37=${snakemake_params[fasta_hg37]}
fasta_hg38=${snakemake_params[fasta_hg38]}

# Select build between GRCh37 or GRCh38
START_BUILD="${snakemake_params[start_build]}"
TARGET_BUILD="${snakemake_params[target_build]}"


########################################
# VALIDATION
########################################

START_BUILD=$(echo "$START_BUILD" | tr '[:lower:]' '[:upper:]')
TARGET_BUILD=$(echo "$TARGET_BUILD" | tr '[:lower:]' '[:upper:]')

if [[ "$START_BUILD" == "$TARGET_BUILD" ]]; then
    echo "ERROR: START_BUILD and TARGET_BUILD are the same."
    exit 1
fi

if [[ ! "$START_BUILD" =~ ^GRCH(37|38)$ ]] || [[ ! "$TARGET_BUILD" =~ ^GRCH(37|38)$ ]]; then
    echo "ERROR: Builds must be GRCh37 or GRCh38."
    exit 1
fi

########################################
# SELECT FILES
########################################

if [[ "$START_BUILD" == "GRCH38" && "$TARGET_BUILD" == "GRCH37" ]]; then
    fasta_src="$fasta_hg38"
    fasta_tgt="$fasta_hg37"
    chain_file="$chain_38to37"
elif [[ "$START_BUILD" == "GRCH37" && "$TARGET_BUILD" == "GRCH38" ]]; then
    fasta_src="$fasta_hg37"
    fasta_tgt="$fasta_hg38"
    chain_file="$chain_37to38"
fi

########################################
# OUTPUT FILES
########################################

base_name=$(basename "$vcf_input" .vcf)
idir=$(dirname "$vcf_input")
odir=$(dirname "$txt_out")

vcf_zip="${idir}/${base_name}.vcf.gz"
vcf_std="${idir}/${base_name}.std.vcf.gz"
vcf_out="${odir}/${base_name}.vcf"


########################################
# LOAD TOOLS
########################################

source /exchange/healthds/singularity_functions

########################################
# LIFTOVER
########################################

echo "Lifting over from $START_BUILD to $TARGET_BUILD"

bgzip -c "$vcf_input" > "$vcf_zip"
tabix -p vcf "$vcf_zip"

bcftools norm \
    -f "$fasta_src" \
    -c s \
    -Oz \
    -o "$vcf_std" \
    "$vcf_zip"

bcftools +liftover \
    --no-version \
    -Ou "$vcf_std" \
    -- \
    -s "$fasta_src" \
    -f "$fasta_tgt" \
    -c "$chain_file" \
    > "$vcf_out"

bcftools query -f '%CHROM\t%POS\t%ID\n' "$vcf_out" > "$txt_out"

rm "$vcf_zip" "$vcf_std" "$vcf_out"

echo "LiftOver completed successfully."
echo "TXT: $txt_out"
