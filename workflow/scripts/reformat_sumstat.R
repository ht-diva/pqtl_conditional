#!/usr/bin/Rscript


# For more info, please look at the issue #40 of 'pqtl_conditional' GitHub repo.

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(Rmpfr))


#----------------------------------------#
#----------       INPUTS        ---------
#----------------------------------------#

# Path to inputs, outputs, and parameters
path_sumstat <- snakemake@input[['gwas']]
path_pos38   <- snakemake@input[['pos38']]
path_freq    <- snakemake@input[['freq']]
path_report  <- snakemake@output[['report']] # report filename
path_ofile   <- snakemake@params[['gwas']] # final gwas filename
seq_name     <- snakemake@params[['seqid']]
loc_name     <- snakemake@params[['locus']]
snp_name     <- snakemake@params[['snp']]
gene_name    <- snakemake@params[['gene']]
uniprot_name <- snakemake@params[['uniprot']]
target_name  <- snakemake@params[['target']]
target_fname <- snakemake@params[['targetf']]
fixed_n      <- snakemake@params[['fixed_n']]


# Return index variant to its original shape
snp_name <- str_replace_all(snp_name, "_", ":")

#----------------------------------------#
#-----      Read GWAS Datasets      -----
#----------------------------------------#

# 1. Read input data sets from paths, including lifted positions 
# 2. Append variants positions in build 38 to marginal sumstat
# 3. Define additional columns requested by MVP in sumstat


# Read marginal GWAS results for the corresponding locus
df_sumstat <- fread(path_sumstat)
df_pos38   <- fread(path_pos38, header = F, col.names = c("CHROM_38", "POS_38", "ID_37"))
df_freq    <- fread(path_freq)


# Check if 'ALT_FREQS' column computed by PLINK2,
# is really the frequency of effect allele (EA) in GWAS
check_freq <- df_freq %>%
  dplyr::mutate(
    a1  = str_extract(ID, "([A-Z])+"),
    a2  = str_extract(ID, "([A-Z])+$"),
    match = ALT == a1
  )


# Number of SNPs in each input
repo <- tibble(
  seqid = seq_name,
  locus = loc_name,
  index = snp_name,
  n_sumstat = nrow(df_sumstat),
  n_liftover = nrow(df_pos38),
  n_freq = nrow(df_freq)
  )

# Save report
fwrite(repo, file = path_report, sep = "\t", row.names = F, quote = F)


# Stop here immediately if there is any allele mismatch
if (any(!check_freq$match)) {
  stop("Execution stopped: at least one variant has allele mismatch in Plink2 output.")
}


#----------------------------------------#
#-----      Recompute P-value       -----
#----------------------------------------#

# Handle NAs when computing p-value from beta & sd
safe_pnorm <- function(b, se, p=FALSE) {
  
  # Ensure the vectors are of the same length
  if(length(b) != length(se)) {
    stop("Beta and SE must be of the same length")
  }
  
  k  <- length(b)
  b  <- as.numeric(b)
  se <- as.numeric(se)
  
  # Initialize result vector with NA values
  result <- rep(NA, k)
  
  # Identify non-missing and non-zero indices
  i <- which(!is.na(b) & !is.na(se) & se != 0)
  
  # compute z-score and take absolute, raise digits with mpfr, apply pnorm for non-missing values
  z_score <- b[i] / se[i]
  z_mpfr <- Rmpfr::mpfr(- abs(z_score), 120)
  p_mpfr <- 2 * pnorm(z_mpfr)
  mlog10p <- - log10(p_mpfr)
  
  # print p-value in character format and mlog10p in numeric
  if(p==TRUE){
    # reformat to mpfr character, then to numeric (don't set digits for MLOG10P)
    mlog10p_mpfr <- Rmpfr::formatMpfr(p_mpfr, scientific = TRUE, digits = 6)
    result[i] <- mlog10p_mpfr
  } else {
    mlog10p_mpfr <- Rmpfr::formatMpfr(mlog10p, scientific = TRUE)
    result[i] <- as.numeric(mlog10p_mpfr)
  }
  
  return(result)
}


#----------------------------------------#
#-----         Merge Datasets       -----
#----------------------------------------#

df_mvp <- df_sumstat %>%
  # append position in build 38
  inner_join(df_pos38, join_by(SNPID == ID_37)) %>%
  # append allele frequencies
  left_join(
    df_freq %>% select(ID, ALT_FREQS),
    join_by(SNPID == ID)
    ) %>%
  # add requested columns by MVP
  dplyr::mutate(
    CHROM_38 = str_remove(CHROM_38, "chr"),
    p = safe_pnorm(BETA, SE, p = TRUE),
    MAF = ifelse(ALT_FREQS < 0.5, ALT_FREQS, 1 - ALT_FREQS), # compute MAF based on EAF from Plink2
    N = fixed_n,                 # fixed value for sample size
    BETA_cond = NA_character_,
    STDERR_cond = NA_character_,
    PVAL_cond = NA_character_,
    MLOG10P_cond = NA_character_,
    SeqID  = seq_name,
    LOCUS_37  = loc_name,
    TISSUE = 'WholeBlood',
    GENE_NAME = gene_name,
    UNIPROT = uniprot_name,
    PROTEIN_NAME = target_name,
    PROTEIN_LONG_NAME = target_fname,
    DATASET = 'INTERVAL_CHRIS_META'
  ) %>%
  # rename columns for consistency with MVP
  dplyr::select(
    CHROM_37 = `##CHR`,
    POS_37 = POS,
    CHROM_38,
    POS_38,
    SNP_37 = SNPID,
    EA,
    NEA,
    EAF = ALT_FREQS,
    MAF,
    BETA_uncond = BETA,
    STDERR_uncond = SE,
    PVAL_uncond = p,
    MLOG10P_uncond = MLOG10P,
    N,
    BETA_cond,
    STDERR_cond,
    PVAL_cond,
    MLOG10P_cond,
    SeqID, LOCUS_37,
    TISSUE, GENE_NAME, 
    UNIPROT, PROTEIN_NAME, PROTEIN_LONG_NAME,
    DATASET
  )


#----------------------#
# Create output directory
odir <- dirname(path_ofile)

if (!dir.exists(odir)) {
  dir.create(odir, recursive = TRUE)
}

# Save reformatted sumstat
data.table::fwrite(
  df_mvp,
  file = path_ofile,
  sep = "\t",
  row.names = F, quote = F
  )

