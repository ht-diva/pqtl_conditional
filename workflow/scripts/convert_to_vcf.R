#!/usr/bin/Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))


#----------------------------------------#
#-----      Inputs and Outputs      -----
#----------------------------------------#

# taking variants file as input
path_gwas <- as.character(snakemake@input[['gwas']])
path_vcf  <- as.character(snakemake@output[['vcf']])


#----------------------------------------#
#-----       VCF Conversion       -----
#----------------------------------------#

# Function to write VCF file
write_vcf <- function(df, out_name) {
  
  # Create a connection to write to a file
  vcf_file <- file(out_name, "w")
  
  # Use tryCatch to ensure the file is closed properly
  tryCatch({
    # Write the VCF header
    writeLines("##fileformat=VCFv4.2", vcf_file)
    writeLines("##source=RScript", vcf_file)
    writeLines("##reference=GRCh38", vcf_file)
    writeLines("##INFO=<ID=EAF,Number=1,Type=Float,Description=Effect Allele Frequency>", vcf_file)
    writeLines("##INFO=<ID=BETA,Number=1,Type=Float,Description=Effect Size Estimate>", vcf_file)
    writeLines("##INFO=<ID=SE,Number=1,Type=Float,Description=Standard Error>", vcf_file)
    writeLines("##INFO=<ID=N,Number=1,Type=Integer,Description=Sample Size>", vcf_file)
    writeLines("##INFO=<ID=MLOG10P,Number=1,Type=Float,Description=Negative Log10 P-value>", vcf_file)
    writeLines("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", vcf_file)
    
    # Write each row as a VCF entry
    for (i in 1:nrow(df)) {
      # Extract chromosome, position, reference allele (NEA), and alternate allele (EA)
      chr <- df$CHR[i]
      pos <- df$POS[i]
      id <- df$SNPID[i]
      ref <- df$NEA[i]
      alt <- df$EA[i]
      qual <- "."
      filter <- "."
      info <- "."
         
      # Write the line to the VCF file
      line <- paste(chr, pos, id, ref, alt, qual, filter, info, sep = "\t")
      
      writeLines(line, vcf_file)
    }
  },
  finally = {
    close(vcf_file) # Ensure the file connection is closed
  })
}

#--------------------#
# Read GWAS subset
df <- fread(path_gwas)

# Remove unwanted columns in GWAS subset
df4vcf <- df %>%
  dplyr::select(
    CHR = `##CHR`,
    POS, SNPID, EA, NEA
    ) %>%
  arrange(CHR, POS) # sort SNPs for converting to VCF


# Create the directory plus all necessary subdirectories
odir <- dirname(path_vcf)

if (!dir.exists(odir)) {
  dir.create(odir, recursive = TRUE)
}

# Convert GWAS subset to VCF
write_vcf(df4vcf, path_vcf)

# return a message showing where the VCF is saved
cat("\nSave VCF here: ", path_vcf)
