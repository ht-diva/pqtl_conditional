#!/usr/bin/Rscript


library(data.table)
library(dplyr)

#----------#
# Get log path from Snakemake, fallback if missing
log_file <- tryCatch(snakemake@log[[1]], error = function(e) "logs/collect_loci/default.log")

# Ensure the directory exists
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

# Redirect stdout and stderr to the log file
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")   # redirect stdout
sink(log_con, type = "message")  # redirect messages / stderr

start_time <- Sys.time()
start_time

#----------#
# taking variants file as input
args <- commandArgs(trailingOnly = TRUE)
loci_path <- snakemake@input
file_path <- snakemake@output[["ofile"]]
nlrp12 <- snakemake@params[["NLRP12"]]
build <- snakemake@params[["build"]]

# convert input to boolean
nlrp12 <- nlrp12 %in% c(TRUE, "yes", "true", "TRUE", "Yes", "1")

#--------------#
# Merge them
loci <- tibble(
  rbindlist(
    fill = TRUE,
    lapply(
      loci_path,
      function(x) fread(x, data.table=F, fill = TRUE)
      )
    )
  )


#-------------------------------------#
#        Filter MHC and NLRP12        #
#-------------------------------------#

# define region depending on the genomic build
if (build == "37") {
  # NLRP12 gene maps to 54,296,995-54,327,657 in GRCh37, but we use suggested positions by Adam () enlarged by +/-20kb.
  nlrp12.start <- 54300000
  nlrp12.end   <- 54360000
  hla.start <- 28477797
  hla.end   <- 33448354
  message("To filter MHC and NLRP12 regions, genomic positions are set in build ", build)
} else if (build == "38") {
  # Using liftover.broadinstitute.org resulted in: chr19:53816370-53836078, then expanded it for 20kb
  nlrp12.start <- 53796000
  nlrp12.end   <- 53856000
  # MHC region maps to chr6:28,510,120-33,480,577 in GRCh38 coordinates.
  hla.start <- 28510120
  hla.end   <- 33480577
  message("To filter MHC and NLRP12 regions, genomic positions are set in build ", build)
}


ex_mhc <- loci %>% 
  arrange(chr) %>%
  #filter(!is.na(chr)) %>%   # remove trait without significant signals
  filter(!(chr == 6 & !(end < hla.start | start > hla.end)))    # remove HLA region

# remove signals overlapping NLRP12 region
if (nlrp12) {
  ex_nlrp12 <- ex_mhc %>% filter(!(chr == 19 & (POS > nlrp12.start & POS < nlrp12.end)))
  message("Removing lead SNPs in NLRP12 region.")
} else {
  ex_nlrp12 <- ex_mhc
  message("Skipping filter on NLRP12 region.")
}

message(nrow(loci), " loci have been identified in total.")
message(nrow(loci) - nrow(ex_mhc), " loci have been removed after filtering MHC region.")
message(nrow(ex_mhc) - nrow(ex_nlrp12), " loci have been removed after filtering NLRP12 region.")
message(nrow(ex_nlrp12), " loci remained in the final list.")


#-------------------------------------#
# compute width of loci and categorize it
loci_final <- ex_nlrp12 %>%
  mutate(
    loci_width = end - start,
    loci_cat = case_when(
      loci_width == 0 ~ "1-SNP",
      loci_width > 0 & loci_width <= 100000 ~ "1bp - 100Kbp",
      loci_width > 100000 & loci_width <= 250000 ~ "100-250Kbp",
      loci_width > 250000 & loci_width <= 500000 ~ "250-500Kbp",
      loci_width > 500000 & loci_width <= 1000000 ~ "500Kbp-1Mbp",
      loci_width > 1000000 & loci_width <= 2000000 ~ "1-2Mbp",
      loci_width > 2000000 & loci_width <= 5000000 ~ "2-5Mbp",
      TRUE ~ ">5Mbp"
      )
    )

#--------------#
# save the joint results
write.csv(loci_final, file = file_path, quote = F, row.names = F)

#--------------#
# Report run time
end_time <- Sys.time()
elapsed_time <- end_time - start_time
end_time

message("Run time: ", round(as.numeric(elapsed_time, units="mins"), 3), " minutes\n")

# Reset sinks at the end
sink(type = "message")
sink(type = "output")
close(log_con)
