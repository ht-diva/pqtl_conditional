from pathlib import Path
import pandas as pd


# read loci list
lb = pd.read_csv(config["loci_file"], sep = '\t')

# SNP id with underscores
lb["snp"] = lb["SNPID"].str.replace(":", "_", regex=False)

# Remove 'chr' from locus
lb["locus"] = lb["locus_START_END_37"].str.replace('chr', '', regex=False)

# Create seqid_locus identifier as wildcard
lb["locuseq"] = lb["SeqID"].astype(str) + "_" + lb["locus"]


# Use only needed columns
my_lb = (
    lb
    .drop_duplicates()
    .set_index("locuseq", drop=False)
    .sort_index()
)


# define the functions generating files' path
def ws_path(file_path):
    return str(Path(config.get("workspace_path"), file_path))

# return features of each locus
def get_column(wildcards, column):
    return str(my_lb.loc[wildcards, column])

def get_gwas(wildcards):
    seqid = my_lb.loc[wildcards, "SeqID"]
    file_path = f"{seqid}/{seqid}.gwaslab.tsv.bgz"
    return str(Path(config.get("path_gwas"), file_path))

def get_ofilename(wildcards):
    seqid = my_lb.loc[wildcards, "SeqID"]
    locus = my_lb.loc[wildcards, "locus"]
    snp   = my_lb.loc[wildcards, "snp"]
    file_path = f"toMVP/marginal_data_{seqid}_locus_{locus}_target_{snp}.tsv"
    return ws_path(file_path)
