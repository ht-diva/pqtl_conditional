## COJO-Conditional analysis Workflow for Protein QTLs 

We started to reimplement this pipeline on protein QTLs in early April, 2024. This workflow was primarily developed in Next-Flow (NF) by the Genomics Team at Human Technopole. We adopted, deployed, and validated it in Snakemake (SMK).

### Locus Breaker

We incorporated **Locus Breaker (LB)** function (see [publication](https://pubmed.ncbi.nlm.nih.gov/37173304/)) in SMK and deployed it on our protein meta-analysis GWAS results in mid April 2024.


### COJO Conditional Analysis
Once running the pipeline, `rule run_cojo` will generate output files below:
- list of independent variants resulted from GCTA cojo-slct (TSV/CSV)
- conditional dataset for each independent signal resulted from GCTA cojo-cond (RDS)
- fine-mapping results using `coloc::coloc.ABF()` function, containing values such as l-ABF, posterior probabilities (PPI) for each variant (RDS)
- colocalization info table containing credible set variants (with cumulative PPI > 0.99) for each independent variant
- regional association plots

These outputs are going to be stored in `workspace_path` provided by the user in `config.yaml` and stored in such directory:
<workspace_path>/results/*/cojo/<seqid>


### Colocalization of Two Proteins
We performed colocalization (Giambartolomei et al., 2014) across the pQTL signals. To meet the fundamental assumption of colocalization of only one causal variant per locus, we used conditional datasets, thus performing one colocalization test per pair of independent SNPs in 2 overlapping loci. For each regional association and each target SNP, we identified a credible set as the set of variants with posterior inclusion probability (PIP) > 0.99 within the region. More precisely, using the conditional dataset, we computed Approximate Bayes Factors (ABF) with the ‘process.dataset’ function in the coloc v5.2.3 R package and calculated posterior probabilities by normalizing ABFs across variants. Variants were ranked, and those with a cumulative posterior probability exceeding 0.99 were included in the credible sets. Among XXX protein pairs with overlapping loci, XXX protein pairs sharing a credible set variant were then tested for colocalization using the ‘coloc.abf’ function. Colocalized pairs were identified when the posterior probability for hypothesis 4 assuming a shared causal variant for two proteins exceeded 0.80.


### New Features on Top of NF pipeline
We also incorporated new features such as exclusion of signals in *HLA* and *NLRP12* regions from the results and follow-up analyses, allowing user to decide through the configuration file.


### NOTE
This SMK pipeline which is designed for pQTLs project **does not** include munging and alignment of input GWAS summary files. Therefore, it is a MUST to have your GWAS results completely harmonized by your genotype data. Eg. variants IDs, reference/alternate (effect/other) alleles should be concordant across your input files. Our GWAS summary stats from REGENIE are already aligned with QC pipeline (adopted by GWASLab) developed by pQTL analysts team at Health Data Science Center.


### How to run the pipeline
You can use the default configuration file in `config/config.yaml`. Otherwise, prepare your configuration in `config/` folder. Then, make sure that `configfile` variable in `workflow/Snakefile` matches with your newly created config file name. Then, run the pipeline by typing below command in bash.

```bash
sbatch submit.sh
```

### Select what to run
Available options:
- locus_breaker
- cojo_conditional
- coloc

If you like to run only locus breaker and **skip running cojo-conditional analyses and colocalization** between your traits, set this parameter in `conf/config.yaml` as follows:
```
run: locus_breaker
```

If you want to run COJO-Conditional analyses and **skip running colocalization** between your traits, set this parameter in `conf/config.yaml` as follows:
```
run: cojo_conditional
```

If you desire to run all the above analyses including colocalization, set this parameter in `conf/config.yaml` as follows:
```
run: coloc
```

### Workflow example

<img src="dag.svg" alt="example workflow">
