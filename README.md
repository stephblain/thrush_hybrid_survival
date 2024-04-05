# thrush_hybrid_survival

R and R markdown scripts associated with manuscript "Reduced hybrid survival in a migratory divide between songbirds"

### 00_functions_survival.R

Functions needed for downstream analysis in "09_phenoytpe_figures.Rmd" and "11_survival_figures.R"

### 01_vcf_to_hiestInput.sh

Process vcf to get divergent SNPs and allele counts for ref panel and target hybrids

### 02_run_hiest.R

Some more data formatting, then run hiest to estimate ancestry and heterozygosity

### 03_phenotypes_cleaning.R

Clean raw phenotype data
Produce metadata for adults and juveniles

### 04_survival_cleaning.R

Clean and visualize raw motus detection and archival tag retrieval data

### 05_CJS_ancestryHeterozygosity.R

Survival analysis - test for selection based on genomic class (ancestry and heterozygosity) with CJS model in juveniles

### 06_MSCJS_ancestryHeterozygosity.R

Survival analysis - test for transitions between migratory states based on genomic class (ancestry and heterozygosity) with multistate CJS model in juveniles

### 07_CJS_phenotypes.R

Phenotype analysis - test for selection against intermediate phenotypes and phenotypic mismatch with CJS model in juveniles

### 08_GLM_adults.random

Survival analysis - test for selection based on genomic class (ancestry and heterozygosity) in adults

### 09_phenoytpes.Rmd

Make phenotype figures (fig 3)

### 10_subspecies_map.R

Map subspecies ranges and motus tracks/towers (fig 1)

### 11_survival_figures.R

Plot survival outcomes based on genomic class from CJS and GLM's for juveniles and adults