# thrush_hybrid_survival

R and R markdown scripts associated with manuscript "Reduced hybrid survival in a migratory divide between songbirds"

### 00_functions_survival.R

Functions needed for downstream analysis in "04_phenoytpes.Rmd" and "05_survival.Rmd"

### 01_phenotypes_cleaning.R

Clean raw phenotype data
Produce metadata for adults and juveniles

### 02_latitude_cleaning.R

Clean and visualize raw motus detection and archival tag retrieval data

### 03_phenotype_network.R

Build correlation network of migratory phenotypes

### 04_phenoytpes.Rmd

Phenotype analysis - test for selection against intermediate phenotypes and phenotypic mismatch

### 05_survival.Rmd

Survival analysis - test for selection based on genomic class (ancestry and heterozygosity)

### 06_simulate_juvies_to2022.R

Simulate expected juvenile genomic classes based on a random mating among birds that survive migration and subspecies migrants

### 07_plot_juvie_sims.Rmd

Visualize output from simulated juvies

### 08_subspecies_map.R

Map subspecies ranges and motus tracks/towers