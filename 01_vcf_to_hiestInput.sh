#!/bin/bash

#Stitch vcf = imputed genotype probabilities for full thrush hybrid dataset
#Pipeline to impute genotypes is here: https://github.com/kdelmore/swth_rnaseq/tree/main/wgs
STITCH=/scratch/user/sblain/stitch_feb2023.vcf.gz

#TAMU-specific modules to load - need bcftools and vcftools
module load  GCC/11.2.0 BCFtools/1.14

#make index file for stitch
tabix $STITCH

#select individuals of interest - adults and juveniles sampled in the hybrid zone plus subspecies ref panel
#only keep biallelic SNPs
#filter with more than 25% missing data and minor allele frequency cut-off of 0.05
bcftools view -S parent_hybrid_sample_list.txt -m2 -M2 -v snps $STITCH | \
	bcftools filter -e 'F_MISSING > 0.25 || MAF <= 0.05' \
	-o stitch_feb2023.f25.vcf

bgzip stitch_feb2023.f25.vcf
tabix stitch_feb2023.f25.vcf.gz

##############LD-pruning with HWE threshold

module load GCC/9.3.0 PLINK/1.9b5 #module changeover to run plink

plink --vcf stitch_feb2023.f25.vcf.gz \
 --allow-extra-chr --geno 0.25 --hwe 0.0001 midp --indep-pairwise 50 10 0.1 \
 --out stitch.f25.r01.hwe \
 --set-missing-var-ids @:# \
 --double-id

sed 's/:/\t/g' stitch.f25.r01.hwe.prune.in > stitch.f25.r01.hwe.prune.regions

module load  GCC/11.2.0 BCFtools/1.14 VCFtools/0.1.16 #module changeover for vcftools

#regions file format chr\tpos
bcftools view -R stitch.f25.r01.hwe.prune.regions stitch_feb2023.f25.vcf.gz -O z -o stitch_feb2023.f25.ldr01.hwe.vcf.gz

# estimate fst
vcftools --gzvcf stitch_feb2023.f25.ldr01.hwe.vcf.gz --weir-fst-pop inland.sample1 --weir-fst-pop coastal.sample1 --out stitch.ldr01.hwe.s1

#remove sites on the Z & W chromosomes
awk '{if($1!="super_scaffold_7")print$0}' stitch.ldr01.hwe.s1.weir.fst | awk '{if($1!="scaffold_26_arrow_ctg1")print$0}' | \
  awk '{if($1!="scaffold_29_arrow_ctg1")print$0}' | awk '{if($1!="super_scaffold_25_w")print$0}' > stitch.ldr01.hwe.s1.weir.fst.noZ

#limit fst to >0.94
awk '{if($3>0.94 && $3!="-nan")print$0}' stitch.ldr01.hwe.s1.weir.fst.noZ | \
  awk -F "\t" '{print $1, $2}' > stitch.ldr01.hwe.s1.fst94.noZ.pos #make positions file - 176 sites if hwe=0.001, 307 if hwe=0.0001

#get genotypes file (0, 1, 2) for hiest input  
vcftools --gzvcf stitch_feb2023.f25.ldr01.hwe.vcf.gz --positions stitch.ldr01.hwe.s1.fst94.noZ.pos --out stitch.ldr01.hwe.s1.fst94.noZ.hiest_format --012
