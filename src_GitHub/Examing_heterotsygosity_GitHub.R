# Script to examine heterozygous markers in dosage (.raw) form data
# 02/2022 Jonna Clancy
#######################
library(tidyverse)
library(data.table)
#######################
# First convert binary genotype data into dosage format 0/1/2 (.raw) in Linux

# #!/bin/bash
# plink \
# --bfile your data \
# --recodeA \
# --out your data

#######################################################################################
# bring in dosage (.raw) form data
snp_data <- fread("your data.raw",
                        data.table = FALSE)

# bring in ID list of individuals of certain HLA haplotype which heterotsygosity will be examined 
ID_list <- read.table(file="your data.txt",
                               header = TRUE, sep="\t", na = NA, dec=".", strip.white=TRUE,
                               stringsAsFactors = FALSE)

# extract the right individuals from the dataset
Heterozygosity_MHC <- right_join(snp_data, ID_list, by=c("IID"="IID"))

# change the order of the data frame
Heterozygosity_MHC <- as.data.frame(t(Heterozygosity_MHC))

# take the amount of heterozygote markers of each individual into a data frame
Heterozygosity_MHC <- as.data.frame(apply(Heterozygosity_MHC, 2, function(x) sum(x==1)))

#change the column name
colnames(FIN_most_common_MHC)<- "1"      # name the column e.g by HLA haplotype

####################################################################
# repeat to all wanted HLA haplotypes
####################################################################

