# Script to examine heterozygous markers in dosage (.raw) form data
# 02/2022 Jonna Clancy
######################################################################################
library(tidyverse)
library(data.table)
######################################################################################

# prerequisite: 
# 1) genotype data in PLINK format your_data_MHC, created in script HLA_imputation.R
# 2) text file containing list of individuals homozygote for certain HLA haplotype
# which heterozygosity you want to examine in column "IID"

#######################################################################################
# prepare your data (prerequisite: PLINK 1.9 installed in Linux)
# from binary to dosage format 0/1/2 (.raw)

system(command = paste0("plink --bfile your_data_MHC --recodeA --out your_data_MHC"))

#######################################################################################
# bring in dosage (.raw) form data
snp_data <- fread("your_data_MHC.raw",
                        data.table = FALSE)

#remove unnecessary columns
snp_data <- snp_data[ ,-c(1,3:6)]

# bring in ID list of individuals of certain HLA haplotype which heterotsygosity will be examined 
ID_list <- read.table(file="your_list_of_samples.txt",
                               header = TRUE, sep="\t", na = NA, dec=".", strip.white=TRUE,
                               stringsAsFactors = FALSE)


# extract the right individuals from the dataset
Heterozygosity_MHC <- right_join(snp_data, ID_list, by=c("IID"="IID"))

# transpose the data frame
Heterozygosity_MHC <- as.data.frame(t(Heterozygosity_MHC))

# find the amount of heterozygote markers of each individual
Heterozygosity_MHC <- as.data.frame(apply(Heterozygosity_MHC, 2, function(x) sum(x==1)))

#change the column name
colnames(Heterozygosity_MHC)<- "Haplotype 1"     # name the column e.g by HLA haplotype

#######################################################################################
# repeat to all wanted HLA haplotypes
# calculate the percentage of heterozygote markers
#######################################################################################