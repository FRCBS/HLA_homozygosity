library(tidyverse)
library(data.table)

#Script to calculate allele frequencies
#Modified from the script created by Mark Christie on 3/1/2012

##########################################################################################
# the file should consists of two columns:
# allele1: chr "01:01"
# allele2: chr "03:01"
##########################################################################################

HLA_A <- read.table(file="your_HLA-A_allele_file",
                        header = TRUE, sep = " ", na = NA, 
                        dec = ".", strip.white=TRUE,
                        stringsAsFactors = FALSE)


# rename the Allele 1 and 2 colums
names(HLA_A)[names(HLA_A) == "allele1"] <- "HLA_A"
names(HLA_A)[names(HLA_A) == "allele2"] <- "HLA_A"

##########################################################
#now calculate the allele frequencies for HLA-A

#how many columns are there?    
L_HLA_A=ncol(HLA_A)

#find the starting column number for each locus     
locus_positions_HLA_A=(2*(unique(round((1:(L_HLA_A-2))/2)))+1)  

#locus names, from the header
lnames_HLA_A=colnames(HLA_A)

#create a null dataset to append allele freqs to
OUT_HLA_A=NULL

for (x in locus_positions_HLA_A) {                     #begin for loop, to calculate frequencies for each locus
  alleles_HLA_A=c(HLA_A[,x],HLA_A[,x+1])        #For example, combine columns 1 and 2 for locus 1 (two columns because they are diploid)
  alleles2_HLA_A=as.data.frame(table(alleles_HLA_A))             #count each allele at locus x
  missing_HLA_A=alleles2_HLA_A[which(alleles2_HLA_A[,1]==0),2]          #count missing data at locus x, entered as '0' in this dataset (not used further for simplicity)
  alleles3_HLA_A=alleles2_HLA_A[-which(alleles2_HLA_A[,1]==0),]          #remove missing data (otherwise 0 would be counted in total number of alleles)
  alleles4_HLA_A=cbind(alleles2_HLA_A,alleles2_HLA_A[,2]/sum(alleles2_HLA_A[,2])) #calculate frequencies
  output_HLA_A=cbind(x,lnames_HLA_A[x],alleles4_HLA_A)                      #combine x, locusname, and frequencies
  
  
  OUT_HLA_A <<- rbind(OUT_HLA_A,output_HLA_A)
  
}


colnames(OUT_HLA_A) <- c("Number","Locus","allele","count","frequency") #add column headers
Allfreqs_HLA_A=OUT_HLA_A[,-1]

write.table(Allfreqs_HLA_A, "your destination/Allelefreqs_HLA_A",
            quote = FALSE, row.names = FALSE)

#############################################################################################
                                #Repeat the above for all loci
#############################################################################################