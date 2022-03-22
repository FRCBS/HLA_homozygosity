# HLA imputation
# script by Kati Hyvärinen and Jonna Clancy 2020
###############################################################################
# load required packages

library(HIBAG)
library(tidyverse)

#########################################
# import imputation reference files
# set variables into HLA-A

Fin_hg38_model_A <- "your HLA-A model in hg38"
model.list_A <- get(load(Fin_hg38_model_A))

# HLA imputation model at HLA-A
hla.id_A <- "A"
model_A <- hlaModelFromObj(model.list_A[[hla.id_A]])

#########################################
# set variables into HLA-B

Fin_hg38_model_B <- "your HLA-B model in hg38"
model.list_B <- get(load(Fin_hg38_model_B))

# HLA imputation model at HLA-B
hla.id_B <- "B"
model_B <- hlaModelFromObj(model.list_B[[hla.id_B]])

########################################
# set variables into HLA-C

Fin_hg38_model_C <- "your HLA-C model in hg38"
model.list_C <- get(load(Fin_hg38_model_C))

hla.id_C <- "C"
model_C <- hlaModelFromObj(model.list_C[[hla.id_C]])
########################################
# set variables into HLA-DRB1

Fin_hg38_model_DRB1 <- "your HLA-DRB1 model in hg38"
model.list_DRB1 <- get(load(Fin_hg38_model_DRB1))

hla.id_DRB1 <- "DRB1"
model_DRB1 <- hlaModelFromObj(model.list_DRB1[[hla.id_DRB1]])
########################################
# set variables into HLA-DQA1

Fin_hg38_model_DQA1 <- "your HLA-DQA1 model in hg38"
model.list_DQA1 <- get(load(Fin_hg38_model_DQA1))

hla.id_DQA1 <- "DQA1"
model_DQA1 <- hlaModelFromObj(model.list_DQA1[[hla.id_DQA1]])
#######################################
# set variables into HLA-DQB1

Fin_hg38_model_DQB1 <- "your HLA-DQB1 model in hg38"
model.list_DQB1 <- get(load(Fin_hg38_model_DQB1))

hla.id_DQB1 <- "DQB1"
model_DQB1 <- hlaModelFromObj(model.list_DQB1[[hla.id_DQB1]])
#######################################
# set variables into HLA-DPB1

Fin_hg38_model_DPB1 <- "your HLA-DPB1 model in hg38"
model.list_DPB1 <- get(load(Fin_hg38_model_DPB1))

hla.id_DPB1 <- "DPB1"
model_DPB1 <- hlaModelFromObj(model.list_DPB1[[hla.id_DPB1]])

#############################################################################
##########################Import your PLINK BED file#########################
#############################################################################

bed.fn <- "your MHC file.bed"
fam.fn <- "your MHC file.fam"
bim.fn <- "your MHC file.bim"

# assign genotypes
imputation_MHC <- hlaBED2Geno(bed.fn = bed.fn, fam.fn = fam.fn, 
                                      bim.fn = bim.fn, assembly = "hg38")
#############################################################################
######################      HLA imputation     ##############################
#############################################################################

# HLA-A prediction

pred.quess_A <- predict(model_A, imputation_MHC, 
                        type = "response+prob", match.type = "Position")

### ### ### ### ### ### ### ### ### ### ###
# inspecting

pred.quess_A$value                       #shows out the value 'data frame'
summary(pred.quess_A)                    #shows the overview of the results

# writing the results out into a table
results_A <- pred.quess_A$value

write.table(results_A, "your destination/imputation_HLA-A",
            quote = FALSE, row.names = FALSE)
#######################################

# HLA-B prediction

pred.quess_B <- predict(model_B, VPU_imputation_ref_MHC,
                        type = "response+prob", match.type = "Position")
summary(pred.quess_B)

# writing the results out into a table

results_B <- pred.quess_B$value

write.table(results_B, file = "your destination/imputation_HLA-B",
            quote = FALSE, row.names = FALSE)

#######################################

# HLA-C prediction

pred.quess_C <- predict(model_C, VPU_imputation_ref_MHC,
                        type = "response+prob", match.type = "Position")
summary(pred.quess_C)

# writing the results out into a table

results_C <- pred.quess_C$value

write.table(results_C, file = "your destination/imputation_HLA-C",
            quote = FALSE, row.names = FALSE)
#######################################

# HLA-DRB1 prediction

pred.quess_DRB1 <- predict(model_DRB1, VPU_imputation_ref_MHC,
                           type = "response+prob", match.type = "Position")
summary(pred.quess_DRB1)

# writing the results out into a table

results_DRB1 <- pred.quess_DRB1$value

write.table(results_DRB1, file = "your destination/imputation_HLA-DRB1",
            quote = FALSE, row.names = FALSE)
########################################

# HLA-DQA1 prediction

pred.quess_DQA1 <- predict(model_DQA1, VPU_imputation_ref_MHC,
                           type = "response+prob", match.type = "Position")

summary(pred.quess_DQA1)

# writing the results out into a table

results_DQA1 <- pred.quess_DQA1$value

write.table(results_DQA1, file = "your destination/imputation_HLA-DQA1",
            quote = FALSE, row.names = FALSE)

########################################

# HLA-DQB1 prediction

pred.quess_DQB1 <- predict(model_DQB1, VPU_imputation_ref_MHC,
                           type = "response+prob", match.type = "Position")

summary(pred.quess_DQB1)

# writing the results out into a table

results_DQB1 <- pred.quess_DQB1$value

write.table(results_DQB1, file = "your destination/imputation_HLA-DQB1",
            quote = FALSE, row.names = FALSE)

########################################

# HLA-DPB1 prediction

pred.quess_DPB1 <- predict(model_DPB1, VPU_imputation_ref_MHC,
                           type = "response+prob", match.type = "Position")
summary(pred.quess_DPB1)

## writing the results out into a table

results_DPB1 <- pred.quess_DPB1$value

write.table(results_DPB1, file = "your destination/imputation_HLA-DPB1",
            quote = FALSE, row.names = FALSE)

########################################