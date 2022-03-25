library(tidyverse)

# script to make a correlation test between the allele frequencies in two datasets
# (imputed and reference)
# 15.3.22 Jonna Clancy
##########################################################################################
# prerequisite: results from HLA_allele_freq.R (e.g. "Allelefreqs_HLA_A")
# and an added column including corresponding reference allele frequency
# name the file as e.g. "your_HLA_A_allele_freq_file.txt"

# example dataset:
# 'data.frame':	27 obs. of  5 variables:
# $ Locus    : chr  "HLA_A" "HLA_A" "HLA_A" "HLA_A" ...
# $ allele   : chr  "01:01" "01:02" "02" "02:01" ...
# $ count    : int  3762 1 17 13551 3 33 40 10562 6 1769 ...
# $ frequency: num  9.07e-02 2.41e-05 4.10e-04 3.27e-01 7.23e-05 ...
# $ reference.freq : num  0.09304 0.000081 NA 0.320013 0.000263 ...

##########################################################################################
#upload data
HLA_A_freq_imp <- read.table(file="your_HLA_A_allele_freq_file.txt",
                             header = TRUE, sep="\t", na = NA, dec=".", strip.white=TRUE,
                             stringsAsFactors = FALSE)

#correlation test
pears_A <- cor.test(HLA_A_freq_imp$frequency, HLA_A_freq_imp$reference.freq)

HLA_A_freq_imp <- HLA_A_freq_imp[c(-3,-7,-9,-12,-16,-19), ]

# pearson stats into a text string
pears.textA <- paste0("\n Pearson's correlation: ", pears_A$estimate %>% signif(4), ' (',
                      pears_A$conf.int[1] %>% signif(4), '-',
                      pears_A$conf.int[2] %>% signif(4), ')\n',
                      ' p-value: ', pears_A$p.value %>% signif(3))

# add stats to the data frame
HLA_A_freq_imp$pearson <- pears.textA
############################################################################################
                         #repeat the same for the rest of the genes
############################################################################################
#upload data
HLA_B_freq_imp <- read.table(file="your_HLA_B_allele_freq_file.txt",
                             header = TRUE, sep="\t", na = NA, dec=".", strip.white=TRUE,
                             stringsAsFactors = FALSE)

pears_B <- cor.test(HLA_B_freq_imp$frequency, HLA_B_freq_imp$reference.freq)

# pearson stats into a text string
pears.textB <- paste0("\n Pearson's correlation: ", pears_B$estimate %>% signif(4), ' (',
                      pears_B$conf.int[1] %>% signif(4), '-',
                      pears_B$conf.int[2] %>% signif(4), ')\n',
                      ' p-value: ', pears_B$p.value %>% signif(3))

# add to the data frame
HLA_B_freq_imp$pearson <- pears.textB
############################################################################################
#upload data
HLA_C_freq_imp <- read.table(file="your_HLA_C_allele_freq_file.txt",
                             header = TRUE, sep="\t", na = NA, dec=".", strip.white=TRUE,
                             stringsAsFactors = FALSE)

#correlation test
pears_C <- cor.test(HLA_C_freq_imp$frequency, HLA_C_freq_imp$reference.freq)

# pearson stats into a text string
pears.textC <- paste0("\n Pearson's correlation: ", pears_C$estimate %>% signif(4), ' (',
                      pears_C$conf.int[1] %>% signif(4), '-',
                      pears_C$conf.int[2] %>% signif(4), ')\n',
                      ' p-value: ', pears_C$p.value %>% signif(3))

# add to the data frame
HLA_C_freq_imp$pearson <- pears.textC
############################################################################################
#upload data
HLA_DRB1_freq_imp <- read.table(file="your_HLA_DRB1_allele_freq_file.txt",
                                header = TRUE, sep="\t", na = NA, dec=".", strip.white=TRUE,
                                stringsAsFactors = FALSE)

#correlation test
pears_DRB1 <- cor.test(HLA_DRB1_freq_imp$frequency, HLA_DRB1_freq_imp$reference.freq)

# pearson stats into a text string
pears.textDRB1 <- paste0("\n Pearson's correlation: ", pears_DRB1$estimate %>% signif(4), ' (',
                         pears_DRB1$conf.int[1] %>% signif(4), '-',
                         pears_DRB1$conf.int[2] %>% signif(4), ')\n',
                         ' p-value: ', pears_DRB1$p.value %>% signif(3))

# add to the data frame
HLA_DRB1_freq_imp$pearson <- pears.textDRB1
############################################################################################
#upload data
HLA_DQA1_freq_imp <- read.table(file="your_HLA_DQA1_allele_freq_file.txt",
                                header = TRUE, sep="\t", na = NA, dec=".", strip.white=TRUE,
                                stringsAsFactors = FALSE)

#correlation test
pears_DQA1 <- cor.test(HLA_DQA1_freq_imp$frequency, HLA_DQA1_freq_imp$reference.freq)

# pearson stats into a text string
pears.textDQA1 <- paste0("\n Pearson's correlation: ", pears_DQA1$estimate %>% signif(4), ' (',
                         pears_DQA1$conf.int[1] %>% signif(4), '-',
                         pears_DQA1$conf.int[2] %>% signif(4), ')\n',
                         ' p-value: ', pears_DQA1$p.value %>% signif(3))

# add to the data frame
HLA_DQA1_freq_imp$pearson <- pears.textDQA1
############################################################################################
#upload data
HLA_DQB1_freq_imp <- read.table(file="your_HLA_DQB1_allele_freq_file.txt",
                                header = TRUE, sep="\t", na = NA, dec=".", strip.white=TRUE,
                                stringsAsFactors = FALSE)

#correlation test
pears_DQB1 <- cor.test(HLA_DQB1_freq_imp$frequency, HLA_DQB1_freq_imp$reference.freq)

# pearson stats into a text string
pears.textDQB1 <- paste0("\n Pearson's correlation:", pears_DQB1$estimate %>% signif(4), ' (',
                         pears_DQB1$conf.int[1] %>% signif(4), '-',
                         pears_DQB1$conf.int[2] %>% signif(4), ')\n',
                         ' p-value: ', pears_DQA1$p.value %>% signif(3))

# add to the data frame
HLA_DQB1_freq_imp$pearson <- pears.textDQB1
############################################################################################
#upload data
HLA_DPB1_freq_imp <- read.table(file="your_HLA_DPB1_allele_freq_file.txt",
                                header = TRUE, sep="\t", na = NA, dec=".", strip.white=TRUE,
                                stringsAsFactors = FALSE)

#correlation test
pears_DPB1 <- cor.test(HLA_DPB1_freq_imp$frequency, HLA_DPB1_freq_imp$reference.freq)

# pearson stats into a text string
pears.textDPB1 <- paste0("\n Pearson's correlation:", pears_DPB1$estimate %>% signif(4), ' (',
                         pears_DPB1$conf.int[1] %>% signif(4), '-',
                         pears_DPB1$conf.int[2] %>% signif(4), ')\n',
                         ' p-value: ', pears_DPB1$p.value %>% signif(3))

# add to the data frame
HLA_DPB1_freq_imp$pearson <- pears.textDPB1
############################################################################################
# combine all the allele freq files into one table
Allele_freq_all <- rbind(HLA_A_freq_imp, HLA_B_freq_imp, HLA_C_freq_imp,
                         HLA_DRB1_freq_imp, HLA_DQA1_freq_imp,
                         HLA_DQB1_freq_imp, HLA_DPB1_freq_imp)

# keep the data in a correct order for plotting
Allele_freq_all$Locus <- factor(Allele_freq_all$Locus,
                                levels = Allele_freq_all$Locus
                                %>% unique)
# plot
ggplot(Allele_freq_all, aes(x=frequency, y=KSR.freq, label=pearson)) +
  geom_point(color="grey30", alpha=0.3, size=3) +
  geom_text(x=-Inf, y=Inf, vjust="inward", hjust="inward", size = 3.5) +
  geom_abline() +
  theme_light() +
  labs(title = "Correlation of allele frequencies",
       x="Imputed allele frequency", y="Reference allele frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 15)) +
  ylim(0, 0.5) +
  facet_wrap(~ Locus)
############################################################################################