---
title: "TexasTech_AMR"
author: "G Diaz"
date: "2023-05-25"
output: html_document
editor_options: 
  chunk_output_type: console
---
# NOTE R Markdowm

```{r IMPORTANT, eval=FALSE}
# to run code, show output, but don't show code. E.g.: To plot figures but not show the code used
  {r chunk-name, echo=FALSE}
# to run code, don't show output, but show code. E.g.: Package load, load data
  {r chunk-name, results="hide"}
# to run code, don't show output, don't show code. E.g.: to run secret package or option that don't want to show
  {r chunk-name, include=FALSE}
# to don't run code, don't show output, but show code. E.g.: Package installation
  {r chunk-name, eval=FALSE}
```

This will set "echo=TRUE" as default option for chunk, which means that all the code will be shown always
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# 0. Load packages

package installation: **DO IT ONLY ONCE**
```{r installation, eval=FALSE}
#BiocManager::install("igraph")
#BiocManager::install("phyloseq")
# install.packages("forcats")
# install.packages("readxl")
# install.packages("dplyr")
# install.packages("tibble")
# BiocManager::install("microbiome")
# BiocManager::install("metagenomeSeq")
# BiocManager::install("limma")
# install.packages("statmod")
# install.packages("psych")
# install.packages("tidyverse")
# install.packages("lme4")
# install.packages("lmerTest")
# install.packages("emmeans")
# install.packages("car")
# install.packages("table1")
# install.packages("lemon")
# install.packages("ggrepel")
```

package loading
```{r packages, results="hide"}
library("phyloseq")
library("forcats")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("microbiome")   # Extra data manipulation
library("metagenomeSeq")  # Data normalization (?) + Statistical testing
# library("limma")          # Statistical testing
# library("statmod")      # Will work with limma for "duplicateCorrelation" function
library("psych")          # Data exploration
# library("tidyverse")    # ggplot2, dplyr, tidyr, forcats, tibble, string, purrr, readr
library("table1")         # for doing tables
library("vegan")
library("lme4")           # Statistical testing
library("lmerTest")       # Statistical testing
library("emmeans")        # Statistical testing
library("car")            # Statistical testing
#library("table1")         # for doing tables
# library("lemon")          # for volcano plot
# library("ggrepel")      # for volcano plot
library("ggthemes")
library("maaslin3")       # Differential abundance testing
library("forester")
library("extrafont")      # To Load fonts from your system (do it only once): font_import()
library("taxize")         # To obtain taxonomic classification using taxid
library("ggpubr")         # For dot plot for contigs
library("RColorBrewer")
```

# 1. Load data
Build phyloseq object
```{r}
  otu_mat<- read_excel("Results_AMRplusplus3_FINAL/dedup_AMR_analytic_matrix_FINAL.xlsx", sheet = "OTU_matrix")
  tax_mat<- read_excel("Results_AMRplusplus3_FINAL/dedup_AMR_analytic_matrix_FINAL.xlsx", sheet = "Taxonomy_table")
  samples_df<- read_excel("Results_AMRplusplus3_FINAL/dedup_AMR_analytic_matrix_FINAL.xlsx", sheet = "Samples")
    
  samples_df$TimePoint= factor(samples_df$TimePoint, levels=c("D0", "D5", "D14"))
  samples_df$SampleID= as.factor(samples_df$SampleID)
  samples_df$TRT= factor(samples_df$TRT, levels = c("0","1"), labels = c("Control", "Treatment"))

  otu_mat <- otu_mat %>%
    tibble::column_to_rownames("ARG") 
  tax_mat <- tax_mat %>% 
    tibble::column_to_rownames("ARG")
   samples_df <- samples_df %>% 
    tibble::column_to_rownames("samples")
   
  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
  
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(samples_df)
  
  AMR <- phyloseq(OTU, TAX, samples)
  # 1708 taxa and 91 samples with 4 taxonomic ranks
  AMR <- subset_samples(AMR, sample_type =="Sample")
  # 1708 taxa and 90 samples with 4 taxonomic ranks
```

# 2. Data preparation
```{r}
#### CLEANING SAMPLE DUE TO LOW-QUALITY SEQUENCING AND ARG CLASSIFICATION #####
  # 84642D5TT23 : Sample ID 84642 Day 5. This sample had 29269 hits for ARGs while the average was 396315.4 and SD= 115573.2 So Mean-2(SD) = 165169 And sample 84642D5TT23 < Mean-2(SD). Demonstration:
a = colSums(otu_table(AMR))
a = sort(a)
a[1] # number of hits for sample 84642D5TT23
b = a[-1]
describe(b) # mean and sd of the ARG hits distributed across all samples.

# Take out the 84642D5TT23 to have the clean ps
AMR_cln <- prune_samples(!(sample_names(AMR) %in% c("84642D5TT23")), AMR)
# NORMALIZATION ALL DAYS 
ps<-AMR_cln
rumen.metaseq <- phyloseq_to_metagenomeSeq(ps) # change here
rumen.metaseq.norm<- cumNorm(rumen.metaseq, p=cumNormStat(rumen.metaseq))
CSS_rumen.metaseq <- MRcounts(rumen.metaseq.norm, norm = TRUE)
CSS_AMR <- merge_phyloseq(otu_table(CSS_rumen.metaseq,
taxa_are_rows=T),sample_data(ps),tax_table(ps))
CSS_AMR
# AGGLOMERATION
  AMR_gene= tax_glom(CSS_AMR, "Class")
# USE FOR BETA-DIV
  AMR_gene
# USE FOR ALPHA-DIV
  AMR_div = tax_glom(AMR_cln, "Gene_group") # Gene level
```
## 2.1. Sup Tab 1
```{r}
# creating table of classified ARGs at any level a
tmp2 = colSums(otu_table(AMR))
# tmp2 = tmp2[-41,] # without 84642D5TT23
tmp2 = as.data.frame(tmp2)
tmp = data.frame(sample_data(AMR)$TRT)
tmp3 = cbind.data.frame(tmp,tmp2)
# generating summary table
Suptab1 = table1(~ tmp2| `sample_data.AMR..TRT`, data=tmp3)
# statistics
wilcox.test(tmp2~`sample_data.AMR..TRT`, data = tmp3) 
# t.test(tmp2~TRT, data = tmp3) 
# RAW READS COMPARISON
df<- as.data.frame(read_excel("Results_AMRplusplus3_FINAL/dedup_AMR_analytic_matrix_FINAL.xlsx", sheet = "Samples"))
# df = df[-49,] # without 84642D5TT23
df<- subset(df, df$sample_type == "Sample")
df_2 <- cbind(df,tmp3) 
# Proportion of on-target reads
df_2 = df_2[-66,]
df_2$prop_on_target <- df_2$tmp2/df_2$ARG_raw_reads * 100
sum(df_2$tmp2)/sum(df_2$ARG_enrich_seq_rawreads) * 100
mean(df_2$prop_on_target); sd(df_2$prop_on_target)
# table
table1(~ ARG_raw_reads| factor(TRT), data=df)
# statistics 
wilcox.test(ARG_raw_reads~TRT, data = df) 

```

# 3. Relative abundance
## 3.1. Figure 3A
```{r}
# Obtaining table for RA
  phy_relative <- transform_sample_counts(AMR_gene, function(x) (x / sum(x))*100 )
  phy_relative_long <- psmelt(phy_relative)
  phy_relative_long <- phy_relative_long %>%
  group_by(Gene_group) %>%
  mutate(mean_relative_abund = mean(Abundance))
  phy_relative_long$Gene_group <- as.character(phy_relative_long$Gene_group)
  phy_relative_long$mean_relative_abund <- as.numeric(phy_relative_long$mean_relative_abund)
  phy_relative_long$Gene_group[phy_relative_long$Abundance < 2] <- "Others (< 2%)"
  
# ordering Relative abundances (Y axis) for the plot ####
# Calculate total abundance for each Gene_group
phy_relative_long <- phy_relative_long %>%
  group_by(Gene_group) %>%
  mutate(TotalAbundance = sum(Abundance)) %>%
  ungroup()
# Reorder the Gene_group factor based on total abundance
phy_relative_long <- phy_relative_long %>%
  mutate(Gene_group = factor(Gene_group, levels = unique(Gene_group[order(TotalAbundance)])))
# Reordering the Samples (X axis) for the plot ####
sample_grouping <- phy_relative_long %>% 
  group_by(Sample) %>% 
  slice_max(order_by = Abundance) %>% 
  select(Gene_group, Sample) %>% 
  rename(peak_class = Gene_group)
phy_relative_long <- phy_relative_long %>% 
  inner_join(sample_grouping, by = "Sample") %>% 
  group_by(peak_class) %>% 
  mutate(rank = rank(Gene_group)) %>%  # rank samples at the level of each peak subtype
  mutate(Sample = reorder(Sample, rank)) %>%   # this reorders samples
  ungroup()


## Plot
  Fig3A <- phy_relative_long %>%
    ggplot(aes(x = Sample, y = Abundance, fill = Gene_group)) +
    geom_bar(stat = "identity",width = 1,color="white") +
    geom_bar(stat = "identity", alpha=1)+theme_bw()+
    facet_wrap(~TimePoint+TRT, nrow=1, scales="free_x")+
    theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
    ggthemes::scale_fill_tableau("Tableau 20")
```

## 3.2. Description of ARGs
```{r}
# Number of unique ARGs overall, at the Mechanism and Gene_group level
AMR_cln
 glom= tax_glom(AMR_cln, "Gene_group", bad_empty = c("NA", NA, "", " ", "\t"))
# Number of reads classified overall, and at the Mechanism and Gene_group level
  sum(taxa_sums(glom))
# Checking mean proportion and SD of ARGs overall
# This way:
  # table1(~ Abundance | Gene_group, data=phy_relative_long)
# Or this way:
summary_stats <- phy_relative_long %>%
  group_by(Gene_group) %>%
  summarise(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    sd_abundance = sd(Abundance, na.rm = TRUE),
    n = n()
  ) %>%
  arrange(-mean_abundance)
```

# 4. Alpha-div
```{r}
# Calculatin at Gene_group level
  div_16S= estimate_richness(AMR_div, measures=c("Observed", "InvSimpson", "Shannon"))
  Alpha_gene= cbind(div_16S,(evenness(AMR_div, 'pielou')))
  names(Alpha_gene) <- paste(names(Alpha_gene), "gene", sep = "_")
# Biding with dataframe
  AMR_df = cbind(data.frame(sample_data(AMR_div)),Alpha_gene)
```

## 4.1. Figure 3.C, D, E
```{r}
# Plot Richness, Shannon's and Pielou's at the gene_group level
Fig3C <- ggplot(AMR_df, aes(TimePoint,Observed_gene, fill=TRT))+
  geom_boxplot(position=position_dodge(0.8),lwd=1, alpha= 0.2)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "Time-point", y = "Richness", fill="Groups")

Fig3D <- ggplot(AMR_df, aes(TimePoint,Shannon_gene, fill=TRT))+
  geom_boxplot(position=position_dodge(0.8),lwd=1, alpha= 0.2)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "Time-point", y = "Shannon's Index", fill="Groups")

Fig3E <- ggplot(AMR_df, aes(TimePoint,pielou_gene, fill=TRT))+
  geom_boxplot(position=position_dodge(0.8),lwd=1, alpha= 0.2)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "Time-point", y = "Eveness Index", fill="Groups")

# STATISTICS: alpha ~ TRT*TimePoint + BCS + (1|SampleID) ####
# formatting BCS variable
AMR_df <- AMR_df %>%
  group_by(SampleID) %>%
  mutate(BCS = BCS_MC[TimePoint == "D0"][1]) %>%
  ungroup()
# formatting classified reads
class_reads = as.data.frame(colSums(otu_table(AMR_cln)))
class_reads <- class_reads %>% rename(classified_ARG = `colSums(otu_table(AMR_cln))`)
AMR_df<- cbind(class_reads,AMR_df)
AMR_df$classified_ARG_resc <- AMR_df$classified_ARG/100
# Model for Richness
AMR_rich = lmer(Observed_gene ~ TRT*TimePoint + BCS + (1|SampleID), data=AMR_df, REML = T)
  summary(AMR_rich)
  Anova(AMR_rich, type = "III")
  emmeans(AMR_rich,pairwise~TRT|TimePoint)
# Model for Shannon's index
AMR_shan = lmer(Shannon_gene ~ TRT*TimePoint + BCS +(1|SampleID), data=AMR_df, REML = T)
  summary(AMR_shan)
  Anova(AMR_shan, type = "III")
  #confint(gene_shan)
  emmeans(AMR_shan,pairwise~TRT|TimePoint)
# Model for Eveness
AMR_even = lmer(pielou_gene ~ TRT*TimePoint + BCS + (1|SampleID), data=AMR_df, REML = T)
  summary(AMR_even)
  Anova(AMR_even, type = "III")
  emmeans(AMR_even,pairwise~TRT|TimePoint)
```

# 5. Beta-div
```{r}
# ALL DAYS
set.seed(143)
genus.ord <- ordinate(AMR_gene, "NMDS", "bray")
```

## 5.1. Fig 3.B
```{r}
# plot
Fig3B<- plot_ordination(AMR_gene, genus.ord, type="sample", color="TRT")+
  geom_point(size=3, alpha=0.8) + theme_classic() + stat_ellipse() + facet_wrap(~TimePoint)+
  labs(color = "Groups")
Fig3B$layers <- Fig3B$layers[-1]
Fig3B

# STATISTICS: beta ~ TRT*TimePoint + BCS + (1|SampleID) ####
dist = vegdist(t(otu_table(AMR_gene)), method ="bray" )
adonis2(dist ~ TimePoint + TRT + BCS, data=AMR_df , by="margin")
# STATISTICS DAY 0
AMR_df_0 = subset(AMR_df, AMR_df$TimePoint == "D0")
AMR_gene_0 = subset_samples(AMR_gene, TimePoint == "D0")
dist_0 = vegdist(t(otu_table(AMR_gene_0)), method ="bray" )
adonis2(dist_0 ~ TRT + BCS, data=AMR_df_0 , by="margin")
# STATISTICS DAY 5
AMR_df_5 = subset(AMR_df, AMR_df$TimePoint == "D5")
AMR_gene_5 = subset_samples(AMR_gene, TimePoint == "D5")
dist_5 = vegdist(t(otu_table(AMR_gene_5)), method ="bray" )
adonis2(dist_5 ~ TRT + BCS, data=AMR_df_5 , by="margin")
# STATISTICS DAY 14
AMR_df_14 = subset(AMR_df, AMR_df$TimePoint == "D14")
AMR_gene_14 = subset_samples(AMR_gene, TimePoint == "D14")
dist_14 = vegdist(t(otu_table(AMR_gene_14)), method ="bray" )
adonis2(dist_14 ~ TRT + BCS, data=AMR_df_14 , by="margin")

# BETADISPERSION
betadisp_0 <- betadisper(dist_0, AMR_df_0$TRT, type = "centroid")
permutest(betadisp_0, permutations = 999)
betadisp_5 <- betadisper(dist_5, AMR_df_5$TRT, type = "centroid")
permutest(betadisp_5, permutations = 999)
betadisp_14 <- betadisper(dist_14, AMR_df_14$TRT, type = "centroid")
permutest(betadisp_14, permutations = 999)

# STATISTICS CONTROL GROUP (all time points)
AMR_df_CTR = subset(AMR_df, AMR_df$TRT == "Control")
AMR_gene_CTR = subset_samples(AMR_gene, TRT == "Control")
AMR_dist_CTR = vegdist(t(otu_table(AMR_gene_CTR)), method ="bray" )
adonis2(AMR_dist_CTR ~ TimePoint + BCS, data=AMR_df_CTR , by="margin")
# STATISTICS TREATMENT GROUP (all time points)
AMR_df_TRT = subset(AMR_df, AMR_df$TRT == "Treatment")
AMR_gene_TRT = subset_samples(AMR_gene, TRT == "Treatment")
AMR_dist_TRT = vegdist(t(otu_table(AMR_gene_TRT)), method ="bray" )
adonis2(AMR_dist_TRT ~ TimePoint + BCS, data=AMR_df_TRT , by="margin")

# MORE PLOTS!!! Plot daily beta diversity
# Day 0
set.seed(143)
Fig3B_1 <- plot_ordination(AMR_gene_0, ordinate(AMR_gene_0, "NMDS", "bray"), color="TRT")+ geom_point(size=3, alpha=0.8) + theme_classic() + stat_ellipse() + labs(color = "Groups")
Fig3B_1$layers <- Fig3B_1$layers[-1]
Fig3B_1 # Stress:     0.05816244 

# Day 5
set.seed(143)
Fig3B_2 <- plot_ordination(AMR_gene_5, ordinate(AMR_gene_5, "NMDS", "bray"), color="TRT")+ geom_point(size=3, alpha=0.8) + theme_classic() + stat_ellipse() + labs(color = "Groups")
Fig3B_2$layers <- Fig3B_2$layers[-1]
Fig3B_2 # Stress:     0.03613829

# Day 14
set.seed(143)
Fig3B_3 <- plot_ordination(AMR_gene_14, ordinate(AMR_gene_14, "NMDS", "bray"), color="TRT")+ geom_point(size=3, alpha=0.8) + theme_classic() + stat_ellipse() + labs(color = "Groups")
Fig3B_3$layers <- Fig3B_3$layers[-1]
Fig3B_3 # Stress:     0.03980002
```
# 6. Diff Abund Test
```{r}
# Get input_metadata from ps sample data
df_input_metadata = AMR_df
# change name of variable to "reads"
df_input_metadata <- df_input_metadata %>%
  rename(reads = ARG_raw_reads)
# Get and format input_data from ps otu table ####
df_input_data = data.frame(otu_table(AMR_div))
  #Create input tax table with Gene_group names from ps object:
  tax_table=tax_table(AMR_div)[,"Gene_group"]
  # Put row names to a column called ARG in tax_table
  tax_table= data.frame(tax_table) %>% rownames_to_column(var='ARG')
  # Put row names to a column called ARG in df_input_data
  tmp = df_input_data %>% rownames_to_column(var='ARG')
  # Merging df_input_data (otu_table) with tax_table so we have genus names instead of ASV1, ASV2, etc.
  tmp2 = left_join(tmp, tax_table, by=c("ARG"))
  # Put column names 'Gene_group' to row names so we have the names of the Genes as row names
  tmp3 = tmp2 %>% column_to_rownames(var='Gene_group')
  # Delete the ARG column as we already have the Gene names as row names. Also Transpose so now the Genes are columns and the samples are rows
  tmp4 = subset(tmp3, select = -c(ARG)) 
  df_input_data=as.data.frame(t(tmp4))
  # Finally (and just for this dataset) I have an "X" at the beginning of the name of the rows, that I need to get rid of
  rownames(df_input_data) <- gsub("^X", "", rownames(df_input_data))


# Statistical model with interactions: ARG ~ TRT * TimePoint + reads + BCS + (1|SampleID) #### 
fit_out3 <- maaslin3(input_data = df_input_data,
                    input_metadata = df_input_metadata,
                    output = 'Results_AMRplusplus3_FINAL/model_interact_ARG',
                    formula = '~ TimePoint * TRT + BCS + reads + (1|SampleID)',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 100,
                    cores = 1)

# Statistical model without interactions: otu ~ TRT + TimePoint + Classified + BCS + (1|SampleID) #### 
fit_out4 <- maaslin3(input_data = df_input_data,
                    input_metadata = df_input_metadata,
                    output = 'Results_AMRplusplus3_FINAL/model_NOinteract_ARG',
                    formula = '~ TimePoint + TRT + BCS + reads + (1|SampleID)',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 100,
                    cores = 1)

# Just to prove that control vs treatment did not have any meaningful differential abundant ARG at DAY 0 (baseline)
df_input_metadata_D0 <- subset(df_input_metadata, TimePoint == "D0")
df_input_data_D0 <- df_input_data[rownames(df_input_metadata_D0), , drop = FALSE]

# Statistical model for DAY0: otu ~ TRT + Classified + BCS #### 
fit_out5 <- maaslin3(input_data = df_input_data,
                    input_metadata = df_input_metadata,
                    output = 'Results_AMRplusplus3_FINAL/model_D0_ARG',
                    formula = '~ TRT + BCS + reads',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 100,
                    cores = 1)
  
```
## 6.1. Figure 4
```{r}
# FORMATTING ####
# After exponentiating in excel and calculating the 95% CI by: coeff + 1.96*Standard error; coeff - 1.96*standard error. For log2 I exponentiated 2 to the coeff and 2 to the lower/upper bound.
# Read table
DAT_gene <- read.csv("Results_AMRplusplus3_FINAL/model_interact_ARG/significant_results.tsv", sep = "\t")
# Formatting p-values FDR
DAT_gene <- subset(DAT_gene, DAT_gene$qval_individual <= 0.1)
DAT_gene$qval_individual[DAT_gene$qval_individual < 0.001] <- "p<0.001"
DAT_gene$qval_individual[DAT_gene$qval_individual < 0.01] <- "p<0.01"
DAT_gene$qval_individual[DAT_gene$qval_individual < 0.05] <- "p<0.05"
DAT_gene$qval_individual[DAT_gene$qval_individual < 0.1] <- "p<0.1"
# calculating the overall prevalence for significant ARGs
DAT_gene$overall_prev <- round((DAT_gene$N.not.zero/DAT_gene$N)*100,1)
# calculating the overall abundance for significant ARGs
  # Calculating the total sum of reads for all Genes
  total_sum_all <- sum(df_input_data)
  # Extract the list of Genes 
  column_names <- DAT_gene$feature
  # Initialize a vector to store proportions
  proportions <- numeric(length(column_names))
  # Calculate proportions for each Genus
  for (i in seq_along(column_names)) {
    col_sum <- sum(df_input_data[[column_names[i]]])
    proportions[i] <- (col_sum / total_sum_all) * 100
  }
  # Round the proportions to 3 decimal places
  proportions <- round(proportions, 3)
  # Create the resulting dataframe
  result_df <- data.frame(
    Column = column_names,
    overall_abund = proportions
  )
  # Merge dataframes
  DAT_gene = cbind(DAT_gene, result_df)

# Getting TABLE 2 ####
# Subset for plot for day 5 and day 14
DAT_gene_5 <- subset(DAT_gene, DAT_gene$name == "TimePointD5:TRTTreatment" & DAT_gene$model == "abundance")
DAT_gene_14 <- subset(DAT_gene, DAT_gene$name == "TimePointD14:TRTTreatment" & DAT_gene$model == "abundance")
# selecting columns for the table
DAT_gene_5 <- DAT_gene_5 %>%
  select('feature2','overall_prev','overall_abund','qval_individual','coef','lower_coef','upper_coef') %>%
  arrange(feature2) %>% filter(feature2 != "Not found_APH2.DPRIME")
DAT_gene_14 <- DAT_gene_14 %>%
  select('feature2','overall_prev','overall_abund','qval_individual','coef','lower_coef','upper_coef') %>%
  arrange(feature2)
# Creating rows for subheadings Day 5 and Day 14
row_5 <- data.frame(feature2= "Day 5", overall_prev = NA, overall_abund =  NA, qval_individual = NA, coef= NA, lower_coef= NA, upper_coef= NA)
row_14 <- data.frame(feature2= "Day 14", overall_prev = NA, overall_abund =  NA, qval_individual = NA, coef= NA, lower_coef= NA, upper_coef= NA)
# Pasting all subsettings with subheadings
DAT_gene_5_14 <- rbind(row_5,DAT_gene_5,row_14,DAT_gene_14)
# Reanaming columns for the table
DAT_gene_5_14 <- DAT_gene_5_14 %>%
  rename(ARG = feature2, `Overall prev. (%)` = overall_prev, `Overall abund. (%)` =  overall_abund, FDR = qval_individual)
# Indenting for Day5 and Day 14
DAT_gene_5_14$ARG <- ifelse(is.na(DAT_gene_5_14$`Overall abund. (%)`), 
                         DAT_gene_5_14$ARG,
                         paste0("   ", DAT_gene_5_14$ARG))
# Loading font for table
loadfonts()
# Printing table 
tab3 <- forester(left_side_data = DAT_gene_5_14[,1:4],
           estimate = DAT_gene_5_14$coef,
           ci_low = DAT_gene_5_14$lower_coef,
           ci_high = DAT_gene_5_14$upper_coef,
           display = TRUE,
           #nudge_y = -.3,
           estimate_precision = 3,
           null_line_at = 0,
           xlim = c(-1, 9),
           #xbreaks = c(-1,0,1,2,3,4,5,6,7,8,9,10),
           #justify = c(0,0.5,0.5,1),
           estimate_col_name = "Log2(Fold change) (95% CI)",
           file_path = here::here("Results_AMRplusplus3_FINAL/FC_ARG_Day5_14.png"),
           #render_as = "rmarkdown",
           font_family = "Arial")


# Reanaming columns for the table
DAT_gene_14 <- DAT_gene_14 %>%
  rename(ARG = feature, `Overall prev. (%)` = overall_prev, `Overall abund. (%)` =  overall_abund, FDR = qval_individual)
# Getting ride of CMY for improving visualization
DAT_gene_14 <- DAT_gene_14 %>% filter(ARG != 'CMY')
# Printing table 
tab2.1 <- forester(left_side_data = DAT_gene_14[,1:4],
           estimate = DAT_gene_14$coef,
           ci_low = DAT_gene_14$lower_coef,
           ci_high = DAT_gene_14$upper_coef,
           display = TRUE,
           #nudge_y = -.05,
           estimate_precision = 3,
           null_line_at = 0,
           xlim = c(-1, 8),
           #xbreaks = c(0,1,2),
           #justify = c(0,0.5,0.5,1),
           estimate_col_name = "Fold Change (95% CI)",
           file_path = here::here("FP_genus_D5.png"),
           render_as = "rmarkdown",
           font_family = "Arial")
```

## 6.2. Sup. Fig 4
```{r}
# DA over time CONTROL & TREATED GROUP #########
# Maaslin 3 DAT:
# IMPORTANT - Change dataframe here: Use AMR_df_CTR or AMR_df_TRT
df_input_metadata = AMR_df_CTR
# Then run!
df_input_metadata <- df_input_metadata %>%
  rename(reads = ARG_raw_reads)
# Get ps object with non-normalized counts for CTR and TRT from ps object used for alpha diversity
AMR_ps_alpha_CTR <- subset_samples(AMR_div, TRT == "Control")  
AMR_ps_alpha_TRT <- subset_samples(AMR_div, TRT == "Treatment")  
# IMPORTANT - Change ps object here: Use AMR_ps_alpha_CTR or AMR_ps_alpha_TRT
df_input_data = data.frame(otu_table(AMR_ps_alpha_TRT))
# Then run!
  tax_table=tax_table(AMR_div)[,"Gene_group"]
  tax_table= data.frame(tax_table) %>% rownames_to_column(var='ARG')
  tmp = df_input_data %>% rownames_to_column(var='ARG')
  tmp2 = left_join(tmp, tax_table, by=c("ARG"))
  tmp3 = tmp2 %>% column_to_rownames(var='Gene_group')
  tmp4 = subset(tmp3, select = -c(ARG)) 
  df_input_data=as.data.frame(t(tmp4))
  rownames(df_input_data) <- gsub("^X", "", rownames(df_input_data))
# Statistical model: otu ~ TimePoint + reads + BCS #### 
fit_out8 <- maaslin3(input_data = df_input_data,
                    input_metadata = df_input_metadata,
                    output = 'Results_AMRplusplus3_FINAL/Timepoint_CONTROL_ARG',
                    formula = '~ TimePoint + BCS + reads',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 100,
                    cores = 1)
# Statistical model with interactions: otu ~ TimePoint + reads + BCS #### 
fit_out9 <- maaslin3(input_data = df_input_data,
                    input_metadata = df_input_metadata,
                    output = 'Results_AMRplusplus3_FINAL/Timepoint_TREATED_ARG',
                    formula = '~ TimePoint + BCS + reads',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 100,
                    cores = 1)

# FOREST PLOT FOR CONTROLS ####
DAT_gene <- read.csv("Results_AMRplusplus3_FINAL/Timepoint_CONTROL_ARG/significant_results.tsv", sep = "\t")
DAT_gene <- subset(DAT_gene, DAT_gene$qval_individual <= 0.1)
DAT_gene$qval_individual[DAT_gene$qval_individual < 0.001] <- "p<0.001"
DAT_gene$qval_individual[DAT_gene$qval_individual < 0.01] <- "p<0.01"
DAT_gene$qval_individual[DAT_gene$qval_individual < 0.05] <- "p<0.05"
DAT_gene$qval_individual[DAT_gene$qval_individual < 0.1] <- "p<0.1"
DAT_gene$overall_prev <- round((DAT_gene$N.not.zero/DAT_gene$N)*100,1)
  total_sum_all <- sum(df_input_data)
  column_names <- DAT_gene$feature
  proportions <- numeric(length(column_names))
  for (i in seq_along(column_names)) {
    col_sum <- sum(df_input_data[[column_names[i]]])
    proportions[i] <- (col_sum / total_sum_all) * 100
  }
  proportions <- round(proportions, 3)
  result_df <- data.frame(
    Column = column_names,
    overall_abund = proportions
  )
  DAT_gene = cbind(DAT_gene, result_df)
DAT_gene_5 <- subset(DAT_gene, DAT_gene$name == "TimePointD5" & DAT_gene$model == "abundance")
DAT_gene_14 <- subset(DAT_gene, DAT_gene$name == "TimePointD14" & DAT_gene$model == "abundance")
DAT_gene_5 <- DAT_gene_5 %>%
  select('feature2','overall_prev','overall_abund','qval_individual','coef','lower_coef','upper_coef') %>%
  arrange(feature2)
DAT_gene_14 <- DAT_gene_14 %>%
  select('feature2','overall_prev','overall_abund','qval_individual','coef','lower_coef','upper_coef') %>%
  arrange(feature2)
row_5 <- data.frame(feature2= "Day 5 vs Day 0", overall_prev = NA, overall_abund =  NA, qval_individual = NA, coef= NA, lower_coef= NA, upper_coef= NA)
row_14 <- data.frame(feature2= "Day 14 vs Day 0", overall_prev = NA, overall_abund =  NA, qval_individual = NA, coef= NA, lower_coef= NA, upper_coef= NA)
DAT_gene_5_14 <- rbind(row_5,DAT_gene_5,row_14,DAT_gene_14)
DAT_gene_5_14 <- DAT_gene_5_14 %>%
  rename(ARG = feature2, `Overall prev. (%)` = overall_prev, `Overall abund. (%)` =  overall_abund, FDR = qval_individual)
DAT_gene_5_14$ARG <- ifelse(is.na(DAT_gene_5_14$`Overall abund. (%)`), 
                         DAT_gene_5_14$ARG,
                         paste0("   ", DAT_gene_5_14$ARG))
loadfonts()
# Printing table 
tab5 <- forester(left_side_data = DAT_gene_5_14[,1:4],
           estimate = DAT_gene_5_14$coef,
           ci_low = DAT_gene_5_14$lower_coef,
           ci_high = DAT_gene_5_14$upper_coef,
           display = TRUE,
           #nudge_y = -.3,
           estimate_precision = 3,
           null_line_at = 0,
           xlim = c(-8, 4.5),
           #xbreaks = c(-1,0,1,2,3,4,5,6,7,8,9,10),
           #justify = c(0,0.5,0.5,1),
           estimate_col_name = "Log2(Fold change) (95% CI)",
           file_path = here::here("Results_AMRplusplus3_FINAL/FP_ARG_CTR.png"),
           #render_as = "rmarkdown",
           font_family = "Arial")


# FOREST PLOT FOR TREATED ####
DAT_gene <- read.csv("Results_AMRplusplus3_FINAL/Timepoint_TREATED_ARG/significant_results.tsv", sep = "\t")
DAT_gene <- subset(DAT_gene, DAT_gene$qval_individual <= 0.1)
DAT_gene$qval_individual[DAT_gene$qval_individual < 0.001] <- "p<0.001"
DAT_gene$qval_individual[DAT_gene$qval_individual < 0.01] <- "p<0.01"
DAT_gene$qval_individual[DAT_gene$qval_individual < 0.05] <- "p<0.05"
DAT_gene$qval_individual[DAT_gene$qval_individual < 0.1] <- "p<0.1"
DAT_gene$overall_prev <- round((DAT_gene$N.not.zero/DAT_gene$N)*100,1)
  total_sum_all <- sum(df_input_data)
  column_names <- DAT_gene$feature
  proportions <- numeric(length(column_names))
  for (i in seq_along(column_names)) {
    col_sum <- sum(df_input_data[[column_names[i]]])
    proportions[i] <- (col_sum / total_sum_all) * 100
  }
  proportions <- round(proportions, 3)
  result_df <- data.frame(
    Column = column_names,
    overall_abund = proportions
  )
  DAT_gene = cbind(DAT_gene, result_df)
DAT_gene_5 <- subset(DAT_gene, DAT_gene$name == "TimePointD5" & DAT_gene$model == "abundance")
DAT_gene_14 <- subset(DAT_gene, DAT_gene$name == "TimePointD14" & DAT_gene$model == "abundance")
DAT_gene_5 <- DAT_gene_5 %>%
  select('feature2','overall_prev','overall_abund','qval_individual','coef','lower_coef','upper_coef') %>%
  arrange(feature2)
DAT_gene_14 <- DAT_gene_14 %>%
  select('feature2','overall_prev','overall_abund','qval_individual','coef','lower_coef','upper_coef') %>%
  arrange(feature2)
row_5 <- data.frame(feature2= "Day 5 vs Day 0", overall_prev = NA, overall_abund =  NA, qval_individual = NA, coef= NA, lower_coef= NA, upper_coef= NA)
row_14 <- data.frame(feature2= "Day 14 vs Day 0", overall_prev = NA, overall_abund =  NA, qval_individual = NA, coef= NA, lower_coef= NA, upper_coef= NA)
DAT_gene_5_14 <- rbind(row_5,DAT_gene_5,row_14,DAT_gene_14)
DAT_gene_5_14 <- DAT_gene_5_14 %>%
  rename(ARG = feature2, `Overall prev. (%)` = overall_prev, `Overall abund. (%)` =  overall_abund, FDR = qval_individual)
DAT_gene_5_14$ARG <- ifelse(is.na(DAT_gene_5_14$`Overall abund. (%)`), 
                         DAT_gene_5_14$ARG,
                         paste0("   ", DAT_gene_5_14$ARG))
loadfonts()
# Printing table 
tab6 <- forester(left_side_data = DAT_gene_5_14[,1:4],
           estimate = DAT_gene_5_14$coef,
           ci_low = DAT_gene_5_14$lower_coef,
           ci_high = DAT_gene_5_14$upper_coef,
           display = TRUE,
           #nudge_y = -.3,
           estimate_precision = 3,
           null_line_at = 0,
           xlim = c(-4, 11),
           #xbreaks = c(-1,0,1,2,3,4,5,6,7,8,9,10),
           #justify = c(0,0.5,0.5,1),
           estimate_col_name = "Log2(Fold change) (95% CI)",
           file_path = here::here("Results_AMRplusplus3_FINAL/FP_ARG_TRT.png"),
           #render_as = "rmarkdown",
           font_family = "Arial")
```

# 7. Contig analysis
## 7.1. Load data CMY
```{r}
# Load data and build phyloseq object
  otu_mat<- read_excel("contig_analysis_CMY.xlsx", sheet = "otu")
  tax_mat<- read_excel("contig_analysis_CMY.xlsx", sheet = "tax")
  samples_df<- read_excel("contig_analysis_CMY.xlsx", sheet = "samples")
    
  samples_df$TimePoint= factor(samples_df$TimePoint, levels=c("D0", "D5", "D14"))
  samples_df$SampleID= as.factor(samples_df$SampleID)
  samples_df$TRT= factor(samples_df$TRT, levels = c("0","1"), labels = c("Control", "Treatment"))

  otu_mat <- otu_mat %>%
    tibble::column_to_rownames("contig") 
  tax_mat <- tax_mat %>% 
    tibble::column_to_rownames("contig")
   samples_df <- samples_df %>% 
    tibble::column_to_rownames("DNA_ID")
   
  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
  
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(samples_df)

  contig <- phyloseq(OTU, TAX, samples)
  contig <- prune_samples(!(sample_names(contig) %in% c("84642D5TT23","Mock")), contig)

contig_raw = tax_glom(contig, "Full_name")
contig_class= tax_glom(contig, "Class")
contig_order= tax_glom(contig, "Order")
contig_family= tax_glom(contig, "Family")
```

## 7.2. Figure 5
```{r}
# Obtaining table for RA
  phy_relative_long <- psmelt(contig_raw)
  phy_relative_long <- phy_relative_long %>%
  group_by(Full_name) %>%
  mutate(mean_relative_abund = mean(Abundance))
  phy_relative_long$Full_name <- as.character(phy_relative_long$Full_name)
  phy_relative_long$mean_relative_abund <- as.numeric(phy_relative_long$mean_relative_abund)

## Plot
phy_relative_long %>%
    ggplot(aes(x = Sample, y = Abundance, fill = Full_name)) +
    geom_bar(stat = "identity",width = 1,color="white") +
    geom_bar(stat = "identity", alpha=1)+theme_bw()+
    scale_y_continuous(breaks = seq(0, 16, by = 2)) +
    facet_wrap(~TimePoint+TRT, nrow=1, scales="free_x")+
    theme(axis.text.x = element_blank())+ylab("N° of contigs harboring CMY gene")+
    ggthemes::scale_fill_tableau(palette = "Classic 10 Medium")
## To know how many contigs per group and time point
sum(taxa_sums(subset_samples(contig, TRT == "Treatment" & TimePoint == "D14")))
```

## 7.3. Load data betalactams
```{r}
# LOAD AND PARSE DATA TO MAKE IT READABLE FOR A PHYLOSEQ OBJECT ####
# Load data column as a list
df <- read_excel("contig_analysis_betalactams.xlsx", sheet = "tax")  # Change 'Sheet1' to your sheet name if different
# Extract the column named 'Bacteria'
tax_ids <- df[["Bacteria"]]
# Get classification for each taxonomic ID
classifications <- classification(tax_ids, db = "ncbi")
# Convert the classification list to a data frame
class_df <- do.call(rbind, lapply(classifications, function(x) {
  setNames(x$name, x$rank)
}))
# Add the taxonomic IDs as row names
rownames(class_df) <- tax_ids
# print the tax table
write.csv(class_df, file = "contig_betalactams_tax_table.csv")

# After doing a lot of manual processing in excel and running a couple of bash scripts I finally got the number of contigs for unknows and only_contig files ####
# Load data and build phyloseq object
  otu_mat<- read_excel("contig_analysis_betalactams.xlsx", sheet = "otu")
  tax_mat<- read_excel("contig_analysis_betalactams.xlsx", sheet = "tax")
  samples_df<- read_excel("contig_analysis_betalactams.xlsx", sheet = "samples")
    
  samples_df$TimePoint= factor(samples_df$TimePoint, levels=c("D0", "D5", "D14"))
  samples_df$SampleID= as.factor(samples_df$SampleID)
  samples_df$TRT= factor(samples_df$TRT, levels = c("0","1"), labels = c("Control", "Treatment"))

  otu_mat <- otu_mat %>%
    tibble::column_to_rownames("Bacteria") 
  tax_mat <- tax_mat %>% 
    tibble::column_to_rownames("Bacteria")
   samples_df <- samples_df %>% 
    tibble::column_to_rownames("DNA_ID")
   
  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
  
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(samples_df)

  contig <- phyloseq(OTU, TAX, samples)
  contig <- prune_samples(!(sample_names(contig) %in% c("84642D5TT23","Mock")), contig)

contig_order= tax_glom(contig, "Order")
contig_family= tax_glom(contig, "Family")
contig_genus= tax_glom(contig, "Genus")
contig_species= tax_glom(contig, "Species")
```

## 7.4. Betalactams contig abundance 
```{r}
# Obtaining table for RA
  phy_relative_long <- psmelt(contig_genus)
  phy_relative_long <- phy_relative_long %>%
  group_by(Genus) %>%
  mutate(mean_relative_abund = mean(Abundance))
  phy_relative_long$Genus <- as.character(phy_relative_long$Genus)
  phy_relative_long$mean_relative_abund <- as.numeric(phy_relative_long$mean_relative_abund)
  phy_relative_long$Family[phy_relative_long$Abundance <= 1] <- "≤ 1 contigs"
# To make a 60 color palette
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
# To filter out taxa that 
filtered_data <- phy_relative_long %>%
  group_by(Genus) %>%
  filter(sum(Abundance > 0) >= 2) %>%
  ungroup()
# PLOT  
ggplot(filtered_data %>% filter(Abundance > 0), 
       aes(x = Sample, y = Genus, color = Genus)) +
  facet_wrap(~TimePoint + TRT, nrow = 1, scales = 'free_x') +
  geom_point(aes(size = Abundance), alpha = 0.9) +
  #geom_point(aes(size = Abundance, fill = Family), alpha = 0.9, shape = 21, stroke = 0.5, color = "black") +  # Black outline for points
  labs(x = "", y = "") +
  theme_minimal() + theme_pubclean() +
  scale_color_manual(values = col_vector) +
  #scale_fill_manual(values = col_vector)+
  scale_size_continuous(breaks = c(1, 50, 100)) +  # Customize point sizes in the legend
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
        axis.text.y = element_text(size = 10)) +
  theme(strip.text.x = element_text(size = 12)) +
    theme(axis.text.x = element_blank())+
theme(legend.position = 'bottom', legend.margin = margin(t = -25)) +
  guides(color = guide_legend(nrow = 2))+
  guides(color = FALSE) + guides(fill = FALSE)+
  theme(panel.grid.major.y = element_blank())

## To know how many contigs per group and time point
sum(taxa_sums(subset_samples(contig, TRT == "Treatment" & TimePoint == "D14")))
```

## 7.6. Figure 6
```{r}
# Getting ARG counts for CMY
CMY = prune_taxa("MEG_1989", AMR_div)
CMY = t(as.data.frame(otu_table(CMY)))

# Fig 6.A = CMY reads
## Load data and build phyloseq object
  otu_mat<- read_excel("MOBfinder_analysis_CMY.xlsx", sheet = "otu_arg")
  tax_mat<- read_excel("MOBfinder_analysis_CMY.xlsx", sheet = "tax_arg")
  samples_df<- read_excel("MOBfinder_analysis_CMY.xlsx", sheet = "samples")
  samples_df$TimePoint= factor(samples_df$TimePoint, levels=c("D0", "D5", "D14"))
  samples_df$SampleID= as.factor(samples_df$SampleID)
  samples_df$TRT= factor(samples_df$TRT, levels = c("0","1"), labels = c("Control", "Treatment"))
  otu_mat <- otu_mat %>%
    tibble::column_to_rownames("contig") 
  tax_mat <- tax_mat %>% 
    tibble::column_to_rownames("contig")
   samples_df <- samples_df %>% 
    tibble::column_to_rownames("DNA_ID")
  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(samples_df)
  CMY_reads <- phyloseq(OTU, TAX, samples)
CMY_reads = tax_glom(CMY_reads, "ARG_counts")
## Plotting 
  phy_relative_long <- psmelt(CMY_reads)
f6a <- phy_relative_long %>%
    ggplot(aes(x = Sample, y = log10(Abundance), fill = ARG_counts)) +
    geom_bar(stat = "identity",width = 1,color="white") +
    geom_bar(stat = "identity", alpha=1)+theme_bw()+
    # scale_y_continuous(breaks = seq(0, 16, by = 2)) +
    facet_wrap(~TimePoint+TRT, nrow=1, scales="free_x")+
    theme(axis.text.x = element_blank())+
  ylab("N° of CMY genes (log10)")+
  scale_fill_manual(values = c("grey"))

# Fig 6.B = CMY-contigs taxonomic ID
## Load data and build phyloseq object
  otu_mat<- read_excel("MOBfinder_analysis_CMY.xlsx", sheet = "otu_bacteria")
  tax_mat<- read_excel("MOBfinder_analysis_CMY.xlsx", sheet = "tax_bacteria")
  samples_df<- read_excel("MOBfinder_analysis_CMY.xlsx", sheet = "samples")
  samples_df$TimePoint= factor(samples_df$TimePoint, levels=c("D0", "D5", "D14"))
  samples_df$SampleID= as.factor(samples_df$SampleID)
  samples_df$TRT= factor(samples_df$TRT, levels = c("0","1"), labels = c("Control", "Treatment"))
  otu_mat <- otu_mat %>%
    tibble::column_to_rownames("contig") 
  tax_mat <- tax_mat %>% 
    tibble::column_to_rownames("contig")
   samples_df <- samples_df %>% 
    tibble::column_to_rownames("DNA_ID")
  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(samples_df)
  CMY_contig <- phyloseq(OTU, TAX, samples)
CMY_contig = tax_glom(CMY_contig, "CMY_contig_tax_ID")
## Plotting 
  phy_relative_long <- psmelt(CMY_contig)
f6b <- phy_relative_long %>%
    ggplot(aes(x = Sample, y = Abundance, fill = CMY_contig_tax_ID)) +
    geom_bar(stat = "identity",width = 1,color="white") +
    geom_bar(stat = "identity", alpha=1)+theme_bw()+
    #scale_y_continuous(breaks = seq(0, 16, by = 2)) +
    facet_wrap(~TimePoint+TRT, nrow=1, scales="free_x")+
    theme(axis.text.x = element_blank())+
  ylab("N° of CMY contigs")+
    ggthemes::scale_fill_tableau(palette = "Classic Green-Orange 6")

# Fig 6.C = MOB types
  otu_mat<- read_excel("MOBfinder_analysis_CMY.xlsx", sheet = "otu_mob")
  tax_mat<- read_excel("MOBfinder_analysis_CMY.xlsx", sheet = "tax_mob")
  samples_df<- read_excel("MOBfinder_analysis_CMY.xlsx", sheet = "samples")
  samples_df$TimePoint= factor(samples_df$TimePoint, levels=c("D0", "D5", "D14"))
  samples_df$SampleID= as.factor(samples_df$SampleID)
  samples_df$TRT= factor(samples_df$TRT, levels = c("0","1"), labels = c("Control", "Treatment"))
  otu_mat <- otu_mat %>%
    tibble::column_to_rownames("contig") 
  tax_mat <- tax_mat %>% 
    tibble::column_to_rownames("contig")
   samples_df <- samples_df %>% 
    tibble::column_to_rownames("DNA_ID")
  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(samples_df)
  CMY_mob <- phyloseq(OTU, TAX, samples)
CMY_mob = tax_glom(CMY_mob, "MOB_type")
## Plotting 
  phy_relative_long <- psmelt(CMY_mob)
f6c <- phy_relative_long %>%
    ggplot(aes(x = Sample, y = Abundance, fill = MOB_type)) +
    geom_bar(stat = "identity",width = 1,color="white") +
    geom_bar(stat = "identity", alpha=1)+theme_bw()+
    #scale_y_continuous(breaks = seq(0, 16, by = 2)) +
    facet_wrap(~TimePoint+TRT, nrow=1, scales="free_x")+
    theme(axis.text.x = element_blank())+
  ylab("N° of MOBs")+
    ggthemes::scale_fill_tableau(palette = "Classic 10 Medium")

f6a/f6b/f6c

```

