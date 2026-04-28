---
title: "TexasTech_Microbiome"
author: "G Diaz"
date: "2024-07-11"
output: html_document
editor_options: 
  chunk_output_type: console
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# 0. Load packages

package installation: **DO IT ONLY ONCE**
```{r installation, eval=FALSE}
# install.packages("BiocManager")
# BiocManager::install("igraph")
# BiocManager::install("phyloseq")
# BiocManager::install("decontam")
# install.packages("forcats")
# install.packages("readxl")
# install.packages("dplyr")
# install.packages("tibble")
# BiocManager::install("microbiome")
# BiocManager::install("metagenomeSeq")
# BiocManager::install("limma")
# library("devtools") # for maaslin3
# devtools::install_github("biobakery/maaslin3") # maaslin3 package
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
# install.packages("ggthemes")
# devtools::install_github("rdboyes/forester")
```

package loading
```{r packages, results="hide"}
library("phyloseq")
library("forcats")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("microbiome")   # Extra data manipulation
library("metagenomeSeq")  # Data normalization (?) + Statistical testing
# library("limma")          # Statistical testing
# library("statmod")      # Will work with limma for "duplicateCorrelation" function
library("psych")          # Data exploration
library("tidyverse")    # ggplot2, dplyr, tidyr, forcats, tibble, string, purrr, readr
library("vegan")
library("lme4")           # Statistical testing
library("lmerTest")       # Statistical testing
library("emmeans")        # Statistical testing
library("car")            # Statistical testing
library("table1")         # for doing tables
# library("lemon")          # for volcano plot
# library("ggrepel")      # for volcano plot
library("ggthemes")
library("decontam")
library("maaslin3")       # Differential abundance testing
library("forester")
library("extrafont")      # To Load fonts from your system (do it only once): font_import()
```

# 1. Load data
```{r}
# PS object was already built and stored as RDS file, so I just need to load it:

ps <- readRDS("Microbiome_SILVA138_2/ps.long.SILVA138_2.rds") # MODIFIED ON 01/17/2026

# factorize categorical data
View(sample_data(ps))
sample_data(ps)$BCS_MC <- as.numeric(sample_data(ps)$BCS_MC)
sample_data(ps)$TimePoint= as.factor(sample_data(ps)$TimePoint)
sample_data(ps)$TimePoint <- factor(sample_data(ps)$TimePoint, levels = c("D0", "D5", "D14"))
sample_data(ps)$TRT= as.factor(sample_data(ps)$TRT)
sample_data(ps)$...15= as.factor(sample_data(ps)$...15)  # wheter the value is a Body condition score OR a Metritis check score
sample_data(ps)$Dystocia= as.factor(sample_data(ps)$Dystocia) # wheter cow experienced or not dystocia
sample_data(ps)$RP= as.factor(sample_data(ps)$RP) # wheter the cow retained placenta or not
sample_data(ps)$Twin= as.factor(sample_data(ps)$Twin) # whether the cow had a twin parturition or not
sample_data(ps)$Stillbirth= as.factor(sample_data(ps)$Stillbirth) # Whether the cow had stillbirth parturirion 
sample_data(ps)$`1stAI outcome` = as.factor(sample_data(ps)$`1stAI outcome`) # Wheter the cow got pregnandt at the first service
sample_data(ps)$`Date of extraction`= as.factor(sample_data(ps)$`Date of extraction`) # to check batch effect later

```

# 2. Decontamination and data cleaning
```{r}
# 1. Decontamination with Decontam with FREQUENCY sample and default threshold (0.1) ##### 
contamdf.freq <- isContaminant(ps, method="frequency", conc="qPCR_copy")
head(contamdf.freq)
# Geting the number of contaminant ASVs
table(contamdf.freq$contaminant)
# knowing which ASVs are the contaminants 
head(which(contamdf.freq$contaminant))
# Plotting the model to confirm that my contaminants were REALLY contaminats 
set.seed(999)
plot_frequency(ps,taxa_names(ps)[sample(which(contamdf.freq$contaminant),4)],conc="qPCR_copy") 
# Plotting the model to compare a NOT-contaminant vs a Contaminant
plot_frequency(ps, taxa_names(ps)[c(1,552)], conc="qPCR_copy")
# Getting your DECONTAMINATED Phyloseq object 
ps_dc <- prune_taxa(!contamdf.freq$contaminant, ps) # you can get the phyloseq object of your contaminants just taking out the "!" symbol
# Getting the number and percentage of total reads removed with the contaminants
a=sum(sample_sums(ps)) # number of total reads for original ps object
b=sum(sample_sums(ps_dc)) # number of total reads for decontaminated ps object
c=a-b # number of reads removed after decomtamination
c/a*100 # proportion of reads removed after decontamination

# 2. getting rid of NAs taxa at kingdom level (unclassified ones) ##### 
tmp <- subset_taxa(ps_dc, !Kingdom %in% c("NA","Eukaryota"))

# 3. Getting rid of low quality sample 426D14TT3 #####
# Exclusion of sample ID  426D14TT3 (sample 426 D14) due to low-quality sequencing
# This sample had 1259 raw 16S sequences sequences, while the average was 171432 (SD=30185). The sample 426D14TT3 is lower than 2 SDs of the mean for denoised sequences. Also 426D14TT3 had 11 hits for ASVs while the average was 86197 (SD=16188). The number of ASVS detected in 426D14TT3 is lower than 2 SDs of the mean.
# Demonstration for number of denoised sequences:
a = data.frame(sample_data(ps)) %>%
  select("SampleID","input")
a = a[-11,]
b = a[-(90:93),]
describe(b$input)

# Demonstration for number of ASVs:
a = colSums(otu_table(ps))
a = a[-(91:94)]
a = sort(a)
a[1] # number of hits for sample 426D14TT3
b = a[-1]
describe(b) 

# Take out the 426D14TT3 to have the clean ps
ps_cln <- prune_samples(!(sample_names(tmp) %in% c("426D14TT3")), tmp)
# Take out controls
ps_samp <- prune_samples(!(sample_names(ps_cln) %in% c("EB1", "EB2","Mock1", "Mock2")), ps_cln)


```

# 3. Table 1 - demographics
```{r}
# Subset only TimePoint 0
tab1 = data.frame(subset(sample_data(ps), sample_data(ps)$TimePoint == "D0"))

# Subset variables to include in table
# TRT
# LACT #number of lactations
# DIM_1 # Days in milk at enrollment
# VLS # Vaginal laceration score at enrollment (0 to 2)
# BCS_MC # Body Condition Score at enrollment or Metritis check at day 5 or 14
# Temp_1 # Rectal temperature at enrollment
# DOPN # Days open (not pregnant)
# Milk4 # Milk yield in fourth month of lactation (kg/d)
# ECM4 # Energy corrected milk yield in fourth month of lactation (kg/d)
tmp = tab1 %>%
  select(`TRT`,`LACT`,`DIM_1`,`VLS`,`BCS_MC`,`Temp_1`,`DOPN`,`Milk4`,`ECM4`)
str(tmp)
# with variables treated as categorical
tab1 = table1(~ factor(LACT) + DIM_1 + factor(VLS) + BCS_MC + Temp_1 + DOPN + Milk4 + ECM4 | TRT, data=tmp)
# with variables treated as continous
tab1.1 = table1(~ LACT + DIM_1 + VLS + BCS_MC + Temp_1 + DOPN + Milk4 + ECM4 | TRT, data=tmp)

# If all they were continous variables, we may use wilcoxon test
wilcox.test(as.numeric(LACT)~TRT, data = tmp) 
wilcox.test(as.numeric(DIM_1)~TRT, data = tmp) 
wilcox.test(as.numeric(VLS)~TRT, data = tmp) 
wilcox.test(as.numeric(BCS_MC)~TRT, data = tmp) 
wilcox.test(as.numeric(Temp_1)~TRT, data = tmp) 
wilcox.test(as.numeric(DOPN)~TRT, data = tmp) 
wilcox.test(as.numeric(Milk4)~TRT, data = tmp) 
wilcox.test(as.numeric(ECM4)~TRT, data = tmp) 
# t.test(ECM4~TRT, data = tmp)  
# NOTE: If they were categorical variables, a 2x2 table would be used and a chisquare test (or fisher exact) can be ran.

# Extra added on 01/28/24 to check differences in temperature and metritis check score:
T_TRTd5 = c(38.9,	38.9,	39.5,	39.2,	39,	38.4,	39.3,	38.7,	39.3,	39.6,	39.2,	38.7,	38.6,	39.1,	38.5)
T_TRTd14 = c(38.2	,38.4,	39,	38.6,	38.6,	38.8,	38.6,	39.1,	39.1,	38.9,	38.9,	38.5,	38.5,	38.3,	38.4)

T_CTRd5 = c(39.2,	39.5,	39,	39.3,	39.1,	39.6,	38.8,	39.2,	39.2,	38.9,	38.8,	38.6,	38.7,	39.8,	38.8)
  
T_CTRd14 = c(38.8,	39.1,	38.3,	38.9,	39.4,	38.7,	39.4,	38.7,	38.7,	39.1,	38.6,	38.9,	39.1,	38.6,	38.9)

M_TRTd5 = c(3,	5,	5,	4,	5	,5,	5	,3,	5	,5	,3,	5	,4,	5,	5)
  
M_TRTd14 = c(4,	4,	3,	4,	4	,4,	5,	3	,4,	4,	3,	4,	4	,4,	4)

M_CTRd5 = c(3	,5,	4	,5,	4	,4,	5	,5,	5,	4	,5,	5	,5,	5	,5)

M_CTRd14 = c(4,	3	,2,	4	,4,	3	,4,	3	,4,	4,	4,	4	,5,	4	,2)
```

# 3.1. Sup. Table 1 - sequencing
```{r}
# subset variables for 2 groups in 3 timepoints to include in table with median + range:
# qPCR_copy
# input
Suptab1 = data.frame(subset(sample_data(ps), sample_data(ps)$FECAL == 1))
# Suptab1 <- Suptab1[rownames(Suptab1) != "426D14TT3", ] # without 426D14TT3
tmp = Suptab1 %>%
  select(`TRT`,`qPCR_copy`,`input`)
str(tmp)
# creating table of classified reads at any level and pasting with previous table to obtain the final one
tmp2 = colSums(otu_table(ps))
tmp2 = tmp2[-(91:94)]
# tmp2 = tmp2[-11] # without 426D14TT3
tmp2 = as.data.frame(tmp2)
Suptab1 = cbind.data.frame(tmp,tmp2)
# generating summary table
table1(~ qPCR_copy + input + tmp2| TRT, data=Suptab1)
# statistics
wilcox.test(qPCR_copy~TRT, data = Suptab1) 
wilcox.test(input~TRT, data = Suptab1) 
wilcox.test(tmp2~TRT, data = Suptab1) 
# t.test(tmp2~TRT, data = Suptab1)  

```

# 4. Controls description
```{r}
# obtaining ps object with only mock and negative controls
ps_neg <- prune_samples((sample_names(ps_cln) %in% c("EB1", "EB2")), ps_cln)
ps_pos <- prune_samples((sample_names(ps_cln) %in% c("Mock1", "Mock2")), ps_cln)
# agglomerating at genus level
glom_neg= tax_glom(ps_neg, "Genus")
glom_pos= tax_glom(ps_pos, "Genus")
# formatting for relative abundance plots
ps = ps_pos
  relative <- transform_sample_counts(ps, function(x) (x / sum(x))*100 )
  relative_long <- psmelt(relative)
  relative_long <- relative_long %>%
  group_by(Genus) %>%
  mutate(mean_relative_abund = mean(Abundance))
  relative_long$Genus <- as.character(relative_long$Genus)
  relative_long$mean_relative_abund <- as.numeric(relative_long$mean_relative_abund)
  relative_long$Genus[relative_long$Abundance < 0.01] <- "Others (< 0.01%)" 

```

## 4.1. Sup. Fig 2 - controls RA
```{r}
# NEGATIVE controls
SupF2A <- relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")
# POSITIVE controls
SupF2B <- relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")
# Getting RA for positive controls
table1(~ Abundance | Genus, data=relative_long)
```

# 5. Relative abundace 
```{r}
# CSS NORMALIZATION 
ps<-ps_samp
ps.metaseq <- phyloseq_to_metagenomeSeq(ps) # change here
ps.norm<- cumNorm(ps.metaseq, p=cumNormStat(ps.metaseq))
css.metaseq <- MRcounts(ps.norm, norm = TRUE)
ps_css <- merge_phyloseq(otu_table(css.metaseq,
taxa_are_rows=T),sample_data(ps),tax_table(ps))
ps_css
# AGGLOMERATION
ps_genus= tax_glom(ps_css, "Genus", bad_empty = c("NA", NA, "", " ", "\t"))
# ps_genus= tax_glom(ps_css, "Genus") # This will keep the NAs

# USE FOR ALPHA-DIV
ps_alpha = tax_glom(ps_samp, "Genus", bad_empty = c("NA", NA, "", " ", "\t")) # Genus level
```

## 5.1. Figure 1.A
```{r}
# Obtaining table for RA
  phy_relative <- transform_sample_counts(ps_genus, function(x) (x / sum(x))*100 )
  phy_relative_long <- psmelt(phy_relative)
  phy_relative_long <- phy_relative_long %>%
  group_by(Genus) %>%
  mutate(mean_relative_abund = mean(Abundance))
  phy_relative_long$Genus <- as.character(phy_relative_long$Genus)
  phy_relative_long$mean_relative_abund <- as.numeric(phy_relative_long$mean_relative_abund)
  phy_relative_long$Genus[phy_relative_long$Abundance < 5] <- "Low-abundance taxa (< 5%)"
  
# ordering Relative abundances (Y axis) for the plot ####
# Calculate total abundance for each Genus
phy_relative_long <- phy_relative_long %>%
  group_by(Genus) %>%
  mutate(TotalAbundance = sum(Abundance)) %>%
  ungroup()
# Reorder the Genus factor based on total abundance
phy_relative_long <- phy_relative_long %>%
  mutate(Genus = factor(Genus, levels = unique(Genus[order(TotalAbundance)])))
# Reordering the Samples (X axis) for the plot ####
sample_grouping <- phy_relative_long %>% 
  group_by(Sample) %>% 
  slice_max(order_by = Abundance) %>% 
  select(Genus, Sample) %>% 
  rename(peak_class = Genus)
phy_relative_long <- phy_relative_long %>% 
  inner_join(sample_grouping, by = "Sample") %>% 
  group_by(peak_class) %>% 
  mutate(rank = rank(Genus)) %>%  # rank samples at the level of each peak subtype
  mutate(Sample = reorder(Sample, -rank)) %>%   # this reorders samples
  ungroup()


## Plot
phy_relative_long$TRT <- factor(phy_relative_long$TRT, levels = c(0, 1), labels = c("Control", "Treatment"))

  Fig1A <- phy_relative_long %>%
    ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity",width = 1,color="white") +
    geom_bar(stat = "identity", alpha=1)+theme_bw()+
    facet_wrap(~TimePoint+TRT, nrow=1, scales="free_x")+
    theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
    ggthemes::scale_fill_tableau("Tableau 20")

```

## 5.2. Description of ASVs
```{r}
# Number of unique ASVs overall, at the phylum and genus level
ps_samp
 glom= tax_glom(ps_samp, "Genus", bad_empty = c("NA", NA, "", " ", "\t"))
# Number of reads classified overall, and at the phylum and genus level
  sum(taxa_sums(glom))
# Checking proportion of genus overall
  table1(~ Abundance | Genus, data=phy_relative_long)
# To get the proportion of low-abundance taxa (<5% abundance)
phy_relative_long %>%
  filter(Abundance < 5) %>%
  group_by(Sample) %>%
  summarise(sum = sum(Abundance, na.rm = TRUE), count = sum(Abundance > 0, na.rm = TRUE)) %>%
  summarise(
    sum_all = mean(sum, na.rm = TRUE),
    sd_all = sd(sum, na.rm = TRUE),
    count_all = mean(count, na.rm = TRUE),
    sdcount_all = sd(count, na.rm = TRUE))

```

# 6. Alpha diversity
```{r}
# Calculating at genus level
  div_16S= estimate_richness(ps_alpha, measures=c("Observed", "InvSimpson", "Shannon"))
  Alpha_genus= cbind(div_16S,(evenness(ps_alpha, 'pielou')))
  names(Alpha_genus) <- paste(names(Alpha_genus), "genus", sep = "_")
# Biding with dataframe
  df = cbind(data.frame(sample_data(ps_alpha)),Alpha_genus)
```

## 6.1. Figure 1.C, D, E
```{r}
# PLOT: Richness and Shannon at gene level
# Pre-processing
df$TRT <- factor(df$TRT, levels = c(0,1), labels = c("Control", "Treatment"))

# Plotting
Fig1C <- ggplot(df, aes(TimePoint,Observed_genus, fill=TRT))+
  geom_boxplot(position=position_dodge(0.8),lwd=1, alpha= 0.2)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "Time-point", y = "Richness", fill="Groups")

Fig1D <- ggplot(df, aes(TimePoint,Shannon_genus, fill=TRT))+
  geom_boxplot(position=position_dodge(0.8),lwd=1, alpha= 0.2)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "Time-point", y = "Shannon's Index", fill="Groups")

Fig1E <- ggplot(df, aes(TimePoint,pielou_genus, fill=TRT))+
  geom_boxplot(position=position_dodge(0.8),lwd=1, alpha= 0.2)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "Time-point", y = "Eveness Index", fill="Groups")

# STATISTICS: alpha ~ TRT*TimePoint + Classified + BCS + (1|SampleID) ####
# formatting BCS variable
df <- df %>%
  group_by(SampleID) %>%
  mutate(BCS = BCS_MC[TimePoint == "D0"][1]) %>%
  ungroup()
# # formatting classified reads
# class_reads = as.data.frame(colSums(otu_table(ps_samp)))
# class_reads <- class_reads %>% rename(classified = `colSums(otu_table(ps_samp))`)
# df<- cbind(class_reads,df)
#df$classified_resc <- df$classified/100 # Re-scaled variable

## Re-scaling sequencing depth variable (input)
df$input_resc <- df$input/100

# Model for Richness
gene_rich = lmer(Observed_genus ~ TRT*TimePoint + input_resc + BCS + (1|SampleID), data=df, REML = T)
  summary(gene_rich)
  Anova(gene_rich, type = "III")
  emmeans(gene_rich,pairwise~TRT|TimePoint)
  #emmeans(gene_rich,pairwise~TimePoint|TRT)
# Model for Shannon's index
gene_shan = lmer(Shannon_genus ~ TRT*TimePoint + classified_resc + BCS +(1|SampleID), data=df, REML = T)
  summary(gene_shan)
  Anova(gene_shan, type = "III")
  #confint(gene_shan)
  emmeans(gene_shan,pairwise~TRT|TimePoint)
  #emmeans(gene_shan,pairwise~TimePoint|TRT)
# Model for Eveness
gene_even = lmer(pielou_genus ~ TRT*TimePoint + classified_resc + BCS + (1|SampleID), data=df, REML = T)
  summary(gene_even)
  Anova(gene_even, type = "III")
  emmeans(gene_even,pairwise~TRT|TimePoint)
  #emmeans(gene_even,pairwise~TimePoint|TRT)

# DESCRIPTIVE STATS ####
describeBy(df$Observed_genus, df$TimePoint)
  describeBy(df$Observed_genus, list(df$TimePoint, df$TRT))
describeBy(df$Shannon_genus, df$TimePoint)
  describeBy(df$Shannon_genus, list(df$TimePoint, df$TRT))
describeBy(df$pielou_genus, df$TimePoint)
  describeBy(df$pielou_genus, list(df$TimePoint, df$TRT))
```

# 7. Beta diversity
```{r}
set.seed(143)
# genus.ord <- ordinate(ps_genus, "NMDS", "bray") # Not used because beta diversity analysis will be done at the ASV level

ASV.ord <- ordinate(ps_css, "NMDS", "bray")

```
## 7.1. Figure 1.B
```{r}
# Pre-processing
sample_data(ps_css)$TRT <- factor(sample_data(ps_css)$TRT, levels = c(0,1), labels = c("Control", "Treatment"))

# plot
Fig1B<- plot_ordination(ps_css, ASV.ord, type="sample", color="TRT")+
  geom_point(size=3, alpha=0.8) + theme_classic() + stat_ellipse() + facet_wrap(~TimePoint)+
  labs(color = "Groups")
Fig1B$layers <- Fig1B$layers[-1]
Fig1B

# STATISTICS: beta ~ TRT*TimePoint + BCS + (1|SampleID) ####
dist = vegdist(t(otu_table(ps_css)), method ="bray" )
adonis2(dist ~ TimePoint + TRT + BCS, data=df , by="margin") # I changed "classified_resc" variable was deleted because we already normalized (accounted for sequencing depth) by CSS
# STATISTICS DAY 0
df_0 = subset(df, df$TimePoint == "D0")
ps_css_0 = subset_samples(ps_css, TimePoint == "D0")
dist_0 = vegdist(t(otu_table(ps_css_0)), method ="bray" )
adonis2(dist_0 ~ TRT + BCS, data=df_0 , by="margin")
# STATISTICS DAY 5
df_5 = subset(df, df$TimePoint == "D5")
ps_css_5 = subset_samples(ps_css, TimePoint == "D5")
dist_5 = vegdist(t(otu_table(ps_css_5)), method ="bray" )
adonis2(dist_5 ~ TRT + BCS, data=df_5 , by="margin")
# STATISTICS DAY 14
df_14 = subset(df, df$TimePoint == "D14")
ps_css_14 = subset_samples(ps_css, TimePoint == "D14")
dist_14 = vegdist(t(otu_table(ps_css_14)), method ="bray" )
adonis2(dist_14 ~ TRT + BCS, data=df_14 , by="margin")

# BETADISPERSION
betadisp_0 <- betadisper(dist_0, df_0$TRT, type = "centroid")
permutest(betadisp_0, permutations = 999)
betadisp_5 <- betadisper(dist_5, df_5$TRT, type = "centroid")
permutest(betadisp_5, permutations = 999)
betadisp_14 <- betadisper(dist_14, df_14$TRT, type = "centroid")
permutest(betadisp_14, permutations = 999)

# STATISTICS CONTROL GROUP (all time points)
df_CTR = subset(df, df$TRT == "Control")
ps_css_CTR = subset_samples(ps_css, TRT == "Control")
dist_CTR = vegdist(t(otu_table(ps_css_CTR)), method ="bray" )
adonis2(dist_CTR ~ TimePoint + BCS, data=df_CTR , by="margin")
# STATISTICS TREATMENT GROUP (all time points)
df_TRT = subset(df, df$TRT == "Treatment")
ps_css_TRT = subset_samples(ps_css, TRT == "Treatment")
dist_TRT = vegdist(t(otu_table(ps_css_TRT)), method ="bray" )
adonis2(dist_TRT ~ TimePoint + BCS, data=df_TRT , by="margin")

# MORE PLOTS!!! Plot daily beta diversity
# Day 0
set.seed(143)
Fig1B_1 <- plot_ordination(ps_css_0, ordinate(ps_css_0, "NMDS", "bray"), color="TRT")+ geom_point(size=3, alpha=0.8) + theme_classic() + stat_ellipse() + labs(color = "Groups")
Fig1B_1$layers <- Fig1B_1$layers[-1]
Fig1B_1 # Stress:     0.1659788

# Day 5
set.seed(143)
# getting rid of sample 83092D5TT11 to improve visuals
ps_css_5_clean <- prune_samples(sample_names(ps_css_5) != "83092D5TT11", ps_css_5)
Fig1B_2 <- plot_ordination(ps_css_5_clean, ordinate(ps_css_5_clean, "NMDS", "bray"), color="TRT")+ geom_point(size=3, alpha=0.8) + theme_classic() + stat_ellipse() + labs(color = "Groups")
Fig1B_2$layers <- Fig1B_2$layers[-1]
Fig1B_2 # Stress:     0.1698644

# Day 14
set.seed(143)
Fig1B_3 <- plot_ordination(ps_css_14, ordinate(ps_css_14, "NMDS", "bray"), color="TRT")+ geom_point(size=3, alpha=0.8) + theme_classic() + stat_ellipse() + labs(color = "Groups")
Fig1B_3$layers <- Fig1B_3$layers[-1]
Fig1B_3 # Stress:     0.1656592 

```

# 8. Diff Abund Test
```{r}
# Get input_metadata from ps sample data
df_input_metadata = df
# change name of variable to "reads"
df_input_metadata <- df_input_metadata %>%
  rename(reads = fastqc.totalsequences)
# New (01/19/26) make row names the same for df_input_metadata & df_input_data so Maaslin3 can run
df_input_metadata$new_id <- paste0(df_input_metadata$SampleID, df_input_metadata$TimePoint, df_input_metadata$Ttrep)
df_input_metadata <- df_input_metadata %>% column_to_rownames("new_id")
# New (01/19/26) make sure the variable BCS is numeric and no string character (text)
df_input_metadata$BCS <- as.numeric(df_input_metadata$BCS)
# Get and format input_data from ps otu table ####
df_input_data = data.frame(otu_table(ps_alpha))
  #Create input tax table file from ps object:
  tax_table=tax_table(ps_alpha)[,"Genus"]
  # Put row names to a column called ASV in tax_table
  tax_table= data.frame(tax_table) %>% rownames_to_column(var='ASV')
  # Put row names to a column called ASV in df_input_data
  tmp = df_input_data %>% rownames_to_column(var='ASV')
  # Merging df_input_data (otu_table) with tax_table so we have genus names instead of ASV1, ASV2, etc.
  tmp2 = left_join(tmp, tax_table, by=c("ASV"))
  # Put column names 'Genus' to row names so we have the names of the Genus as row names
  #tmp = tmp2 %>% column_to_rownames(var='Genus')
  # I had non unique genus names so I did the following to fix it out and then I could move the names in the column to row names
  # checking which genus/ASVs were repeated 
  non_unique_rows <- tmp2 %>%
    group_by(Genus) %>%
    filter(n() > 1) %>%
    ungroup()
  # I had "UCG-001" genus repeated (ASV_202, ASV_294)
  # I just added 1, 2 and 3 at the end of the repated name
  tmp2$Genus <- ave(tmp2$Genus, tmp2$Genus, FUN = function(x) {
    if(length(x) > 1) {
      paste0(x, seq_along(x))
    } else {
      x
    }
  })
  # Put column names 'Genus' to row names so we have the names of the Genus as row names
  tmp3 = tmp2 %>% column_to_rownames(var='Genus')
  # Delete the ASV column as we already have the genus names as row names. Also Transpose so now the Genus are columns and the samples are rows
  tmp4 = subset(tmp3, select = -c(ASV)) 
  df_input_data=as.data.frame(t(tmp4))
  # Finally (and just for this dataset) I have an "X" at the beginning of the name of the rows, that I need to get rid of
  rownames(df_input_data) <- gsub("^X", "", rownames(df_input_data))



# Statistical model with interactions: otu ~ TRT * TimePoint + reads + BCS + (1|SampleID) #### 
fit_out <- maaslin3(input_data = df_input_data,
                    input_metadata = df_input_metadata,
                    output = 'Microbiome_SILVA138_2/model_interact',
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

# Statistical model without interactions: otu ~ TRT + TimePoint + reads + BCS + (1|SampleID) #### 
fit_out2 <- maaslin3(input_data = df_input_data,
                    input_metadata = df_input_metadata,
                    output = 'Microbiome_SILVA138_2/model_NOinteract',
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
```
## 8.1. Figure 2
```{r}
# FORMATTING ####
# After modyfing the names of the genera (because Maaslin3 adds points to the names where there are special characters like space, dashes, etc), exponentiating in excel and calculating the 95% CI by: coeff + 1.96*Standard error; coeff - 1.96*standard error. For log2 I exponentiated 2 to the coeff or 2 to the lower/upper bound.
# Read table
DAT_genus <- read.csv("Microbiome_SILVA138_2/model_interact/significant_results.tsv", sep = "\t")
# Formatting p-values FDR
DAT_genus$qval_individual[DAT_genus$qval_individual < 0.001] <- "p<0.001"
DAT_genus$qval_individual[DAT_genus$qval_individual < 0.01] <- "p<0.01"
DAT_genus$qval_individual[DAT_genus$qval_individual < 0.05] <- "p<0.05"
DAT_genus$qval_individual[DAT_genus$qval_individual < 0.1] <- "p<0.1"
# calculating the overall prevalence for significant Genera
DAT_genus$overall_prev <- round((DAT_genus$N.not.zero/DAT_genus$N)*100,1)
# calculating the overall abundance for significant Genera
  # Calculating the total sum of reads for all genera
  total_sum_all <- sum(df_input_data)
  # Extract the list of Genera 
  column_names <- DAT_genus$feature
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
  # Extra code because of the problem of duplication of Incertae Sedis. The overall abundance was calculated by hand (2955/5698965*100 = 0.052)
  # result_df$overall_abund[result_df$overall_abund == 0] <- 0.052
  # Merge dataframes
  DAT_genus = cbind(DAT_genus, result_df)

# Getting TABLE 2 ####
# Subset for plot for day 5 and 14
DAT_genus_5 <- subset(DAT_genus, DAT_genus$name == "TimePointD5:TRTTreatment" & DAT_genus$model == "abundance") # Changed "TimePointD5:TRTYes" for "TimePointD5:TRTTreatment" to align with new naming of the Treatment group
DAT_genus_14 <- subset(DAT_genus, DAT_genus$name == "TimePointD14:TRTTreatment" & DAT_genus$model == "abundance") # Changed "TimePointD5:TRTYes" for "TimePointD5:TRTTreatment" to align with new naming of the Treatment group
# selecting columns for the table
DAT_genus_5 <- DAT_genus_5 %>%
  select('feature','overall_prev','overall_abund','qval_individual','coef','lower_coef','upper_coef') %>%
  arrange(coef)
DAT_genus_14 <- DAT_genus_14 %>%
  select('feature','overall_prev','overall_abund','qval_individual','coef','lower_coef','upper_coef') %>%
  arrange(coef)
# Creating rows for subheadings Day 5 and Day 14
row_5 <- data.frame(feature= "Day 5", overall_prev = NA, overall_abund =  NA, qval_individual = NA, coef= NA, lower_coef= NA, upper_coef= NA)
row_14 <- data.frame(feature= "Day 14", overall_prev = NA, overall_abund =  NA, qval_individual = NA, coef= NA, lower_coef= NA, upper_coef= NA)
# Pasting all subsettings with subheadings
DAT_genus_5_14 <- rbind(row_5,DAT_genus_5,row_14,DAT_genus_14)
# Reanaming columns for the table
DAT_genus_5_14 <- DAT_genus_5_14 %>%
  rename(Genus = feature, `Overall prev. (%)` = overall_prev, `Overall abund. (%)` =  overall_abund, FDR = qval_individual)
# Indenting for Day5 and Day 14
DAT_genus_5_14$Genus <- ifelse(is.na(DAT_genus_5_14$`Overall abund. (%)`), 
                         DAT_genus_5_14$Genus,
                         paste0("   ", DAT_genus_5_14$Genus))
# Eliminating Eubacterium siraeum group on DAY 14 because it has qval_individual (or FDR) of 1...I don't know how it got there because it was suposed to filter FDR < 0.1
DAT_genus_5_14 <- DAT_genus_5_14[-13,]
# Loading font for table
loadfonts()
# Printing table 
tab2 <- forester(left_side_data = DAT_genus_5_14[,1:4],
           estimate = as.numeric(DAT_genus_5_14$coef),
           ci_low = as.numeric(DAT_genus_5_14$lower_coef),
           ci_high = as.numeric(DAT_genus_5_14$upper_coef),
           display = TRUE,
           #nudge_y = -.3,
           estimate_precision = 3,
           null_line_at = 0,
           xlim = c(-7, 4),
           #xbreaks = c(-1,0,1,2,3,4,5,6,7,8,9,10),
           #justify = c(0,0.5,0.5,1),
           estimate_col_name = "Log2(Fold change) (95% CI)",
           file_path = here::here("Microbiome_SILVA138_2/FP_genus_D5.png"),
           #render_as = "rmarkdown",
           font_family = "Arial")
```

## 8.2. Sup. Fig 3
```{r}
# DA over time CONTROL & TREATED GROUP #########
# Maaslin 3 DAT:
# IMPORTANT - Change dataframe here: Use df_CTR or df_TRT
df_input_metadata = df_TRT
# Then run!
df_input_metadata <- df_input_metadata %>%
  rename(reads = fastqc.totalsequences)
df_input_metadata$new_id <- paste0(df_input_metadata$SampleID, df_input_metadata$TimePoint, df_input_metadata$Ttrep)
df_input_metadata <- df_input_metadata %>% column_to_rownames("new_id")
df_input_metadata$BCS <- as.numeric(df_input_metadata$BCS)
# Get ps object with non-normalized counts for CTR and TRT from ps object used for alpha diversity
ps_alpha_CTR <- subset_samples(ps_alpha, TRT == 0)  
ps_alpha_TRT <- subset_samples(ps_alpha, TRT == 1)  
# IMPORTANT - Change ps object here: Use ps_alpha_CTR or ps_alpha_TRT
df_input_data = data.frame(otu_table(ps_alpha_TRT))
# Then run!
  tax_table=tax_table(ps_alpha)[,"Genus"]
  tax_table= data.frame(tax_table) %>% rownames_to_column(var='ASV')
  tmp = df_input_data %>% rownames_to_column(var='ASV')
  tmp2 = left_join(tmp, tax_table, by=c("ASV"))
  non_unique_rows <- tmp2 %>%
    group_by(Genus) %>%
    filter(n() > 1) %>%
    ungroup()
  tmp2$Genus <- ave(tmp2$Genus, tmp2$Genus, FUN = function(x) {
    if(length(x) > 1) {
      paste0(x, seq_along(x))
    } else {
      x
    }
  })
  tmp3 = tmp2 %>% column_to_rownames(var='Genus')
  tmp4 = subset(tmp3, select = -c(ASV)) 
  df_input_data=as.data.frame(t(tmp4))
  rownames(df_input_data) <- gsub("^X", "", rownames(df_input_data))
# Statistical model: otu ~ TimePoint + reads + BCS #### 
fit_out6 <- maaslin3(input_data = df_input_data,
                    input_metadata = df_input_metadata,
                    output = 'Microbiome_SILVA138_2/Timepoint_CONTROL',
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
fit_out7 <- maaslin3(input_data = df_input_data,
                    input_metadata = df_input_metadata,
                    output = 'Microbiome_SILVA138_2/Timepoint_TREATED',
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

# FOREST PLOT ####
DAT_genus <- read.csv("Microbiome_SILVA138_2/Timepoint_TREATED/significant_results.tsv", sep = "\t")
DAT_genus <- subset(DAT_genus, DAT_genus$qval_individual <= 0.1)
DAT_genus$qval_individual[DAT_genus$qval_individual < 0.001] <- "p<0.001"
DAT_genus$qval_individual[DAT_genus$qval_individual < 0.01] <- "p<0.01"
DAT_genus$qval_individual[DAT_genus$qval_individual < 0.05] <- "p<0.05"
DAT_genus$qval_individual[DAT_genus$qval_individual < 0.1] <- "p<0.1"
DAT_genus$overall_prev <- round((DAT_genus$N.not.zero/DAT_genus$N)*100,1)
  total_sum_all <- sum(df_input_data)
  column_names <- DAT_genus$feature
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
  DAT_genus = cbind(DAT_genus, result_df)
DAT_genus_5 <- subset(DAT_genus, DAT_genus$name == "TimePointD5" & DAT_genus$model == "abundance")
DAT_genus_14 <- subset(DAT_genus, DAT_genus$name == "TimePointD14" & DAT_genus$model == "abundance") 
DAT_genus_5 <- DAT_genus_5 %>%
  select('feature','overall_prev','overall_abund','qval_individual','coef','lower_coef','upper_coef') %>%
  arrange(coef)
DAT_genus_14 <- DAT_genus_14 %>%
  select('feature','overall_prev','overall_abund','qval_individual','coef','lower_coef','upper_coef') %>%
  arrange(coef)
row_5 <- data.frame(feature= "Day 5 vs Day 0", overall_prev = NA, overall_abund =  NA, qval_individual = NA, coef= NA, lower_coef= NA, upper_coef= NA)
row_14 <- data.frame(feature= "Day 14 vs Day 0", overall_prev = NA, overall_abund =  NA, qval_individual = NA, coef= NA, lower_coef= NA, upper_coef= NA)
DAT_genus_5_14 <- rbind(row_5,DAT_genus_5,row_14,DAT_genus_14)
DAT_genus_5_14 <- DAT_genus_5_14 %>%
  rename(Genus = feature, `Overall prev. (%)` = overall_prev, `Overall abund. (%)` =  overall_abund, FDR = qval_individual)
DAT_genus_5_14$Genus <- ifelse(is.na(DAT_genus_5_14$`Overall abund. (%)`), 
                         DAT_genus_5_14$Genus,
                         paste0("   ", DAT_genus_5_14$Genus))
loadfonts()
# Printing table 
tab4 <- forester(left_side_data = DAT_genus_5_14[,1:4],
           estimate = as.numeric(DAT_genus_5_14$coef),
           ci_low = as.numeric(DAT_genus_5_14$lower_coef),
           ci_high = as.numeric(DAT_genus_5_14$upper_coef),
           display = TRUE,
           #nudge_y = -.3,
           estimate_precision = 3,
           null_line_at = 0,
           #xlim = c(-7, 4),
           #xbreaks = c(-1,0,1,2,3,4,5,6,7,8,9,10),
           #justify = c(0,0.5,0.5,1),
           estimate_col_name = "Log2(Fold change) (95% CI)",
           file_path = here::here("Microbiome_SILVA138_2/FP_genus_TRT.png"),
           #render_as = "rmarkdown",
           font_family = "Arial")
```
