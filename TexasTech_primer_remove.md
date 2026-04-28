---
title: "TexasTech_16S"
author: "G Diaz"
date: "2023-05-24"
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
# For R version 4.3 and greater you need Bioconductor 3.17 or greater
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")
# will download latest dada2 version, in this case 1.28 and is compatible with arm64 architecture of M1 chips of Apple
BiocManager::install("dada2")
# will download "ShortRead" latest version
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ShortRead")
# will download "Biostrings" latest version
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
# installation of cutadapt for macOS with M1 processor. RUN IN TERMINAL
  sudo python3 -m pip install cutadapt # will install in "/usr/local/bin" directory
```

package loading
```{r packages, results="hide"}
library(dada2); packageVersion("dada2")           # General analysis 
library(ShortRead); packageVersion("ShortRead")   # For primer removal
library(Biostrings); packageVersion("Biostrings") # For primer removal
# give path for cutadapt
cutadapt <- "/usr/local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R
```

# 1. Load fastq files
```{r fastq files}
path <- "/Users/gerardodiaz/Desktop/Gerardo/UMN/Side_projects/TexasTech_resistome/Noyes_Project_051_ReSeq"  ## CHANGE ME to the directory containing the fastq files.
list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))
# plot quality for 16 first raw FORWARD files
plotQualityProfile(fnFs[1:16])
# plot quality for 16 first raw REVERSE files
plotQualityProfile(fnRs[1:16])
```

# 2. Remove ambiguous bases and make sure that your primers are in the sequences 
```{r abiguous bases}
# Remove ambiguous bases (N)
	fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
	fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
	filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
# Getting the 4 orientations for each of your 2 primers 
  allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
  }
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  FWD.orients
# Looking for how many sequences have your primers forward and reverse (and their reverse complements)
  primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
```

# 3. Remove primers
```{r primer removal}
# Identify your primers (asign your primers to an object) 
  FWD <- "CCTACGGGAGGCAGCAG"  ## CHANGE ME to your forward primer sequence
	REV <- "GGACTACHVGGGTWTCTAAT"  ## CHANGE ME to your reverse primer sequence
# Remove them 
	path.cut <- file.path(path, "cutadapt")
	if(!dir.exists(path.cut)) dir.create(path.cut)
	fnFs.cut <- file.path(path.cut, basename(fnFs))
	fnRs.cut <- file.path(path.cut, basename(fnRs))
	FWD.RC <- dada2:::rc(FWD)
	REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
	R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
	R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
		for(i in seq_along(fnFs)) {
		system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,"-m", 1:1, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
		}
# Sanity check to confirm primer removal
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
```

