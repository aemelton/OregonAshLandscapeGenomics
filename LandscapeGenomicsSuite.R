#### Header ####
# Anthony E. Melton
# Postdoctoral Scholar
# Hamilton Lab, Penn State University
# Landscape genomics script suite, 2023
### ~~~ ###

#### Description #### 
# The goal of this set of scripts is to analyse a set of SNP data derived from
# *Fraxinus latifolia* samples across the entire species range to inform 
# conservation and restoration efforts. This script contains code for analyses
# to estimate ploidy levels, genetic diversity, population genetic structure, 
# gene-environment interactions, genetic offset under climate change, and 
# spatial visualizations of estimated migration rates. This script does require
# some software outside of R, though commands are available within this script
# for programs that can be called from a command-line interface. Please note,
# the settings used in the analyses in this script should be adjusted for each
# research case. Please follow best practices for each analysis. Links to
# vignettes have been provided, when possible.
### ~~~ ###

#### Other required software #### 
# Plink https://www.cog-genomics.org/plink/
# Plink2  https://www.cog-genomics.org/plink/2.0/
# fastStructure https://rajanil.github.io/fastStructure/ or https://gitlab.com/StuntsPT/Structure_threader
# EEMS  https://github.com/dipetkov/eems
# vcftools  https://vcftools.sourceforge.net/

#### Public data downloads #### 
# climateNA https://climatena.ca/ for climatic data
# Fraxinus latifolia species distribution shapefile https://databasin.org/datasets/50fbdf842f584e4ca628433ff5f9cd86/

#### Test data sets #### 
# To aid in testing this script and learning how to set up the analyses, I have
# added a test data set that can be downloaded from
# https://github.com/aemelton/TestData. It contains 25 samples and 550 loci. 
# There is also a csv file containing latitude andlongitude coordinates. Please 
# note, THE VCF AND CSV DATA ARE NOT LINKED IN THE REAL WORLD. They are just 
# here for test purposes and do not represent real linked data.

#### Setting a color palette ####
# Color palette w/'force vector'  http://medialab.github.io/iwanthue/

# Color palette for the test vcf. Plots currently call on the "my.cols.test"
# object. Please check your palette before plotting!
my.cols.test <- c("#78a477","#b738bd","#49cb4a","#1653d8","#d0c90b")

# Color palettes used in Melton et al., In Prep, for Fraxinus latifolia analyses.
my.cols <- c("#78a477","#b738bd","#49cb4a","#1653d8","#d0c90b","#742ea9","#5da600",
             "#901493","#009c25","#f62893","#00b05e","#f82363","#42deae","#ce005f",
             "#00a77a","#c7091c","#00bbe4","#c43200","#2c9fff","#e0ae00","#9295ff",
             "#c7ce2e","#304e9d","#ff9b1a","#7fbeff","#f44739","#009aa9","#ffad39",
             "#bcacff","#a6a600","#a60076","#a9d46b","#ff77c7","#006f2e","#ff6f97",
             "#476d00","#fdabf3","#9e9300","#823573","#c9cc60","#ff616e","#018667",
             "#ff9d45","#02978f","#933d00","#83d7ae","#943029","#aad37b","#853d3b",
             "#ebc152","#255e2d","#ff9396","#846600","#af6871","#ae6800","#d1c38b",
             "#7b4526","#f8ba78","#77703c","#ff9074","#ffac8d")

# A subset of the color palette for analyses with outlier populations removed
my.cols.north <- c("#78a477","#49cb4a","#1653d8","#d0c90b","#742ea9","#5da600",
                   "#901493","#009c25","#f62893","#00b05e","#f82363","#42deae","#ce005f",
                   "#00a77a","#c7091c","#00bbe4","#c43200","#2c9fff","#e0ae00","#9295ff",
                   "#c7ce2e","#304e9d","#ff9b1a","#7fbeff","#f44739","#009aa9","#ffad39",
                   "#bcacff","#a6a600","#a60076","#a9d46b","#006f2e","#ff6f97",
                   "#476d00","#fdabf3","#9e9300","#823573","#ff616e","#018667",
                   "#ff9d45","#02978f","#933d00","#83d7ae","#943029","#aad37b","#853d3b",
                   "#ebc152","#846600","#af6871","#d1c38b",
                   "#7b4526","#f8ba78","#77703c","#ff9074","#ffac8d")


#### Install required libraries that are not on CRAN ####
# Apex, required by strataG https://cran.r-project.org/src/contrib/Archive/apex/

if (any(installed_packages == FALSE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  # BiocManager::install("biomaRt")
  # BiocManager::install("Biostrings")
  # BiocManager::install("DECIPHER")
  # BiocManager::install("LEA")
  # BiocManager::install("msa")
  BiocManager::install("SNPRelate")
  # BiocManager::install("snpStats")
  
  if (!require('devtools')) install.packages('devtools')
  #devtools::install_github("thomasp85/patchwork")
  devtools::install_github('ericarcher/strataG', build_vignettes = FALSE)
  
  remotes::install_github('royfrancis/pophelper')
  
}

packages <- c("adegenet",
              "biomartr",
              "Biostrings",
              "boot",
              "colorRamps",
              "dartR",
              "DECIPHER",
              "dplyr",
              "elevatr",
              "gbs2ploidy",
              "gdm",
              "ggplot2",
              "ggpubr",
              "gtools",
              "LEA",
              "mapdata",
              "maps",
              "msa",
              "MultiPhen",
              "patchwork",
              "pegas",
              "pinfsc50",
              "poppr",
              "rEEMSplots",
              "rentrez",
              "rgdal",
              # "RgoogleMaps",
              "rworldmap",
              "related", # Must install related using a tar.gz file
              "simpleboot",
              "SNPRelate",
              "snpStats",
              "strataG",
              "stringr",
              "svMisc",
              "usmap",
              "vcfR")

#### Load libraries ####

library(colorRamps)
library(data.table)
library(elevatr)
library(gbs2ploidy)
library(gdm)
library(geosphere)
# library(ggcorrplot)
library(ggforce)
library(ggmap)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggsci)
library(grid)
library(gridExtra)
library(gridGraphics)
library(hierfstat)
library(mapdata)
library(maps)
library(patchwork)
library(pophelper)
library(poppr)
library(reshape2)
library(scales)
library(scatterpie)
library(SNPRelate)
library(strataG)
library(tidyverse)
library(vcfR)
library(vegan)
library(wesanderson)

#### Working environment set-up ####
# These scripts are set-up to work over apre-determined folder hierarchy. Please
# read over the script and change folder paths as needed. All data and output
# folders should be within one shared "project" folder.

project.folder <- getwd()
setwd(project.folder)

##### Check folders! #####
if (file.exists("0_Data/")){ 
  print("The data folder already exists!")
} else { 
  # create a new sub directory inside 
  dir.create("0_Data/") 
}

if (file.exists("1_gbs2ploidy/")){ 
  print("The gbs2ploidy folder already exists!")
} else { 
  # create a new sub directory inside 
  dir.create("1_gbs2ploidy/") 
}

if (file.exists("2_SNP_PCA/")){ 
  print("The snp pca folder already exists!")
} else { 
  # create a new sub directory inside 
  dir.create("2_SNP_PCA/") 
}

if (file.exists("3_fastStructure/")){ 
  print("The fastStructure folder already exists!")
} else { 
  # create a new sub directory inside 
  dir.create("3_fastStructure/") 
}

if (file.exists("4_GeneticDiversity/")){ 
  print("The genetic diversity folder already exists!")
} else { 
  # create a new sub directory inside 
  dir.create("4_GeneticDiversity/") 
}

if (file.exists("5_Kinship/")){ 
  print("The kinship folder already exists!")
} else { 
  # create a new sub directory inside 
  dir.create("5_Kinship/") 
}

if (file.exists("6_EffectivePopulationSize/")){ 
  print("The effective population size folder already exists!")
} else { 
  # create a new sub directory inside 
  dir.create("6_EffectivePopulationSize/") 
}

if (file.exists("7_GEA/")){ 
  print("The gea folder already exists!")
} else { 
  # create a new sub directory inside 
  dir.create("7_GEA/") 
}

if (file.exists("8_GDM/")){ 
  print("The gdm folder already exists!")
} else { 
  # create a new sub directory inside 
  dir.create("8_GDM/") 
}

if (file.exists("9_EEMS/")){ 
  print("The eems folder already exists!")
} else { 
  # create a new sub directory inside 
  dir.create("9_EEMS/") 
}

#### Download test data ####
# This is a test data set for running through this script to help learn the 
# functions and outputs. Please reach out (https://github.com/aemelton) if you
# run into any hurdles.
setwd("0_Data/")
download.file(url = "https://github.com/aemelton/TestData/blob/main/SampleData_5pop_25samp_550snp.csv",
              destfile = "SampleData_5pop_25samp_550snp.csv")
download.file("https://github.com/aemelton/TestData/blob/main/TestVCF_5pop_25samp_550snp.vcf.gz",
              destfile = "TestVCF_5pop_25samp_550snp.vcf.gz")

#### Load population data ####
# "all" or "a" pops are the total 61 populations which were initially included for the study.
# "north" or "n" are pops that were included in the final round of analyses, excluding those
# at the southern terminus of the species range which were likely polyploid hybrids.

# If you do not have a population summary file, this section will take long/lat
# means and count samples per population for you
setwd(project.folder)
sampleData <- read.csv("0_Data/SampleData_5pop_25samp_550snp.csv")
longitude <- aggregate(sampleData$Longitude, list(sampleData$Population), FUN=mean)
latitude <- aggregate(sampleData$Latitude, list(sampleData$Population), FUN=mean)
sample.counts <- data.frame(table(sampleData$Population))

Pop_Sum <- data.frame("Population" = unique(sampleData$Population),
                      "State" = c("Q", "Q", "Q", "R", "R"),
                      "Longitude" = longitude$x,
                      "Latitude" = latitude$x,
                      "n" = sample.counts$Freq)
head(Pop_Sum)
write.csv(x = Pop_Sum, file = "0_Data/Pop_Sum.csv", row.names = F)

# Pop_Sum == Population summary
Pop_Sum <- read.csv("0_Data/Pop_Sum.csv")

sampleData <- read.csv("0_Data/SampleData_5pop_25samp_550snp.csv", header = T)
nrow(sampleData) # Total number of samples

# Requires a CSV file with population ID and sample size per population.
setwd(project.folder)
Pop_Sum <- read.csv("1_FRLA_data/all/data/allPop_Sum.csv", header = T) 
nrow(Pop_Sum) # Total number of populations
sort(Pop_Sum$Pop) # List of population names

sampleData <- read.csv("1_FRLA_data/all/data/allPop_ID_Sum.csv", header = T)
nrow(sampleData) # Total number of samples

#### Climatic data processing ####
Pop_Sum <- read.csv("1_FRLA_data/north/data/nPop_Sum.csv") # load pop data

# Set working directory to where the climate data are stored. I generally have
# a separate folder for these data, rather than a sub-folder within a project
# folder, becuase I use environmental data for a lot of different projects.

longlat.df.samp <- sampleData[,c("Longitude","Latitude")]
colnames(longlat.df.samp) <- c("x", "y")
ele.df.samp <- get_elev_point(locations = longlat.df.samp, prj = "EPSG:4326")

longlat.df.pop <- Pop_Sum[,c("Longitude","Latitude")]
colnames(longlat.df.pop) <- c("x", "y")
ele.df.pop <- get_elev_point(locations = longlat.df.pop, prj = "EPSG:4326")

setwd("~/Dropbox/Environmental_Data/ClimateNA/") # Load the environmental data
file.list <- list.files()
envStack <- raster::stack(file.list)
envDF.samp <- raster::extract(x = envStack, y = sampleData[,c("Longitude","Latitude")])
envDF.samp <- cbind(sampleData, envDF.samp)
Elevation <- data.frame("Elevation" = ele.df.samp$elevation)
envDF.samp <- cbind(envDF.samp, Elevation)
head(envDF.samp)

envDF.pop <- raster::extract(x = envStack, y = Pop_Sum[,c("Longitude","Latitude")])
envDF.pop <- cbind(Pop_Sum[, 1:2], envDF.pop)
Elevation <- data.frame("Elevation" = ele.df.pop$elevation)
envDF.pop <- cbind(envDF.pop, Elevation)
head(envDF.pop)

setwd(project.folder)
write.csv(envDF.samp, "0_Data/env_samp.csv", row.names = F) # Save table of env data
write.csv(envDF.pop, "0_Data/env_pop.csv", row.names = F) # Save table of env data

envDF.samp <- read.csv("0_Data/env_samp.csv")
envDF.samp <- left_join(sampleData, envDF.samp)

envDF.pop <- read.csv("0_Data/env_pop.csv")
envDF.pop <- left_join(Pop_Sum, envDF.pop)

##### Correlation tests #####
env_keep <- as.data.frame(scale(envDF.pop[c("Latitude", "Longitude",
                                           "Elevation", "FFP", "MAP", "MAT", "RH", "PAS", "TD")]))

corr <- cor(env_keep, use = "complete.obs", method = "pearson")
# ggcorrplot(corr)

corr_NA <- round(corr, 2)
corr_NA[which(abs(corr) > 0.75)] <- NA
cat("\n\nvariables kept:\n")

abs(corr_NA)

env_final <- data.frame(State = envDF.pop$State,
                        Pop = envDF.pop$Pop,
                        env_keep)

write.csv(env_final, "0_Data/env_final.csv", row.names = F)

##### merge with Pop_Sum #####
# envDF.samp <- read.csv("0_Data/env_samp.csv")
# env_indv_final <- envDF.samp[c("State", "Population", "sample", "Latitude",
#                              "Longitude", "Elevation", "FFP", "MAP", "MAT", "RH", "PAS",
#                              "TD")]

Pop_Sum <- read.csv("0_Data/Pop_Sum.csv")
Pop_env <- left_join(Pop_Sum, envDF.pop)

write.csv(Pop_env, "0_Data/Pop_Sum_env.csv", row.names = F)
##### Make environmental distance matrix #####

env_final <- read.csv("0_Data/env_final.csv")
EnvDist <- as.matrix(dist(env_final[, -c(1:4)]))
write.table(EnvDist, "0_Data/EnvDist.txt", row.names = F, col.names = F,
            quote = F)
####

#### Generate a color-coded table of sampling; Table S1 ####
# Set working directory and load data
setwd(project.folder)
Pop_Sum <- read.csv("0_Data/Pop_Sum.csv")

# Sort the populations by decreasing latitude
Pop.df <- Pop_Sum[order(Pop_Sum$Population, decreasing = F),]
Pop.df$Col <- my.cols.test

# Load population data
pop.tab <- read.csv("0_Data/Pop_Sum_env.csv", header = T) # This table has
# environmental data within it. These data are downloaded from 'climateNA' and
# extracted at each population using the raster R package. Please see 'Climate
# data processing' above.

###### If pop file has a *notes* column ######
# pop.tab <- pop.tab[,-ncol(pop.tab)] # Drop the last column, which has notes and not data

colnames(pop.tab) <- c("State",
                       "Population",
                       "Longitude",
                       "Latitude",
                       "n",
                       "FFP",
                       "MAP",
                       "MAT",
                       "PAS",
                       "RH",
                       "TD",
                       "Elevation (m)")

pop.tab <- pop.tab[order(pop.tab$Pop, decreasing = F),]

pop.tab$'Color code' <- my.cols.test # add color palette to the data frame to make
# keeping all the pieces together easier

# Order by latitude
pop.tab <- pop.tab[order(pop.tab$Latitude, decreasing = T),]
#Then turn it back into a factor with the levels in the correct order
pop.tab$Population <- factor(pop.tab$Population, levels=unique(pop.tab$Population))

pop.tab$`Elevation (m)` <- signif(x = pop.tab$`Elevation (m)`, digits = 2)
pop.tab$Latitude <- signif(x = pop.tab$Latitude, digits = 5)
pop.tab$Longitude <- signif(x = pop.tab$Longitude, digits = 5)

#
tt <- ttheme_minimal(
  core=list(bg_params = list(fill = pop.tab$`Color code`, col=NA),
            fg_params=list(fontface=1)),
  colhead=list(fg_params=list(col="black", fontface=1L)),
  rowhead=list(fg_params=list(col="black", fontface=1L)))

##### Plot; Table S1 #####
pdf("0_Data/Population_and_ColorCode.pdf", height = 18, width = 12) # Check dimensions for your data!
grid.arrange(tableGrob(pop.tab, theme=tt, rows = NULL))
dev.off()
####


#### Get summary statistics and plots for sexed samples used in analyses ####
# This section requires a csv with sample data, including the sex of the samples
# for dioecious species.

# Test samples currently are not dioecious, but evolve to be as such in the future.
# Load sample data
final.pops <- read.csv("1_FRLA_data/north/data/nPop_Sum.csv")
final.samps <- read.csv("1_FRLA_data/north/data/nPop_ID_Sum.csv")
all.samps <- read.csv("1_FRLA_data/aem2022_FRALAT_NaturalPopulations_DNAExtractions.csv")

nrow(all.samps)

# Subset the data to include sex info for only the samples used in analyses
all.samps.sub <- all.samps[which(gsub(pattern = "-", replacement = "_", x = all.samps$Individual_ID) %in% final.samps$All),]
all.samps.sub.sub <- all.samps.sub[which(all.samps.sub$Pop.Code %in% final.pops$Pop),]
nrow(all.samps.sub)
nrow(all.samps.sub.sub)

all.samps.sub.sexd <- all.samps.sub[which(all.samps.sub$Sex.M_or_F. != "<NA>"),]
tail(all.samps.sub.sexd)
nrow(all.samps.sub.sexd)
unique(all.samps.sub.sexd$Pop.Code)

all.samps.sub.sexd.sub <- all.samps.sub.sexd[which(order(all.samps.sub.sexd$Individual_ID) %in% order(final.samps$All)),]
all.samps.sub.sexd.sub <- all.samps.sub.sexd.sub[order(all.samps.sub.sexd.sub$Pop.Code, decreasing = T),]
all.samps.sub.sexd.sub
nrow(all.samps.sub.sexd.sub)

# Build a table with samples, their sex, and population
sex.table <- data.frame(all.samps.sub.sexd.sub %>%
  group_by(Pop.Code) %>%
  summarise(
    "Sexed" = n(),
    across(Sex.M_or_F., list(M = ~ sum(. == "M"),
                          F  = ~ sum(. == "F")))))

sex.table

# Save file
write.csv(sex.table, file = "SexTable.csv", row.names = F)

# Load file
sex.table <- read.csv("SexTable.csv")
head(sex.table)

final.pops.sub <- final.pops[which(final.pops$Pop %in% sex.table$Pop.Code),]
final.pops.sub <- final.pops.sub[order(final.pops.sub$Pop, decreasing = F),]
sex.table$TotalSampled <- final.pops.sub$n
colnames(sex.table) <- c("Population", "n sexed", "male", "female", "n sampled")

Pop.df.sub <- Pop.df[which(Pop.df$Pop %in% sex.table$Population),]

sex.table$PercentSexed <- sex.table$`n sexed`/sex.table$`n sampled`

sex.table$Latitude <- Pop.df.sub$Latitude
sex.table$Col <- Pop.df.sub$Col

sex.table <- sex.table[order(sex.table$Latitude, decreasing = T),]
#Then turn it back into a factor with the levels in the correct order
sex.table$Population <- factor(sex.table$Population, levels=unique(sex.table$Population))

dfm <- melt(sex.table[,c('Population','female','male')], id.vars = 1)

# Plot the percent of samples for each population that had sex information
percent.plot <- ggplot(sex.table, aes(x = Population, y = PercentSexed*100)) + 
  geom_bar(stat = 'identity', fill = sex.table$Col, col = sex.table$Col) + 
  geom_hline(yintercept = mean(sex.table$PercentSexed)*100, color="grey", linetype="dashed") +
  ylab("Percent of samples sexed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none")

# Plot F:M for each population
sex.plot <-  ggplot(dfm,aes(x = Population, y = value)) + 
  geom_bar(aes(fill = variable), stat = "identity", position = "dodge") + 
  ylab("Number of samples") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill=guide_legend(title="n sampled per sex")) +
  theme(legend.position = c(0.15, 0.85)) # c(0.2, 0.8)

##### Merged plot; Figure S7 #####
pdf("3_Plots/SexedPlantPlot.pdf")
ggarrange(percent.plot, sex.plot,
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()
####
 
#### Ploidy level estimation ####
# gbs2ploidy #
##### Modify function to write out the stdev from estploidy PCA ####
estploidy2 <- function (alphas = NA, het = NA, depth = NA, train = FALSE, pl = NA, 
          set = NA, nclasses = 2, ids = NA, pcs = 1:2) 
{
  outH <- lm(het ~ depth)
  H <- outH$residuals
  Nind <- length(alphas)
  Nprops <- dim(alphas[[1]])[1]
  ad <- matrix(NA, nrow = Nind, ncol = Nprops)
  for (i in 1:Nind) {
    ad[i, ] <- alphas[[i]][, 3]
  }
  datamatrix <- cbind(H, ad)
  pcout <- prcomp(datamatrix, center = TRUE, scale = TRUE)
  if (train == FALSE) {
    kout <- kmeans(x = pcout$x[, pcs], centers = nclasses)
    lout <- lda(kout$cluster ~ pcout$x[, pcs], CV = TRUE, 
                prior = rep(1/nclasses, nclasses))
    pp <- cbind(ids, lout$posterior)
  }
  if (train == TRUE) {
    df <- data.frame(pl = pl, pcout$x[, pcs])
    ltrn <- lda(pl ~ ., df, subset = set, prior = rep(1/nclasses, 
                                                      nclasses))
    lout <- predict(object = ltrn, newdata = df[-set, ])
    pp <- lout$posterior
  }
  out <- list(pp = pp, pcwghts = pcout$rotation, pcsdev = pcout$sdev, pcscrs = pcout$x)
  out
}

setwd(project.folder)

# load vcf file
# This takes a long time. The "test data" set has been included in the current
# code for demonstation purposes.
# vcf_file <- "1_FRLA_data/all/FRLA_all_NoDGR.recode.vcf.gz"

vcf_file <- "0_Data/TestVCF_5pop_25samp_550snp.vcf.gz"

#SampleData <- read.csv("1_FRLA_data/all/data/all_pca_df_AEM.csv")
SampleData <- read.csv("0_Data/SampleData_5pop_25samp_550snp.csv")
nrow(sampleData)

##### Read data into memory #####
vcf <- read.vcfR(vcf_file, verbose = FALSE)
vcf

##### data prep #####
ad <- extract.gt(vcf, element = 'AD')
gt <- extract.gt(vcf, element = 'GT')
hets <- is_het(gt)
is.na( ad[ !hets ] ) <- TRUE

allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)

ad1 <- allele1 / (allele1 + allele2)
ad2 <- allele2 / (allele1 + allele2)

length(ad1[!is.na(ad1)])

ad1 <- matrix(na.omit(ad1))
ad2 <- matrix(na.omit(ad2))

##### estprops #####
# Parameters used in research; plugging in fewer reps for testing.
# props <- estprops(cov1=allele1,
#                   cov2=allele2,
#                   mcmc.steps=10000,
#                   mcmc.burnin=1000,
#                   mcmc.thin=100)

props <- estprops(cov1=allele1,
                  cov2=allele2,
                  mcmc.steps=100,
                  mcmc.burnin=10,
                  mcmc.thin=10)

head(props)

# calculate observed heterozygosity and depth of coverage from the allele count
# data
hx <- apply(is.na(t(allele1)+t(allele2))==FALSE,1,mean)
dx <- apply(t(allele1)+t(allele2),1,mean,na.rm=TRUE)

##### estploidy #####
pl <- estploidy2(alphas=props,
                het=hx,
                depth=dx,
                train=FALSE,
                pl=NA,
                set=NA,
                nclasses=2,
                ids=sampleData$Population,
                pcs=1:2)

# Inspect outputs
head(pl$pp)
pl$pcsdev[1]/sum(pl$pcsdev)

# Make a list of the putative polyploids
BigPloidy <- pl$pp[which(as.numeric(pl$pp[,"2"]) > 0.5),]
unique(BigPloidy[,"ids"])

##### Save the R data #####
save(props, file = "1_gbs2ploidy/AlleleProps.RDA")
save(pl, file = "1_gbs2ploidy/PloidyProbs.RDA")
write.csv(noquote(pl$pp), file = "1_gbs2ploidy/GroupProbabilities.csv")

##### Load saved R data #####
load("gbs2ploidy/AlleleProps.RDA")
load("gbs2ploidy/PloidyProbs.RDA")

# options(scipen = 999) # options(scipen = 0) for default

# Add group to the data frame to make manipulation easier
# Check your data to see what are the likely groups. Fraxinus latifolia is 
# generally diploid, so deviations from the norm here are likely polyploids.
# There was also a geographic signal, with the outliers being from a region
# where hybridization, and possibly polyploidy, was thought to occur.

Group <- NULL
for (r in 1:nrow(pl$pp)) {
  if (as.numeric(pl$pp[r,"1"]) > 0.5) {
    Group[r] <- 21
  } else {
    Group[r] <- 24
  }
}

Group # print groups to check 

pca.df <- data.frame(pl$pcscrs)
pca.df$Pop <- sampleData$Population
tail(pca.df)

pca.df$Group <- Group
tail(pca.df)

# Add color to the data frame to make manipulation easier
Color <- NULL
pop.tab$Population # From TableS1
for (p in 1:nrow(pop.tab)) {
  for (r in 1:nrow(pca.df)) {
    if (pca.df$Pop[r] == pop.tab$Population[p]) {
      Color[r] <- pop.tab$`Color code`[p]
    }
  }
}
Color

pca.df$Color <- Color

pca.df[which(pca.df$Group == 24),]

##### Generate plot #####
main.p <- ggplot(pca.df, aes(x=PC1,
                     y=PC2)) + 
  geom_point(data = pca.df,
             aes(fill=Pop,
                 shape=as.factor(Group)),
             size = 3) +
  scale_shape_manual(name = "Estimated ploidy level",
                     labels = c("2x", "4x"),
                     values = c(21,24),
                     guide = NULL) +
  scale_fill_manual(name = "Population",
                    values = my.cols.test,
                    guide = NULL) +
  geom_vline(xintercept=0,linetype="dashed",col="grey") +
  geom_hline(yintercept=0,linetype="dashed",col="grey") +
  xlab("PC1 (37.8%)") + #paste0("PC1 (", round(pl$pcsdev[1]/sum(pl$pcsdev), 3)*100, "%)")
  ylab("PC2 (26.1%)") + #paste0("PC2 (", round(pl$pcsdev[2]/sum(pl$pcsdev), 3)*100, "%)")
  theme_classic()

legend.p <- ggplot(pca.df, aes(x=PC1,
                               y=PC2)) + 
  geom_point(data = pca.df,
             aes(color=Pop,
                 shape=as.factor(Group)),
             size = 3) +
  scale_shape_manual(name = "Estimated ploidy level",
                     labels = c("2x", "4x"),
                     values = c(21,24)) +
  scale_color_manual(name = "Population",
                    values = my.cols.test)

# Extract the legend. Returns a gtable
leg <- get_legend(legend.p)

# Convert to a ggplot and print
legend.plot <- as_ggplot(leg)

pdf("1_gbs2ploidy/gbs2ploidy.pdf") # Fig. S1
ggarrange(main.p, legend.plot, ncol = 2, labels = NULL)
dev.off()

####


#### Spatial genetic structure ####
Pop_Sum <- read.csv("1_FRLA_data/all/data/allPop_Sum.csv")

sort(Pop_Sum$Pop)
sort(unique(SampleData$Pop)) 

#### Create map with ALL populations ####
# No shapefile for the test data as this species has an unknown range.
fralat <- shapefile("~/Downloads/Oregon ash (Fraxinus latifolia) extent, North America/data/commondata/data0/fraxlati.shp")

df.sp <- spTransform(fralat, CRS("+proj=longlat +init=epsg:3857"))

# Set up a google maps account and register your KEY
register_google(key = "KEY", write = TRUE)

# map <- get_googlemap(center = c(-123,42.25),
#                      zoom = 5,
#                      size = c(600,600),
#                      maptype = "satellite")

map <- get_googlemap(center = c(-85,35),
                     zoom = 5,
                     size = c(600,600),
                     maptype = "satellite")

mapplot <- ggmap(map) +
  # geom_polygon(data = df.sp,
  #              aes(x = long, y = lat, group = group),
  #              fill = "grey", colour = "darkgrey", alpha = 0.33) +
  geom_point(data = Pop_Sum, aes(x = Longitude,
     y = Latitude),
     size = 2,
     col = "black",
     fill = "white",
     pch = 21) +
  geom_label_repel(data = Pop_Sum,
     aes(x = Longitude,
      y = Latitude,
      label = Population,
      fill = Population),
     box.padding = 0,
     size = 1.5,
     max.overlaps = 20,
     colour = "black",
     fontface = "bold") +
  scale_fill_manual(name = "Population:", values = my.cols.test) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  
  # scale_x_continuous(limits = c(-125, -117.5), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(33.5, 49), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-90, -75), expand = c(0, 0)) +
  scale_y_continuous(limits = c(25, 45), expand = c(0, 0)) +
  theme_bw() + 
  #scale_size_area(c(3,20)) +
  theme(axis.text = element_text(size = 12),
    axis.title = element_text(size = 15, colour = "black", face = "bold"),
    panel.border = element_rect(linewidth = 1, colour = "black"),
    legend.text = element_text(size = 10), legend.position = "none",
    legend.title = element_text(size = 12, face = "bold"), panel.grid = element_blank())

mapplot

mapplot.g <- ggarrange(mapplot, ncol = 1, nrow = 1, labels = "A")

mapplot.g

ggsave("0_Data/Map_pop.pdf", height = 16, width = 12)

#### Perform PCA for "All" ####

# snpgdsVCF2GDS(vcf.fn = "1_FRLA_data/all/FRLA_all_NoDGR.recode.vcf.gz", out.fn = "1_FRLA_data/all/all_filt.gds")
# genofile <- snpgdsOpen("1_FRLA_data/all/all_filt.gds")

snpgdsVCF2GDS(vcf.fn = "0_Data/TestVCF_5pop_25samp_550snp.vcf.gz", out.fn = "0_Data/TestVCF.gds")
genofile <- snpgdsOpen("0_Data/TestVCF.gds")

a_pca_out <- snpgdsPCA(gdsobj = genofile, autosome.only = F, remove.monosnp = T)
a_pca_df <- data.frame(a_pca_out$eigenvect[,1:5])
a_pca_df$pop <- sampleData$Population
colnames(a_pca_df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "Population")
a_pve <- round(a_pca_out$varprop[1:5], 3)

write.csv(x = a_pca_df, file = "2_SNP_PCA/PCA_df.csv", row.names = F)
write.csv(a_pve, file = "2_SNP_PCA/a_pve.csv", row.names = F)

snpgdsClose(genofile)

##### PCA plots for "All" #####
a_pca_df <- read.csv("2_SNP_PCA/PCA_df.csv")
a_pve <- unlist(read.csv("2_SNP_PCA/a_pve.csv"))

State <- NULL
for (p in 1:nrow(Pop_Sum)) {
  for (r in 1:nrow(a_pca_df)) {
    if (a_pca_df$Population[r] == Pop_Sum$Population[p]) {
      State[r] <- Pop_Sum$State[p]
    }
  }
}
State

a_pca_df$State <- State

Apca12_pop <- ggplot(data = a_pca_df) +
  geom_point(aes(x=PC1, y=PC2, fill=Population, shape=State), size = 3) + 
  scale_shape_manual(values = c(21,22,24,25), name = 'State') +
  scale_fill_manual(values = my.cols.test, name = 'Population') +
  xlab(paste("PC",1," (",a_pve[1]*100,"%)",sep="")) + 
  xlim(c(round(min(a_pca_df$PC1),2)-0.01, round(max(a_pca_df$PC1),2)+0.01)) + 
  ylab(paste("PC",2," (",a_pve[2]*100,"%)",sep="")) +
  theme_bw() + 
  theme(legend.position = 'none',
        axis.text = element_text(size=10),
        axis.title = element_text(size = 16, colour="black",
                                  face = "bold",vjust = 1),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Apca23_pop <- ggplot(data = a_pca_df) +
  geom_point(aes(x=PC3, y=PC2, fill=Population, shape=State), size = 3) + 
  scale_shape_manual(values = c(21,22,24,25), name = 'State') +
  scale_fill_manual(values = my.cols.test, name = 'Population') +
  xlab(paste("PC",3," (",a_pve[3]*100,"%)",sep="")) + 
  xlim(c(round(min(a_pca_df$PC3),2)-0.01, round(max(a_pca_df$PC3),2)+0.01)) + 
  ylab(paste("PC",2," (",a_pve[2]*100,"%)",sep="")) +
  #scale_fill_manual(name = 'Population:', values = my.cols) +
  theme_bw() + 
  theme(legend.position = 'none',
        axis.text = element_text(size=10),
        axis.title = element_text(size = 16, colour="black",
                                  face = "bold",vjust = 1),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Apca34_pop <- ggplot(data = a_pca_df) +
  geom_point(aes(x=PC4, y=PC3, fill=Population, shape=State), size = 3) + 
  scale_shape_manual(values = c(21,22,24,25), name = 'State') +
  scale_fill_manual(values = my.cols.test, name = 'Population') +
  xlab(paste("PC",4," (",a_pve[4]*100,"%)",sep="")) + 
  xlim(c(round(min(a_pca_df$PC4),2)-0.01, round(max(a_pca_df$PC4),2)+0.01)) + 
  ylab(paste("PC",3," (",a_pve[3]*100,"%)",sep="")) +
  #scale_fill_manual(name = 'Population:', values = my.cols) +
  theme_bw() + 
  theme(legend.position = 'none',
        axis.text = element_text(size=10),
        axis.title = element_text(size = 16, colour="black",
                                  face = "bold",vjust = 1),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Apca12_pop.legend <- 
  ggplot(data = a_pca_df) +
  geom_point(aes(x=PC1, y=PC2, color=Population, shape=State), size = 3) + 
  scale_shape_manual(values = c(21,22,24,25), name = 'State') +
  scale_color_manual(values = my.cols.test, name = 'Population') +
  xlab(paste("PC",1," (",a_pve[1]*100,"%)",sep="")) + 
  xlim(c(round(min(a_pca_df$PC1),2)-0.01, round(max(a_pca_df$PC1),2)+0.01)) + 
  ylab(paste("PC",2," (",a_pve[2]*100,"%)",sep="")) +
  theme_bw() + 
  theme(#legend.position = 'none',
        axis.text = element_text(size=10),
        axis.title = element_text(size = 16, colour="black",
                                  face = "bold",vjust = 1),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# # Extract the legend. Returns a gtable
leg <- get_legend(Apca12_pop.legend)
# 
# # Convert to a ggplot and print
legend.plot <- as_ggplot(leg)
legend.plot

##### Perform PCA for "North" #####
# No subsets for test data
sampleData <- read.csv("1_FRLA_data/north/data/nPop_ID_Sum.csv", header = T)
nrow(sampleData)

nPop_Sum <- read.csv("1_FRLA_data/north/data/nPop_Sum.csv")
nrow(nPop_Sum)
# sampleData$site <- paste0(sampleData$State, "_", sampleData$Pop) 

snpgdsVCF2GDS(vcf.fn = "1_FRLA_data/north/FRLA_north_NoDGR.recode.vcf.gz", out.fn = "1_FRLA_data/north/north_filt.gds")
genofile <- snpgdsOpen("1_FRLA_data/north/north_filt.gds")

n_pca_out <- snpgdsPCA(gdsobj = genofile, autosome.only = F, remove.monosnp = T)
n_pca_df <- data.frame(n_pca_out$eigenvect[,1:5])
n_pca_df$pop <- sampleData$Pop
colnames(n_pca_df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "Population")
n_pve <- round(n_pca_out$varprop[1:5], 3)

write.csv(x = n_pca_df, file = "1_FRLA_data/north/data/NorthPCA_df_AEM.csv", row.names = F)
write.csv(n_pve, file = "1_FRLA_data/north/data/n_pve.csv", row.names = F)

snpgdsClose(genofile)

##### PCA plot, all samples #####
pca.grid <- ggarrange(Apca12_pop,
                      Apca23_pop,
                      nrow = 2,
                      ncol = 1,
                      labels = c("B", "C")) #, "D", "E"

pdf("2_SNP_PCA//Map_and_PCAs.pdf", width = 12) # Fig. 1
grid.arrange(mapplot.g, pca.grid, legend.plot,
             ncol=3, 
             widths=c(1,2,1))
dev.off()

##### PCA plots for "north" #####

n_pca_df <- read.csv("1_FRLA_data/north/data/NorthPCA_df_AEM.csv")
n_pve <- unlist(read.csv("1_FRLA_data/North/data/n_pve.csv"))

State <- NULL
for (p in 1:nrow(nPop_Sum)) {
  for (r in 1:nrow(n_pca_df)) {
    if (n_pca_df$Population[r] == nPop_Sum$Pop[p]) {
      State[r] <- nPop_Sum$State[p]
    }
  }
}
State

n_pca_df$State <- State

Npca12_pop <- ggplot(data = n_pca_df) +
  geom_point(aes(x=PC2, y=PC1, fill=Population, shape=State), size = 3) + 
  scale_shape_manual(values = c(21,22,24,25), name = 'State') +
  scale_fill_manual(values = my.cols.north, name = 'Population') +
  xlab(paste("PC",2," (",n_pve[2]*100,"%)",sep="")) + 
  xlim(c(round(min(n_pca_df$PC2),2)-0.01, round(max(n_pca_df$PC2),2)+0.01)) + 
  ylab(paste("PC",1," (",n_pve[1]*100,"%)",sep="")) +
  #scale_fill_manual(name = 'Population:',values = my.cols.north) +
  theme_bw() + 
  theme(legend.position = 'none',
        axis.text = element_text(size=10),
        axis.title = element_text(size = 16, colour="black",
                                  face = "bold",vjust = 1),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Npca23_pop <-   ggplot(data = n_pca_df) +
  geom_point(aes(x=PC3, y=PC2, fill=Population, shape=State), size = 3) + 
  scale_shape_manual(values = c(21,22,24,25), name = 'State') +
  scale_fill_manual(values = my.cols.north, name = 'Population') +
  xlab(paste("PC",3," (",n_pve[3]*100,"%)",sep="")) + 
  xlim(c(round(min(n_pca_df$PC2),2)-0.01, round(max(n_pca_df$PC2),2)+0.01)) + 
  ylab(paste("PC",2," (",n_pve[2]*100,"%)",sep="")) +
  #scale_fill_manual(name = 'Population:', values = my.cols.north) +
  theme_bw() + 
  theme(legend.position = 'none',
        axis.text = element_text(size=10),
        axis.title = element_text(size = 16, colour="black",
                                  face = "bold",vjust = 1),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

##### Generate gridded figure #####

pca.grid <- ggarrange(Apca12_pop,
                         #Apca23_pop,
                         Npca12_pop,
                         #Npca23_pop,
                         nrow = 2,
                         ncol = 1,
                         labels = c("B", "C")) #, "D", "E"

pdf("3_Plots/Map_and_PCAsV4.pdf", width = 12) # Fig. 1
grid.arrange(mapplot.g, pca.grid, legend.plot,
             ncol=3, 
             widths=c(1,2,1))
dev.off()


# pca.grid <- ggarrange(Apca12_pop,
#                       Apca23_pop,
#                       Apca34_pop,
#                       nrow = 3,
#                       ncol = 1,
#                       labels = c("A", "B", "C")) #, "D", "E"
# 
# pdf("3_Plots/AllPCAs1-4.pdf", width = 12)
# grid.arrange(mapplot.g, pca.grid, legend.plot,
#              ncol=3,
#              widths=c(1,2,1))
# dev.off()


pca.grid <- ggarrange(Apca23_pop,
                      Npca23_pop,
                      nrow = 2,
                      ncol = 1,
                      labels = c("B", "C")) #, "D", "E"

pdf("3_Plots/SUPPLEMENTALMap_and_PCAsV3.pdf", width = 12) # Fig. S2
grid.arrange(mapplot.g, pca.grid, legend.plot,
             ncol=3, 
             widths=c(1,2,1))
dev.off()
####

#### Spatial structure using 'fastStructure', as implemented in 'structure-threader' ####
# https://www.royfrancis.com/pophelper/articles/index.html

# Plink2 command to generate the bed file required for fastStructure.
# structure_threader commands to run fastStructure analysis and plot. Plots were
# remade using the 'pophelper' R package below.
write.table(Pop_Sum[c("Population", "n")], file = "3_fastStructure/popfile.txt", row.names = F, col.names = F, sep = "\t")
system("~/Desktop/plink2 --vcf 0_Data/TestVCF_5pop_25samp_550snp.vcf.gz --make-bed --allow-extra-chr --out 0_Data/Test")

# structure-threader seems to only run from CLI as of right now; will troubleshoot later.
# system("structure_threader run -Klist 1 2 3 4 5 6 7 8 9 10 -i 0_Data/Test.bed -o 3_fastStructure -t 24 -fs ~/Desktop/fastStructure --pop 3_fastStructure/popfile.txt")
# system("structure_threader plot -i 3_fastStructure/. -f faststructure -K 2 3 4 -o 3_fastStructure/2_4_plots --pop 3_fastStructure/popfile.txt")

setwd(project.folder)

samp.df <- read.csv("0_Data/SampleData_5pop_25samp_550snp.csv")
pop.order.df <- read.csv("0_Data/Pop_Sum.csv")

samp.order.df <- NULL
for (p in 1:nrow(pop.order.df)) {
  for (s in 1:nrow(samp.df)) {
    if (samp.df$Pop[s] == pop.order.df$Pop[p]) {
      tmp.df <- samp.df[s,]
      samp.order.df <- rbind(samp.order.df, tmp.df)
    }
  }
}
samp.order.df

setwd("3_fastStructure/")
q2 <- read.csv("fS_run_K.2.meanQ", header = F, sep = "")
q3 <- read.csv("fS_run_K.3.meanQ", header = F, sep = "")
q4 <- read.csv("fS_run_K.4.meanQ", header = F, sep = "")

q2 <- cbind(q2, samp.df$sample)
q3 <- cbind(q3, samp.df$sample)
q4 <- cbind(q4, samp.df$sample)

q2.ordered <- q2[match(samp.order.df$sample, q2$`samp.df$sample`),]
q3.ordered <- q3[match(samp.order.df$sample, q3$`samp.df$sample`),]
q4.ordered <- q4[match(samp.order.df$sample, q4$`samp.df$sample`),]

# Save the re-ordered Q matrices
# write.table(x = q2.ordered[,1:2], file = "k2.meanQ.ORDERED", row.names = F, col.names = F)
# write.table(x = q3.ordered[,1:3], file = "k3.meanQ.ORDERED", row.names = F, col.names = F)
# write.table(x = q4.ordered[,1:4], file = "k4.meanQ.ORDERED", row.names = F, col.names = F)

sfiles <- list.files(pattern = "*meanQ", full.names=T)
# basic usage
slist <- readQ(files=sfiles)

tr1 <- tabulateQ(qlist=slist[3:5]) # Which runs do you want to use?
sr1 <- summariseQ(tr1)
summariseQ(tr1, writetable=TRUE, exportpath=getwd())

GroupLabels <- data.frame("Population" = samp.order.df[,1], "State" = samp.order.df[,2])
q3.ordered.test <- cbind(q3.ordered[,1:3], GroupLabels) # check the table

# Create a data frame for the group labels, or call on another data frame with
# the relevent data
p1 <- plotQ(slist[3:5],
            imgoutput = "join",
            returnplot = T,
            exportplot = F,
            basesize = 11,
            clustercol = c("#ebc152", "#f44739", "#1653d8", "#c43200"),
            splab = c("K=2", "K=3", "K=4"),
            grplabangle = 90,
            grplab = GroupLabels,
            grplabsize = 2.5,
            grplabheight = 0.1,
            grplabpos = 0.5)

# Save the plot output
pdf("fastStructurePlot.pdf", width = 14, height = 7)
grid.arrange(p1$plot[[1]])
dev.off()

##### Get likelihood scores; generate likelihood curve plot #####
log.file.list <- list.files(pattern = ".log")
log.df <- NULL # Generate empty object to fill with the marginal likelihoods

for (f in 1:length(log.file.list)) {
  tmp <- readLines(con = log.file.list[f])
  tmp.marg.like <- grep(pattern = "Marginal Likelihood = ", x = tmp, value = T)
  tmp.marg.like <- gsub(pattern = "Marginal Likelihood = ", replacement = "", x = tmp.marg.like)
  log.df <- data.frame(rbind(log.df, as.numeric(tmp.marg.like)))
} 

log.df <- cbind(log.df, as.numeric(c("2", "3", "4", "5", "6")))
colnames(log.df) <- c("Marginal Likelihood", "K")
 
pdf("fastStructure_marglike_plot.pdf")
plot(log.df$`Marginal Likelihood` ~ log.df$K,
     col = "#f44739",
     pch = 19,
     ylab = "Marginal Likelihood",
     xlab = expression(italic("k")))
dev.off()

#### Population genetic diversity analyses ####
##### Pi, W, Tajimas D #####
# angsd and angsd-wrapper were used to estimate genetic diversity metrics. This section
# of code does not computer nucleotide diversity; just makes plots.
# Test data have no ANSGD
setwd(project.folder)
Pop_div <- read.csv("1_FRLA_data/all/data/aemPop_angsd.csv")

Pop_div <- Pop_div[order(Pop_div$Pop, decreasing = F),]

Pop_div$col <- my.cols.north

Pop_div_lat <- Pop_div[order(Pop_div$Latitude, decreasing = T),]
Pop_div_lat$Pop <- factor(Pop_div_lat$Pop, levels = unique(Pop_div_lat$Pop))

cor.test(Pop_div$Pi, Pop_div$Latitude)
cor.test(Pop_div$W, Pop_div$Latitude)
cor.test(Pop_div$TajD, Pop_div$Latitude)

##### Pi #####
Pi_plot <- ggplot(data = Pop_div_lat, aes(x = Pop, y = Pi, fill = Pop)) +
  geom_errorbar(aes(ymin = Pi - Pi_sd, ymax = Pi + Pi_sd), width = 1, colour = "black", size = 1) +
  geom_point(size = 4,pch = 21, colour = "black") +
  ylab("\u03C0") + # Diversity (Pi)
  xlab("Population") +
  scale_fill_manual(values = Pop_div_lat$col) +
  theme_bw() +
  theme(legend.position = "None",
        plot.title = element_text(size = 26, colour = "black", face = "bold"),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.title = element_text(size = 15, colour = "black", face = "bold"),
        panel.border = element_rect(linewidth = 1.5,colour = "black"),
        legend.title = element_text(size = 20, colour = "black", face = "bold", vjust = 1),
        legend.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background = element_rect(size = 1.5, colour = "#333333", fill = "#CCCCCC"))
        
Pi_plot
##### W #####                                                                                                                                                                                                                                                                                               
W_plot <- ggplot(data = Pop_div_lat, aes(x = Pop, y = W, fill = Pop)) +
  geom_errorbar(aes(ymin = W - W_sd, ymax = W + W_sd), width = 1, colour = "black", size = 1) +
  geom_point(size = 4, pch = 21,colour = "black") +
  ylab("\u0398 W") + # Diversity (Watterson's)
  xlab("Population") +
  scale_fill_manual(values = Pop_div_lat$col) +
  theme_bw() +
  theme(legend.position = "None",
        plot.title = element_text(size = 26, colour = "black", face = "bold"),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.title = element_text(size = 15, colour = "black", face = "bold"),
        panel.border = element_rect(linewidth = 1.5,colour = "black"),
        legend.title = element_text(size = 20, colour = "black", face = "bold", vjust = 1),
        legend.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background = element_rect(size = 1.5, colour = "#333333", fill = "#CCCCCC"))                                                                                                                                           
W_plot

##### TajD #####
TajD_plot <- TajD_plot <- ggplot(data = Pop_div_lat, aes(x = Pop, y = TajD, fill = Pop)) +
  geom_errorbar(aes(ymin = TajD - TajD_sd, ymax = TajD + TajD_sd), width = 1, colour = "black", size = 1) +
  geom_point(size = 4, pch = 21, colour = "black") +
  ylab("D") + # Tajima's D
  xlab("Population") +
  scale_fill_manual(values = Pop_div_lat$col) +
  theme_bw() +
  theme(legend.position = "None",
    plot.title = element_text(size = 26, colour = "black", face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 90),
    axis.title = element_text(size = 15, colour = "black", face = "bold"),
    panel.border = element_rect(linewidth = 1.5, colour = "black"),
    legend.title = element_text(size = 20, colour = "black", face = "bold", vjust = 1),
    legend.text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 22, face = "bold"),
    strip.background = element_rect(size = 1.5, colour = "#333333", fill = "#CCCCCC"))
TajD_plot

# pdf("3_Plots/DivPlots_AEM.pdf")
ggarrange(Pi_plot,W_plot,TajD_plot, ncol = 1, labels = c("A", "B", "C"))

ggsave(filename = "3_Plots/DivPlots_AEMv2.pdf", width = 7.5, height = 10, dpi = 600)
ggsave(filename = "3_Plots/DivPlots_AEMv2.png", width = 7.5, height = 10, dpi = 600)


#### Nei's D and Fst ####
# the 012 file is an output of 'vcftools'
setwd(project.folder)
system("vcftools --gzvcf 0_Data/TestVCF_5pop_25samp_550snp.vcf.gz --maf 0.02 --out 0_Data/TestVCF_maf02 --012")

df012 <- fread("0_Data/TestVCF_maf02.012", sep = "\t", data.table = F) 
dim(df012)

df012[,1]
df012 <- df012[,-1]
dim(df012)

Pop_Sum <- read.csv("0_Data/Pop_Sum.csv")
nrow(Pop_Sum)
pca_df <- read.csv("2_SNP_PCA/PCA_df.csv")
nrow(pca_df)

nrow(Pop_Sum)
nrow(pca_df)

##### remove missing data #####
df012 <- apply(df012, 2, function(d) gsub(-1, NA, d, fixed = TRUE))
df012 <- apply(df012, 2, function(d) as.numeric(d))

Pop_Sum$St_Pop <- paste0(Pop_Sum$State, "_", Pop_Sum$Population)

##### Calc allele freq by population ######
g <- t(df012)
uniq <- unique(as.character(Pop_Sum$St_Pop))

allele_df <- matrix(ncol = length(uniq), nrow = nrow(g))
for (i in 1:length(uniq)) {
  index <- which(Pop_Sum$St_Pop == uniq[i])
  df <- g[, index]
  if (length(index) == 1) {
    allele_df[, i] <- df/2
  } else {
    af <- apply(df, 1, function(d) sum(d, na.rm = TRUE)/(2 *
                                                           length(which(!is.na(d)))))
    allele_df[, i] <- af
  }
}
afreq <- t(allele_df)

hist(afreq[, 1:100])

write.table(round(afreq, digits = 5), file = "4_GeneticDiversity/afreq.txt",
            quote = F, row.names = F, col.names = F, sep = ",")


##### Nei's D #####
neisD <- function(P) {
  ## rows = pops, cols = loci
  N <- dim(P)[1]
  L <- dim(P)[2]
  D <- matrix(0, nrow = N, ncol = N)
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      Pi <- as.numeric(P[i, ])
      Pj <- as.numeric(P[j, ])
      D[i, j] <- 1 - (sum(sqrt((1 - Pi) * (1 - Pj)), na.rm = TRUE) +
                        sum(sqrt(Pi * Pj), na.rm = TRUE))/L
      D[j, i] <- D[i, j]
    }
  }
  D <- as.dist(D)
  return(D)
}

neiDcommon <- neisD(afreq)

neiD <- as.matrix(neiDcommon)
rownames(neiD) <- unique(Pop_Sum$St_Pop)
write.table(neiD, file = "4_GeneticDiversity/maf02_neiD.txt", col.names = F, quote = F)

##### Fst and D / Fst triangle #####
##### read in afreq #####
afreq <- fread("4_GeneticDiversity/afreq.txt", sep = ",", data.table = F)
afreq <- as.matrix(afreq)

##### pairwise hudsons Fst #####

hudsonFst2 <- function(p1 = NA, p2 = NA, n1 = NA, n2 = NA) {
  numerator <- p1 * (1 - p1) + p2 * (1 - p2)
  denominator <- p1 * (1 - p2) + p2 * (1 - p1)
  fst <- 1 - numerator/denominator
  out <- cbind(numerator, denominator, fst)
  return(out)
}

Pop_Sum <- read.csv("0_Data/Pop_Sum.csv")

Pop_Sum$St_Pop <- paste0(Pop_Sum$State, "_", Pop_Sum$Population)

pops <- unique(Pop_Sum$St_Pop)
pops <- as.matrix(pops)

pop.fst.mean.med <- matrix(data = NA, nrow = length(unique(Pop_Sum$St_Pop)),
                           ncol = length(unique(Pop_Sum$St_Pop)))

for (i in 1:nrow(afreq) - 1) {
  for (j in i:nrow(afreq)) {
    fst <- hudsonFst2(p1 = as.numeric(afreq[i, ]), p2 = as.numeric(afreq[j,
    ]))
    pop.fst.mean.med[i, j] <- mean(fst[, 3], na.rm = T)
    pop.fst.mean.med[j, i] <- quantile(fst[, 3], probs = c(0.5),
                                       na.rm = T)
  }
}


colnames(pop.fst.mean.med) <- pops
write.csv(cbind(pops, pop.fst.mean.med), file = "4_GeneticDiversity/pop_fst_mean_med.csv",
          col.names = F, quote = F, row.names = F)


summary(pop.fst.mean.med[upper.tri(pop.fst.mean.med)])
summary(pop.fst.mean.med[lower.tri(pop.fst.mean.med)])

##### Create a data frame with half Fst and half Neis D #####
###### Sorted by Lat ######
pop_lat <- as.character(unique(Pop_Sum$St_Pop[order(Pop_Sum$Latitude,
                                                       decreasing = T)]))

# Fst
pwFst <- read.csv("4_GeneticDiversity/pop_fst_mean_med.csv")
rownames(pwFst) <- pwFst$X
pwFst <- pwFst[, -1]
colnames(pwFst) <- rownames(pwFst)
pwFst[nrow(pwFst), ncol(pwFst)] <- 0  #last element always NA, idk why

# keep only upper tri
for (i in 1:nrow(pwFst)) {
  for (j in 1:nrow(pwFst)) {
    pwFst[j, i] <- pwFst[i, j]
  }
}

pwFst <- pwFst[sapply(pop_lat, function(s) which(rownames(pwFst) ==
                                                   s)), sapply(pop_lat, function(s) which(colnames(pwFst) ==
                                                                                            s))]


# neiD
neiD <- as.data.frame(read_table2("4_GeneticDiversity/maf02_neiD.txt", col_names = FALSE))
rownames(neiD) <- neiD$X1
neiD <- neiD[, -1]
colnames(neiD) <- rownames(neiD)
neiD <- neiD[sapply(pop_lat, function(s) which(rownames(neiD) ==
                                                 s)), sapply(pop_lat, function(s) which(colnames(neiD) ==
                                                                                          s))]

pwFst_neiD <- pwFst
pwFst_neiD[lower.tri(pwFst_neiD)] <- neiD[lower.tri(neiD)]

write.csv(pwFst_neiD, "4_GeneticDiversity/pwFst_neiD.csv", row.names = F)

###### make df for plotting ######
pwFst_neiD_df <- data.frame(var1 = NA, var2 = NA, value = NA,
                            stat = NA)

pw_mat <- lower.tri(neiD)
pw_mat[lower.tri(neiD)] <- "neiD"
pw_mat[upper.tri(neiD)] <- "Fst"
pw_mat[which(pw_mat == "FALSE")] <- "diag"

for (i in 1:ncol(pwFst_neiD)) {
  for (j in 1:nrow(pwFst_neiD)) {
    var1 <- rownames(pwFst_neiD)[j]
    var2 <- colnames(pwFst_neiD)[i]
    value <- pwFst_neiD[j, i]
    stat <- pw_mat[j, i]
    new_df <- data.frame(var1, var2, value, stat)
    pwFst_neiD_df <- rbind(pwFst_neiD_df, new_df)
  }
}
pwFst_neiD_df <- pwFst_neiD_df[-1, ]

pwFst_neiD_df$value_scale <- NA
pwFst_neiD_df$value_scale[which(pwFst_neiD_df$stat == "neiD")] <- scale(pwFst_neiD_df$value[which(pwFst_neiD_df$stat ==
                                                                                                    "neiD")])
pwFst_neiD_df$value_scale[which(pwFst_neiD_df$stat == "Fst")] <- scale(pwFst_neiD_df$value[which(pwFst_neiD_df$stat ==
                                                                                                   "Fst")])
write.csv(pwFst_neiD_df, "4_GeneticDiversity/pwFst_neiD_df.csv", row.names = F)

##### Nei D and Fst plot #####
col <- wes_palette("Zissou1", 1000, type = "continuous")

pw_plot <- ggplot(data=pwFst_neiD_df,aes(#x=var1,
  x=reorder(var1, desc(var1)),
  #y=reorder(var2, desc(var2)),
  y=var2,
  fill=value_scale)) +
  #ggplot(data=pwFst_neiD_df,aes(x=var1_label,y=var2_label)) +
  #geom_tile(aes(fill=value_scale),color='white') + 
  geom_tile(color='white') + 
  #geom_text(aes(x=var2_num,y=var1_num,
  #              label = as.character(round(value,3))),fontface='bold',size=5) +
  geom_text(aes(label = round(value,3)),fontface='bold',size=5) +
  scale_fill_gradientn(colours=col) +
  coord_flip() + 
  scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0),position = 'right',guide = guide_axis(angle = 90)) +
  theme_bw() + theme(
    legend.position = 'none',
    #axis.text.y = element_blank(),
    axis.text.y = element_text(size=16,face='bold',colour='black'),
    axis.text.x = element_text(size=16,face='bold',colour='black'),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank(),
    #plot.margin = margin(10, 30, 10, 10),
    panel.border = element_rect(size = 1.5, colour = "black"))
pw_plot

ggsave('4_GeneticDiversity/Fst_NeiD_heatmap.pdf',height=26,width=30)
####


#### Estimating kinship and effective population size ####
##### SNPRelate relatedness metrics #####
# Vignette / tutorial: https://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html

# Set working directory
setwd(project.folder)

#### File conversion ####

# snpgdsVCF2GDS(vcf.fn = "1_FRLA_data/north/FRLA_north_NoDGR.recode.vcf.gz", 
#               out.fn = "north.ccm.gds") # vcf.fn = vcf_file
genofile <- openfn.gds("0_Data/TestVCF.gds")

# Estimate IBD coefficients MoM method
Pops <- unique(sampleData$Population)
for (p in 1:length(Pops)) {
  
  sample_ids <- sampleData$sample[sampleData$Population == Pops[p]]
  
  ibd.MoM <- snpgdsIBDMoM(genofile,
                            sample.id=sample_ids,
                            num.thread=10, 
                            autosome.only = F)
  ibd.MoM.coeff <- snpgdsIBDSelection(ibd.MoM)
  write.csv(x = ibd.MoM.coeff, file = paste0("5_Kinship/MoM_", Pops[p], ".csv"), row.names = F)
  
  ibd.MLE <- snpgdsIBDMLE(genofile,
                      sample.id=sample_ids,
                      num.thread=10,
                      autosome.only = F)
  ibd.MLE.coeff <- snpgdsIBDSelection(ibd.MLE)
  write.csv(x = ibd.MLE.coeff, file = paste0("5_Kinship/MLE_", Pops[p], ".csv"), row.names = F)
  
}

closefn.gds(genofile)

##### Get averages #####
setwd("5_Kinship/")

ibd.MoM.files <- list.files(pattern = "^MoM_*{3}")
ibd.MoM.all <- NULL
for (f in 1:length(ibd.MoM.files)) {
  tmp.csv <- read.csv(ibd.MoM.files[f])
  tmp.pop <- gsub(pattern = "MoM_", replacement = "", x = ibd.MoM.files[f])
  tmp.pop <- gsub(pattern = ".csv", replacement = "", x = tmp.pop)
  
  k0.mean <- mean(tmp.csv$k0)
  k1.mean <- mean(tmp.csv$k1)
  kinship.mean <- mean(tmp.csv$kinship)
  
  tmp.df <- data.frame(tmp.pop, k0.mean, k1.mean, kinship.mean)
  
  ibd.MoM.all <- rbind(ibd.MoM.all, tmp.df)
  
}
write.csv(x = ibd.MoM.all, file = "ibd_MoM_Means.csv", row.names = F)

ibd.MLE.files <- list.files(pattern = "^MLE_*{3}")
ibd.MLE.all <- NULL
for (f in 1:length(ibd.MLE.files)) {
  tmp.csv <- read.csv(ibd.MLE.files[f])
  tmp.pop <- gsub(pattern = "MLE_", replacement = "", x = ibd.MLE.files[f])
  tmp.pop <- gsub(pattern = ".csv", replacement = "", x = tmp.pop)
  
  k0.mean <- mean(tmp.csv$k0)
  k1.mean <- mean(tmp.csv$k1)
  kinship.mean <- mean(tmp.csv$kinship)
  
  tmp.df <- data.frame(tmp.pop, k0.mean, k1.mean, kinship.mean)
  
  ibd.MLE.all <- rbind(ibd.MLE.all, tmp.df)
  
}
write.csv(x = ibd.MLE.all, file = "ibd_MLE_Means.csv", row.names = F)

##### All samples, not just means #####
ibd.MoM.files <- list.files(pattern = "^MoM_*{3}")
ibd.MoM.all.samps <- NULL
for (f in 1:length(ibd.MoM.files)) {
  tmp.csv <- read.csv(ibd.MoM.files[f])
  tmp.pop <- gsub(pattern = "MoM_", replacement = "", x = ibd.MoM.files[f])
  tmp.pop <- gsub(pattern = ".csv", replacement = "", x = tmp.pop)
  
  tmp.csv$Population <- tmp.pop
  ibd.MoM.all.samps <- rbind(tmp.csv, ibd.MoM.all.samps)
  
}
write.csv(x = ibd.MoM.all.samps, file = "ibd_MoM_AllSamps.csv", row.names = F)

ibd.MLE.files <- list.files(pattern = "^MLE_*{3}")
ibd.MLE.all.samps <- NULL
for (f in 1:length(ibd.MLE.files)) {
  tmp.csv <- read.csv(ibd.MLE.files[f])
  tmp.pop <- gsub(pattern = "MLE_", replacement = "", x = ibd.MLE.files[f])
  tmp.pop <- gsub(pattern = ".csv", replacement = "", x = tmp.pop)
  
  tmp.csv$Population <- tmp.pop
  ibd.MLE.all.samps <- rbind(tmp.csv, ibd.MLE.all.samps)
  
}
write.csv(x = ibd.MLE.all.samps, file = "ibd_MLE_AllSamps.csv", row.names = F)
#### Generate plots ####
pdf(file = "KinshipEstimate.pdf")
par(mfrow=c(2,1))
plot(x = as.factor(ibd.MoM.all$tmp.pop),
     y = ibd.MoM.all$kinship.mean,
     xlab = "Population",
     ylab = "Mean kinship",
     main = "Method of moments",
     cex.axis=0.5,
     las = 2)
arrows()

plot(x = as.factor(ibd.MLE.all$tmp.pop),
     y = ibd.MLE.all$kinship.mean,
     xlab = "Population",
     ylab = "Mean kinship",
     main = "Maximum likelihood estimation",
     cex.axis=0.5,
     las = 2)
arrows()

dev.off()

##### Read in data, if needed #####
setwd(paste0(project.folder, "/5_Kinship/"))
Pop_Sum <- read.csv("../0_Data/Pop_Sum.csv")
ibd.MoM.all.samps <- read.csv("ibd_MoM_AllSamps.csv")
ibd.MLE.all.samps <- read.csv("ibd_MLE_AllSamps.csv")

Pop.df <- Pop_Sum[order(Pop_Sum$Population, decreasing = F),]
Pop.df$Col <- my.cols.test # Check your color palette!

##### MoM sort #####
#Turn your 'treatment' column into a character vector
ibd.MoM.all.samps$Population <- as.character(ibd.MoM.all.samps$Population)

ibd.MoM.all.samps$Latitude <- rep(0, nrow(ibd.MoM.all.samps))
for (r in 1:nrow(ibd.MoM.all.samps)) {
  for (q in 1:nrow(Pop.df)) {
    if (ibd.MoM.all.samps$Population[r] == Pop.df$Pop[q]) {
      ibd.MoM.all.samps$Latitude[r] <- Pop.df$Latitude[q]
      ibd.MoM.all.samps$Col[r] <- Pop.df$Col[q]
    }
  } 
}

ibd.MoM.all.samps <- ibd.MoM.all.samps[order(ibd.MoM.all.samps$Latitude, decreasing = T),]
#Then turn it back into a factor with the levels in the correct order
ibd.MoM.all.samps$Population <- factor(ibd.MoM.all.samps$Population, levels=unique(ibd.MoM.all.samps$Population))

##### MLE sort #####
#Turn your 'treatment' column into a character vector
ibd.MLE.all.samps$Population <- as.character(ibd.MLE.all.samps$Population)

ibd.MLE.all.samps$Latitude <- rep(0, nrow(ibd.MLE.all.samps))
for (r in 1:nrow(ibd.MLE.all.samps)) {
  for (q in 1:nrow(Pop.df)) {
    if (ibd.MLE.all.samps$Population[r] == Pop.df$Pop[q]) {
      ibd.MLE.all.samps$Latitude[r] <- Pop.df$Latitude[q]
      ibd.MLE.all.samps$Col[r] <- Pop.df$Col[q]
    }
  } 
}

ibd.MLE.all.samps <- ibd.MLE.all.samps[order(ibd.MLE.all.samps$Latitude, decreasing = T),]
#Then turn it back into a factor with the levels in the correct order
ibd.MLE.all.samps$Population <- factor(ibd.MLE.all.samps$Population, levels=unique(ibd.MLE.all.samps$Population))


##### Plot all samples by population #####
# Scatter plot by group

MoM.plot.LEGENDARY <- ggplot(data = ibd.MoM.all.samps, aes(x = Population, y = kinship, fill = Population)) +
  geom_point(pch=21,color='black',cex=2) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10)) +
  theme(legend.direction="horizontal", legend.position = "bottom") +
  guides(fill=guide_legend(ncol=14)) +
  scale_fill_manual(name = 'Population:', values = unique(ibd.MoM.all.samps$Col)) 

MoM.plot <- ggplot(data = ibd.MoM.all.samps, aes(x = Population, y = kinship, fill = Population)) +
  geom_hline(yintercept = 0.25, color="grey", linetype="dashed") +
  # geom_violin(show.legend = FALSE) +
  geom_point(pch=21,color='black',cex=2, show.legend = FALSE) + 
  ylab("kinship (MoM)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10)) +
  scale_fill_manual(name = 'Population:', values = unique(ibd.MoM.all.samps$Col))

MLE.plot <- ggplot(data = ibd.MLE.all.samps, aes(x = Population, y = kinship, fill = Population)) +
  geom_hline(yintercept = 0.25, color="grey", linetype="dashed") +
  # geom_violin(show.legend = FALSE) +
  geom_point(pch=21,color='black',cex=2, show.legend = FALSE) + 
  theme_classic() +
  ylab("kinship (MLE)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10)) +
  scale_fill_manual(name = 'Population:', values = unique(ibd.MLE.all.samps$Col)) 

# Extract the legend. Returns a gtable
leg <- get_legend(MoM.plot.LEGENDARY)

# Convert to a ggplot and print
legend.plot <- as_ggplot(leg)

###### Save pdf ######
pdf(file = "KinshipEstimates.pdf")
ggarrange(MoM.plot, MLE.plot,
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
# (MoM.plot / MLE.plot / legend.plot)
dev.off()

##### Summary stats per pop #####
setwd(paste0(project.folder, "/5_Kinship/"))
ibd.MoM.all.samps <- read.csv("ibd_MoM_AllSamps.csv")
ibd.MLE.all.samps <- read.csv("ibd_MLE_AllSamps.csv")

big.kin <- ibd.MLE.all.samps[which(round(ibd.MLE.all.samps$kinship, 2) >= 0.25),]
lil.kin <- ibd.MLE.all.samps[which(round(ibd.MLE.all.samps$kinship, 2) < 0.25),]

(nrow(big.kin)/nrow(ibd.MLE.all.samps))*100
nrow(big.kin)
length(unique(big.kin$Population))

nrow(lil.kin)

PopList <- unique(ibd.MoM.all.samps$Population)

MoM.sum <- NULL
for (p in 1:length(PopList)) {
  tmp.subset <- subset(ibd.MoM.all.samps, ibd.MoM.all.samps$Population == PopList[p])
  tmp.sum <- summary(tmp.subset$kinship)
  tmp.sum$Population <- PopList[p]
  MoM.sum <- rbind(MoM.sum, tmp.sum)
}
MoM.sum

MLE.sum <- NULL
for (p in 1:length(PopList)) {
  tmp.subset <- subset(ibd.MLE.all.samps, ibd.MLE.all.samps$Population == PopList[p])
  tmp.sum <- summary(tmp.subset$kinship)
  tmp.sum$Population <- PopList[p]
  MLE.sum <- rbind(MLE.sum, tmp.sum)
}
MLE.sum

####

#### LD-Ne estimation using StrataG ####
# Vignette / tutorial: https://github.com/EricArcher/strataG

##### StrataG LDNE #####
setwd(paste0(project.folder, "/6_EffectivePopulationSize/"))
vcf <- vcfR::read.vcfR(file = "../0_Data/TestVCF_5pop_25samp_550snp.vcf.gz") # load vcf
genind.obj <- vcfR::vcfR2genind(x = vcf) # Convert vcf format to genind format
sampleData <- read.csv("../0_Data/SampleData_5pop_25samp_550snp.csv", header = T)
nrow(sampleData)

genind.obj@pop <- as.factor(sampleData$Population) # Add pop ID to genind object

df <- genind2df(genind.obj, sep = "/", usepop = T, oneColPerAll = F)
for(j in 2:ncol(df)){
  sel <- which(df[,j]=="NA")
  df[sel,j] <- "-9/-9"
  sel <- which(is.na(df[,j])==T)
  df[sel,j] <- "-9/-9"
}

indname <- rownames(df)
popname <- df[,1]
locname <- colnames(df)[2:ncol(df)]
df <- df[,-1]
  
gen <- df2genind(df, sep = "/", ncode = NULL, ind.names = indname,  pop=popname,
                   loc.names = locname, NA.char = "-9", ploidy = 2,
                   type = "codom")

gtypes <- genind2gtypes(x = genind.obj)
eff <- ldNe(gtypes,
            maf.threshold = 0.02,
            by.strata = F,
            ci = 0.95,
            num.cores = NULL,
            drop.missing = F)

eff
write.csv(eff, file = "StrataG.csv", row.names = F)
##### Plot results #####
eff <- read.csv(file = "StrataG.csv")

# Check you color palette!
eff$Col <- my.cols.test

Pop_Sum <- read.csv("../0_Data/Pop_Sum.csv")

Pop.df <- Pop_Sum[order(Pop_Sum$Population, decreasing = F),]

all.equal(Pop.df$Population, eff$stratum)

eff$Latitude <- Pop.df$Latitude

# Order by latitude
eff$stratum <- as.character(eff$stratum)
eff <- eff[order(eff$Latitude, decreasing = T),]
eff$stratum <- factor(eff$stratum, levels=eff$stratum)

# plot
pdf(file = "StrataG.pdf")
ggplot(data = eff, aes(x = stratum, y = Ne, fill = stratum)) +
  geom_point(pch=21,color='black',cex=2, show.legend = F) +
  xlab("Population") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10)) +
  geom_errorbar(data = eff, 
                aes(x=stratum, 
                    ymin=param.lci, 
                    ymax=param.uci),
                width=0.5,
                # colour=my.cols,
                alpha=0.5,
                size=0.75) +
  scale_fill_manual(name = 'Population:', values = eff$Col)
dev.off()

# Send to object to combine with kinship plots
StrataG <- ggplot(data = eff, aes(x = stratum, y = Ne, fill = stratum)) +
  geom_hline(yintercept = 50, color="grey", linetype="dashed") +
  geom_hline(yintercept = 500, color="grey", linetype="dashed") +
  geom_point(pch=21,color='black',cex=2, show.legend = F) + 
  xlab("Population") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10)) +
  geom_errorbar(data = eff, 
                aes(x=stratum, 
                    ymin=param.lci, 
                    ymax=param.uci),
                width=0.5,
                # colour=my.cols,
                alpha=0.5,
                size=0.75) +
  scale_fill_manual(name = 'Population', values = eff$Col)
  
# Pull the "StrataG" object from LDNe script
pdf("Kinship_Ne.pdf") # Fig. 3
ggarrange(MoM.plot, MLE.plot, StrataG,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)
dev.off()

####
#### Estimating Effective Migration Surfaces (EEMS) ####
### ~~~ Write SNP data to plink format for EEMS analysis.
### ~~~ Plot results of EEMS analysis
# Vignette / tutorial: https://github.com/dipetkov/eems

#### Set working directory ####
setwd(paste0(project.folder, "/9_EEMS/"))

#### Functions
# Compute the diffs matrix using the "mean allele frequency"
# imputation method
bed2diffs_v2 <- function(genotypes) {
  
  nIndiv <- nrow(genotypes)
  nSites <- ncol(genotypes)
  missing <- is.na(genotypes)
  
  ## Impute NAs with the column means (= twice the allele frequencies)
  geno_means <- colMeans(genotypes, na.rm = TRUE)
  # nIndiv rows of genotype means
  geno_means <- matrix(geno_means, nrow = nIndiv, ncol = nSites, byrow = TRUE) 
  
  ## Set the means which correspond to observed genotypes to 0
  geno_means[missing == FALSE] <- 0
  ## Set the missing genotypes to 0 (used to be NA) 
  genotypes[missing == TRUE] <- 0
  genotypes <- genotypes + geno_means
  
  similarities <- genotypes %*% t(genotypes) / nSites
  self_similarities <- diag(similarities)
  vector1s <- rep(1, nIndiv)
  
  diffs <- 
    self_similarities %*% t(vector1s) + 
    vector1s %*% t(self_similarities) - 2 * similarities
  diffs
}

### ~~~ 

##### Plink #####
# Download plink2 binary from https://www.cog-genomics.org/plink/2.0/

system("~/Desktop/plink2 --vcf ../0_Data/TestVCF_5pop_25samp_550snp.vcf.gz --make-bed --allow-extra-chr --out Test")

##### bed2diff #####
# Run bed2diff script
fam <- read.table("Test.fam")
fam <- fam[,-1]
fam$V1 <- fam$V2
fam_df <- data.frame(V1=fam$V1, V2=fam$V2, V3=fam$V3, V4=fam$V4, V5=fam$V5, V6=fam$V6)
write.table(fam_df,"Test.fam", sep="\t", quote=F, col.names = F, row.names = F )
plink_dataset <- "Test"
genotypes <- MultiPhen::read.plink(plink_dataset, indiv = NULL)
genotypes

# genotypes
dim(genotypes) #  

##### Convert format and write table #####
diffs2 <- bed2diffs_v2(genotypes)
View(diffs2[1:300,1:300])
write.table(diffs2, "Test.diffs",sep="\t",quote=F, col.names = F, row.names = F )

#### Plotting EEMS results ####
# ~ Run eems on the HPC. Move output files to local directories for plotting.

# plots
# Download EEMS outputs to local folders called eemsR#_output, if run on HPC
##### Loop over the output folders to make plots that go in the eemsR#_plots folders #####
for (i in 1:12) {
 
  mcmcpath <- paste0("R", i, "_output/")
  plotpath <- paste0("R", i, "_plots/")
  
  # eems.plots(mcmcpath, plotpath, longlat = TRUE, out.png = FALSE)
  
  # Better plotting
  
  map_world <- getMap()
  map_NA <- map_world[which(map_world@data$continent=="North America"),]
  
  # eems.plots(mcmcpath, plotpath, longlat = TRUE ,  m.plot.xy={plot(map_NA, col=NA, add=TRUE)},q.plot.xy={plot(map_NA, col=NA, add=TRUE)}, out.png=FALSE,add.grid =FALSE, add.demes = TRUE, add.outline = TRUE, col.outline = "black", lwd.outline = 2, col.demes = "black", pch.demes =18, min.cex.demes = 0.5, max.cex.demes =2.0)
  eems.plots(mcmcpath,
             plotpath,
             longlat = TRUE,
             m.plot.xy={plot(map_NA, col=NA, add=TRUE)},
             q.plot.xy={plot(map_NA, col=NA, add=TRUE)},
             out.png=FALSE,
             add.grid =TRUE,
             add.demes = TRUE,
             add.outline = TRUE,
             col.outline = "black",
             lwd.outline = 2,
             col.demes = "black",
             pch.demes =18,
             min.cex.demes = 0.5,
             max.cex.demes =2.0)
  
}
#### reemsplots2 ####

# library("devtools")
# install_github("dipetkov/reemsplots2")
library(reemsplots2)

occs <- read.csv("../1_FRLA_data/north/data/pca_df_NoDGR.csv")[,c("Longitude","Latitude")]
head(occs)
nrow(occs)

mcmcpath <- paste0("NorthOnly/R5_output/")
plotpath <- paste0("../3_Plots/EEMS/")

# eems.plots(mcmcpath, plotpath, longlat = TRUE, out.png = FALSE)

# Better plotting
###
library(raster)
states    <- c('California','Oregon', 'Washington')
provinces <- c("British Columbia")

us <- getData("GADM",country="USA",level=1)
canada <- getData("GADM",country="CAN",level=1)

us.states <- us[us$NAME_1 %in% states,]
ca.provinces <- canada[canada$NAME_1 %in% provinces,]

us.bbox <- bbox(us.states)
ca.bbox <- bbox(ca.provinces)
xlim <- c(min(us.bbox[1,1],ca.bbox[1,1]),max(us.bbox[1,2],ca.bbox[1,2]))
ylim <- c(min(us.bbox[2,1],ca.bbox[2,1]),max(us.bbox[2,2],ca.bbox[2,2]))

map <- get_googlemap(center = c(-123,42.25),
                     zoom = 5,
                     size = c(600,600), 
                     maptype = "satellite") # %>% ggmap()

mapplot <- ggmap(map) +
  geom_point(data = nPop_Sum, aes(x = Longitude,
                                 y = Latitude),
             size = 2,
             col = "black",
             fill = "white",
             pch = 21) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  scale_x_continuous(limits = c(-125, -120), breaks = c(-125, -122.5, -120), expand = c(0, 0)) +
  scale_y_continuous(limits = c(37.5, 50), expand = c(0, 0)) +
  theme_bw() # + 
  # scale_size_area(c(3,20)) +
  # theme(axis.text = element_text(size = 12),
  #       axis.title = element_text(size = 15, colour = "black", face = "bold"),
  #       panel.border = element_rect(linewidth = 1, colour = "black"),
  #       legend.text = element_text(size = 10), legend.position = "none",
  #       legend.title = element_text(size = 12, face = "bold"), panel.grid = element_blank())

mapplot
ggsave("gMap.pdf", plot = mapplot)

#plot(ca.provinces, xlim=xlim, ylim=ylim)
#par(new=TRUE)
#plot(us.states, xlim=xlim, ylim=ylim)
plot(mapplot)
par(new=TRUE)
eems.plot$mrates01
points(x = occs$Longitude, y = occs$Latitude)

#### the real command, commented out ####
# eems.plots(mcmcpath, plotpath, longlat = TRUE ,  m.plot.xy={plot(map_NA, col=NA, add=TRUE)},q.plot.xy={plot(map_NA, col=NA, add=TRUE)}, out.png=FALSE,add.grid =FALSE, add.demes = TRUE, add.outline = TRUE, col.outline = "black", lwd.outline = 2, col.demes = "black", pch.demes =18, min.cex.demes = 0.5, max.cex.demes =2.0)
#### other tries ####
eems.plot <- make_eems_plots(mcmcpath,
                             longlat = T,
                             #add_grid = T,
                             add_demes = T)

eems.plot$mrates01

###
eems.plot$mrates01 +
  ggplot(us.states,aes(x=long,y=lat,group=group))+
  geom_path()+
  geom_path(data=ca.provinces)+
  coord_map()
####
##### Parameter testing for EEMS analyses #####
#### Set working directory ####
setwd(paste0(project.folder, "/9_EEMS/"))

# Read in the log files
FileList <- list.files(path = "LogFiles/")
head(FileList)

# Loop over the log files, extract relevant info, pick "best" run based on
# logLikelihood score. Adding in nDemes vs observed demes per KRL.

eemsDF <- NULL # Create empty object to add to in loop
for (l in 1:length(FileList)) {
  tmpFile <- readLines(con = paste0("LogFiles/", FileList[l]))
  Run <- gsub(pattern = "aem_eems", replacement = "", x = FileList[l])
  Run <- gsub(pattern = ".pbs.*.", replacement = "", x = Run)
  # head(tmpFile)
  n <- noquote(gsub(pattern = "[[:space:]]*.There are ", replacement = "", x = tmpFile[grep(pattern = "[[:space:]]+There are \\d+ samples", x = tmpFile)]))
  n <- as.numeric(str_extract(string = n, pattern = "\\d+")) # noquote(gsub(pattern = " observed demes (out of \\[1-9]. demes)", replacement = "", x = ObsDemes))
  
  ObsDemes <- noquote(gsub(pattern = "[[:space:]]*.There are ", replacement = "", x = tmpFile[grep(pattern = "[[:space:]]*.There are \\d+ observed demes", x = tmpFile)]))
  ObsDemes <- as.numeric(str_extract(string = ObsDemes, pattern = "\\d+")) # noquote(gsub(pattern = " observed demes (out of \\[1-9]. demes)", replacement = "", x = ObsDemes))
  
  nDemesGrid <- noquote(gsub(pattern = "[[:space:]]*.The population grid has ", replacement = "", x = tmpFile[grep(pattern = "[[:space:]]*.The population grid has ", x = tmpFile)]))
  nDemesGrid <- as.numeric(str_extract(string = nDemesGrid, pattern = "\\d+")) # noquote(gsub(pattern = " observed demes (out of \\[1-9]. demes)", replacement = "", x = ObsDemes))
  
  nDemes <- as.numeric(gsub(pattern = "[[:space:]]*.nDemes = ", replacement = "", x = tmpFile[grep(pattern = "nDemes = ", x = tmpFile)]))
  
  numMCMCIter <- as.numeric(gsub(pattern = "[[:space:]]*.numMCMCIter = ", replacement = "", x = tmpFile[grep(pattern = "numMCMCIter = ", x = tmpFile)]))
  
  numBurnIter <-  as.numeric(gsub(pattern = "[[:space:]]*.numBurnIter = ", replacement = "", x = tmpFile[grep(pattern = "numBurnIter = ", x = tmpFile)]))
  
  numThinIter <-   as.numeric(gsub(pattern = "[[:space:]]*.numThinIter = ", replacement = "", x = tmpFile[grep(pattern = "numThinIter = ", x = tmpFile)]))
  
  DoF <- as.numeric(gsub(pattern = " and effective degrees of freedom = ", replacement = "", x = tail(tmpFile[grep(pattern = " and effective degrees of freedom = ", x = tmpFile)], 1)))

  final.logPrior <- as.numeric(gsub(pattern = "Final log prior: ", replacement = "", x = tmpFile[grep(pattern = "Final log prior: ", x = tmpFile)]))
  
  final.logLike <- as.numeric(gsub(pattern = "Final log llike: ", replacement = "", x = tmpFile[grep(pattern = "Final log llike: ", x = tmpFile)]))
  
  tmp.df <- data.frame(Run, n, nDemes, ObsDemes, nDemesGrid, numMCMCIter, numBurnIter, numThinIter, DoF, final.logPrior, final.logLike)
  eemsDF <- rbind(eemsDF, tmp.df)
}

eemsDF
eemsDF[which(eemsDF$final.logLike == max(eemsDF$final.logLike)),]

###### logLikelihood by nDemes ######
plot(eemsDF$final.logLike ~ eemsDF$nDemes,
     ylab = "logLikelihood",
     xlab = "nDemes",
     main = "Effects of nDemes on logLikelihood",
     pch = 19,
     col = "red")

###### logLikelihood by degrees of freedom ######
plot(eemsDF$final.logLike ~ eemsDF$DoF,
     ylab = "logLikelihood",
     xlab = "Degrees of Freedom",
     main = "Effects of Degrees on Freedom on logLikelihood",
     pch = 19,
     col = "red")

###### Pick best by changes in slope by nDemes ######
eemsDF <- eemsDF[order(eemsDF$nDemes),]
eemsDF$slope <- c(NA, diff(eemsDF$ObsDemes)/diff(eemsDF$nDemes))
eemsDF$slope_chg <- c(NA, round(diff(eemsDF$slope),5))
eemsDF$change <- ifelse(eemsDF$slope_chg != 0, "change","")
eemsDF 

# Plot changes in observed demes per change in nDemes parameter. Add vertical line
# where change in slope is == 0
plot(eemsDF$ObsDemes ~ eemsDF$nDemes,
     ylab = "Observed Demes",
     xlab = "nDemes",
     main = "Effects of nDemes on Observed Demes",
     pch = 19,
     col = "red")
abline(v = eemsDF$nDemes[which(eemsDF$slope_chg <= 0)][1]) # Plot line at first point
# that has slope change == 0

###### Calculate AIC ######

for (m in 1:nrow(eemsDF)) {
  eemsDF$AIC[m] <- -2*(eemsDF$final.logLike[m])+2*eemsDF$DoF[m]
}

eemsDF[which(eemsDF$AIC == max(eemsDF$AIC)),]
eemsDF[which(eemsDF$AIC == min(eemsDF$AIC)),]

# Calculate delta AIC scores
for (m in 1:nrow(eemsDF)) {
  eemsDF$dAIC[m] <- eemsDF$AIC[m] - min(eemsDF$AIC)
}
eemsDF

# Calculate AIC weights
for (m in 1:nrow(eemsDF)) {
eemsDF$AICw <- round(-0.5 * eemsDF$dAIC)/sum((-0.5 * eemsDF$dAIC)) # IS THIS RIGHT?!?
}
eemsDF

####

#### Genetic-Environment Associations ####
# Useful workshop and tutorial links for RDA and variance partitioning analyses:
# https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
# https://r.qcbs.ca/workshop10/book-en/partial-redundancy-analysis.html
# https://r.qcbs.ca/workshop10/book-en/variation-partitioning.html
# https://www.mooreecology.com/uploads/2/4/2/1/24213970/constrained_ordination.pdf
# https://www.davidzeleny.net/anadat-r/doku.php/en:rda_cca_examples

##### RDA #####
setwd(project.folder)

# genetic data

# Use system() function to run the vcftools command
# system("vcftools --gzvcf 0_Data/TestVCF_5pop_25samp_550snp.vcf.gz --maf 0.02 --out 0_Data/TestVCF_maf02 --012")

df012 <- fread("0_Data/TestVCF_maf02.012", sep = "\t", data.table = F)
df012 <- df012[,-1]
dim(df012)

Pop_Sum <- read.csv("0_Data/Pop_Sum.csv")
pca_df <- read.csv("2_SNP_PCA/PCA_df.csv")

##### Normalize df012 using patterson et al 2006 standardization #####
df012 <- apply(df012, 2, function(d) gsub(-1, NA, d, fixed = TRUE))
df012 <- apply(df012, 2, function(d) as.numeric(d))

colmean <- apply(df012, 2, mean, na.rm = TRUE)
normalize <- matrix(nrow = nrow(df012), ncol = ncol(df012))
af <- colmean/2
for (m in 1:length(af)) {
  nr <- df012[, m] - colmean[m]
  dn <- sqrt(af[m] * (1 - af[m]))
  normalize[, m] <- nr/dn
}

# if NA, make equal to the mean (which is 0)
normalize[is.na(normalize)] <- 0

##### load environmental variables #####
env <- read.csv("0_Data/env_samp.csv")
names(env)
nrow(env)

env <- env[, c(4:ncol(env))]

# run rda with no conditioning
m <- rda(formula = normalize ~ ., scale = FALSE, data = env)
summary(m)

saveRDS(m, file = "7_GEA/RDA_latlong.RDS") # Save the RDA model

m <- readRDS("7_GEA/RDA_latlong.RDS") # reload data, if needed

summary(m)

##### Run ordiplot command #####
ordiplot(m, scaling = 2, type = "text")

##### Run ordiR2step command to get variable importance ####
fwd.sel <- ordiR2step(rda(formula = normalize ~ 1, scale = FALSE, data = env), m) #rda(formula = normalize ~ ., scale = FALSE, data = env))
fwd.sel
fwd.sel$anova

write.csv(fwd.sel$anova, "7_GEA/ordiR2step_out.csv")

##### Run anova.cca command to get variable importance ####
cca.term <- anova.cca(m, step = 1000, by = "term")
cca.term
cca.term[order(cca.term$F, decreasing = T),]
fwd.sel$anova[order(fwd.sel$anova$F, decreasing = T),]
fwd.sel

write.csv(cca.term[order(cca.term$F, decreasing = T),], "7_GEA/anova_cca.csv")

perc <- round(100*(summary(m)$cont$importance[2, 1:2]), 2)

RsquareAdj(m)

summary(eigenvals(m, model = "constrained"))[, 1:6]

# RDA1
RDA1imp <- sort(abs(summary(m)$biplot[, 1]), decreasing = TRUE)
RDA1imp

# weighted total loadings
RDAimp <- sort(rowSums(abs(summary(m)$biplot) * summary(eigenvals(m, model = "constrained"))[2, ]), decreasing = TRUE)
RDAimp

R.sum <- summary(m)
R.sum$cont   # Prints the "Importance of components" table
R.sum$cont$importance[2, "RDA1"]
# 0.74785
R.sum$cont$importance[2, "RDA2"]
# [1] 0.19804

##### RDA plot #####
# genetic data 
m <- readRDS("7_GEA/RDA_latlong.RDS")
sampleData <- read.csv("0_Data/SampleData_5pop_25samp_550snp.csv")

pca_df <- read.csv("2_SNP_PCA/PCA_df.csv")

rda_df <- data.frame(State=as.character(sampleData$State),
                     Population = as.character(sampleData$Population),
                     ID=as.character(sampleData$sample))

m_sum_pve <- summary(m)
RDA1_pve <- paste("RDA1 (",round((m_sum_pve$concont$importance[2,1]*100),digits=2),"%)", sep="")
RDA2_pve <- paste("RDA2 (",round((m_sum_pve$concont$importance[2,2]*100),digits=2),"%)", sep="")
RDA3_pve <- paste("RDA3 (",round((m_sum_pve$concont$importance[2,3]*100),digits=2),"%)", sep="")
RDA4_pve <- paste("RDA4 (",round((m_sum_pve$concont$importance[2,4]*100),digits=2),"%)", sep="")

m_sum <- summary(m)
m_sum <- scores(m,choices = c(1,2,3,4),display=c("sp","sites", "bp")) 
rda_snp <- as.data.frame(m_sum$species)
rda_indv <- as.data.frame(m_sum$sites)
rda_indv <- cbind(rda_df,rda_indv)
rda_biplot <- as.data.frame(m_sum$biplot)
rda_biplot$var <- row.names(rda_biplot)

# Add color to the data frame to make manipulation easier
Color <- NULL
rda_indv$Population # From TableS1
for (p in 1:nrow(pop.tab)) {
  for (r in 1:nrow(rda_indv)) {
    if (rda_indv$Population[r] == pop.tab$Population[p]) {
      Color[r] <- pop.tab$`Color code`[p]
    }
  }
}
Color

rda_indv <- cbind(rda_indv, Color)

rda_indv[which(rda_indv$Population == "C"),]


###### function to get what to multiple scores by for plotting ######
getMult <- function(m,choices=c(1,2),fill=1){
  site_score <- scores(m, choices = choices, display = "sites")
  bp_score <- scores(m, choices = choices, display = "bp")
  u <- c(range(site_score[, 1], na.rm = TRUE), 
         range(site_score[, 2], na.rm = TRUE))
  r <- c(range(bp_score[, 1], na.rm = TRUE), 
         range(bp_score[, 2], na.rm = TRUE))
  u <- u/r
  u <- u[is.finite(u) & u > 0]
  return(fill * min(u))
}

mult12 <- getMult(m)
mult23 <- getMult(m,choices = c(2,3))
mult34 <- getMult(m,choices = c(3,4))


# individual plots

rda12_plot <- ggplot(data=rda_indv) + 
  geom_point(data=rda_indv, aes(x=RDA2, y=RDA1,shape=State,fill=Population),size=3) +
  scale_shape_manual(values = c(21,22,24,25), name = 'State') +
  scale_fill_manual(values = my.cols.test, name = 'Population') +
  geom_segment(data = rda_biplot,
               aes(x = 0, xend = mult12 * RDA2,y = 0, yend = mult12 * RDA1),
               arrow = arrow(length = unit(0.25, "cm")), colour = "darkgrey") + #grid is required for arrow to work.
  geom_label_repel(data = rda_biplot,
                   aes(x= (mult12 + mult12/8) * RDA2, y = (mult12 + mult12/8) * RDA1, #we add 10% to the text to push it slightly out from arrows
                       label = var), #otherwise you could use hjust and vjust. I prefer this option
                   size = 4,fontface = "bold") + 
  xlab(RDA2_pve) +
  ylab(RDA1_pve) +
  theme_bw()  + 
  geom_vline(xintercept=0,linetype="dashed",col="grey") +
  geom_hline(yintercept=0,linetype="dashed",col="grey") +
  theme(legend.position = "none",
        axis.text = element_text(size=13), 
        axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
rda12_plot  

rda34_plot <- ggplot(data=rda_indv) + 
  geom_point(data=rda_indv, aes(x=RDA4, y=RDA3*-1,shape=State,fill=Population),size=3) +
  scale_shape_manual(values = c(21,22,24,25), name = 'State') +
  scale_fill_manual(values = my.cols.test, name = 'Population') +
  geom_segment(data = rda_biplot,
               aes(x = 0, xend = mult34 * RDA4,y = 0, yend = mult34 * RDA3*-1),
               arrow = arrow(length = unit(0.25, "cm")), colour = "darkgrey") + #grid is required for arrow to work.
  geom_label_repel(data = rda_biplot,
                   aes(x= (mult34 + mult34/8) * RDA4, y = (mult34 + mult34/8) * RDA3*-1, #we add 10% to the text to push it slightly out from arrows
                       label = var), #otherwise you could use hjust and vjust. I prefer this option
                   size = 4,fontface = "bold") + 
  xlab(RDA4_pve) + 
  ylab(RDA3_pve) +
  theme_bw()  + 
  geom_vline(xintercept=0,linetype="dashed",col="grey") +
  geom_hline(yintercept=0,linetype="dashed",col="grey") +
  theme(legend.position = "none",
        axis.text = element_text(size=13), 
        axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
rda34_plot 

rda.plots <- ggarrange(rda12_plot,rda34_plot, ncol = 1, labels = c("B", "C"))

Pop_Sum <- read.csv("0_Data/Pop_Sum.csv", header = T)
nrow(Pop_Sum)
sort(Pop_Sum$Population)

# map <- get_stamenmap(bbox = c(left = -125, bottom = 38, right = -120,
#                                top = 49), zoom = 7, maptype = "terrain")

map <- get_googlemap(c(-123,43.75), zoom = 5, size = c(400,400), maptype = "satellite") # %>% ggmap()

mapplot <- ggmap(map) +
  geom_point(data = nPop_Sum, aes(x = Longitude,
                                  y = Latitude),
             size = 2,
             col = "black",
             fill = "white",
             pch = 21) +
  geom_label_repel(data = nPop_Sum,
                   aes(x = Longitude,
                       y = Latitude,
                       label = Pop,
                       fill = Pop),
                   box.padding = 0,
                   size = 1.5,
                   max.overlaps = 25,
                   colour = "black",
                   fontface = "bold") +
  scale_fill_manual(name = "Population:", values = my.cols.north) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme_bw() + 
  scale_x_continuous(limits = c(-125, -120), expand = c(0, 0)) +
  scale_y_continuous(limits = c(38, 49), expand = c(0, 0)) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15, colour = "black", face = "bold"),
        panel.border = element_rect(linewidth = 1, colour = "black"),
        legend.text = element_text(size = 10), legend.position = "none",
        legend.title = element_text(size = 12, face = "bold"), panel.grid = element_blank())


mapplot.N <- ggarrange(mapplot, ncol = 1, nrow = 1, labels = "A")

mapplot.N


#pdf("RDA1-3v2.pdf") # Fig. 2
ggarrange(mapplot.N, rda.plots, legend.plot, ncol = 3, labels = NULL)
#dev.off()

ggsave("=RDA1-3.pdf", height=6,width=11.5)


#### Pairwise Geographical distance ####

Pop_Sum <- read.csv("0_Data/Pop_Sum.csv")

GeoDist <- matrix(data = NA, nrow = nrow(Pop_Sum), ncol = nrow(Pop_Sum))
for (i in 1:nrow(GeoDist)) {
  for (j in 1:nrow(GeoDist)) {
    dH <- distHaversine(Pop_Sum[i, c(4, 3)], Pop_Sum[j, c(4,
                                                          3)])  #long,lat 
    GeoDist[i, j] <- dH
  }
}

##### convert geographic distance to km #####
GeoDist <- GeoDist/1000

write.table(GeoDist, file = "0_Data/GeoDist.txt", row.names = F,
            col.names = F, quote = F)

##### Load data for isolation analyses #####
GeoDist <- read.table('0_Data/GeoDist.txt')
neiD <- read.table('4_GeneticDiversity/maf02_neiD.txt')

PopOrder <- neiD[,1]

neiD <- neiD[,-1]

Env <- read.csv("0_Data/env_final.csv")
Env <- Env[,5:ncol(Env)] # Check which columns hace the environmental data you want!
EnvDist <- as.matrix(dist(Env))

##### IBD, IBE #####
# vegan::mantel(as.matrix(neiD),as.matrix(GeoDist))

# vegan::mantel(as.matrix(neiD),as.matrix(EnvDist))

vegan::mantel(as.matrix(EnvDist),as.matrix(GeoDist))

vmp.d <- vegan::mantel.partial(as.matrix(neiD), as.matrix(GeoDist), as.matrix(EnvDist)) # IBD w/ IBE removed
vmp.d

vmp.e <- vegan::mantel.partial(as.matrix(neiD), as.matrix(EnvDist), as.matrix(GeoDist)) # IBE w/ IBD removed
vmp.e

##### IB plot #####
neiD <- as.vector(as.dist(neiD))
geo <- as.vector(as.dist(GeoDist))
env <- as.vector(as.dist(EnvDist))

IBD_df <- data.frame(neiD,geo,env)

# neiD 
IBD_plot <- ggplot(data=IBD_df,aes(x=geo,y=neiD)) + 
  geom_point(pch=21,size=3,alpha=.8,fill='#83B9B9',colour='black') +
  #geom_smooth(method = "lm",se=FALSE) +
  stat_smooth(method = "lm",alpha=.5,fill='#83B9B9',linetype=0) +
  stat_smooth(method = "lm",colour='grey30',se=FALSE,linetype='dashed') +
  theme_bw() + xlab("Geographic Distance (km)") + ylab("Genetic Distance (Nei's D)") + 
  theme(axis.text = element_text(size=13), 
        axis.title = element_text(size = 16, colour="black",face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2.5, colour = "black"),
        legend.text = element_text(size = 11),
        #legend.position = "bottom",
        legend.position = 'none',
        legend.title = element_text(size = 13, face = "bold"))

#neiD env 
IBE_plot <- ggplot(data=IBD_df,aes(x=env/100,y=neiD)) + 
  geom_point(pch=21,size=3,alpha=.8,fill='#F7D28A',colour='black') +
  #geom_smooth(method = "lm",se=FALSE) +
  stat_smooth(method = "lm",alpha=.5,fill='#F7D28A',linetype=0) +
  stat_smooth(method = "lm",colour='grey30',se=FALSE,linetype='dashed') +
  theme_bw() + xlab("Environmental Distance") + ylab("Genetic Distance (Nei's D)") + 
  theme(axis.text = element_text(size=13), 
        axis.title = element_text(size = 16, colour="black",face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2.5, colour = "black"),
        legend.text = element_text(size = 11),
        #legend.position = "bottom",
        legend.position = 'none',
        legend.title = element_text(size = 13, face = "bold"))

### Geo Env ###
EBD_plot <- ggplot(data=IBD_df,aes(x=geo,y=env/100)) + 
  geom_point(pch=21,size=3,alpha=.8,fill='grey65',colour='black') +
  #geom_smooth(method = "lm",se=FALSE) +
  stat_smooth(method = "lm",alpha=.5,fill='grey65',linetype=0) +
  stat_smooth(method = "lm",colour='black',se=FALSE,linetype='dashed') +
  theme_bw() + ylab("Environmental Distance") + xlab("Geographic Distance (km)") + 
  theme(axis.text = element_text(size=13), 
        axis.title = element_text(size = 16, colour="black",face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2.5, colour = "black"),
        legend.text = element_text(size = 11),
        #legend.position = "bottom",
        legend.position = 'none',
        legend.title = element_text(size = 13, face = "bold"))


##### plot all #####

IBD_IBE_plot <- ggarrange(IBD_plot,IBE_plot,EBD_plot,align='hv',nrow = 1, labels = "AUTO")
IBD_IBE_plot
ggsave('7_GEA/IBD_IBE_plot.pdf',IBD_IBE_plot,height = 6,width=12)

##### varpart #####
GeoDist <- read.table('0_Data/GeoDist.txt')
neiD <- read.table('4_GeneticDiversity/maf02_neiD.txt')

PopOrder <- neiD[,1]

neiD <- neiD[,-1]

Env <- read.csv("0_Data/env_final.csv")
LatLon <- Env[,3:4]
Env <- Env[,-c(1:2)]
EnvDist <- as.matrix(dist(Env))

vp <- vegan::varpart(Y = neiD, ~ Latitude+Longitude, ~ Elevation+FFP+MAP+MAT+RH+PAS, data = Env)
summary(vp)

#pdf("3_Plots/VarParts.pdf")
par(mar = c(0.5, 0.5, 0.5, 0.5)) 
plot(vp,
     Xnames = c("Geography", "Environment"), # name the partitions
     bg = c("#83B9B9", "#F7D28A"), 
     alpha = 100, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1,
     xaxs="i", yaxs="i")
VarPlot <- recordPlot()

VarPlotplot <- ggarrange(VarPlot, nrow = 1, hjust = -1, vjust = 2, labels = "D") +
  theme(plot.margin = margin(0,0,0,0, 'lines'))
VarPlotplot

ggsave(plot = VarPlotplot, filename = "7_GEA/varplotplot.pdf")

dev.off()

##### IB and VarPart plots together ####
bigplot <- ggarrange(IBD_IBE_plot, VarPlotplot, nrow = 2, labels = NULL, align = "hv") 

bigplot

ggsave(filename = "7_GEA/IB_and_VarPlot.pdf", plot = bigplot)
####

#### Predicting genetic offset ####
### ~~~ gdm: Generalized Dissimilarity Modeling
# Vignette / tutorial: https://github.com/fitzLab-AL/gdm

##### Set working directory #####
setwd(paste0(project.folder, "/8_GDM/"))

##### Load sample data #####
Pop_Sum <- read.csv("../0_Data/Pop_Sum.csv")
nrow(Pop_Sum)

envDF.pop <- read.csv("../0_Data/env_pop.csv")

##### Set up the GDM data sets #####
# biological data; get columns with xy, site ID, and species data
sppTab <- Pop_Sum[, c("State", "Population", "Longitude", "Latitude")]
colnames(sppTab) <- c("State", "site", "Longitude", "Latitude")

# get columns with site ID, env. data, and xy-coordinates
envTab <- data.frame(site = Pop_Sum$Population,
                     Latitude = Pop_Sum$Latitude,
                     Longitude = Pop_Sum$Longitude,
                     envDF.pop[,3:8])  # Check which variables get added!
envTab

# x-y species list example
gdmTab <- formatsitepair(bioData = sppTab, 
                         bioFormat = 2, #x-y spp list
                         XColumn = "Longitude", 
                         YColumn = "Latitude",
                         sppColumn = "State", 
                         siteColumn = "site", 
                         predData = envTab)

##### Generate gdmDissim object #####
sampleData <- read.csv("../0_Data/SampleData_5pop_25samp_550snp.csv")
nrow(sampleData)

# Subset by Pop
sampleData$site <- paste0(sampleData$State, "_", sampleData$Population)
nrow(sampleData)

##### Load vcf #####
setwd(project.folder)
vcf.file <- "../0_Data/TestVCF_5pop_25samp_550snp.vcf.gz" 
vcf.for.snp <- read.vcfR( vcf.file, verbose = FALSE )
snp <- vcfR2genind(vcf.for.snp)
# You can access the data using the "@" sign:
snp
snp.genclone <- as.genclone(snp)
snp.genclone@pop <- as.factor(sampleData$site)

##### Calculate genetic distance #####
snp@pop <- as.factor(sampleData$site)
gen.pop <- genind2genpop(snp)
gen.dist <- dist.genpop(x = gen.pop, upper = F, method = 4)
gdmDissim <- as.matrix(gen.dist)

# Make sure data are in same order
# Biological distance matrix example
dim(gdmDissim)
#> [1] 94 94
gdmDissim[1:5, 1:5]

# get the site column from sppTab
site <- unique(sppTab$site)
# bind to gdmDissim
gdmDissim <- cbind(site = sort(site), gdmDissim)
gdmDissim[1:5, 1:5]

##### Set up for genetic distance gdm #####
gdmTab.dis <- formatsitepair(bioData = gdmDissim, 
                             bioFormat = 3, #diss matrix 
                             XColumn = "Longitude", 
                             YColumn = "Latitude", 
                             predData = envTab, 
                             siteColumn = "site")

##### Crop layers to match distribution #####
# No cropping for test species; uncomment to add cropping.
# fralat <- sf::st_read("~/Dropbox/PennStateU/Research/1_Fraxinus_latifolia_PopGen/11_ENM/BIEN_ranges/Fraxinus_latifolia.shp")
# fralat

swBioclims <- raster::stack(list.files("~/Dropbox/Environmental_Data/ClimateNA/", full.names = T))
# swBioclims.c <- raster::crop(x = swBioclims.full, fralat[1])
# swBioclims.c <- raster::mask(x = swBioclims.c, mask = fralat[1])

# swBioclims <- raster::crop(x = swBioclims.c, c(-125.0625, -114.9375, 38, 49.64033))
plot(swBioclims[[1]])
points(x = SampleData$Longitude, y = SampleData$Latitude)

#### Fitting GDM ####

gdm.1 <- gdm(data = gdmTab.dis, geo = F)
summary(gdm.1)

##### Plots #####
length(gdm.1$predictors) # get ideal of number of panels

gdm.1.splineDat <- isplineExtract(gdm.1)
str(gdm.1.splineDat)

gdm.1.pred <- predict(object = gdm.1, data = gdmTab.dis)

head(gdm.1.pred)

#
pdf("gdmTab_EucDist.pdf")
plot(gdmTab$distance, 
     gdm.1.pred, 
     xlab = "Observed dissimilarity", 
     ylab = "Predicted dissimilarity", 
     xlim = c(0,1), 
     ylim = c(0,1), 
     pch = 20, 
     col = rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))
dev.off()

##### Project gdm to future #####
###### Crop future layers ######
# fralat <- sf::st_read("~/Dropbox/PennStateU/Research/1_Fraxinus_latifolia_PopGen/11_ENM/BIEN_ranges/Fraxinus_latifolia.shp")
futRasts <- raster::stack(list.files("~/Dropbox/Environmental_Data/FutureClimateNA_370/", full.names = T))
# futRasts.c <- raster::crop(x = futRasts.full, fralat[1])
# futRasts.c <- raster::mask(x = futRasts.c, mask = fralat[1])
# futRasts <- raster::crop(x = futRasts.c, c(-125.0625, -114.9375, 38, 49.64033))
plot(futRasts[[1]])
points(x = sampleData$Longitude, y = sampleData$Latitude)

# Check the different model types and predictions and their effects on turnover
timePred <- predict(object = gdm.1, data = swBioclims, time = T, predRasts = futRasts)
raster::writeRaster(x = timePred, filename = "timePred_370.tif", overwrite = T)

# Plot future prediction
pdf("gdm_EucDist_FuturePrediction_370.pdf")
raster::plot(timePred,
             col = rgb.tables(1000),
             main = " Ensemble ssp585 2041-2070")
points(SampleData[c("Longitude", "Latitude")],
       pch = 19,
       cex = 0.5)
dev.off()

# Extract predicted offset values
timePred <- raster::raster("timePred_245.tif")
PredOff <- raster::extract(x = timePred, y = SampleData[,c("Longitude","Latitude")])
range(PredOff, na.rm = T)
mean(PredOff, na.rm = T)

summary(PredOff)
transRasts <- gdm.transform(model = gdm.1, data = swBioclims)

pdf("north_eucDistTransformed_layersv2_370.pdf")
raster::plot(transRasts, col=rgb.tables(1000), main = "Transformed layers")
dev.off()

###### Visualizing multi-dimensional biological patterns ######
# Get the data from the gdm transformed rasters as a table
rastDat <- na.omit(raster::getValues(transRasts))

# The PCA can be fit on a sample of grid cells if the rasters are large
rastDat <- raster::sampleRandom(transRasts, 50000) 

# perform the principle components analysis
pcaSamp <- prcomp(rastDat)
summary(pcaSamp)
pcaSamp$rotation

a1 <- pcaSamp$x[,1]
a2 <- pcaSamp$x[,2]
a3 <- pcaSamp$x[,3]
r <- a1+a2
b <- a3+a2-a1
g <- -a2
r <- (r-min(r)) / (max(r)-min(r)) * 255
g <- (g-min(g)) / (max(g)-min(g)) * 255
b <- (b-min(b)) / (max(b)-min(b)) * 255
cols <- rgb(r,g,b,max=255)
cols2 <- col2rgb(cols)
cols3 <- t(cols2)
gradients <- cbind(pcaSamp$x[,c("PC1","PC2")],cols3)

##Biplot specifics
nvs <- dim(pcaSamp$rotation)[1]
vec <- rownames(pcaSamp$rotation)
lv <- length(vec)
vind <- rownames(pcaSamp$rotation) %in% vec
scal <-100
xrng <- range(pcaSamp$x[, 1], pcaSamp$rotation[, 1]/scal) *1
yrng <- range(pcaSamp$x[, 2], pcaSamp$rotation[, 2]/scal) *1

# PCsites <- predict(transRasts, pcaSamp, index = 1:3) %>% # s.data.frame(raster::predict(transRasts, pcaSamp, index = 1:3)) # 
#   as.data.frame() %>%
#   add_column(Group_ID = sample.sites$Group_ID)

pt.cols <- rgb(r, g, b, maxColorValue = 255)
pcs.df <- as.data.frame(pcaSamp$x) %>%
  mutate(group = paste0("Group_", PC1, "_", PC2))

line.scale <- 10
lines.env <- as.data.frame(pcaSamp$rotation)/line.scale

# From https://github.com/mgdesaix/bcrf-climate/tree/main
p.biplot <- ggplot() +
  geom_point(data = pcs.df, aes(x = PC1, y = PC2, color = group),
             show.legend = F,
             size = 0.5, shape = 15) +
  # geom_point(data = PCsites, aes(x = layer.1, y=layer.2),
  #            shape = 19, size = 1, stroke = 1) +
  scale_color_manual(values = pt.cols) +
  geom_segment(aes(x = 0, xend = PC1, y = 0, yend = PC2),
               linewidth = 1,
               data = lines.env,
               arrow = grid::arrow(length = unit(0.1, "cm")),
               inherit.aes = F) +
  geom_label(aes(x = PC1,
                 y = PC2,
                 label = rownames(lines.env),
                 #angle = 0.5,
                 #hjust = -0.5,
                 fontface = "bold"), 
             label.size = NA,
             data = lines.env, 
             fill = '#dddddd80') +
  xlab("PC1") +
  ylab("PC2") +
  theme_bw() +
  coord_equal() +
  theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    # plot.background = element_rect(fill = "transparent",colour = NA),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )

p.biplot

ggsave(plot = p.biplot, filename = "biplotv3.pdf", width = 8, height = 6)

# Predict the first three principle components for every cell in the rasters
# note the use of the 'index' argument
pcaRast <- raster::predict(transRasts, pcaSamp, index = 1:3)

# scale the PCA rasters to make full use of the colour spectrum
pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255
pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255
pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255

raster::writeRaster(x = pcaRast, filename = "pcaRastv2_370.tif", overwrite = T)
# Plot the three PCA rasters simultaneously, each representing a different colour 
#  (red, green, blue)


pcaRast.full <- raster::stack("pcaRastv2_370.tif")
plot(pcaRast.full)
pcaRast.c <- raster::crop(x = pcaRast.full, fralat[1])
pcaRast.c <- raster::mask(x = pcaRast.c, mask = fralat[1])
pcaRast.c <- raster::crop(x = pcaRast.c, c(-125.0625, -114.9375, 38, 49.64033))
plot(pcaRast.c[[1]])

##### Future layer contrubition to offset by area #####
pdf("CROPnorth_eucDistFutureTransformedLayersv2_370.pdf")
raster::plotRGB(pcaRast.c,
                r = 1, g = 2, b = 3,
                main = "Variables contributing to turnover",
                axes = T,
                ext = raster::extent(-130,-110,38,50))
dev.off()

raster::writeRaster(x = pcaRast.c, filename = "CROPpcaRastv2_370.tif", overwrite = T)
####

##### Plot multiple GDM predictions ####
pred.245 <- raster::raster("timePred_245.tif")
pred.370 <- raster::raster("timePred_370.tif")
pred.585 <- raster::raster("timePred_585.tif")
pcaRast.c <- raster::brick("CROPpcaRastv2_585.tif")

#raster::extent(pcaRast.c) <- c(-130,-110,38,50)

line = 1
cex = 1
side = 3
adj=-0.5

pdf("../3_Plots/GDM_MultiPlot.pdf")
par( mfrow=c(2,2), mar=c(5,5,5,6) )
raster::plot(pred.245,
             col = rgb.tables(1000),
             main = " Ensemble ssp245 2041-2070")#,
             #ext = raster::extent(-130,-110,38,50))
  # points(SampleData[c("Longitude", "Latitude")],
  #        pch = 19,
  #        cex = 0.5)
mtext("A", side=side, line=line, cex=cex, adj=adj)

raster::plot(pred.370,
             col = rgb.tables(1000),
             main = " Ensemble ssp370 2041-2070")#,
             #ext = raster::extent(-130,-110,38,50))
# points(SampleData[c("Longitude", "Latitude")],
#        pch = 19,
#        cex = 0.5)
mtext("B", side=side, line=line, cex=cex, adj=adj)

raster::plot(pred.585,
             col = rgb.tables(1000),
             main = " Ensemble ssp585 2041-2070")#,
             #ext = raster::extent(-130,-110,38,50))
# points(SampleData[c("Longitude", "Latitude")],
#        pch = 19,
#        cex = 0.5)
mtext("C", side=side, line=line, cex=cex, adj=adj)

raster::plotRGB(pcaRast.c,
                r = 1, g = 2, b = 3,
                main = "Variables contributing to offset",
                axes = T,
                ext = raster::extent(-140,-100,38,50))
mtext("D", side=side, line=line, cex=cex, adj=adj)

dev.off()
####

#### Check the gdm predictions; error on test data :( #####