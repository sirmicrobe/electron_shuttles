#Ed shuttles project
sessionInfo()
theme_set(theme_bw())

#
###################### DEMULTIPLEXING #####################
#need to demultiplex all files
#srun --nodes=1 --ntasks=16 #to run interactive job on slurm
#MiSeq_20150914
source(idemp -b /mmfs1/home/marshallc/projects/ed_shuttles/OneDrive_1_3-4-2021/dmux_MiSeq_20150914.txt -I1 /mmfs1/home/marshallc/projects/ed_shuttles/OneDrive_1_3-4-2021/MiSeq_20150914/Undetermined_S0_L001_I1_001.fastq.gz -R1 /mmfs1/home/marshallc/projects/ed_shuttles/OneDrive_1_3-4-2021/MiSeq_20150914/Undetermined_S0_L001_R1_001.fastq.gz -R2 /mmfs1/home/marshallc/projects/ed_shuttles/OneDrive_1_3-4-2021/MiSeq_20150914/Undetermined_S0_L001_R2_001.fastq.gz -m 1 -o /mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20150914_dmux)
#MiSeq_Run_062016
source(idemp -b /mmfs1/home/marshallc/projects/ed_shuttles/OneDrive_1_3-4-2021/dmux_MiSeq_Run_062016.txt -I1 /mmfs1/home/marshallc/projects/ed_shuttles/OneDrive_1_3-4-2021/Misc_MiSeq_Run_062016/Undetermined_S0_L001_I1_001.fastq.gz -R1 /mmfs1/home/marshallc/projects/ed_shuttles/OneDrive_1_3-4-2021/Misc_MiSeq_Run_062016/Undetermined_S0_L001_R1_001.fastq.gz -R2 /mmfs1/home/marshallc/projects/ed_shuttles/OneDrive_1_3-4-2021/Misc_MiSeq_Run_062016/Undetermined_S0_L001_R2_001.fastq.gz -m 1 -o /mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_Run_062016_dmux)
#MiSeq_20170407
source(idemp -b /mmfs1/home/marshallc/projects/ed_shuttles/OneDrive_1_3-4-2021/dmux_March2017.txt -I1 /mmfs1/home/marshallc/projects/ed_shuttles/OneDrive_1_3-4-2021/MiSeq_20170407/Undetermined_S0_L001_I1_001.fastq.gz -R1 /mmfs1/home/marshallc/projects/ed_shuttles/OneDrive_1_3-4-2021/MiSeq_20170407/Undetermined_S0_L001_R1_001.fastq.gz -R2 /mmfs1/home/marshallc/projects/ed_shuttles/OneDrive_1_3-4-2021/MiSeq_20170407/Undetermined_S0_L001_R2_001.fastq.gz -m 1 -o /mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20170407_dmux)

source(mkdir FWD)
source(mkdir REV)
source(mv *_R1_001.fastq* FWD/)
source(mv *_R2_001.fastq* REV/)

ex_file <- "Undetermined_S0_L001_R2_001.fastq.gz_NQL.C.05.fastq.gz"
ex_file <- sapply(strsplit(basename(ex_file), "*001.fastq.gz_"), `[`, 2)
ex_file <- sapply(strsplit(basename(ex_file), ".fastq.gz"), `[`, 1)

##################################################################################
# 2015 #
##################################################################################

#Filtering
#R 3.6.1
library(dada2); packageVersion("dada2") #1.14.1
# File parsing
pathF <- "/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20150914_dmux/FWD" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20150914_dmux//REV" # CHANGE ME ...
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS (151 bp reads)
filtered_out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),truncLen=c(145,140), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,compress=TRUE, verbose=TRUE, minLen=130, multithread=TRUE)
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),truncLen=c(145,140), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,compress=TRUE, verbose=TRUE, minLen=130, multithread=TRUE)

save.image("/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20150914_dmux/miseq2015.RData")
load("/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20150914_dmux/miseq2015.RData")

#Infer sequence variants

library(dada2); packageVersion("dada2") #1.18
# File parsing
filtpathF <- "/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20150914_dmux/FWD/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20150914_dmux/REV/filtered" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
#sample names parsing
sample.names <- sapply(strsplit(basename(filtFs), "*001.fastq.gz_"), `[`, 2)
sample.names <- sapply(strsplit(basename(sample.names), ".fastq.gz"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "*001.fastq.gz_"), `[`, 2)
sample.namesR <- sapply(strsplit(basename(sample.namesR), ".fastq.gz"), `[`, 1)


if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.namesR
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "/mmfs1/home/marshallc/projects/ed_shuttles/2015_seqtab.rds") # CHANGE ME to where you want sequence table saved



##################################################################################
# 2016 #
##################################################################################

setwd("/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_Run_062016_dmux")

source(mkdir FWD)
source(mkdir REV)
source(mv *_R1_001.fastq* FWD/)
source(mv *_R2_001.fastq* REV/)

#Filtering
library("dada2")
# File parsing
pathF <- "/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_Run_062016_dmux/FWD" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_Run_062016_dmux/REV" # CHANGE ME ...
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS (151 bp reads)
filtered_out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),truncLen=c(145,140), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,compress=TRUE, verbose=TRUE, minLen=130, multithread=TRUE)
#The filter removed all reads: /mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_Run_062016_dmux/FWD/filtered/Undetermined_S0_L001_R1_001.fastq.gz_AQS.100.uM.Batch2.C.07.fastq.gz and /mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_Run_062016_dmux/REV/filtered/Undetermined_S0_L001_R2_001.fastq.gz_AQS.100.uM.Batch2.C.07.fastq.gz not written.

save.image("/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_Run_062016_dmux/miseq2016.RData")
load("/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_Run_062016_dmux/miseq2016.RData")


#Infer sequence variants

# File parsing
filtpathF <- "/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_Run_062016_dmux/FWD/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_Run_062016_dmux/REV/filtered" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
#sample names parsing
sample.names <- sapply(strsplit(basename(filtFs), "*001.fastq.gz_"), `[`, 2)
sample.names <- sapply(strsplit(basename(sample.names), ".fastq.gz"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "*001.fastq.gz_"), `[`, 2)
sample.namesR <- sapply(strsplit(basename(sample.namesR), ".fastq.gz"), `[`, 1)


if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.namesR
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "/mmfs1/home/marshallc/projects/ed_shuttles/2016_seqtab.rds") # CHANGE ME to where you want sequence table saved

save.image("/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_Run_062016_dmux/miseq2016.RData")


################################################################################
# 2017 #
##################################################################################
source(rm Undetermined_S0_L001_R2_001.fastq.gz_unsigned.fastq.gz)
setwd("/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20170407_dmux")

source(mkdir FWD)
source(mkdir REV)
source(mv *_R1_001.fastq* FWD/)
source(mv *_R2_001.fastq* REV/)

#Filtering
#R 3.6.1
library(dada2); packageVersion("dada2") #1.14.1
# File parsing
pathF <- "/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20170407_dmux/FWD" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20170407_dmux/REV" # CHANGE ME ...
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS (151 bp reads)
filtered_out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),truncLen=c(145,140), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,compress=TRUE, verbose=TRUE, minLen=130, multithread=TRUE)

save.image("/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20170407_dmux/miseq2017.RData")
load("/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20170407_dmuxx/miseq2017.RData")


#Infer sequence variants

# File parsing
filtpathF <- "/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20170407_dmux/FWD/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20170407_dmux/REV/filtered" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
#sample names parsing
sample.names <- sapply(strsplit(basename(filtFs), "*001.fastq.gz_"), `[`, 2)
sample.names <- sapply(strsplit(basename(sample.names), ".fastq.gz"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "*001.fastq.gz_"), `[`, 2)
sample.namesR <- sapply(strsplit(basename(sample.namesR), ".fastq.gz"), `[`, 1)


if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.namesR
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "/mmfs1/home/marshallc/projects/ed_shuttles/2017_seqtab.rds") # CHANGE ME to where you want sequence table saved


save.image("/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20170407_dmux/miseq2017.RData")



##################################################################################
# Merge runs, Remove Chemeras, Assign taxonomy
##################################################################################

library(dada2); packageVersion("dada2")
# Merge multiple runs (if necessary)
st1 <- readRDS("/mmfs1/home/marshallc/projects/ed_shuttles/2015_seqtab.rds")
st2 <- readRDS("/mmfs1/home/marshallc/projects/ed_shuttles/2016_seqtab.rds")
st3 <- readRDS("/mmfs1/home/marshallc/projects/ed_shuttles/2017_seqtab.rds")
st.all <- mergeSequenceTables(st1, st2, st3)
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy - SILVA from Mar 11, 2021
tax <- assignTaxonomy(seqtab, "/mmfs1/home/marshallc/database/silva_db/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
tax <- addSpecies(tax, "/mmfs1/home/marshallc/database/silva_db/silva_species_assignment_v138.fa.gz")
# Write to disk
saveRDS(seqtab, "/mmfs1/home/marshallc/projects/ed_shuttles/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
saveRDS(tax, "/mmfs1/home/marshallc/projects/ed_shuttles/tax_final.rds") # CHANGE ME ...

#constructing a phylogenic tree 

library(DECIPHER)
seqs <- getSequences(seqtab)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

# phangorn tree building
library("phangorn")
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

#negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
#this has been running for 5 days and still not finished. could switch rearrangement = "stochastic" to rearrangement = "NNI" and it will run 100x faster with maybe not much difference in tree quality

saveRDS(fitGTR, "/mmfs1/home/marshallc/projects/ed_shuttles/tree.rds") 

#Tidying Up Sample Data
samples.out <- data.frame(Sample=rownames(seqtab))

save.image("/mmfs1/home/marshallc/projects/ed_shuttles/MiSeq_20170407_dmux/assign_tax_all.RData")


######################## bring down locally ####################
setwd("/Users/chrismarshall/Documents/1Marquette/Projects/Ed_shuttles")
#read in sequence and taxa tables
seqtab <- readRDS("/Users/chrismarshall/Documents/1Marquette/Projects/Ed_shuttles/seqtab_final.rds")
taxa <- readRDS("/Users/chrismarshall/Documents/1Marquette/Projects/Ed_shuttles/tax_final.rds")
fitGTR_HPC <- readRDS("/Users/chrismarshall/Documents/1Marquette/Projects/Ed_shuttles/tree.rds")

#read in mapping files/ sample metadata
sam2015 <- read.delim("/Users/chrismarshall/Documents/1Marquette/Projects/Ed_shuttles/MappingFile_eShuttles_MiSeq_20150914.txt",sep="\t")
sam2016 <- read.delim("/Users/chrismarshall/Documents/1Marquette/Projects/Ed_shuttles/MappingFile_eShuttles_Misc_MiSeq_Run_062016.txt",sep="\t")
sam2017 <- read.delim("/Users/chrismarshall/Documents/1Marquette/Projects/Ed_shuttles/Flynn_March2017Samples_MappingFile_JCK04062017.txt",sep="\t")

samdf <- rbind(sam2015,sam2016,sam2017)

#################
#sample data sheet ins and outs
#################
#reinoc samples were transferred from AQC batch 1 and into fresh media. (Rep A broke after event 6)
#microcosm sampels have higher sediment-to-media ratio so they are more diverse

#are the sample names the same?
rownames(seqtab)[rownames(seqtab) %in% samdf$X.SampleID ]
rownames(seqtab)[!(rownames(seqtab) %in% samdf$X.SampleID)] 
samdf$X.SampleID[!(samdf$X.SampleID%in% rownames(seqtab))] #this file is missing AQS.100.uM.Batch2.C.07

#remove this sample from samdf AQS.100.uM.Batch2.C.07
samdf <- subset(samdf, X.SampleID != "AQS.100.uM.Batch2.C.07")
samdf <- droplevels(samdf)
samdf$X.SampleID <- as.character(samdf$X.SampleID)
rownames(samdf) <- samdf$X.SampleID
which(samdf$ElectronShuttle=="inoculum") #77 201
samdf$ElectronShuttle[77] <- "Inoculum"
samdf$ElectronShuttle[201] <- "Inoculum"
samdf$eShuttle.Conc[387] <- "inoc"
samdf <- droplevels(samdf)
write.csv(samdf,"/Users/chrismarshall/Documents/1Marquette/Projects/Ed_shuttles/3map_samples.csv")
write.csv(data.frame(sort(as.character(samdf$X.SampleID)),sort(rownames(seqtab))),"samples.csv")


#R v 4.0.4
library(phyloseq); packageVersion("phyloseq") #1.34.0
library(ggplot2); packageVersion("ggplot2") #3.3.3
library(tidyverse); packageVersion("tidyverse") #1.3.1

####################### PHYLOSEQ OBJECT ################
theme_set(theme_bw())
setwd("/Users/chrismarshall/Documents/1Marquette/Projects/Ed_shuttles")
#reread data
samdf_read <- read.csv("3map_samples.csv", header = T)
rownames(samdf_read) <- samdf_read$X.SampleID

#are the sample names the same?
rownames(seqtab)[rownames(seqtab) %in% samdf_read$X.SampleID ]
rownames(seqtab)[!(rownames(seqtab) %in% samdf_read$X.SampleID)] 
samdf_read$X.SampleID[!(samdf_read$X.SampleID%in% rownames(seqtab))] # no, changed AQC.A.10 to AQC.A.07


#Create phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf_read), 
               tax_table(taxa),
               phy_tree(fitGTR$tree))
#19237 taxa and 387 samples

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps_orig <- ps

#AQC.100.um.Batch2.A.07 is major outlier on ordination
ps1 <- subset_samples(ps,X.SampleID!="AQC.100.uM.Batch2.A.07")
ps1 = prune_taxa(taxa_sums(ps1) > 0, ps1)

#Remove microcosm experiment from dataset
ps1_paper <- subset_samples(ps1,Experiment !="Microcosm.AQC.100.uM")
ps1_paper <- subset_samples(ps1_paper,Experiment !="Microcosm.AQC.1000.uM")
ps1_paper <- subset_samples(ps1_paper,Experiment !="AQC.Reinoc")

ps1_paper = prune_taxa(taxa_sums(ps1_paper) > 0, ps1_paper)

ps1_paper_noconc <- subset_samples(ps1_paper,eShuttle.Conc %in% c("100","inoc","none"))
ps1_paper_noconc = prune_taxa(taxa_sums(ps1_paper_noconc) > 0, ps1_paper_noconc)
#remove h2:co2 data
ps1_paper_noconc_noh2 <- subset_samples(ps1_paper_noconc, Experiment != "AQC.H2.CO2.100.uM")
ps1_paper_noconc_noh2 = prune_taxa(taxa_sums(ps1_paper_noconc_noh2) > 0, ps1_paper_noconc_noh2)

#remove batch 2 of aqc data
ps1_paper_noconc_noh2 <- subset_samples(ps1_paper_noconc_noh2, Experiment != "AQC.100.uM")
ps1_paper_noconc_noh2 = prune_taxa(taxa_sums(ps1_paper_noconc_noh2) > 0, ps1_paper_noconc_noh2) #11878 taxa and 217 samples 

sample_data(ps1_paper_noconc_noh2)$sample.rep.batch.days <- paste(sample_data(ps1_paper_noconc_noh2)$ElectronShuttle,sample_data(ps1_paper_noconc_noh2)$Replicate,sample_data(ps1_paper_noconc_noh2)$Batch,sample_data(ps1_paper_noconc_noh2)$days_smooth,sep=".")

sum(sample_sums(ps1_paper_noconc_noh2)) #6913275
median(sample_sums(ps1_paper_noconc_noh2)) #31242
mean(sample_sums(ps1_paper_noconc_noh2)) #31858.41
sd(sample_sums(ps1_paper_noconc_noh2)) #20617.82

sort(get_taxa_unique(ps1_paper_noconc_noh2, "Genus"))

#re-order factor levels for shuttles
sample_data(ps1_paper_noconc_noh2)$ElectronShuttle<- (as.factor(sample_data(ps1_paper_noconc_noh2)$ElectronShuttle))

sample_data(ps1_paper_noconc_noh2)$ElectronShuttle <- factor(sample_data(ps1_paper_noconc_noh2)$ElectronShuttle, levels=c('Inoculum', 'no.shuttle', 'AQC', 'AQDS', 'AQS', 'AQZ','NQJ', "NQL", "NQS", "Riboflavin"))
levels(sample_data(ps1_paper_noconc_noh2)$ElectronShuttle)

#reaorder sample.rep.batch.days

sample.rep.batch.days.levels <- data.frame(sample_data(ps1_paper_noconc_noh2))
sample.rep.batch.days.levels <- sample.rep.batch.days.levels[,c("sample.rep.batch.days","days_smooth","ElectronShuttle")]
sample.rep.batch.days.levels %>% 
  group_by(sample.rep.batch.days) %>%
  mutate(days_smooth=sort(as.numeric(days_smooth),decreasing=F)) %>%
  ungroup()

df.sample.rep.batch.days.levels <- sample.rep.batch.days.levels %>% 
  group_by(sample.rep.batch.days) %>% 
  arrange(ElectronShuttle,as.numeric(days_smooth)) %>%
  ungroup()

sample_data(ps1_paper_noconc_noh2)$sample.rep.batch.days <- factor(sample_data(ps1_paper_noconc_noh2)$sample.rep.batch.days, levels=df.sample.rep.batch.days.levels$sample.rep.batch.days)
levels(sample_data(ps1_paper_noconc_noh2)$sample.rep.batch.days)



##################################### Alpha Diversity ########################
#plotting Shannon Alpha diversity 


plot_richness(ps1_paper_noconc_noh2, x="ElectronShuttle", measures=c("Shannon")) + geom_boxplot() + geom_point() 
#with this command the graph will just be black and white

plot_richness(ps1_paper_noconc_noh2, x="ElectronShuttle", measures=c("Shannon"),color="ElectronShuttle") + geom_boxplot() + geom_point() + scale_color_manual(values=c("grey26","red","royalblue", "chartreuse3","darkorange","cyan2", "darkgreen","blue4", "yellow1","deepskyblue", "mediumorchid3","orangered3","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770"))
plot_richness(ps1_paper_noconc_noh2, x="ElectronShuttle", measures=c("Shannon"), color="Replicate") + geom_boxplot() + geom_point(position = position_dodge(width = 0.75))
#timepoints
plot_richness(ps1_paper_noconc_noh2, x="ElectronShuttle", measures=c("Shannon"), color="timepoint") + geom_boxplot() + geom_point(position = position_dodge(width = 0.75))

d_rich <- estimate_richness(ps1_paper_noconc_noh2,measures="Shannon")

d_rich$sample_names <- sample_data(ps1_paper_noconc_noh2)$X.SampleID
d_rich$ElectronShuttle <- sample_data(ps1_paper_noconc_noh2)$ElectronShuttle
d_rich$sample.rep.batch.days <- sample_data(ps1_paper_noconc_noh2)$sample.rep.batch.days

#Kruskal wallis test
kruskal.test(Shannon ~ ElectronShuttle, data = d_rich) #Kruskal-Wallis chi-squared = 122.9, df = 9, p-value < 2.2e-16
library("FSA") #v0.8.30 Ogle, D.H., P. Wheeler, and A. Dinno. 2020. FSA: Fisheries Stock Analysis. R package version 0.8.30,https://github.com/droglenc/FSA.
dunnTest(Shannon ~ as.factor(ElectronShuttle),
         data=d_rich,
         method="bh") #multiple comparison test 
#                   Comparison         Z      P.unadj        P.adj
1               AQC - AQDS -5.9075282 3.472788e-09 3.906887e-08
2                AQC - AQS  1.7434164 8.126090e-02 1.108103e-01
3               AQDS - AQS  7.2942170 3.004010e-13 1.351804e-11
4                AQC - AQZ -0.4438640 6.571409e-01 6.877056e-01
5               AQDS - AQZ  5.3216782 1.028143e-07 7.711075e-07
6                AQS - AQZ -2.1139236 3.452179e-02 5.356830e-02
7           AQC - Inoculum -4.0497548 5.127132e-05 1.538140e-04
8          AQDS - Inoculum -1.0230162 3.063002e-01 3.445877e-01
9           AQS - Inoculum -4.8580413 1.185527e-06 5.927633e-06
10          AQZ - Inoculum -3.8031068 1.428927e-04 4.018858e-04
11        AQC - no.shuttle -5.0787105 3.800053e-07 2.137530e-06
12       AQDS - no.shuttle  2.2252150 2.606681e-02 4.344468e-02
13        AQS - no.shuttle -6.7829730 1.177276e-11 1.765914e-10
14        AQZ - no.shuttle -4.3240389 1.531982e-05 5.303016e-05
15   Inoculum - no.shuttle  2.0939115 3.626786e-02 5.440179e-02
16               AQC - NQJ -0.6418420 5.209758e-01 5.581883e-01
17              AQDS - NQJ  5.1375071 2.784069e-07 1.789759e-06
18               AQS - NQJ -2.3032628 2.126406e-02 3.680318e-02
19               AQZ - NQJ -0.1916913 8.479840e-01 8.479840e-01
20          Inoculum - NQJ  3.7072611 2.095129e-04 5.545929e-04
21        no.shuttle - NQJ  4.0928505 4.261025e-05 1.369615e-04
22               AQC - NQL -1.9715409 4.866204e-02 6.843100e-02
23              AQDS - NQL  3.6817790 2.316121e-04 5.790303e-04
24               AQS - NQL -3.5167914 4.367969e-04 1.034519e-03
25               AQZ - NQL -1.5009155 1.333774e-01 1.765289e-01
26          Inoculum - NQL  2.9910099 2.780565e-03 5.687520e-03
27        no.shuttle - NQL  2.3042913 2.120630e-02 3.817134e-02
28               NQJ - NQL -1.3167445 1.879243e-01 2.168357e-01
29               AQC - NQS -3.4344473 5.937635e-04 1.335968e-03
30              AQDS - NQS  2.0353963 4.181101e-02 6.069340e-02
31               AQS - NQS -4.8548897 1.204538e-06 5.420419e-06
32               AQZ - NQS -2.9514886 3.162462e-03 6.187427e-03
33          Inoculum - NQS  2.1338152 3.285791e-02 5.280736e-02
34        no.shuttle - NQS  0.3737611 7.085821e-01 7.246862e-01
35               NQJ - NQS -2.7764991 5.494778e-03 1.030271e-02
36               NQL - NQS -1.4750422 1.402012e-01 1.752515e-01
37        AQC - Riboflavin -5.4313855 5.591818e-08 5.032636e-07
38       AQDS - Riboflavin  0.6819839 4.952491e-01 5.435661e-01
39        AQS - Riboflavin -6.8838149 5.827056e-12 1.311088e-10
40        AQZ - Riboflavin -4.8291469 1.371192e-06 5.609422e-06
41   Inoculum - Riboflavin  1.3885333 1.649747e-01 1.953648e-01
42 no.shuttle - Riboflavin -1.5001314 1.335804e-01 1.717462e-01
43        NQJ - Riboflavin -4.6374556 3.527244e-06 1.322716e-05
44        NQL - Riboflavin -3.1387788 1.696535e-03 3.635432e-03
45        NQS - Riboflavin -1.4568993 1.451442e-01 1.765267e-01

rich_anova <- aov(Shannon ~ ElectronShuttle, data = d_rich) 
summary(rich_anova)                 
                #Df Sum Sq Mean Sq F value Pr(>F)    
#ElectronShuttle   9 112.47  12.497   40.81 <2e-16 ***
TukeyHSD(rich_anova)
#                             diff         lwr         upr     p adj
no.shuttle-Inoculum   -2.73366948 -3.78258570 -1.68475325 0.0000000
AQC-Inoculum          -4.00792276 -5.09181241 -2.92403310 0.0000000
AQDS-Inoculum         -2.44692977 -3.55070897 -1.34315057 0.0000000
AQS-Inoculum          -4.32055059 -5.41641730 -3.22468389 0.0000000
AQZ-Inoculum          -3.64809274 -4.74055083 -2.55563465 0.0000000
NQJ-Inoculum          -3.60767661 -4.70013470 -2.51521852 0.0000000
NQL-Inoculum          -3.18928799 -4.29306719 -2.08550879 0.0000000
NQS-Inoculum          -2.80074029 -3.92017672 -1.68130387 0.0000000
Riboflavin-Inoculum   -2.54949820 -3.64195629 -1.45704011 0.0000000
AQC-no.shuttle        -1.27425328 -1.70608526 -0.84242130 0.0000000
AQDS-no.shuttle        0.28673970 -0.19283364  0.76631304 0.6606335
AQS-no.shuttle        -1.58688112 -2.04795150 -1.12581073 0.0000000
AQZ-no.shuttle        -0.91442327 -1.36733246 -0.46151407 0.0000000
NQJ-no.shuttle        -0.87400713 -1.32691633 -0.42109794 0.0000002
NQL-no.shuttle        -0.45561852 -0.93519186  0.02395482 0.0781527
NQS-no.shuttle        -0.06707082 -0.58165800  0.44751637 0.9999936
Riboflavin-no.shuttle  0.18417128 -0.26873792  0.63708047 0.9525792
AQDS-AQC               1.56099298  1.00910338  2.11288258 0.0000000
AQS-AQC               -0.31262784 -0.84851723  0.22326155 0.6918623
AQZ-AQC                0.35983001 -0.16905399  0.88871401 0.4765192
NQJ-AQC                0.40024614 -0.12863785  0.92913014 0.3191485
NQL-AQC                0.81863476  0.26674516  1.37052436 0.0001659
NQS-AQC                1.20718246  0.62460928  1.78975565 0.0000000
Riboflavin-AQC         1.45842456  0.92954056  1.98730855 0.0000000
AQS-AQDS              -1.87362082 -2.44867662 -1.29856502 0.0000000
AQZ-AQDS              -1.20116297 -1.76969619 -0.63262975 0.0000000
NQJ-AQDS              -1.16074684 -1.72928005 -0.59221362 0.0000000
NQL-AQDS              -0.74235822 -1.33235302 -0.15236342 0.0031416
NQS-AQDS              -0.35381052 -0.97260229  0.26498125 0.7161525
Riboflavin-AQDS       -0.10256843 -0.67110164  0.46596479 0.9998971
AQZ-AQS                0.67245785  0.11944309  1.22547261 0.0051842
NQJ-AQS                0.71287398  0.15985922  1.26588874 0.0021567
NQL-AQS                1.13126260  0.55620680  1.70631840 0.0000001
NQS-AQS                1.51981030  0.91524553  2.12437507 0.0000000
Riboflavin-AQS         1.77105239  1.21803763  2.32406715 0.0000000
NQJ-AQZ                0.04041613 -0.50581291  0.58664518 1.0000000
NQL-AQZ                0.45880475 -0.10972847  1.02733797 0.2333064
NQS-AQZ                0.84735245  0.24898851  1.44571639 0.0004179
Riboflavin-AQZ         1.09859454  0.55236550  1.64482359 0.0000000
NQL-NQJ                0.41838862 -0.15014460  0.98692183 0.3593644
NQS-NQJ                0.80693632  0.20857238  1.40530026 0.0010178
Riboflavin-NQJ         1.05817841  0.51194937  1.60440746 0.0000001
NQS-NQL                0.38854770 -0.23024407  1.00733947 0.5941102
Riboflavin-NQL         0.63978979  0.07125658  1.20832301 0.0142857
Riboflavin-NQS         0.25124209 -0.34712185  0.84960603 0.9422825

########################### BETA DIVERSITY #########################
# preparing data as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps1_paper_noconc_noh2, function(otu) otu/sum(otu))

ord.nmds.bray <- ordinate(ps1_paper_noconc_noh2, method="NMDS", distance="bray") 
#Run 6 stress 0.1762635 
#... New best solution
#... Procrustes: rmse 0.03885477  max resid 0.1786448 
plot_ordination(ps.prop, ord.nmds.bray, color="ElectronShuttle", shape="Batch",title="Bray NMDS") + scale_shape_manual(values=c(15,17,8,19,19,19)) + scale_color_manual(values=c("black","grey","darkgreen", "chartreuse3","darkorange","cyan2","red","blue4", "yellow1","deepskyblue", "mediumorchid3","royalblue"))
#shuttle+batch
plot_ordination(ps.prop, ord.nmds.bray, color="ElectronShuttle", shape="Batch",title="Bray NMDS") + scale_shape_manual(values=c(15,17,8)) + scale_color_manual(values=c("red", "chartreuse3","darkorange","cyan2","black","grey","darkgreen","blue4", "yellow1","deepskyblue", "mediumorchid3","royalblue")) + geom_point(size=3)
#experiment
plot_ordination(ps.prop, ord.nmds.bray, color="Experiment", shape="eShuttle.Conc",title="Bray NMDS") + scale_shape_manual(values=c(15,17,8,19,19,19)) + scale_color_manual(values=c("darkgreen","darkgreen","darkolivegreen","darkolivegreen4","darkolivegreen3", "chartreuse1","chartreuse1","chartreuse1","darkorange","cyan2","black", "olivedrab3","olivedrab3","grey","grey60","grey80","red","blue4", "yellow1","deepskyblue"))#+stat_ellipse()

#calculate distance
library("vegan") #vegan 2.5-7
bray_ps.prop <- distance(ps.prop, method="bray")
bray_dist_matrix<- as.matrix(bray_ps.prop)


ps.prop.batch1 <- subset_samples(ps.prop, Batch=="Batch1")
bray_ps.prop.b1 <- distance(ps.prop.batch1 , method="bray")
bray_dist_b1_matrix<- as.matrix(bray_ps.prop.b1)

ps.prop.batch2 <- subset_samples(ps.prop, Batch=="Batch2")
bray_ps.prop.b2 <- distance(ps.prop.batch2, method="bray")
bray_dist_b2_matrix<- as.matrix(bray_ps.prop.b2)

ps.prop.batch3 <- subset_samples(ps.prop, Batch=="Batch3")
bray_ps.prop.b3 <- distance(ps.prop.batch3, method="bray")
bray_dist_b3_matrix<- as.matrix(bray_ps.prop.b3)



library("reshape2")
df_bray_dist <- melt(as.matrix(bray_ps.prop), varnames = c("row", "col"))
colnames(df_bray_dist) <- c("row","X.SampleID", "Bray")

df_bray_dist_samples <-left_join(df_bray_dist,data.frame(sample_data(ps1_paper_noconc_noh2))%>%select(X.SampleID,Experiment,ElectronShuttle,Batch),by="X.SampleID")
df_bray_dist_samples <-left_join(df_bray_dist_samples,data.frame(sample_data(ps1_paper_noconc_noh2))%>%select(X.SampleID,Experiment,ElectronShuttle,Batch),by=c("row" ="X.SampleID"),suffix=c(".col",".row"))
levels(as.factor(df_bray_dist_samples$Experiment.col))

set.seed(12345)
all_test <- adonis(bray_dist_matrix~ElectronShuttle, data=data.frame(sample_data(ps1_paper_noconc_noh2)),permutations=10000)
#Permutation: free
#Number of permutations: 10000

#Terms added sequentially (first to last)

#    Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)    
# ElectronShuttle   9    38.808  4.3120   35.67 0.60798 9.999e-05 ***
#  Residuals       207    25.023  0.1209         0.39202             
# Total           216    63.831                 1.00000           
all_test$aov.tab$`Pr(>F)`[1] #gives just p-value

adonis(bray_dist_matrix~ElectronShuttle*Batch, data=data.frame(sample_data(ps1_paper_noconc_noh2)),permutations=10000)
#Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
#ElectronShuttle         9    38.808  4.3120  43.135 0.60798 9.999e-05 ***
#  Batch                   2     4.042  2.0208  20.215 0.06332 9.999e-05 ***
#  ElectronShuttle:Batch   2     0.689  0.3444   3.445 0.01079 9.999e-05 ***
#  Residuals             203    20.293  0.1000         0.31791              
#Total                 216    63.831                 1.00000              

#bray_dist_matrix_md <- as_tibble(bray_dist_matrix)
#bray_dist_matrix_md$X.SampleID <- rownames(bray_dist_matrix)
#meta_dist_matrix <- inner_join(data.frame(sample_data(ps1_paper_noconc_noh2)),bray_dist_matrix_md,by="X.SampleID")
#subset(meta_dist_matrix, Batch=="Batch1" )

adonis(bray_dist_b1_matrix ~ ElectronShuttle, data=data.frame(sample_data(ps.prop.batch1)),permutations=10000)
#             Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
#ElectronShuttle  4   12.6770  3.1693  33.921 0.65019 9.999e-05 ***
#  Residuals       73    6.8205  0.0934         0.34981              
#Total           77   19.4975                 1.00000        


bd.b1 <- betadisper(bray_ps.prop.b1,data.frame(sample_data(ps.prop.batch1))$ElectronShuttle)
anova(bd.b1)
permutest(bd.b1,pairwise=T) #<0.05 means signficant different variation in samples
#nothing significantly different

bd.b2 <- betadisper(bray_ps.prop.b2,data.frame(sample_data(ps.prop.batch2))$ElectronShuttle)
anova(bd.b2)
permutest(bd.b2,pairwise=T) # AQS and AQZ significantly different variation from no shuttle

bd.b3 <- betadisper(bray_ps.prop.b3,data.frame(sample_data(ps.prop.batch3))$ElectronShuttle)
anova(bd.b3)
permutest(bd.b3,pairwise=T) # NQJ only significant difference with no shuttle

######### doing all pairwise adonis ######


# 'Inoculum', 'no.shuttle', 'AQC', 'AQDS', 'AQS', 'AQZ','NQJ', "NQL", "NQS", "Riboflavin"
shuttles <- c('AQC', 'AQDS', 'AQS', 'AQZ','NQJ', "NQL", "NQS", "Riboflavin")
shuttles.b1 <- c('AQC', 'AQDS',  "NQL")
shuttles.b2 <- c('AQS', 'AQZ')
shuttles.b3 <- c('NQJ', "NQS", "Riboflavin")


pairwise_p1 <- numeric()
for(i in 1:length(shuttles.b1)){
  singleshuttle <- shuttles.b1[i]
  ps.prop.b1.singleshuttle <- subset_samples(ps.prop.batch1, ElectronShuttle %in% c("no.shuttle",singleshuttle))
  bray_ps.prop.b1.shuttle <- distance(ps.prop.b1.singleshuttle , method="bray")
  bray_dist_b1_matrix.shuttle<- as.matrix(bray_ps.prop.b1.shuttle)
  set.seed(12345)
  permanova.pair.b1 <- adonis(bray_dist_b1_matrix.shuttle ~ ElectronShuttle, data=data.frame(sample_data(ps.prop.b1.singleshuttle)),permutations=10000)
  pairwise_p1[shuttles.b1[i]]<- permanova.pair.b1[["aov.tab"]][["Pr(>F)"]][[1]]
  print(shuttles.b1[i])
  #print(permanova.pair.b1)
  print(pairwise_p1)

}
  

pairwise_p2 <- numeric()
for(i in 1:length(shuttles.b2)){
  singleshuttle <- shuttles.b2[i]
  ps.prop.b2.singleshuttle <- subset_samples(ps.prop.batch2, ElectronShuttle %in% c("no.shuttle",singleshuttle))
  bray_ps.prop.b2.shuttle <- distance(ps.prop.b2.singleshuttle , method="bray")
  bray_dist_b2_matrix.shuttle<- as.matrix(bray_ps.prop.b2.shuttle)
  set.seed(12345)
  permanova.pair.b2 <- adonis(bray_dist_b2_matrix.shuttle ~ ElectronShuttle, data=data.frame(sample_data(ps.prop.b2.singleshuttle)),permutations=10000)
  pairwise_p2[shuttles.b2[i]]<- permanova.pair.b2[["aov.tab"]][["Pr(>F)"]][[1]]
  print(shuttles.b2[i])
  #print(permanova.pair.b1)
  print(pairwise_p2)
  
}



pairwise_p3 <- numeric()
for(i in 1:length(shuttles.b3)){
  singleshuttle <- shuttles.b3[i]
  ps.prop.b3.singleshuttle <- subset_samples(ps.prop.batch3, ElectronShuttle %in% c("no.shuttle",singleshuttle))
  bray_ps.prop.b3.shuttle <- distance(ps.prop.b3.singleshuttle , method="bray")
  bray_dist_b3_matrix.shuttle<- as.matrix(bray_ps.prop.b3.shuttle)
  set.seed(12345)
  permanova.pair.b3 <- adonis(bray_dist_b3_matrix.shuttle ~ ElectronShuttle, data=data.frame(sample_data(ps.prop.b3.singleshuttle)),permutations=10000)
  pairwise_p3[shuttles.b3[i]]<- permanova.pair.b3[["aov.tab"]][["Pr(>F)"]][[1]]
  print(shuttles.b3[i])
  #print(permanova.pair.b1)
  print(pairwise_p3)
  
}

all.pairwise <- c(pairwise_p1, pairwise_p2, pairwise_p3)
p.adjust(all.pairwise,method="BH")
#all pairwise adonis significantly different from no shuttle

ns.allbatch <- subset_samples(ps.prop, ElectronShuttle == "no.shuttle")
bray_ps.prop.no.shuttle <- distance(ns.allbatch , method="bray")
dist.bray_ps.prop.no.shuttle<- as.matrix(bray_ps.prop.no.shuttle)
set.seed(12345)
permanova.pair.ns <- adonis(dist.bray_ps.prop.no.shuttle ~ Batch, data=data.frame(sample_data(ns.allbatch)),permutations=10000) #p=9.99e-5***


nmds.batch1 <- ordinate(ps.prop.batch1,"NMDS","bray")
nmds.batch2 <- ordinate(ps.prop.batch2,"NMDS","bray")
nmds.batch3 <- ordinate(ps.prop.batch3,"NMDS","bray")

plot_ordination(ps.prop.batch1, nmds.batch1, color="ElectronShuttle",title="Bray NMDS")  + scale_color_manual(values=c("black", "gray","darkgreen","green","darkblue"))#+stat_ellipse()

plot_ordination(ps.prop.batch2, nmds.batch2, color="ElectronShuttle",title="Bray NMDS")  + scale_color_manual(values=c("black", "gray","darkorange","cyan2","red","blue4", "yellow1","deepskyblue"))#+stat_ellipse()

plot_ordination(ps.prop.batch3, nmds.batch3, color="ElectronShuttle",title="Bray NMDS")  + scale_color_manual(values=c("black", "gray","red","yellow1","deepskyblue"))#+stat_ellipse()



#
##################### transitioning to barplot ########################



p_all <- plot_bar(ps.prop, x="X.SampleID", fill="Phylum") 
p_all + geom_bar(stat="identity", position="stack") + scale_fill_manual(values = c("grey26","red","royalblue", "chartreuse3","darkorange","cyan2", "darkgreen","blue4", "yellow1","deepskyblue", "mediumorchid3","orangered3","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030",   "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861","brown2","darkcyan","orchid1","orange1","mediumpurple1","skyblue4","salmon1","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4","grey26", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7","#ABC6DF","#272617","#D6D676","#8C378E","#56133D","#C9734D","#6E7424","#A0276E","#6D6944","#DC6025","#4B2A30","#AB7187","#4E66B8","#E6B554","#6670EA","#A62E4C","#5EA568","#DE3E36","#7A5A1E","#6890AA","#BEE82C","#8E3631","#9294C6","#6D4C67","#412457","#B2E967","#E0937A","#A9895A","#406B3C","#4B2012","#E4C738","#605453","#91A787","#D7BFE7","#D35BEA","#326D5B","#2A412B","#274146","#DDE3DE","#392039","#A69695","#4E871C","#8A4E45","#904366")) + facet_wrap(~ElectronShuttle, scales="free_y") +coord_flip() +theme(legend.position = "none")

top100 <- names(sort(taxa_sums(ps.prop), decreasing=TRUE))[1:100]
ps.top100 <- prune_taxa(top100, ps.prop)

plot_bar(ps.top100, x="X.SampleID", fill="Family") + scale_fill_manual(values = c("grey26","red","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4")) +coord_flip() + facet_wrap(~ElectronShuttle, scales="free_y")
plot_bar(ps.top100, x="sample.rep.batch.days", fill="Family") + scale_fill_manual(values = c("grey26","red","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4")) +coord_flip() + facet_wrap(~ElectronShuttle, scales="free_y")

plot_bar(subset_taxa(ps.top100, Phylum %in% c("Euryarchaeota","Halobacterota")) , x="sample.rep.batch.days", fill="Genus")
plot_bar(subset_taxa(ps.prop, Genus=="Methanosarcina") , x="sample.rep.batch.days", fill="Genus") +facet_grid(~ElectronShuttle, scales="free_x")


sample_data(ps.prop)$sample.rep.batch.days 

plot_bar(subset_taxa(ps.prop, Family=="Geobacteraceae") , x="sample.rep.batch.days", fill="Genus") + facet_wrap(~ElectronShuttle,scales="free_x") +scale_fill_manual(values=c("darkgreen","yellowgreen","mediumpurple1","cyan2"))
geobacter_relab <- subset_taxa(ps.prop, Family=="Geobacteraceae") 

geo_all_relab <- subset_taxa(ps.prop, Family %in% c("Geobacteraceae", "Geothermobacter")) #"Prolixibacteraceae"

plot_richness(geo_all_relab, x="ElectronShuttle", measures=c("Shannon")) + geom_boxplot() + geom_point() 



library("tidyverse")
otu_geobacter_relab <- as.data.frame(otu_table(geobacter_relab))
otu_geobacter_relab$samples <- rownames(otu_geobacter_relab)
tib_otu_geobacter_relab <- otu_geobacter_relab %>% pivot_longer(cols=c(0:162), values_to="Frequency",names_to="Geobacteraceae")

top20_geobacter_names <- names(sort(taxa_sums(geobacter_relab ), decreasing=TRUE))[1:20]
geobacter_colors <- c("mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","royalblue","red","yellow1","salmon1","violetred", "#89C5DA","#d1972c" ,"#8A4E45", "#004461","#C84248", "#673770","#508578","skyblue4","yellowgreen")
names(geobacter_colors) <- top20_geobacter_names

all_geo_asv <- levels(as.factor(tib_otu_geobacter_relab$Geobacteraceae))
df_geo_colors <- as.data.frame(geobacter_colors)
df_geo_colors$Geobacteraceae <- rownames(df_geo_colors)
df_all_geo_asv <- as.data.frame(all_geo_asv)
colnames(df_all_geo_asv) <- "Geobacteraceae"
df_geo_colors_all <- right_join(df_geo_colors,df_all_geo_asv,by="Geobacteraceae")
df_geo_colors_all$geobacter_colors <- df_geo_colors_all$geobacter_colors %>% replace_na("grey")

tib_otu_geobacter_relab_color <- right_join(tib_otu_geobacter_relab,df_geo_colors_all,by="Geobacteraceae")
tib_otu_geobacter_relab_color$sample_name <- tib_otu_geobacter_relab_color$samples
tib_otu_geobacter_relab_color <- left_join(tib_otu_geobacter_relab_color,data.frame(sample_data(ps1_paper_noconc_noh2))%>%select(X.SampleID,ElectronShuttle,sample.rep.batch.days),by=c("samples"="X.SampleID"))

tib_otu_geobacter_relab_color$Geobacteraceae <- as.factor(tib_otu_geobacter_relab_color$Geobacteraceae)
all_geobacter_colors<- tib_otu_geobacter_relab_color$geobacter_colors
names(all_geobacter_colors) <- as.character(tib_otu_geobacter_relab_color$Geobacteraceae)


tib_otu_geobacter_relab_color$shuttle<- (as.factor(tib_otu_geobacter_relab_color$ElectronShuttle))


p_geo_only <- ggplot(data=tib_otu_geobacter_relab_color, aes(x=sample.rep.batch.days,y=Frequency,fill=Geobacteraceae)) + 
  geom_bar(stat="identity") +
  ylim(0,1) + 
  scale_fill_manual(values=all_geobacter_colors) + 
  coord_flip() 
p_geo_only + facet_wrap(~shuttle,scales="free") +theme(legend.position = "none")


################# AQC only ################
ps_AQC <- subset_samples(ps.prop, ElectronShuttle=="AQC")
top100_aqc <- names(sort(taxa_sums(ps_AQC), decreasing=TRUE))[1:100]
AQC.top100 <- prune_taxa(top100_aqc, ps_AQC)
plot_bar(AQC.top100 , x="sample.rep.batch.days", fill="Family") + scale_fill_manual(values = c("grey26","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#8A7C64","red","#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "yellow1", "#599861","brown2","#904366","darkcyan","orchid1","orange1","steelblue2","skyblue4","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45")) +coord_flip() 
plot_bar(AQC.top100 , x="sample.rep.batch.days", fill="Genus") + scale_fill_manual(values = c("grey26","red","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()

#ps_AQC_paper <- subset_samples(ps_AQC, Experiment !="AQC.H2.CO2.100.uM")

top10_aqc <- names(sort(taxa_sums(ps_AQC), decreasing=TRUE))[1:10]
AQC.top10 <- prune_taxa(top10_aqc, ps_AQC)
p10_AQC <- plot_bar(AQC.top10 , x="sample.rep.batch.days", fill="Family")+theme(legend.position = "none") + custom_colors +ylim(0,1)+coord_flip() 


#AQC communities
aqc.ord.nmds.bray <- ordinate(ps_AQC, method="NMDS", distance="bray") #Run 20 stress 0.2003899 ; Procrustes: rmse 0.0003951823  max resid 0.00377728 
plot_ordination( ps_AQC, aqc.ord.nmds.bray, color="Experiment", shape="Replicate",title="Bray NMDS") + scale_shape_manual(values=c(15,17,19,19,19)) + scale_color_manual(values=c("red","blue4", "mediumorchid3","grey","black","darkgreen","darkorange","darkcyan")) +geom_point(size=3)+geom_text(aes(label=timepoint),position=position_jitter(width=.05,height=.05))



########################### No Shuttle #################
#no shuttle only
ps_noshuttle <- subset_samples(ps.prop, ElectronShuttle %in% c("Inoculum","no.shuttle"))
top100_noshuttle <- names(sort(taxa_sums(ps_noshuttle), decreasing=TRUE))[1:100]
noshuttle.top100 <- prune_taxa(top100_noshuttle, ps_noshuttle)
plot_bar(noshuttle.top100 , x="sample.rep.batch.days", fill="Family") + scale_fill_manual(values = c("grey26","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "red","#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "steelblue2","#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()
plot_bar(noshuttle.top100  , x="sample.rep.batch.days", fill="Genus") + scale_fill_manual(values = c("grey26","red","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()

top10_noshuttle <- names(sort(taxa_sums(ps_noshuttle), decreasing=TRUE))[1:10]
noshuttle.top10 <- prune_taxa(top10_noshuttle, ps_noshuttle)
p10_noshuttle <- plot_bar(noshuttle.top10 , x="sample.rep.batch.days", fill="Family")+theme(legend.position = "none") + custom_colors +ylim(0,1)+coord_flip()
p10_noshuttle_legend <- plot_bar(noshuttle.top10 , x="sample.rep.batch.days", fill="Family") + custom_colors +ylim(0,1)+coord_flip()


#noshuttle communities
noshuttle.ord.nmds.bray <- ordinate(ps_noshuttle, method="NMDS", distance="bray") #Run 20 stress 0.1105012 ; Procrustes: rmse 0.001435094  max resid 0.00516266  
plot_ordination( ps_noshuttle, noshuttle.ord.nmds.bray, color="Experiment", shape="Replicate",title="Bray NMDS") + scale_shape_manual(values=c(15,17,19,19,19)) + scale_color_manual(values=c("red","blue4", "mediumorchid3","grey","black","darkgreen","darkorange","darkcyan")) +geom_point(size=3)+geom_text(aes(label=timepoint),position=position_jitter(width=.05,height=.05))

########################### AQZ #################

ps_aqz <- subset_samples(ps.prop, ElectronShuttle == "AQZ")
top100_aqz <- names(sort(taxa_sums(ps_aqz), decreasing=TRUE))[1:100]
aqz.top100 <- prune_taxa(top100_aqz, ps_aqz)
plot_bar(aqz.top100 , x="sample.rep.batch.days", fill="Family") + scale_fill_manual(values = c("grey26","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "red","#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "steelblue2","#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()
plot_bar(aqz.top100  , x="sample.rep.batch.days", fill="Genus") + scale_fill_manual(values = c("grey26","red","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()

top10_aqz <- names(sort(taxa_sums(ps_aqz), decreasing=TRUE))[1:10]
aqz.top10 <- prune_taxa(top10_aqz, ps_aqz)
p10_aqz <- plot_bar(aqz.top10 , x="sample.rep.batch.days", fill="Family") +theme(legend.position = "none")+ custom_colors +ylim(0,1) +coord_flip()


#aqz communities
aqz.ord.nmds.bray <- ordinate(ps_aqz, method="NMDS", distance="bray") #stress too close to 0
plot_ordination( ps_aqz, aqz.ord.nmds.bray, color="Experiment", shape="Replicate",title="Bray NMDS") + scale_shape_manual(values=c(15,17,19,19,19)) + scale_color_manual(values=c("red","blue4", "mediumorchid3","grey","black","darkgreen","darkorange","darkcyan")) +geom_point(size=3)#+geom_text(aes(label=timepoint),position=position_jitter(width=.05,height=.05))


######################### AQS ####################

ps_aqs <- subset_samples(ps.prop, ElectronShuttle == "AQS")
top100_aqs <- names(sort(taxa_sums(ps_aqs), decreasing=TRUE))[1:100]
aqs.top100 <- prune_taxa(top100_aqs, ps_aqs)
plot_bar(aqs.top100 , x="sample.rep.batch.days", fill="Family") + scale_fill_manual(values = c("grey26","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "red","#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "steelblue2","#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()
plot_bar(aqs.top100  , x="sample.rep.batch.days", fill="Genus") + scale_fill_manual(values = c("grey26","red","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()

top10_aqs <- names(sort(taxa_sums(ps_aqs), decreasing=TRUE))[1:10]
aqs.top10 <- prune_taxa(top10_aqs, ps_aqs)
p10_aqs <- plot_bar(aqs.top10 , x="sample.rep.batch.days", fill="Family") +theme(legend.position = "none")+ custom_colors +ylim(0,1) +coord_flip()


#aqz communities
aqs.ord.nmds.bray <- ordinate(ps_aqs, method="NMDS", distance="bray") #stress too close to 0
plot_ordination( ps_aqs, aqs.ord.nmds.bray, color="Experiment", shape="Replicate",title="Bray NMDS") + scale_shape_manual(values=c(15,17,19,19,19)) + scale_color_manual(values=c("red","blue4", "mediumorchid3","grey","black","darkgreen","darkorange","darkcyan")) +geom_point(size=3)#+geom_text(aes(label=timepoint),position=position_jitter(width=.05,height=.05))


########################### AQDS #################

ps_aqds <- subset_samples(ps.prop, ElectronShuttle == "AQDS")
top100_aqds <- names(sort(taxa_sums(ps_aqds), decreasing=TRUE))[1:100]
aqds.top100 <- prune_taxa(top100_aqds, ps_aqds)
plot_bar(aqds.top100 , x="sample.rep.batch.days", fill="Family") + scale_fill_manual(values = c("grey26","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "red","#508578","#74D944", "#C84248", "#673770", "#D3D93E",  "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "steelblue2","#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()
plot_bar(aqds.top100  , x="sample.rep.batch.days", fill="Genus") + scale_fill_manual(values = c("grey26","red","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()
plot_bar(aqds.top100  , x="sample.rep.batch.days", fill="Class") + scale_fill_manual(values = c("grey26","red","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()

#ps_aqds_paper <- subset_samples(ps_aqds, eShuttle.Conc == "100")
#sum(taxa_sums(ps_aqds_paper) == 0)
#ps_aqds_paper = prune_taxa(taxa_sums(ps_aqds_paper) > 0, ps_aqds_paper) #[1593 taxa and 18 samples ]

top10_aqds <- names(sort(taxa_sums(ps_aqds), decreasing=TRUE))[1:10]
aqds.top10 <- prune_taxa(top10_aqds, ps_aqds)
p10_aqds <- plot_bar(aqds.top10 , x="sample.rep.batch.days", fill="Family")+theme(legend.position = "none") + custom_colors+ylim(0,1) +coord_flip()


#aqds communities
aqds.ord.nmds.bray <- ordinate(ps_aqds, method="NMDS", distance="bray") #Run 20 stress 0.1025913  ; Procrustes: rmse 3.544963e-06  max resid 1.125111e-05 
plot_ordination( ps_aqds, aqds.ord.nmds.bray, color="Experiment", shape="Replicate",title="Bray NMDS") + scale_shape_manual(values=c(15,17,19,19,19)) + scale_color_manual(values=c("chartreuse2","green3","forestgreen")) +geom_point(size=3)+geom_text(aes(label=timepoint),position=position_jitter(width=.05,height=.05))

########################### NQL #################

ps_nql <- subset_samples(ps.prop, ElectronShuttle == "NQL")
top100_nql <- names(sort(taxa_sums(ps_nql), decreasing=TRUE))[1:100]
nql.top100 <- prune_taxa(top100_nql, ps_nql)
plot_bar(nql.top100 , x="sample.rep.batch.days", fill="Family") + scale_fill_manual(values = c("grey26","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724","#74D944", "#C84248", "#673770","red","#508578", "#D3D93E",  "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "steelblue2","#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()
plot_bar(nql.top100  , x="sample.rep.batch.days", fill="Genus") + scale_fill_manual(values = c("grey26","red","blue4","royalblue","mediumpurple1", "#6DDE88","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()

top10_nql <- names(sort(taxa_sums(ps_nql), decreasing=TRUE))[1:10]
nql.top10 <- prune_taxa(top10_nql, ps_nql)
p10_nql <- plot_bar(nql.top10 , x="sample.rep.batch.days", fill="Family")+theme(legend.position = "none") +custom_colors+ylim(0,1) +coord_flip()


#nql communities
nql.ord.nmds.bray <- ordinate(ps_nql, method="NMDS", distance="bray") #insufficient data
plot_ordination( ps_nql, nql.ord.nmds.bray, color="Experiment", shape="Replicate",title="Bray NMDS") + scale_shape_manual(values=c(15,17,19,19,19)) + scale_color_manual(values=c("blue4"))# +geom_point(size=3)+geom_text(aes(label=timepoint),position=position_jitter(width=.05,height=.05))

########################### NQJ #################

ps_nqj <- subset_samples(ps.prop, ElectronShuttle == "NQJ")
top100_nqj <- names(sort(taxa_sums(ps_nqj), decreasing=TRUE))[1:100]
nqj.top100 <- prune_taxa(top100_nqj, ps_nqj)
plot_bar(nqj.top100 , x="sample.rep.batch.days", fill="Family") + scale_fill_manual(values = c("grey26","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724","#74D944", "#C84248", "#673770","red","#508578", "#D3D93E",  "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "steelblue2","#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()
plot_bar(nqj.top100  , x="sample.rep.batch.days", fill="Genus") + scale_fill_manual(values = c("grey26","red","blue4","royalblue","mediumpurple1", "#6DDE88","cyan2", "darkgreen","darkorange","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()

top10_nqj <- names(sort(taxa_sums(ps_nqj), decreasing=TRUE))[1:10]
nqj.top10 <- prune_taxa(top10_nqj, ps_nqj)
p10_nqj<-plot_bar(nqj.top10 , x="sample.rep.batch.days", fill="Family") +theme(legend.position = "none")+ custom_colors +ylim(0,1) +coord_flip()



#nqj communities
nqj.ord.nmds.bray <- ordinate(ps_nqj, method="NMDS", distance="bray") #insufficient data; Procrustes: rmse 2.750087e-06  max resid 7.474323e-06 
plot_ordination( ps_nqj, nqj.ord.nmds.bray, color="Experiment", shape="Replicate",title="Bray NMDS") + scale_shape_manual(values=c(15,17,19,19,19)) + scale_color_manual(values=c("red","blue4", "mediumorchid3","grey","black","darkgreen","darkorange","darkcyan")) +geom_point(size=3)+geom_text(aes(label=timepoint),position=position_jitter(width=.05,height=.05))


########################### NQS #################

ps_nqs <- subset_samples(ps.prop, ElectronShuttle == "NQS")
top100_nqs <- names(sort(taxa_sums(ps_nqs), decreasing=TRUE))[1:100]
nqs.top100 <- prune_taxa(top100_nqs, ps_nqs)
plot_bar(nqs.top100 , x="sample.rep.batch.days", fill="Family") + scale_fill_manual(values = c("grey26","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724","#74D944", "#C84248", "#673770","red","#508578", "#D3D93E",  "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "steelblue2","#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()
plot_bar(nqs.top100  , x="sample.rep.batch.days", fill="Genus") + scale_fill_manual(values = c("grey26","red","blue4","royalblue","mediumpurple1", "#6DDE88","cyan2", "darkgreen","darkorange","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()

top10_nqs <- names(sort(taxa_sums(ps_nqs), decreasing=TRUE))[1:10]
nqs.top10 <- prune_taxa(top10_nqs, ps_nqs)
p10_nqs<-plot_bar(nqs.top10 , x="sample.rep.batch.days", fill="Family") +theme(legend.position = "none")+ custom_colors +ylim(0,1) +coord_flip()



#nqs communities
nqs.ord.nmds.bray <- ordinate(ps_nqs, method="NMDS", distance="bray") #insufficient data; Procrustes: rmse 2.750087e-06  max resid 7.474323e-06 
plot_ordination( ps_nqs, nqs.ord.nmds.bray, color="Experiment", shape="Replicate",title="Bray NMDS") + scale_shape_manual(values=c(15,17,19,19,19)) + scale_color_manual(values=c("red","blue4", "mediumorchid3","grey","black","darkgreen","darkorange","darkcyan")) +geom_point(size=3)+geom_text(aes(label=timepoint),position=position_jitter(width=.05,height=.05))


########################### Riboflavin #################

ps_flavin <- subset_samples(ps.prop, ElectronShuttle == "Riboflavin")
top100_flavin <- names(sort(taxa_sums(ps_flavin), decreasing=TRUE))[1:100]
flavin.top100 <- prune_taxa(top100_flavin, ps_flavin)
plot_bar(flavin.top100 , x="sample.rep.batch.days", fill="Family") + scale_fill_manual(values = c("grey26","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724","#74D944", "#C84248", "#673770","red","#508578", "#D3D93E",  "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0","steelblue2", "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()
plot_bar(flavin.top100  , x="sample.rep.batch.days", fill="Genus") + scale_fill_manual(values = c("grey26","red","blue4","royalblue","mediumpurple1", "#6DDE88","cyan2","darkorange", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C","#8A4E45","#904366")) +coord_flip()
top10_flavin <- names(sort(taxa_sums(ps_flavin), decreasing=TRUE))[1:10]
flavin.top10 <- prune_taxa(top10_flavin, ps_flavin)
p10_flavin <- plot_bar(flavin.top10 , x="sample.rep.batch.days", fill="Family") + custom_colors +ylim(0,1)+theme(legend.position = "none")+coord_flip()
p10_flavin

#flavin communities
flavin.ord.nmds.bray <- ordinate(ps_flavin, method="NMDS", distance="bray") #Run 20 stress 0.08824323 ; Procrustes: rmse 2.770151e-06  max resid 6.333081e-06 
plot_ordination( ps_flavin, flavin.ord.nmds.bray, color="Experiment", shape="Replicate",title="Bray NMDS") + scale_shape_manual(values=c(15,17,19,19,19)) + scale_color_manual(values=c("deepskyblue")) +geom_point(size=3)+geom_text(aes(label=timepoint),position=position_jitter(width=.05,height=.05))

###############################

#define custom color scale
myColors <-  c("mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","royalblue","red","yellow1","salmon1","violetred", "#89C5DA","#d1972c" ,"#8A4E45", "#004461","#C84248", "#673770","#508578")
names(myColors) <- c("Acholeplasmataceae","Aeromonadaceae","Burkholderiaceae","Chromobacteriaceae","Comamonadaceae","Erysipelotrichaceae","Geobacteraceae","Geothermobacter","Lachnospiraceae","Methanosarcinaceae","Paludibacteraceae", "Prolixibacteraceae" , "Pseudomonadaceae","Rhodocyclaceae","Rikenellaceae","Spirochaetaceae", "Sulfurospirillaceae","Trichloromonas" )
custom_colors <- scale_fill_manual(name = "Family", values = myColors)



#D3D93E",  "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", ", "#6DDE88", "#652926", "#7FDCC0","steelblue2", "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861","darkcyan","orchid1","orange1","skyblue4","yellowgreen", "lightslategrey","peachpuff","turquoise4","#4E871C",,"#904366"

#AQDS family NAs are all "Desulfuromonadia" Class
library("patchwork")
p10_noshuttle_legend
( p10_AQC + p10_aqds + p10_aqz + p10_aqs)|( p10_flavin + p10_nqj + p10_nql + p10_nqs)

 ############################# corncob #################
install.packages("corncob")
library("corncob")
library(phyloseq)
library(magrittr)
ps1_paper_noconc_noh2
set.seed(123)
da_analysis <- differentialTest(formula= ~ ElectronShuttle,
                                phi.formula = ~ ElectronShuttle,
                                formula_null = ~1,
                                phi.formula_null = ~ ElectronShuttle,
                                test = "Wald", boot=F,
                                data=ps1_paper_noconc_noh2,
                                fdr_cutoff = 0.05)
da_analysis$significant_taxa
plot(da_analysis)

#looking at reduced dataset - pairwise vs. no shuttle control
#AQDS
View(sample_data(ps1_paper_noconc_noh2))
ps1_aqds_v_no <- subset_samples(ps1_paper,ElectronShuttle %in% c("AQDS","no.shuttle"))
ps1_aqds_v_no = prune_taxa(taxa_sums(ps1_aqds_v_no) > 0, ps1_aqds_v_no)

da_analysis_aqds <- differentialTest(formula= ~ ElectronShuttle,
                                phi.formula = ~ ElectronShuttle,
                                formula_null = ~1,
                                phi.formula_null = ~ ElectronShuttle,
                                test = "Wald", boot=F,
                                data=ps1_aqds_v_no,
                                fdr_cutoff = 0.05)
da_analysis_aqds $significant_taxa
plot(da_analysis_aqds)


#Flavin

#AQZ

#AQC

ps1_aqc_v_no <- subset_samples(ps1_paper_noconc_noh2,ElectronShuttle %in% c("AQC","no.shuttle"))
ps1_aqc_v_no = prune_taxa(taxa_sums(ps1_aqc_v_no) > 0, ps1_aqc_v_no)

da_analysis_aqc <- differentialTest(formula= ~ ElectronShuttle,
                                     phi.formula = ~ ElectronShuttle,
                                     formula_null = ~1,
                                     phi.formula_null = ~ ElectronShuttle,
                                     test = "Wald", boot=F,
                                     data=ps1_aqc_v_no,
                                     fdr_cutoff = 0.05)
da_analysis_aqc$significant_taxa
plot(da_analysis_aqc)

