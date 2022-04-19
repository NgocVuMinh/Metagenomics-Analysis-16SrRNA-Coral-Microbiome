
# PAPER: 
# Nhung, D. T., & Van Ngoc, B. (2020). Bioinformatic approaches for analysis of coral-associated bacteria using R programming language. Vietnam Journal of Biotechnology, 18(4), 733-743.

# In-detail intruction for the DADA2 package: https://benjjneb.github.io/dada2/tutorial_1_8.html
# 
# Download files according to instructions of NCBI: https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
# After downloading, all files are in .sra format. 
# To convert them to fastq, use these commands: 
#      dir > list.txt # append all file names into a text file
#      for i in $(cat list.txt); do echo $i; date; fasterq-dump -S $i; done
# I only selected 20 out of 78 files downloaded for this script to run.


# ____________ INSTALLING AND LOADING PACKAGES ____________

# installing the dada2 package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

#  installing the phyloseq package
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

#  packages for visualization
install.packages("viridis")
install.packages("RColorBrewer")
install.packages("treemap")

library(phyloseq)
library(dada2)
library(Rcpp)
library(ggplot2)
library(dplyr)
library(treemap)
library(viridis)
library(RColorBrewer)


# ____________ FILES MANUPILATION ____________

# Create a path for writing and downloading data
path = "~/Desktop/fastq_20files" 

# List all file names 
list.files(path) 

# The forward and reverse sequences are grouped into two different groups: fnFs and fnRs.
fnFs = sort(list.files(path, pattern="_1.fastq", full.names = T))
fnRs = sort(list.files(path, pattern="_2.fastq", full.names = T))

# Visualize the quality of the first two reads of the sequence files
plotQualityProfile(fnFs[1: 2]) 
plotQualityProfile(fnRs[1: 2])

# Extract sample names
basename(fnFs)
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Put filtered items in the filtered sub-directory
filt_path = file.path(path, "filtered")

# Place filtered files in filtered/subdirectory
filtFs = file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs = file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Filter and trim low-quality sequences
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(220,220),
                     maxN = 0, maxEE = c(2,2), truncQ = 2, trimLeft = c(25,20), 
                     rm.phix = TRUE, compress = TRUE, multithread = T) 
out1 = as.data.frame(out)
View(out1) # TABLE 1


# ______________ LEARN THE ERROR RATES _________________

# Calculating the sequence error rates of reads
errF = learnErrors(filtFs, multithread = T)  
errR = learnErrors(filtRs, multithread = T)


# ______________ DEREPLICATION _________________

# Combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance” equal to the number of reads with that unique sequence
derepFs = derepFastq(filtFs, verbose=TRUE)
derepRs = derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names   
names(derepRs) <- sample.names


# _______ APPLYING THE DADA2 ALGORITHM TO THE DEREPLICATED DATA _______

dadaFs <- dada(derepFs, err=errF, multithread=F)
dadaRs <- dada(derepRs, err=errR, multithread=F)
dadaFs[[1]]   # Inspecting the returned dada-class object


# ______________ MERGING PAIRED READS _________________

# merge together the inferred forward and reverse sequences and eliminated the residual errors (non-overlap paired sequences) to construct the ASVs table (ASVs and their abundances).
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])


# ________________ CONSTRUCTING ASV TABLE _______________

seqtab <- makeSequenceTable(mergers)
# remove chimera
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)  

# calculate the percentage of remained sequences after removing chimera
sum(seqtab.nochim)/sum(seqtab) 


# ____________ ASSIGNING TAXONOMY ____________

# Download training data sets at https://zenodo.org/record/4587955#.YlecbchBy00

# Assign taxonomy to genus level
taxaRC <- assignTaxonomy(seqtab.nochim, "~/Desktop/silva_nr99_v138.1_train_set.fa", tryRC=TRUE, multithread = T) 

# Assign taxonomy to species level
taxaSp <- addSpecies(taxaRC, "~/Desktop/silva_species_assignment_v138.1.fa") 

# Inspecting taxonomy assignment
taxaRC_print = taxaRC
rownames(taxaRC_print) = NULL # remove sequence row names for display
class(taxaRC) # matrix, array
df_taxaRC = as.data.frame(taxaRC_print) # convert to data frame
View(df_taxaRC) 

taxaSp_print = taxaSp
rownames(taxaSp_print) = NULL
df_taxaSp = as.data.frame(taxaSp_print)
View(df_taxaSp) # TABLE 2


# ____________ MICROBIOM ANALYSIS ____________

# Contruct a sample data frame (data table) of the information above
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

# Combine the ASVs table, taxonomic table and sample data table into phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim,taxa_are_rows=FALSE), tax_table(taxaSp), sample_data(samdf))

# Filtering phyloseq objects at the "Phylum" level
Bacphy = ps%>%
  tax_glom(taxrank = "Phylum")%>% # agglomerate data of the same phylum
  transform_sample_counts(function(x) {x/sum(x)})%>% # convert to a diversity analysis object
  psmelt()%>% # combine format(?), melt and merge all into a data frame
  filter(Abundance > 0.02)%>% # filter out low diversity taxa
  arrange(Phylum) # arrange data frame alphabetically by phylum
class(Bacphy) # data frame

# Filtering phyloseq objects at the "Class" level
Bacclass = ps%>%
  tax_glom(taxrank = "Class")%>%
  transform_sample_counts(function(x) {x/sum(x)})%>%
  psmelt()%>%
  filter(Abundance > 0.02)%>%
  arrange(Class)
class(Bacclass) 

# Filtering phyloseq objects at the "Genus" level
Bacgenus = ps%>%
  tax_glom(taxrank = "Genus")%>%
  transform_sample_counts(function(x) {x/sum(x)})%>%
  psmelt()%>%
  filter(Abundance > 0.02)%>%
  arrange(Genus) 


# ____________ VISUALIZATION ____________

# Phylum Composition - FIGURE 1A
ggplot(data = Bacphy, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = viridis(10)) + # fill color
  ylab("Relative Abundance") + # assign y axis name
  scale_y_continuous(expand = c(0,0)) + # clear space below 0 on y-axis
  labs(title = "Phylum Composition of Microbiota", caption = "Visual by Ngoc") +
  theme(axis.text.x = element_text(angle = 45))

# Class Composition - FIGURE 1B
ggplot(data = Bacclass, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = viridis(10)) +
  ylab("Relative Abundance") + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(title = "Class Composition of Microbiota", caption = "Visual by Ngoc") +
  theme(axis.text.x = element_text(angle = 45))

# Genus Composition - FIGURE 2
pal = c(brewer.pal(10, "Paired"), brewer.pal(5, "RdBu"), 
        brewer.pal(5, "PuOr"), brewer.pal(5, "PiYG"), 
        brewer.pal(5, "BrBG")) # this is apparently not the best color palatte choice
mycolors <- colorRampPalette(pal)(34)
ggplot(data = Bacgenus, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = mycolors) +
  ylab('Relative Abundance') +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "Genus Composition of Microbiota", caption = "Visual by ngoc")+
  theme(axis.text.x = element_text(angle = 45))


# TREEMAP

# Classification levels of Class and Order - FIGURE 3A
treemap(Bacgenus, index = c("Class", "Order"),
        vSize = "Abundance", type = "index",
        fontsize.labels = c(15,12), 
        fontcolor.labels = c("white", "black"),
        fontface.labels = c(2,1), #1,2,3,4: lower, bold, italic, bold_italic
        bg.labels = "transparent", # background color of labels
        align.labels = list(c("center", "center"), c("left", "bottom")),
        overlap.labels = 0.5, # tolerance for overlap between labels
        inflate.labels = F, # if true, the label is larger if rectangle is larger,
        fontsize.title = 12,
        palette = "Blues")

# Classification levels of Family and Genus - FIGURE 3B
treemap(Bacgenus, index = c("Family", "Genus"),
        vSize = "Abundance", type = "index",
        fontsize.labels = c(15,12), 
        fontcolor.labels = c("white", "black"),
        fontface.labels = c(2,1), #1,2,3,4: lower, bold, italic, bold_italic
        bg.labels = "transparent", # background color of labels
        align.labels = list(c("center", "center"), c("left", "bottom")),
        overlap.labels = 0.5, # tolerance for overlap between labels
        inflate.labels = F, # if true, the label is larger if rectangle is larger,
        fontsize.title = 12,
        palette = "Blues")



