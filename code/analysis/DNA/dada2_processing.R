# Use this script to process sequencing data (filter trimmed reads by quality, look 
# for chimeras, assign taxonomy with dada2, and remove contaminants)

# this is setting the working directory to the location we were just working in
#setwd(paste0("/export/data1/projects/daniela/SR2113/June2022_16S/"))
setwd(paste0("/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/"))
#run this to confirm the working directory is correct
getwd()
#run this to get a list of your files
list.files()
#run this to get a list of the trimmed files in the folder "primer-trimmed-reads/"
list.files("primer-trimmed-reads/")
#Load the dada2 package and print its version
library(dada2)
packageVersion("dada2")

#Now we'll create some variables that will make running the following commands easier
# reading in our sample IDs into an R object
samples <- scan("unique-sample-IDs.txt", what = "character")
samples # will print them out

# creating a variable that holds the paths to all the forward reads
forward_reads <- paste0("primer-trimmed-reads/", samples, "-trimmed-R1.fastq.gz")
forward_reads # to peek at it

# and one for all the reverse reads
reverse_reads <- paste0("primer-trimmed-reads/", samples, "-trimmed-R2.fastq.gz")

# here we are making a directory to hold the output quality-filtered reads
dir.create("quality-filtered-reads")

# And here we are creating variables holding what WILL be the paths to
# all the quality-filtered reads when we actually make those files.
# These files don't exist yet at this point, but these variables will be the output arguments
# we give a command in a bit to make them
filtered_forward_reads <- paste0("quality-filtered-reads/", samples, "-trimmed-filtered-R1.fastq.gz")
filtered_reverse_reads <- paste0("quality-filtered-reads/", samples, "-trimmed-filtered-R2.fastq.gz")

#Quality trimming/filtering
#We can get an overview summary image of the quality of our reads with the
#following function, to start we are just going to look at one forward and one reverse read file

#Forward reads of one sample
plotQualityProfile(forward_reads[1])

#Forward reads of all samples
plotQualityProfile(forward_reads)

#Reverse reads of one sample
plotQualityProfile(reverse_reads[1])

#Reverse reads of all samples
#plotQualityProfile(reverse_reads)

#### OPTIONAL: quantitatively figure out trim lengths
#Skip this if uninterested in optimizing read lengths
```{r message=F, eval=T, error=F, warning=F}
get_stats_for_fastq <- function(raw_fastq, qthresh=30, n_subset=1e5) {
  quals <- ShortRead::qa(dirPath = raw_fastq, n=1e5)[['perCycle']]$quality
  
  quals <- quals %>% group_by(lane, Cycle) %>% summarise(long_score = rep(Score, Count)) %>% 
    summarise(mean_score = mean(long_score), sd_score=sd(long_score), med_score = median(long_score),
              mad_score=stats::mad(long_score, constant = 1))
  
  quals$quality <- 'good_quality'
  quals$quality[quals$mean_score < qthresh] <- 'bad_quality'
  
  #quals$lane <- factor(quals$lane, levels = unique(c(grep("_R1_", quals$lane, value=T), grep("_R2_", quals$lane, value=T))))
  quals$lane <- factor(quals$lane, levels = unique(c(grep("R1.fast", quals$lane, value=T), grep("R2.fast", quals$lane, value=T))))
  
  return(quals)
}

diving_deriv <- function(df, thresh=-0.06) {
  x <- predict(loess(mean_score~Cycle, df))  # loess smooth out the jitters
  y <- df$Cycle
  
  fderiv <- diff(x) / diff(y)
  
  possibles <- which(fderiv < thresh)
  
  return(min(possibles))
}
```
#Pick a quality threshold (30 is conventional), and label each position as being good or bad quality, based on each position's mean quality score being at least or below the threshold. Then do some metrics to figure out where quality starts to dip: 
# * first position with bad quality 
# * last position with bad quality 
# * middle between first and last, 
# * Estimate based on the position where 1st derivative starts decreasing quickly (very cautious)

```{r message=F, error=F, warning=F, eval=T}
quality_thresh <- 30
options(dplyr.summarise.groups=FALSE)

stats_quality <- get_stats_for_fastq(c(forward_reads[1:3], reverse_reads[1:3]), qthresh = quality_thresh)

first_bad <- sapply(levels(stats_quality$lane), function(x) 
  min(stats_quality$Cycle[stats_quality$quality == 'bad_quality' & (stats_quality$lane == x)]))
last_good <- sapply(levels(stats_quality$lane), function(x) 
  max(stats_quality$Cycle[stats_quality$quality == 'good_quality' & (stats_quality$lane == x)]))
middle <- first_bad + floor((last_good - first_bad)/2)
estimated <- sapply(levels(stats_quality$lane), 
                    function(x) diving_deriv(stats_quality[stats_quality$lane == x,], thresh = -0.06))

summary_df <- data.frame(lane=names(estimated), first_bad=first_bad, middle=middle, last_good=last_good, estimated=estimated)
#summary_df$lane <- factor(summary_df$lane, 
#                          levels = unique(c(grep("_R1_", summary_df$lane, value=T), 
#                                            grep("_R2_", summary_df$lane, value=T))))
summary_df$lane <- factor(summary_df$lane, 
                          levels = unique(c(grep("R1", summary_df$lane, value=T), 
                                            grep("R2", summary_df$lane, value=T))))

print(summary_df, row.names = F)

ggplot(stats_quality, aes(x=Cycle, y= mean_score)) + 
  geom_linerange(aes(ymin=mean_score - sd_score, ymax=mean_score + sd_score)) + 
  geom_point(pch="_", aes(color=quality)) +
  facet_wrap(~lane, ncol=3) +
  geom_hline(yintercept = 30, color='grey20', lty=2) +
  geom_vline(data=summary_df, aes(xintercept=first_bad), color='blue', lty=3) +
  geom_vline(data=summary_df, aes(xintercept=middle), color='grey50', lty=3) +
  geom_vline(data=summary_df, aes(xintercept=last_good), color='red', lty=3) +
  geom_vline(data=summary_df, aes(xintercept=estimated), color='cyan', lty=3) +
  theme_classic()
```

#Apply trimming of sequences based on the quality plots
#Since the upper limit for trimming for 16S (411 bp) is of 170, 
#to have at least 12 bp of overlap
#Here we left 

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                              rm.phix=TRUE, minLen=150, truncLen=c(240,200))

# We can see the filtered reads in this file
list.files("quality-filtered-reads/")
# This command allows to see how many sequences were filtered out
filtered_out
# Here we can see what % of reads they kept
filtered_out[, 2] / filtered_out[, 1] * 100

#To look at the quality-filtered profiles:
#For forward samples
plotQualityProfile(filtered_forward_reads[1:9])
#For one reverse sample
plotQualityProfile(filtered_reverse_reads[1:9])

#Generating an error model of the data
#Specific-error signature of the forward reads
err_forward_reads <- learnErrors(filtered_forward_reads, multithread = 4)
#Specific-error signature of the reverse reads
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread = 4)

# Plot to visualize how the error rates match with the observed
plotErrors(err_forward_reads[1:9])
# Plot to visualize how the error rates match with the observed
plotErrors(err_reverse_reads[1:9])

#Dereplicate sequences to find unique sequences
derep_forward <- derepFastq(filtered_forward_reads, verbose=FALSE)
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=FALSE)

#Inferring amplicon sequencing variants (ASVs)

#Inference of biological sequences from the forward filtered reads
dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread=4)

#Inference of biological sequences from the reverse filtered reads
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread=4)

#Merging forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, 
                               derep_forward, 
                               dada_reverse,
                               derep_reverse,
                               trimOverhang=TRUE,
                               minOverlap=12)

#Generate a count table
seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab) # see number of unique sequences across our dataset

# Chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose = TRUE, method = 'consensus')
# Identified 9834 bimeras out of 14856 input sequences.

#See the fraction of merged sequences
print("Merged fraction:")
rowSums(seqtab) / rowSums(makeSequenceTable(dada_forward))
# Although we lost almost 90% of our unique sequences, 
# we don’t know if those particular sequences held a lot in terms of abundance yet. Here is one quick way we can look at that
sum(seqtab.nochim) / sum(seqtab)
# In terms of abundance, we are retaining about 63% of reads following chimeral removal.

# Making an overview of counts through the process
# This makes a good supplemental figure
# Making a little helper function
getN <- function(x) sum(getUniques(x))

# Making a summary table
summary_tab <- data.frame(row.names = samples,
                          dada2_input = filtered_out[,1],
                          filtered = filtered_out[,2],
                          dada_f = sapply(dada_forward, getN),
                          dada_r = sapply(dada_reverse, getN), 
                          merged = sapply(merged_amplicons, getN),
                          nonchim = rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))


#To see the summary table                                               
summary_tab

#Write out the table to see it outside of R
write.table(summary_tab, "read-count-tracking.tsv", 
            quote = FALSE, sep = "\t", 
            col.names = NA)

library(reshape)
library(ggplot2)
x <- melt(summary_tab[,-ncol(summary_tab)], variable.name='Step', value.name='Reads')
x
ggplot(melt(summary_tab[,-ncol(summary_tab)], variable.name='Step', value.name='Reads'), aes(x = variable, y = value)) +
  #geom_bar(aes(fill = step), stat = 'identity') + theme_bw() +
  geom_boxplot(aes(color = value), outlier.shape = NA) + theme_bw() +
  geom_jitter(height=0, width=0.2, color='black', alpha=0.6) +
  theme(axis.text.x = element_text(angle = 90, color='black', hjust=1, vjust=0.5), panel.background = element_blank(),
        axis.line.x = element_line(), axis.line.y = element_line())

s <- read.csv("read-count-tracking.tsv", sep="\t", header=T)
print(s)
ggplot(melt(s[,-ncol(s)], variable.name='Step', value.name='Reads'), aes(x = variable, y = value)) +
  geom_boxplot(aes(color = value), outlier.shape = NA) + theme_bw() +
  geom_jitter(height=0, width=0.2, color='black', alpha=0.6) +
  theme(axis.text.x = element_text(angle = 90, color='black', hjust=1, vjust=0.5), panel.background = element_blank(),
        axis.line.x = element_line(), axis.line.y = element_line())
ggsave(file='/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/read-count-tracking_percent.pdf', width=8, height=6)

#Assigning taxonomy
# Load in the latest SILVA dataset
load("/export/data1/db/16S_tag_processing_db/SILVA_SSU_r138_2019_orphanlab_k8-210919.RData")  # SILVA + orphan seqs

# loading DECIPHER
library(DECIPHER)
packageVersion("DECIPHER") # v2.20.0

# creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

# and classifying
tax_info <- IdTaxa(test = dna, 
                   trainingSet = trainingSet, 
                   strand = "both", 
                   processors = 4)
# here is how we can save the R object if we want
# and next is read it in
save(tax_info, file = "decipher-tax-obj.RData")

#Make a directory to save files to export
dir.create("export")
newDirPath <- "export"
newFilePath <- "/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/decipher-tax-obj.RData"
file.copy(newFilePath, newDirPath)

#Extracting the standard goods from dada2
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
#Make an empty vector that is going to hold the new names
asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}
#Check that the headers were changed
asv_headers

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

a <-t(seqtab.nochim)

# Create count table:
asv_tab <- t(seqtab.nochim)
colnames(asv_tab) <- samples
#Find the > symbol and replace it with nothing in asv_headers and assign them to row.names
row.names(asv_tab) <- sub(">", "", asv_headers)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
# this part def gets busy, and i would need to break it down to remember
# what everything is doing/how it works
# so don't worry about it too much right now, it's more R-heavy than
# we're getting into currently

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
#asv_tax <- t(sapply(tax_info, function(x) {
#  m <- match(ranks, x$rank)
#  taxa <- x$taxon[m]
#  taxa[startsWith(taxa, "unclassified_")] <- NA
#  taxa
#}))
#colnames(asv_tax) <- ranks
#rownames(asv_tax) <- gsub(pattern = ">", replacement = "", x = asv_headers)

asv_tax <- lapply(tax_info, function(x) {
  taxa <- unlist(strsplit(x$taxon, ";"))[-1]  # drop "Root"
  taxa <- c(taxa, rep(taxa[length(taxa)], max(0,length(ranks)-length(taxa))))
  taxa <- taxa[1:length(ranks)]  # pesky euks
})
asv_tax <- matrix(unlist(asv_tax), ncol=length(ranks), byrow = T, dimnames = list(names(asv_tax), ranks))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)
head(asv_tax)
# Export the table outside of R
write.table(asv_tax, "ASVs-taxonomy.tsv", sep = "\t", quote = FALSE, col.names = NA)
#Read table
tax_tab <- read.table("ASVs-taxonomy.tsv", header = TRUE, row.names = 1,
                      check.names = FALSE, sep = "\t")
tax_tab$ASV <- rownames(tax_tab)

#Get ASVs from contaminants in the taxonomy table
ASV_contams<-subset(tax_tab, genus=="Corynebacterium" | genus=="Curtobacterium" | genus=="Cutibacterium" | genus=="Cutibacterium_1" |
                      genus=="Erysipelotrichaceae UCG-003" | genus=='Adhaeribacter' | genus=="Skermanella" | genus=="Staphylococcus_1" |
                      genus=="Staphylococcus" | genus=='Lactobacillus' | genus=='Lactococcus' | genus=='Negativicoccus' |
                      genus=="Turicella" | genus=="Hymenobacter" | genus=='Anaerococcus' | genus=='Atopobium' | genus=='Prevotella_7' |
                      genus=='Veillonella' | genus=='Streptococcus' | genus=='Stenotrophomonas' | genus=='Taibaiella' |
                      genus=='Auritidibacter' | genus=='Neisseria' | genus=='Gemella' | genus=='Finegoldia' | genus=='Rothia' |
                      genus=='Parvimonas' | genus=='Alloiococcus' | genus=='Prevotella' | genus=='Empedobacter' | genus=='Lawsonella' |
                      genus=='Capnocytophaga' | genus=='Pseudonocardia' | genus=='Lactobacillus_1' |
                      genus=='unclassified_Staphylococcaceae' | genus=='unclassified_Streptococcaceae' | genus=='Atopobium_2' |
                      genus=='Corynebacterium_1' | genus=='unclassified_Corynebacteriaceae' | genus=='unclassified_Pseudomonadaceae' |
                      genus=='unclassified_Lactobacillales' | genus=='Anaerococcus_1' | 
                      genus=='unclassified_Propionibacteriaceae' | genus=='Pseudomonas_1' | genus=='unclassified_Enterobacterales' |
                      genus=='unclassified_Enterobacteriaceae_1' | genus=='unclassified_Staphylococcales' | genus=='Acinetobacter' | genus=='unclassified_Staphylococcaceae' |
                      genus=='unclassified_Capnodiales' | genus=='unclassified_Intrasporangiaceae_1' | genus=='Campylobacter' |
                      genus=='Megasphaera_1' | genus=='Lautropia' | genus=='Prevotella_6' | genus=='Helcococcus' | genus=='Prevotella_5' | genus=='Empedobacter_1' |
                      genus=='Murdochiella' | genus=='Prevotella_9'| genus=='Prevotella_4' | genus=="Cutibacterium_2" | genus=='Prevotella_2'|
                      genus=='Finegoldia_1' | genus=='Corynebacteriaceae' | genus=='Afipia_2' ,select = ASV)

#Remove contaminant ASVs from counts_tab
write.table(asv_tab, "ASVs-counts.tsv", sep = "\t", quote = FALSE, col.names = NA)
counts_tab <- read.table("ASVs-counts.tsv", header = TRUE, row.names = 1,
                         check.names = FALSE, sep = "\t")
counts_tab$ASV <- rownames(counts_tab)
asv_tab <- counts_tab[ ! counts_tab$ASV %in% ASV_contams$ASV, ]
write.table(asv_tab, "ASVs-counts.tsv", sep = "\t", quote = FALSE, col.names = NA)

#Filter contamination out of the taxonomy table
asv_tax <- subset(tax_tab, genus!="Corynebacterium" & genus!="Curtobacterium" & genus!="Cutibacterium" & genus!="Cutibacterium_1" &
                    genus!="Erysipelotrichaceae UCG-003" & genus!='Adhaeribacter' & genus!="Skermanella" & genus!="Staphylococcus_1" &
                    genus!="Staphylococcus" & genus!='Lactobacillus' & genus!='Lactococcus' & genus!='Negativicoccus' &
                    genus!="Turicella" & genus!="Hymenobacter" & genus!='Anaerococcus' & genus!='Atopobium' & genus!='Prevotella_7' &
                    genus!='Veillonella' & genus!='Streptococcus' & genus!='Stenotrophomonas' & genus!='Taibaiella' &
                    genus!='Auritidibacter' & genus!='Neisseria' & genus!='Gemella' & genus!='Finegoldia' & genus!='Rothia' &
                    genus!='Parvimonas' & genus!='Alloiococcus' & genus!='Prevotella' & genus!='Empedobacter' & genus!='Lawsonella' &
                    genus!='Capnocytophaga' & genus!='Pseudonocardia' & genus!='Lactobacillus_1' &
                    genus!='unclassified_Staphylococcaceae' & genus!='unclassified_Streptococcaceae' & genus!='Atopobium_2' &
                    genus!='Corynebacterium_1' & genus!='unclassified_Corynebacteriaceae' & genus!='unclassified_Pseudomonadaceae' &
                    genus!='unclassified_Lactobacillales' & genus!='Anaerococcus_1' &
                    genus!='unclassified_Propionibacteriaceae' & genus!='Pseudomonas_1' & genus!='unclassified_Enterobacterales' &
                    genus!='unclassified_Enterobacteriaceae_1' & genus!='unclassified_Staphylococcales' & genus!='Acinetobacter' & genus!='unclassified_Staphylococcaceae' &
                    genus!='unclassified_Capnodiales' & genus!='unclassified_Intrasporangiaceae_1' & genus!='Campylobacter' &
                    genus!='Megasphaera_1' & genus!='Lautropia' & genus!='Prevotella_6' & genus!='Prevotella_4' & genus!='Helcococcus' & genus!='Prevotella_5' & genus!='Empedobacter_1' &
                    genus!='Murdochiella' & genus!='Prevotella_9' & genus!="Cutibacterium_2" & genus!='Prevotella_2' & genus!='Finegoldia_1'
                  & genus!='Corynebacteriaceae' & genus!='Afipia_2' )

# Export the table outside of R
write.table(asv_tax, "ASVs-taxonomy.tsv", sep = "\t", quote = FALSE, col.names = NA)

#Create a count table only for sulfur utilizing clades
ASV_proteobacteria<-subset(tax_tab, phylum=="Proteobacteria"
                           | phylum=="Bacteroidota"
                           | phylum=="Chloroflexi"
                           | phylum=="Dadabacteria"
                           | phylum=="Desulfobacterota_1"
                           | phylum=="Desulfobacterota_2"
                           | phylum=="Nitrospirota"
                           | phylum=="Halobacterota"
                           | phylum=="Actinobacteriota"
                           | phylum=="Planctomycetota"
                           | phylum=="Crenarchaeota"
                           | phylum=="Euryarchaeota")
proteobacteria_counts_tab <- asv_tab[asv_tab$ASV %in% ASV_proteobacteria$ASV, ]

#Create a count table only for Firmicutes
ASV_proteobacteria<-subset(tax_tab, phylum=="Proteobacteria"
                           | phylum=="Bacteroidota"
                           | phylum=="Chloroflexi"
                           | phylum=="Dadabacteria"
                           | phylum=="Desulfobacterota_1"
                           | phylum=="Desulfobacterota_2"
                           | phylum=="Nitrospirota"
                           | phylum=="Halobacterota"
                           | phylum=="Actinobacteriota"
                           | phylum=="Planctomycetota"
                           | phylum=="Crenarchaeota"
                           | phylum=="Euryarchaeota")
proteobacteria_counts_tab <- asv_tab[asv_tab$ASV %in% ASV_proteobacteria$ASV, ]
# Export the table outside of R
write.table(proteobacteria_counts_tab, "ASVs-counts_proteobacteria.tsv", sep = "\t", quote = FALSE, col.names = NA)

# Confirm that the table file was created
list.files()

head(counts_tab)
head(tax_tab)
head(asv_tab)

