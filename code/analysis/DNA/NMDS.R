#### VEGAN NMDS TEST
#### https://www.rpubs.com/RGrieger/545184
```{r setup, include=FALSE}
### RUN THIS LINE IN THE CONSOLE AGAIN:
setwd(paste0("[your path]/data/processed/16S_rRNA_seqs")) 
```
```{r eval=F}
#load libraries
library(vegan)
library(tidyr)
library(dplyr)
library(ggplot2)
```
```{r eval=F}
We will start by loading the tables of general and taxonomy level-specific ASV normalized counts. 
#MUST USE RELATIVE ABUNDANCES OR NORMALIZED COUNTS, not raw counts
counts_tab<- read.csv("ASVs-counts_norm.tsv",header=T, row.names=1,check.names=F, sep="\t")
colnames(counts_tab)
head(counts_tab)
counts_tab_phylum<- read.table("SR2113_phylum.tsv",header=T, row.names=1,check.names=F, sep="\t")
counts_tab_order<- read.table("SR2113_order.tsv",header=T, row.names=1,check.names=F, sep="\t")
counts_tab_genus<- read.table("SR2113_genus.tsv",header=T, row.names=1,check.names=F, sep="\t")
head(counts_tab_phylum)

#Make replacements in count table
new_colnames <- c("GC04_136",
                  "GC06_111",
                  "GC06_81",
                  "GC06_156",
                  "GC04_156",
                  "GC06_166",
                  "GC06_171",
                  "GC06_146",
                  "GC06_71",
                  "GC06_121",
                  "GC06_41",
                  "GC02_63.5",
                  "GC02_90.5",
                  "GC04_206",
                  "GC04_171",
                  "GC04_191",
                  "GC06_186",
                  "GC06_176",
                  "MC01_23",
                  "MC01_18",
                  "MC01_12",
                  "MC01_8",
                  "MC01_2",
                  "GC06_101",
                  "MC01_4",
                  "GC06_51",
                  "GC06_91",
                  "GC06_21",
                  "GC02_114.5",
                  "GC04_201",
                  "GC02_99.5",
                  "GC06_131",
                  "GC06_61",
                  "GC02_84.5",
                  "GC02_105.5",
                  "GC02_126.5",
                  "GC04_86",
                  "GC02_111.5",
                  "GC02_45.5",
                  "GC02_57.5",
                  "GC02_72.5",
                  "GC06_31",
                  "GC02_81.5",
                  "GC02_18.5",
                  "GC02_51.5",
                  "GC06_181",
                  "GC02_6.5",
                  "GC02_39.5",
                  "GC04_186",
                  "GC04_166",
                  "GC04_196",
                  "GC06_16",
                  "GC02_21.5",
                  "GC06_161",
                  "GC04_146",
                  "GC04_71",
                  "GC06_141",
                  "GC04_106",
                  "GC04_96",
                  "GC06_126",
                  "GC04_26",
                  "GC02_15.5",
                  "MC01_27",
                  "GC02_60.5",
                  "GC02_93.5",
                  "GC04_176",
                  "GC04_16",
                  "GC02_117.5",
                  "MC01_31",
                  "MC01_36",
                  "GC06_56",
                  "GC04_76",
                  "GC02_36.5",
                  "GC02_96.5",
                  "GC04_41",
                  "GC06_11",
                  "GC02_87.5",
                  "GC04_51",
                  "GC04_46",
                  "GC04_61",
                  "GC02_9.5",
                  "GC04_66",
                  "GC04_11",
                  "GC04_126",
                  "GC02_12.5",
                  "GC04_31",
                  "GC04_116")

colnames (counts_tab) <- new_colnames
colnames (counts_tab_phylum) <- c("Taxa", new_colnames)
colnames (counts_tab_order) <- c("Taxa", new_colnames)
colnames (counts_tab_genus) <- c("Taxa", new_colnames)

#Dropping the blank samples
#counts_tab <- counts_tab [!grepl("np|ne", colnames(counts_tab))]
colnames(counts_tab)

# Create relative abundance dataframes 
rownames (counts_tab_phylum) <- counts_tab_phylum$Taxa
counts_tab_phylum$Taxa<- NULL
relative_tab_phylum <- apply(counts_tab_phylum, 2, function(x) 100*(x/sum(x)))
rownames (counts_tab_order) <- counts_tab_order$Taxa
counts_tab_order$Taxa<- NULL
relative_tab_order <- apply(counts_tab_order, 2, function(x) 100*(x/sum(x)))
rownames (counts_tab_genus) <- counts_tab_genus$Taxa
counts_tab_genus$Taxa<- NULL
relative_tab_genus <- apply(counts_tab_genus, 2, function(x) 100*(x/sum(x)))

#Remove samples from the counts tab that have less than 600 AVSs 
#Filter out depths with less than 600 total ASVs
depths_to_keep <- subset(df_sum, as.numeric(Total) <= 900)
filtered_absolute_asvs <- phylum_tab_for_plot.g2
#filtered_absolute_asvs <- filter(phylum_tab_for_plot.g2,
#                                      phylum_tab_for_plot.g2$Depth %in% depths_to_keep$Depth)

#Drop unnecesary columns
filtered_absolute_asvs <- subset(filtered_absolute_asvs, select = -c(Core, Depth, core_colors, Depth_bin))

# now we'll transform the table with AVSs only for depths with more than 600 total, to broad format
filtered_absolute_asvs_for_manip <- 
  filtered_absolute_asvs %>% 
  pivot_wider(names_from = Sample, values_from = Proportion) 

colnames(filtered_absolute_asvs_for_manip)

counts_tab_filtered<-counts_tab[,colnames(counts_tab) %in% colnames(filtered_absolute_asvs_for_manip)]
colnames(counts_tab_filtered)

#sample_info_tab is created from a metadata file
sample_info_tab <- read.table("SR2113_metadata_final.txt", sep = "\t", header = TRUE)
#create dataframe of ASV's and sample names in order needed for vegan package
data1 <- data.frame(t(counts_tab_filtered[])) #96 rows
datap <- data.frame(t(relative_tab_phylum[]))
datao <- data.frame(t(relative_tab_order[]))
datag <- data.frame(t(relative_tab_genus[]))
#data1 <- setNames(data.frame(t(counts_tab_filtered[,-1])), counts_tab_filtered[,1])
data1[data1 < 0] <- 0 
datap[datap < 0] <- 0 
datao[datao < 0] <- 0 
datag[datag < 0] <- 0 
#convert sample info / environmental variables into format needed for vegan
#samp1 <- sample_info_tab[,-1]
samp1 <- sample_info_tab
rownames(samp1) <- sample_info_tab[,1] #ID is in column N 1
#Drop values in the samp1 df that are not in the data1 df
samp2<-samp1[rownames(samp1) %in% rownames(data1) ,]
data2<-data1[rownames(data1) %in% rownames(samp2) ,]
datap2<-datap[rownames(datap) %in% rownames(samp2) ,]
datao2<-datao[rownames(datao) %in% rownames(samp2) ,]
datag2<-datag[rownames(datag) %in% rownames(samp2) ,]
#remove unneeded columns for ordination fitting
#THIS WILL BE UNIQUE TO YOUR DATA SET
#I am removing ALL columns except for the environmental data I want to test to see if it correlates with the NMDS plot
samp2 <- samp2 %>% select(-c(Core, Silica))
#change values to numeric value type
samp2$Calcium = as.numeric(samp2$Calcium)
samp2$Carbon_in_OM = as.numeric(samp2$Carbon_in_OM)
samp2$Iron = as.numeric(samp2$Iron)
samp2$Magnesium = as.numeric(samp2$Magnesium)
samp2$Manganese = as.numeric(samp2$Manganese)
samp2$Nitrate = as.numeric(samp2$Nitrate)
samp2$O_in_sulfate = as.numeric(samp2$O_in_sulfate)
samp2$S_in_Sulfate = as.numeric(samp2$S_in_Sulfate)
samp2$Strontium = as.numeric(samp2$Strontium)
samp2$Sulfate = as.numeric(samp2$Sulfate)
samp2$Sulfide = as.numeric(samp2$Sulfide)
samp2$Sulfite = as.numeric(samp2$Sulfite)
samp2$Thiosulfate = as.numeric(samp2$Thiosulfate)
samp2$TOC_percent = as.numeric(samp2$TOC_percent)
#create ordination
#CHECK THAT STRESS VALUES ARE GOOD
stations.mds <- metaMDS(data2, distance = "bray", autotransform = FALSE)
#create environmental vectors fit
stations.envfit <- envfit(stations.mds, samp2, permutations = 999, na.rm = TRUE)
#create species fit - will run for a long time
#will tell you what species are driving the changes seen, have not fully tested
stations.spp.fit <- envfit(stations.mds, data1, permutations = 999)
names(stations.mds)
#Extract sample and variable (species) scores
variableScores <- stations.mds$species
sampleScores <- stations.mds$points
stations.mds$stress
```

```{r eval=F}
#create ordination for taxonomic categories
#CHECK THAT STRESS VALUES ARE GOOD
stations.mds.p <- metaMDS(datap2, distance = "bray", autotransform = FALSE)
stations.mds.o <- metaMDS(datao2, distance = "bray", autotransform = FALSE)
stations.mds.g <- metaMDS(datag2, distance = "bray", autotransform = FALSE)
#create environmental vectors fit
stations.envfit.p <- envfit(stations.mds.p, samp2, permutations = 999, na.rm = TRUE)
stations.envfit.o <- envfit(stations.mds.o, samp2, permutations = 999, na.rm = TRUE)
stations.envfit.g <- envfit(stations.mds.g, samp2, permutations = 999, na.rm = TRUE)
#create species fit 
names(stations.mds.p)
names(stations.mds.o)
names(stations.mds.g)
#Extract sample and variable (species) scores
variableScores.p <- stations.mds.p$species
variableScores.o <- stations.mds.o$species
variableScores.g <- stations.mds.g$species

sampleScores.p <- stations.mds.p$points
sampleScores.o <- stations.mds.o$points
sampleScores.g <- stations.mds.g$points

stations.mds.p$stress
stations.mds.o$stress
stations.mds.g$stress
```
#### PLOT VEGAN NMDS
```{r eval=F}
#basic plot, this is a version of an NMDS, but probably not the one you want to publish as it shows all the ASVs and makes it look messy
#red dots are ASVs
#black and white circles are sites
#RUN THIS WHOLE CODE CHUNK TOGETHER
plot(stations.mds)
points(stations.mds, display = "species", col = "red")
points(stations.mds, display = "sites", col = "black")
```
#To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. 
# We will do this by calling the scores of the MDS.
```{r eval=F}
site.scrs <- as.data.frame(scores(stations.mds, display = "sites")) #save NMDS results into dataframe.
site.scrs <- cbind(site.scrs, Station = samp2$Station) #add grouping variable "station" to dataframe
site.scrs <- cbind(site.scrs, Calcium =samp2$Calcium)
site.scrs <- cbind(site.scrs, Carbon_in_OM = samp2$Carbon_in_OM)
site.scrs <- cbind(site.scrs, Iron = samp2$Iron)
site.scrs <- cbind(site.scrs, Magnesium = samp2$Magnesium)
site.scrs <- cbind(site.scrs, Manganese = samp2$Manganese)
site.scrs <- cbind(site.scrs, Nitrate = samp2$Nitrate)
site.scrs <- cbind(site.scrs, O_in_sulfate = samp2$O_in_sulfate)
site.scrs <- cbind(site.scrs, S_in_Sulfate =samp2$S_in_Sulfate)
site.scrs <- cbind(site.scrs, Strontium = samp2$Strontium)
site.scrs <- cbind(site.scrs, Sulfate = samp2$Sulfate)
site.scrs <- cbind(site.scrs, Sulfide = samp2$Sulfide)
site.scrs <- cbind(site.scrs, Sulfite = samp2$Sulfite)
site.scrs <- cbind(site.scrs, Thiosulfate = samp2$Thiosulfate)
site.scrs <- cbind(site.scrs, TOC_percent = samp2$TOC_percent)
site.scrs <- cbind(site.scrs, Real_depth_cm = samp2$Real_depth_cm)
site.scrs <- cbind(site.scrs, ID = samp2$ID)
site.scrs$ID <-gsub("^.*_", "", site.scrs$ID)

site.scrs
```

```{r eval=F}
#For phyla
site.scrs.p <- as.data.frame(scores(stations.mds.p, display = "sites")) #save NMDS results into dataframe.
site.scrs.p <- cbind(site.scrs.p, Station = samp2$Station) #add grouping variable "station" to dataframe
site.scrs.p <- cbind(site.scrs.p, Calcium =samp2$Calcium)
site.scrs.p <- cbind(site.scrs.p, Carbon_in_OM = samp2$Carbon_in_OM)
site.scrs.p <- cbind(site.scrs.p, Iron = samp2$Iron)
site.scrs.p <- cbind(site.scrs.p, Magnesium = samp2$Magnesium)
site.scrs.p <- cbind(site.scrs.p, Manganese = samp2$Manganese)
site.scrs.p <- cbind(site.scrs.p, Nitrate = samp2$Nitrate)
site.scrs.p <- cbind(site.scrs.p, O_in_sulfate = samp2$O_in_sulfate)
site.scrs.p <- cbind(site.scrs.p, S_in_Sulfate =samp2$S_in_Sulfate)
site.scrs.p <- cbind(site.scrs.p, Strontium = samp2$Strontium)
site.scrs.p <- cbind(site.scrs.p, Sulfate = samp2$Sulfate)
site.scrs.p <- cbind(site.scrs.p, Sulfide = samp2$Sulfide)
site.scrs.p <- cbind(site.scrs.p, Sulfite = samp2$Sulfite)
site.scrs.p <- cbind(site.scrs.p, Thiosulfate = samp2$Thiosulfate)
site.scrs.p <- cbind(site.scrs.p, TOC_percent = samp2$TOC_percent)
site.scrs.p <- cbind(site.scrs.p, Real_depth_cm = samp2$Real_depth_cm)
site.scrs.p <- cbind(site.scrs.p, ID = samp2$ID)
site.scrs.p$ID <-gsub("^.*_", "", site.scrs.p$ID)

site.scrs.p
```
```{r eval=F}
#For order
site.scrs.o <- as.data.frame(scores(stations.mds.o, display = "sites")) #save NMDS results into dataframe.
site.scrs.o <- cbind(site.scrs.o, Station = samp2$Station) #add grouping variable "station" to dataframe
site.scrs.o <- cbind(site.scrs.o, Calcium =samp2$Calcium)
site.scrs.o <- cbind(site.scrs.o, Carbon_in_OM = samp2$Carbon_in_OM)
site.scrs.o <- cbind(site.scrs.o, Iron = samp2$Iron)
site.scrs.o <- cbind(site.scrs.o, Magnesium = samp2$Magnesium)
site.scrs.o <- cbind(site.scrs.o, Manganese = samp2$Manganese)
site.scrs.o <- cbind(site.scrs.o, Nitrate = samp2$Nitrate)
site.scrs.o <- cbind(site.scrs.o, O_in_sulfate = samp2$O_in_sulfate)
site.scrs.o <- cbind(site.scrs.o, S_in_Sulfate =samp2$S_in_Sulfate)
site.scrs.o <- cbind(site.scrs.o, Strontium = samp2$Strontium)
site.scrs.o <- cbind(site.scrs.o, Sulfate = samp2$Sulfate)
site.scrs.o <- cbind(site.scrs.o, Sulfide = samp2$Sulfide)
site.scrs.o <- cbind(site.scrs.o, Sulfite = samp2$Sulfite)
site.scrs.o <- cbind(site.scrs.o, Thiosulfate = samp2$Thiosulfate)
site.scrs.o <- cbind(site.scrs.o, TOC_percent = samp2$TOC_percent)
site.scrs.o <- cbind(site.scrs.o, Real_depth_cm = samp2$Real_depth_cm)
site.scrs.o <- cbind(site.scrs.o, ID = samp2$ID)
site.scrs.o$ID <-gsub("^.*_", "", site.scrs.o$ID)

site.scrs.o
```
```{r eval=F}
#For genus
site.scrs.g <- as.data.frame(scores(stations.mds.g, display = "sites")) #save NMDS results into dataframe.
site.scrs.g <- cbind(site.scrs.g, Station = samp2$Station) #add grouping variable "station" to dataframe
site.scrs.g <- cbind(site.scrs.g, Calcium =samp2$Calcium)
site.scrs.g <- cbind(site.scrs.g, Carbon_in_OM = samp2$Carbon_in_OM)
site.scrs.g <- cbind(site.scrs.g, Iron = samp2$Iron)
site.scrs.g <- cbind(site.scrs.g, Magnesium = samp2$Magnesium)
site.scrs.g <- cbind(site.scrs.g, Manganese = samp2$Manganese)
site.scrs.g <- cbind(site.scrs.g, Nitrate = samp2$Nitrate)
site.scrs.g <- cbind(site.scrs.g, O_in_sulfate = samp2$O_in_sulfate)
site.scrs.g <- cbind(site.scrs.g, S_in_Sulfate =samp2$S_in_Sulfate)
site.scrs.g <- cbind(site.scrs.g, Strontium = samp2$Strontium)
site.scrs.g <- cbind(site.scrs.g, Sulfate = samp2$Sulfate)
site.scrs.g <- cbind(site.scrs.g, Sulfide = samp2$Sulfide)
site.scrs.g <- cbind(site.scrs.g, Sulfite = samp2$Sulfite)
site.scrs.g <- cbind(site.scrs.g, Thiosulfate = samp2$Thiosulfate)
site.scrs.g <- cbind(site.scrs.g, TOC_percent = samp2$TOC_percent)
site.scrs.g <- cbind(site.scrs.g, Real_depth_cm = samp2$Real_depth_cm)
site.scrs.g <- cbind(site.scrs.g, ID = samp2$ID)
site.scrs.g$ID <-gsub("^.*_", "", site.scrs.g$ID)

site.scrs.g
```

To show environmental extrinsic variables another datasheet needs to be created
```{r eval=F}
env.scores.stations <- as.data.frame(scores(stations.envfit, display = "vectors")) #extracts relevant scores from envifit
env.scores.stations <- cbind(env.scores.stations, env.variables = rownames(env.scores.stations)) #and then gives them their names
env.scores.stations <- cbind(env.scores.stations, pval = stations.envfit$vectors$pvals) # add pvalues to dataframe
#sig.env.scrs <- subset(env.scores.stations, pval<=0.05) #subset data to show variables significant at 0.05
#you can look at the table to see the p values associated with each environmental data. The significant ones are saved separated in sig.env.scrs and are used below in the next code chunk.
include_list <- c("Carbon_in_OM", "Sulfide","Sulfate","Nitrate","Manganese","O_in_Sulfate")
sig.env.scrs <-subset(env.scores.stations, rownames(env.scores.stations) %in% include_list)
env.scores.stations
write.table(env.scores.stations,file="MDS_env_vars.txt")

```
```{r eval=F}
#For Phylum
env.scores.stations.p <- as.data.frame(scores(stations.envfit.p, display = "vectors")) #extracts relevant scores from envifit
env.scores.stations.p <- cbind(env.scores.stations.p, env.variables = rownames(env.scores.stations.p)) #and then gives them their names
env.scores.stations.p <- cbind(env.scores.stations.p, pval = stations.envfit.p$vectors$pvals) # add pvalues to dataframe
#sig.env.scrs <- subset(env.scores.stations, pval<=0.05) #subset data to show variables significant at 0.05
#you can look at the table to see the p values associated with each environmental data. The significant ones are saved separated in sig.env.scrs and are used below in the next code chunk.
include_list <- c("Carbon_in_OM", "Sulfide","Sulfate","Nitrate","Manganese","O_in_Sulfate")
sig.env.scrs.p <-subset(env.scores.stations.p, rownames(env.scores.stations.p) %in% include_list)
env.scores.stations.p
write.table(env.scores.stations.p,file="MDS_env_vars_p.txt")

```
```{r eval=F}
#For Order
env.scores.stations.o <- as.data.frame(scores(stations.envfit.o, display = "vectors")) #extracts relevant scores from envifit
env.scores.stations.o <- cbind(env.scores.stations.o, env.variables = rownames(env.scores.stations.o)) #and then gives them their names
env.scores.stations.o <- cbind(env.scores.stations.o, pval = stations.envfit.o$vectors$pvals) # add pvalues to dataframe
#sig.env.scrs <- subset(env.scores.stations, pval<=0.05) #subset data to show variables significant at 0.05
#you can look at the table to see the p values associated with each environmental data. The significant ones are saved separated in sig.env.scrs and are used below in the next code chunk.
include_list <- c("Carbon_in_OM", "Sulfide","Sulfate","Nitrate","Manganese","O_in_Sulfate")
sig.env.scrs.o <-subset(env.scores.stations.o, rownames(env.scores.stations.o) %in% include_list)
env.scores.stations.o
write.table(env.scores.stations.o,file="MDS_env_vars_o.txt")
```
```{r eval=F}
#For Genus
env.scores.stations.g <- as.data.frame(scores(stations.envfit.g, display = "vectors")) #extracts relevant scores from envifit
env.scores.stations.g <- cbind(env.scores.stations.g, env.variables = rownames(env.scores.stations.g)) #and then gives them their names
env.scores.stations.g <- cbind(env.scores.stations.g, pval = stations.envfit.g$vectors$pvals) # add pvalues to dataframe
#sig.env.scrs <- subset(env.scores.stations, pval<=0.05) #subset data to show variables significant at 0.05
#you can look at the table to see the p values associated with each environmental data. The significant ones are saved separated in sig.env.scrs and are used below in the next code chunk.
include_list <- c("Carbon_in_OM", "Sulfide","Sulfate","Nitrate","Manganese","O_in_Sulfate")
sig.env.scrs.g <-subset(env.scores.stations.g, rownames(env.scores.stations.g) %in% include_list)
env.scores.stations.g
write.table(env.scores.stations.g,file="MDS_env_vars_g.txt")
```
#### PLOT
```{r eval=F}
only_samples <- site.scrs
#filter(locale != "neg")

#MDS plot for samples by core
nmds.plot.stations <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, shape = Station, color = Station), size = 2)+ #adds site points to plot
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Station", shape = "Station")+ # add legend labels 
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot
nmds.plot.stations + labs(title = "Basic ordination plot") #displays plot

#Declare plot for MDS of samples with arrows and environmental variables
nmds.plot.stations_only_samples <- ggplot(only_samples, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1,NMDS2, shape = Station,color = as.numeric(Real_depth_cm)),size = 2)+ 
  coord_fixed()+
  #scale_color_viridis_d() +
  scale_colour_gradient(high = "#132B43", low = "#56B1F7") +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", linetype = "solid"))+
  #=textNudge <- 1.2
  #text(variableScores[, 1]*textNudge, variableScores[, 2]*textNudge, rownames(variableScores), #cex=0.7)
  geom_segment(data = sig.env.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.2, "cm")), colour = "grey10", lwd=0.3) +
  ggrepel::geom_text_repel(data = sig.env.scrs, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 3, direction = "both", segment.size = 1)+ #add labels for env variables
  ggtitle(paste0('stress: ', round(stations.mds$stress,3)))+
  geom_text(aes(label = ID), hjust = -0.2) #add labels to points

# labs(colour = "locale", shape = "station")+ # add legend labels for Management and Landuse
theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot
#scale_shape_manual(values=c(15, 16, 17, 8))

#Display plot with arrows and environmental variables
nmds.plot.stations_only_samples 

```
```{r eval=F}
#Plot for Phylum
only_samples.p <- site.scrs.p
#filter(locale != "neg")

#Declare plot for MDS of samples with arrows and environmental variables
nmds.plot.stations_only_samples.p <- ggplot(only_samples.p, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1,NMDS2, shape = Station,color = as.numeric(Real_depth_cm)),size = 2)+ 
  coord_fixed()+
  #scale_color_viridis_d() +
  scale_colour_gradient(high = "#132B43", low = "#56B1F7") +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", linetype = "solid"))+
  #=textNudge <- 1.2
  #text(variableScores[, 1]*textNudge, variableScores[, 2]*textNudge, rownames(variableScores), #cex=0.7)
  geom_segment(data = sig.env.scrs.p, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.2, "cm")), colour = "grey10", lwd=0.3) +
  ggrepel::geom_text_repel(data = sig.env.scrs.p, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 3, direction = "both", segment.size = 1)+ #add labels for env variables
  ggtitle(paste0('stress: ', round(stations.mds.p$stress,3)))+
  geom_text(aes(label = ID), hjust = -0.2) #add labels to points

# labs(colour = "locale", shape = "station")+ # add legend labels 
theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot
#scale_shape_manual(values=c(15, 16, 17, 8))

#Display plot with arrows and environmental variables
nmds.plot.stations_only_samples.p 

```

```{r eval=F}
#Plot for Order
only_samples.o <- site.scrs.o
#filter(locale != "neg")

#Declare plot for MDS of samples with arrows and environmental variables
nmds.plot.stations_only_samples.o <- ggplot(only_samples.o, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1,NMDS2, shape = Station,color = as.numeric(Real_depth_cm)),size = 2)+ 
  coord_fixed()+
  #scale_color_viridis_d() +
  scale_colour_gradient(high = "#132B43", low = "#56B1F7") +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", linetype = "solid"))+
  #=textNudge <- 1.2
  #text(variableScores[, 1]*textNudge, variableScores[, 2]*textNudge, rownames(variableScores), #cex=0.7)
  geom_segment(data = sig.env.scrs.o, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.2, "cm")), colour = "grey10", lwd=0.3) +
  ggrepel::geom_text_repel(data = sig.env.scrs.o, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 3, direction = "both", segment.size = 1)+ #add labels for env variables
  ggtitle(paste0('stress: ', round(stations.mds.o$stress,3)))+
  geom_text(aes(label = ID), hjust = -0.2) #add labels to points

# labs(colour = "locale", shape = "station")+ # add legend labels for Management and Landuse
theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot
#scale_shape_manual(values=c(15, 16, 17, 8))

#Display plot with arrows and environmental variables
nmds.plot.stations_only_samples.o 

```

```{r eval=F}
#Plot for Genus
only_samples.g <- site.scrs.g
#filter(locale != "neg")

#Declare plot for MDS of samples with arrows and environmental variables
nmds.plot.stations_only_samples.g <- ggplot(only_samples.g, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1,NMDS2, shape = Station,color = as.numeric(Real_depth_cm)),size = 2)+ 
  coord_fixed()+
  #scale_color_viridis_d() +
  scale_colour_gradient(high = "#132B43", low = "#56B1F7") +
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", linetype = "solid"))+
  #=textNudge <- 1.2
  #text(variableScores[, 1]*textNudge, variableScores[, 2]*textNudge, rownames(variableScores), #cex=0.7)
  geom_segment(data = sig.env.scrs.g, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.2, "cm")), colour = "grey10", lwd=0.3) +
  ggrepel::geom_text_repel(data = sig.env.scrs.g, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 3, direction = "both", segment.size = 1)+ #add labels for env variables
  ggtitle(paste0('stress: ', round(stations.mds.g$stress,3)))+
  geom_text(aes(label = ID), hjust = -0.2) #add labels to points

# labs(colour = "locale", shape = "station")+ # add legend labels for Management and Landuse
#theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot
#scale_shape_manual(values=c(15, 16, 17, 8))

#Display plot with arrows and environmental variables
nmds.plot.stations_only_samples.g 

```
#CONTOUR LINES ON NMDS
```{r eval=F}
#shows contour lines of environmental data and how it is associated with the plotted NMDS points
plot(stations.mds, type = "n") 
points(stations.mds, display = "sites") 
ordisurf(stations.mds,samp2$Sulfide, add = TRUE)
plot(stations.mds, type = "n") 
points(stations.mds, display = "sites") 
ordisurf(stations.mds,samp2$Carbon_in_OM, add = TRUE)
plot(stations.mds, type = "n") 
points(stations.mds, display = "sites")
ordisurf(stations.mds,samp2$Manganese, add = TRUE)
plot(stations.mds, type = "n") 
points(stations.mds, display = "sites")
ordisurf(stations.mds,samp2$Thiosulfate, add = TRUE)
plot(stations.mds, type = "n") 
points(stations.mds, display = "sites")
```
```{r eval=F}
# function for ellipsess 
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

library(vegan)#data for ellipse, in this case using the station factor
df_ell.station <- data.frame() #sets up a data frame before running the function.
for(g in levels(as.factor(site.scrs$Station))){
  df_ell.station <- rbind(df_ell.station, cbind(as.data.frame(with(site.scrs [site.scrs$Station==g,],
                                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,Station=g))
}

# data for labelling the ellipse
NMDS.mean.station=aggregate(site.scrs[ ,c("NMDS1", "NMDS2")], 
                            list(group = site.scrs$Station), mean)

# data for labelling the ellipse
NMDS.mean=aggregate(site.scrs[,c("NMDS1", "NMDS2")], 
                    list(group = site.scrs$Station), mean)
nmds.plot.stations+ 
  geom_path(data = df_ell.station, aes(x = NMDS1, y = NMDS2, group = Station)) #this is the ellipse, separate ones by Site.