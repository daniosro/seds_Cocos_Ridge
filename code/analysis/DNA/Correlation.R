---
  title: "Correlation"
output: html_document
---
  
  #### Correlation maps
  ```{r setup, include=FALSE}
### RUN THIS LINE IN THE CONSOLE AGAIN:
setwd(paste0("[your path]/data/processed/16S_rRNA_seqs")) 
setwd(paste0("/Users/dosorior/OneDrive - University of Southern California/git/seds_Cocos_Ridge/data/processed/16S_rRNA_seqs"))
```
```{r eval=F}
#load libraries
library(vegan)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
install.packages("Hmisc")
library(Hmisc)
install.packages("PerformanceAnalytics", dependencies = TRUE)
library(PerformanceAnalytics)
```
```{r eval=F}

#Load the taxonomy table
phylum_tab<- read.csv("SR2113_phylum.tsv",header=T, row.names=1,check.names=F, sep="\t")
head(phylum_tab)

#Load the table with sample information
mod_sample_info_tab <- read.table("mod_sample_info_tab.tsv",header=T, row.names=1,check.names=F, sep="\t")
mod_sample_info_tab$Sample <-gsub("^[^_]*_", "", mod_sample_info_tab$Sample)

#Make replacements in count table
colnames (phylum_tab) <- gsub("\\.", "-", colnames(phylum_tab))

#Create a copy of the table for manipulation
phylum_tab_for_plot <- data.frame(phylum_tab)
head(phylum_tab_for_plot)
# now we'll transform the table into narrow, or long, format
phylum_tab_for_plot.g <- 
  phylum_tab_for_plot %>% 
  pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
#Add depth and core to table
phylum_tab_for_plot.g2 <- 
  phylum_tab_for_plot.g %>% left_join(mod_sample_info_tab)
head(phylum_tab_for_plot.g2)
#Replace dots by dashes in the samples names
phylum_tab_for_plot.g2$Sample <- gsub("^[^_]*_", "", phylum_tab_for_plot.g2$Sample)
head(phylum_tab_for_plot.g2)
#Get total Phylum ASVs for depth
df_sum <-
  phylum_tab_for_plot.g2 %>%
  group_by(Depth, Core) %>%
  summarise(Total = sum(Proportion))


#Filter out depths with less than 600 total ASVs
depths_to_keep <- subset(df_sum, as.numeric(Total) <= 900)
filtered_absolute_asvs <- phylum_tab_for_plot.g2
#filtered_absolute_asvs <- filter(phylum_tab_for_plot.g2,
#                                      phylum_tab_for_plot.g2$Depth %in% depths_to_keep$Depth)

#Drop unnecesary columns
filtered_absolute_asvs <- subset(filtered_absolute_asvs, select = -c(Core, Depth, core_colors, Depth_bin,depth_colors))

# now we'll transform the table with AVSs only for depths with more than 600 total, to broad format
filtered_absolute_asvs_for_manip <- 
  filtered_absolute_asvs %>% 
  pivot_wider(names_from = Sample, values_from = Proportion) 

colnames(filtered_absolute_asvs_for_manip)


#Create abundance table
abund_table <-t(filtered_absolute_asvs_for_manip)
abund_table <- abund_table[grepl("NA", rownames(abund_table))==F,]
#rownames(abund_table) <-gsub("^[^_]*_", "", rownames(abund_table))
rownames(abund_table)
#Set Phylum as row names
colnames(abund_table) <- lapply(abund_table[1, ], as.character)
abund_table <- abund_table[-1,] 

#env_vars_tab is created from a metadata file
meta_table <- read.table("/Users/dosorior/Library/CloudStorage/OneDrive-UniversityofSouthernCalifornia/git/seds_Cocos_Ridge/data/raw/SR2113_metadata_final.txt", sep = "\t", header = TRUE)
meta_table <- read.table("/[PATH]/git/seds_Cocos_Ridge/data/raw/SR2113_metadata_final.txt", sep = "\t", header = TRUE)
rownames(meta_table) <- meta_table$ID
rownames(meta_table)
rownames(abund_table)

#Extract the corresponding meta_table for the samples in abund_table
meta_table<-meta_table[rownames(abund_table),]
meta_table <- meta_table[grepl("NA", rownames(meta_table))==F,]

colnames(meta_table)

#Filter by core
#meta_table_gc02 <- meta_table %>% filter(grepl("GC02", Core))
#head(meta_table_gc02)
#meta_table_mc01 <- meta_table %>% filter(grepl("MC01", Core))
#meta_table_gc04 <- meta_table %>% filter(grepl("GC04", Core))
#meta_table_gc06 <- meta_table %>% filter(grepl("GC06", Core))

head(phylum_tab_for_plot.g2)
```


```{r eval=F}
# Use sel_env to specify the variables you want to use and sel_env_label to specify the labes for the panel
sel_env<-c("Real_depth_cm","Iron","Manganese","Nitrate","Sulfate","Sulfide","Sulfite","Thiosulfate", "Carbon_in_OM", "TOC_percent", "O_in_sulfate", "S_in_Sulfate")
sel_env_label <- list(
  'Real_depth_cm'="Depth (cm)",
  'Iron'="Iron",
  'Manganese'="Manganese",
  'Nitrate'="Nitrate",
  'Sulfate'="Sulfate",
  'Sulfide'='Sulfide',
  'Sulfite'='Sulfite',
  "Thiosulfate"="Thiosulfate",
  #"ASVs"="ASVs",
  "Carbon_in_OM"="d13C OM",
  "TOC_percent" = "% TOC",
  'O_in_sulfate'="d18O Sulfate",
  'S_in_Sulfate'="d34S Sulfate"
)

sel_env_label<-t(as.data.frame(sel_env_label))
sel_env_label<-as.data.frame(sel_env_label)
colnames(sel_env_label)<-c("Trans")
sel_env_label$Trans<-as.character(sel_env_label$Trans)


#Now get a filtered table based on sel_env
meta_table_filtered<-meta_table[,sel_env]
meta_table_filtered <- meta_table_filtered[grepl("^NA", rownames(meta_table_filtered))==F,]
rownames(meta_table_filtered)
abund_table_filtered<-abund_table[rownames(meta_table_filtered),]
rownames(abund_table)
#Convert abundance table to numeric
abund_table_filtered <- apply(abund_table_filtered, 2, as.numeric)
rownames(abund_table_filtered) <- rownames(meta_table_filtered)

#Add core and depth to abundance table

abund_table_cd <- data.frame(abund_table_filtered) 
abund_table_cd$Sample <- rownames(abund_table_cd)
abund_table_cd_1 <- abund_table_cd %>% left_join(mod_sample_info_tab)

#Apply normalisation (either use relative or log-relative transformation)
#x<-abund_table_filtered/rowSums(abund_table_filtered)
x<-log((abund_table_filtered+1)/(rowSums(abund_table_filtered)+dim(abund_table_filtered)[2]))
x<-x[,order(colSums(x),decreasing=TRUE)]
rownames(x)
#Extract list of top N Taxa
#N<-51
#taxa_list<-colnames(x)[1:N]
taxa_list<-colnames(x)
#remove "__Unknown__" and add it to others
#taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
N<-length(taxa_list)
x<-data.frame(x[,colnames(x) %in% taxa_list])
y<-meta_table_filtered
rownames(y)
#Get grouping information
grouping_info<-meta_table$Station
#grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
#grouping_info = grouping_info[-1,]
grouping_info
```

```{r eval=F}
#Create lists of categories

#Merge taxa and env. variables dataframes
#Merge taxa and env. variables dataframes
data_frame_merge <- merge(x, y,
                          by = 'row.names', all = TRUE)
#Creating core and Depth columns in the df 
x_1 <- abs(x)
x_1 <- data.frame(x_1,do.call(rbind,str_split(rownames(x_1),"_")))
# Changing column names
colnames(x_1)[colnames(x_1) == 'X1'] <- 'Core'
colnames(x_1)[colnames(x_1) == 'X2'] <- 'Depth'
x_1$Depth<-as.double(x_1$Depth)
#Add total ASVs
x_1 <- x_1 %>% left_join(df_sum)
colnames(x_1)[colnames(x_1) == 'Total'] <- 'Total_ASVs'

#Create category based on total ASVs
ASV_groups <- c()

for (entry in x_1 %>% pull(Total_ASVs)) {
  
  if (as.numeric(entry) %>% between(0,1000)) {
    
    ASV_groups <- c(ASV_groups, "below 1000")
    
  } else if (as.numeric(entry) %>% between(1001,4999)) {
    
    ASV_groups <- c(ASV_groups, "1000-5000")
    
  } else {
    
    ASV_groups <- c(ASV_groups, "above 5000")
    
  }
  
}
#Append it to table
x_1$ASV_groups <- ASV_groups

#Subset by core
x_1_mc01 <- subset(x_1, grepl('^MC01', rownames(x)))
x_1_gc02 <- subset(x_1, grepl('^GC02', rownames(x)))
x_1_gc04 <- subset(x_1, grepl('^GC04', rownames(x)))
x_1_gc06 <- subset(x_1, grepl('^GC06', rownames(x)))

#Order
x_1 <- x_1[order(x_1$Core,x_1$Depth ),]


#Get absolute value and convert to matrix
x.matrix<-as.matrix(abs(x))
x_mc01.matrix<-as.matrix(abs(x_mc01))
x_gc02.matrix<-as.matrix(abs(x_gc02))
x_gc04.matrix<-as.matrix(abs(x_gc04))
x_gc06.matrix<-as.matrix(abs(x_gc06))

# Create distance matrices
x.dist<-vegdist(x.matrix, method='bray')
x_mc01.dist<-vegdist(x_mc01.matrix, method='bray')
x_gc02.dist<-vegdist(x_gc02.matrix, method='bray')
x_gc04.dist<-vegdist(x_gc04.matrix, method='bray')
x_gc06.dist<-vegdist(x_gc06.matrix, method='bray')

#Perform permANOVA to determine if there are differences in taxonomy (phylum) between cores
x.div<-adonis2(x.dist~Core, data=x_1, permutations = 999, method="bray")
x.div

#Define pairwise adonis function
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni',reduce=NULL,perm=999)
{
  
  co <- combn(unique(as.character(factors)),2)
  pairs <- c()
  Df <- c()
  SumsOfSqs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  
  for(elem in 1:ncol(co)){
    if(inherits(x, 'dist')){
      x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                      factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
    }
    
    else  (
      if (sim.function == 'daisy'){
        x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
      }
      else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    )
    
    x2 = data.frame(Fac = factors[factors %in% c(co[1,elem],co[2,elem])])
    
    ad <- adonis2(x1 ~ Fac, data = x2,
                  permutations = perm);
    pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    Df <- c(Df,ad$Df[1])
    SumsOfSqs <- c(SumsOfSqs,ad$SumOfSqs[1])
    F.Model <- c(F.Model,ad$F[1]);
    R2 <- c(R2,ad$R2[1]);
    p.value <- c(p.value,ad$`Pr(>F)`[1])
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
  
  if(!is.null(reduce)){
    pairw.res <- subset (pairw.res, grepl(reduce,pairs))
    pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
    
    sig = c(rep('',length(pairw.res$p.adjusted)))
    sig[pairw.res$p.adjusted <= 0.1] <-'.'
    sig[pairw.res$p.adjusted <= 0.05] <-'*'
    sig[pairw.res$p.adjusted <= 0.01] <-'**'
    sig[pairw.res$p.adjusted <= 0.001] <-'***'
    pairw.res <- data.frame(pairw.res[,1:7],sig)
  }
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
}

### Method summary
summary.pwadonis = function(object, ...) {
  cat("Result of pairwise.adonis:\n")
  cat("\n")
  print(object, ...)
  cat("\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}

#Run pairwise permANOVA
b <- pairwise.adonis(x, x_1$ASV_groups, sim.method = "euclidean",
                     p.adjust.m = "bonferroni")
b

#There does not seem to be a difference in taxonomic groups between samples 
#with different numbers of ASVs, so there's no need to filter out samples with low number of ASVs.
#Now we will make an example plot of correlation betwwen a particular group and an environmental variable.

#Filter by cores
mc01_phyla_df <- dplyr::filter(data_frame_merge, grepl('MC01', Row.names))
gc02_phyla_df <- dplyr::filter(data_frame_merge, grepl('GC02', Row.names))
gc04_phyla_df <- dplyr::filter(data_frame_merge, grepl('GC04', Row.names))
gc06_phyla_df <- dplyr::filter(data_frame_merge, grepl('GC06', Row.names))
head(gc02_phyla_df)

`%nin%` = Negate(`%in%`)
updated_proteo_gc02 <- subset(gc02_phyla_df, Row.names!= "GC02_6.5")
updated_proteo_gc04 <- subset(gc04_phyla_df, Row.names %nin% c("GC04_11", "GC04_6.5"))
#Make scatter plots to visualize correlations
p <- ggplot(updated_proteo_gc04) +
  aes(x = Proteobacteria, y = Manganese) +
  geom_point(colour = "#0c4c8a") +
  geom_text(aes(label=Row.names), size=3) +
  theme_minimal()
print (p)

```

```{r eval=F}
#Let us group on cores
#groups<-grouping_info[,1]
groups<-grouping_info

#You can use kendall, spearman, or pearson below:
method<-"kendall"

#We will remove columns we will not use from y because they don't have enough data for all the cores
y_sub<- subset(y, select = -c(Carbon_in_OM,TOC_percent,O_in_sulfate,S_in_Sulfate))
#Now we will drop NA   
y_sub <- drop_na(y_sub)
head(y_sub)

#Now calculate the correlation between individual Taxa and the environmental data
df<-NULL
for(i in colnames(x)){
  for(j in colnames(y_sub)){
    for(k in unique(groups)){
      a<-x[groups==k,i,drop=F]
      b<-y[groups==k,j,drop=F]
      tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method,exact=FALSE)$p.value,k)
      if(is.null(df)){
        df<-tmp  
      }
      else{
        df<-rbind(df,tmp)
      }    
    }
  }
}

df<-data.frame(row.names=NULL,df)
colnames(df)<-c("Taxa","Env","Correlation","Pvalue","Type")
df$Pvalue<-as.numeric(as.character(df$Pvalue))
df$AdjPvalue<-rep(0,dim(df)[1])
df$Correlation<-as.numeric(as.character(df$Correlation))

#You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
# 1 -> donot adjust
# 2 -> adjust Env + Type (column on the correlation plot)
# 3 -> adjust Taxa + Type (row on the correlation plot for each type)
# 4 -> adjust Taxa (row on the correlation plot)
# 5 -> adjust Env (panel on the correlation plot)
adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")
adjustment<-5

if(adjustment==1){
  df$AdjPvalue<-df$Pvalue
} else if (adjustment==2){
  for(i in unique(df$Env)){
    for(j in unique(df$Type)){
      sel<-df$Env==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==3){
  for(i in unique(df$Taxa)){
    for(j in unique(df$Type)){
      sel<-df$Taxa==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==4){
  for(i in unique(df$Taxa)){
    sel<-df$Taxa==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
} else if (adjustment==5){
  for(i in unique(df$Env)){
    sel<-df$Env==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
}
```


```{r eval=F}
#Now we generate the labels for significant values
df$Significance<-cut(df$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

#We ignore NAs
df<-df[complete.cases(df),]

#We want to reorganize the Env data based on how they appear
#df$Env<-factor(df$Env,as.character(df$Env))

#We use the function to change the labels for facet_grid in ggplot2
Env_labeller <- function(variable,value){
  return(sel_env_label[as.character(value),"Trans"])
}

#Make plot, organizing taxa in alphabetical order
p <- ggplot(aes(x=Type, y=factor(Taxa, 
                                 levels = rev(levels(factor(Taxa)))), fill=Correlation), data=df)
p <- p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") 
p<-p+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=0.5))
p<-p+geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL, fill=method)
p<-p+facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x")#),labeller=Env_labeller)
#pdf(paste("plots/Correlation_low_ASVs_",adjustment_label[adjustment],".pdf",sep=""),height=8,width=22)
print(p)
#dev.off()
head(df)
```

```{r eval=F}
#Correlation analysis just among environmental variables
#Load metatable
meta_table <- read.table("SR2113_metadata_final.txt", sep = "\t", header = TRUE)
head(meta_table)
#Filter by core
meta_table_mc01 <- filter(meta_table, Core == "MC01")
meta_table_gc02 <- filter(meta_table, Core == "GC02")
meta_table_gc04 <- filter(meta_table, Core == "GC04")
meta_table_gc06 <- filter(meta_table, Core == "GC06")

#Remove absent environmental variables from cores SCB01 and GC06
sel_env_mc01<-c("Real_depth_cm","Iron","Manganese","Nitrate","Sulfate","Sulfide","Sulfite","Thiosulfate")
sel_env_gc06<-c("Real_depth_cm","Iron","Manganese","Nitrate","Sulfate","Sulfide","Sulfite","Thiosulfate")

#keep only relevant variables
meta_table_mc_01_filtered<-meta_table_mc01[,sel_env_mc01]
meta_table_gc_02_filtered<-meta_table_gc02[,sel_env]
meta_table_gc_04_filtered<-meta_table_gc04[,sel_env]
meta_table_gc_06_filtered<-meta_table_gc06[,sel_env_gc06]

#Remove outliers
#meta_table_mc_01_filtered_2 <- subset(meta_table_mc_01_filtered, Real_depth_cm %nin% c("151", "91", "146"))
#meta_table_gc_02_filtered_2 <- subset(meta_table_gc_02_filtered, Real_depth_cm %nin% c())
meta_table_gc_04_filtered_2 <- subset(meta_table_gc_04_filtered, Real_depth_cm %nin% c("26", "21"))
meta_table_gc_06_filtered_2 <- subset(meta_table_gc_06_filtered, Real_depth_cm %nin% c("151", "91", "146"))

head(meta_table_gc_06_filtered_2)
```
```{r eval=F}
updated_meta_table_gc04 <- subset(meta_table_gc04, ID %nin% c("GC04_26", "GC04_21"))
#Make scatter plots to visualize correlations
p <- ggplot(updated_meta_table_gc04) +
  aes(x = Sulfide, y = Sulfite) +
  geom_point(colour = "#0c4c8a") +
  geom_text(aes(label=ID), size=3) +
  theme_minimal()
print (p)
updated_meta_table_gc06 <- subset(meta_table_gc06, ID %nin% c("GC06_151", "GC06_91", "GC06_146"))
p <- ggplot(updated_meta_table_gc06) +
  aes(x = Sulfide, y = Sulfite) +
  geom_point(colour = "#0c4c8a") +
  geom_text(aes(label=ID), size=3, hjust=-0.5) +
  theme_minimal()
print (p)
```

```{r eval=F}
#Function to flatten correlation matrix
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

#Create correlation matrices for each core
res_mc_01 <- rcorr(as.matrix(meta_table_mc_01_filtered_1), type = "pearson")
res_gc_02 <- rcorr(as.matrix(meta_table_gc_02_filtered_1), type = "pearson")
res_gc_04 <- rcorr(as.matrix(meta_table_gc_04_filtered_2), type = "pearson")
res_gc_06 <- rcorr(as.matrix(meta_table_gc_06_filtered_2), type = "pearson")

#Flatten corr matrix for MC 01
flat_corr_mc_01 <- flattenCorrMatrix(res_mc_01$r, res_mc_01$P)
#Sort by row name
df_corr_mc_01 <- with(flat_corr_mc_01,  flat_corr_mc_01[order(row) , ])
#Round numeric columns to 3 digits
df_corr_mc_01 <- df_corr_mc_01 %>% 
  mutate(across(where(is.numeric), round, digits=3))
#See dataframe
df_corr_mc_01
#Build correlation plot
chart.Correlation(meta_table_mc_01_filtered_1)
```

```{r eval=F}
#Flatten corr matrix for GC 02
flat_corr_gc_02 <- flattenCorrMatrix(res_gc_02$r, res_gc_02$P)
#Sort by row name
df_corr_gc_02 <- with(flat_corr_gc_02,  flat_corr_gc_02[order(row) , ])
#Round numeric columns to 3 digits
df_corr_gc_02 <- df_corr_gc_02 %>% 
  mutate(across(where(is.numeric), round, digits=3))
#See dataframe
df_corr_gc_02
#Build correlation plot
chart.Correlation(meta_table_gc_02_filtered_1)
```

```{r eval=F}
#Flatten corr matrix for GC 04
flat_corr_gc_04 <- flattenCorrMatrix(res_gc_04$r, res_gc_04$P)
#Sort by row name
df_corr_gc_04 <- with(flat_corr_gc_04,  flat_corr_gc_04[order(row) , ])
#Round numeric columns to 3 digits
df_corr_gc_04 <- df_corr_gc_04 %>% 
  mutate(across(where(is.numeric), round, digits=3))
#See dataframe
df_corr_gc_04
#Build correlation plot
chart.Correlation(meta_table_gc_04_filtered_2)
```

```{r eval=F}
#Flatten corr matrix for GC 06
flat_corr_gc_06 <- flattenCorrMatrix(res_gc_06$r, res_gc_06$P)
#Sort by row name
df_corr_gc_06 <- with(flat_corr_gc_06,  flat_corr_gc_06[order(row) , ])
#Round numeric columns to 3 digits
df_corr_gc_06 <- df_corr_gc_06 %>% 
  mutate(across(where(is.numeric), round, digits=3))
#Build correlation plot
chart.Correlation(meta_table_gc_06_filtered_2)

```
