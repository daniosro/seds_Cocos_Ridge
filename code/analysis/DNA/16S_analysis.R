#We'll first change the location we are working in

# this is just getting our unique user ID
user_id <- system("echo $USER", intern = TRUE)

# this is setting the working directory to the location we were just working in
#setwd(paste0("/export/data1/projects/daniela/SR2113/June2022_16S/"))
setwd(paste0("/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/"))
#run this to confirm the working directory is correct
getwd()

#Load the packages we are going to use
library(tidyverse)
library(phyloseq)
packageVersion("phyloseq")
library(vegan)
library(dendextend)
library(RColorBrewer)
library(DESeq2)
library(DataCombine)
library(scales)
library(reshape2)
library(tidyr)
library(viridis)
library(Biostrings)
library(DECIPHER)

#Read in our processed files
## count table
counts_tab <- read.table("ASVs-counts.tsv", header = TRUE, row.names = 1,
                         check.names = FALSE, sep = "\t")
counts_tab$ASV <- NULL
proteobacteria_counts_tab <- read.table("ASVs-counts_proteobacteria.tsv", header = TRUE, row.names = 1,
                                        check.names = FALSE, sep = "\t")
proteobacteria_counts_tab$ASV <- NULL
colnames(proteobacteria_counts_tab)
# since we already used our "blank", we're going to 
# remove it from our data table here
# the sign - removes that column from the variable
#counts_tab <- counts_tab %>% select(-`MLW-DNAexBlank-50cyc_S193`)

relative_tab <- apply(counts_tab, 2, function(x) 100*(x/sum(x)))

## taxonomy table

tax_tab <- read.table("ASVs-taxonomy.tsv", header = TRUE, row.names = 1,
                      check.names = FALSE, sep = "\t")
tax_tab$ASV <- rownames(tax_tab)
asv_seqs <- readDNAStringSet("ASVs.fa")

# Make a sample information table
#First let’s read in our unique sample names again:
samples <- scan("unique-sample-IDs.txt", what = "character")
samples # will print them out
#Replace the names of the samples in all the files to make them suitable for binning by categories
#Create replacements dataframe
Replaces <-data.frame (find = c("DO-101_S223",
                                "DO-102_S255",
                                "DO-10_S245",
                                "DO-110_S207",
                                "DO-113_S242",
                                "DO-114_S208",
                                "DO-115_S209",
                                "DO-116_S206",
                                "DO-11_S254",
                                "DO-12_S204",
                                "DO-13_S253",
                                "DO-153_S233",
                                "DO-154_S218",
                                "DO-156_S235",
                                "DO-159_S221",
                                "DO-15_S251",
                                "DO-161_S243",
                                "DO-162_S219",
                                "DO-17_S203",
                                "DO-19_S249",
                                "DO-200_S225",
                                "DO-202_S224",
                                "DO-203_S195",
                                "DO-206_S194",
                                "DO-209_S192",
                                "DO-20_S202",
                                "DO-211_S193",
                                "DO-21_S200",
                                "DO-23_S246",
                                "DO-2_S252",
                                "DO-30_S197",
                                "DO-31_S198",
                                "DO-35_S229",
                                "DO-37_S236",
                                "DO-39_S205",
                                "DO-3_S201",
                                "DO-46_S234",
                                "DO-48_S237",
                                "DO-49_S239",
                                "DO-50_S213",
                                "DO-51_S238",
                                "DO-53_S230",
                                "DO-54_S232",
                                "DO-55_S244",
                                "DO-58_S196",
                                "DO-5_S241",
                                "DO-60_S228",
                                "DO-61_S226",
                                "DO-63_S231",
                                "DO-6_S250",
                                "DO-70_S210",
                                "DO-72_S211",
                                "DO-73_S199",
                                "DO-74_S217",
                                "DO-77_S220",
                                "DO-7_S222",
                                "DO-81_S227",
                                "DO-83_S248",
                                "DO-86_S216",
                                "DO-88_S212",
                                "DO-8_S247",
                                "DO-98_S215",
                                "DO-99_S240",
                                "DO-a_S214",
                                "DO-01_S234",
                                "DO-105_S217",
                                "DO-107_S233",
                                "DO-14_S213",
                                "DO-152_S231",
                                "DO-152-star_S236",
                                "DO-154_S211",
                                "DO-157_S216",
                                "DO-158_S230",
                                "DO-15_S214",
                                "DO-16_S215",
                                "DO-18_S235",
                                "DO-28_S227",
                                "DO-35_S209",
                                "DO-44_S232",
                                "DO-47_S218",
                                "DO-55_S205",
                                "DO-59_S229",
                                "DO-68_S219",
                                "DO-69_S226",
                                "DO-71_S220",
                                "DO-78_S208",
                                "DO-79_S228",
                                "DO-84_S221",
                                "DO-85_S225",
                                "DO-87_S224",
                                "DO-93_S223",
                                "DO-98_S210",
                                "DO-ne-35_S212",
                                "DO-ne-45_S206",
                                "DO-np-35_S222",
                                "DO-np-45_S207"), 
                       replace = c("DO-101_GC04_136",
                                   "DO-102_GC06_111",
                                   "DO-10_GC06_81",
                                   "DO-110_GC06_156",
                                   "DO-113_GC04_156",
                                   "DO-114_GC06_166",
                                   "DO-115_GC06_171",
                                   "DO-116_GC06_146",
                                   "DO-11_GC06_71",
                                   "DO-12_GC06_121",
                                   "DO-13_GC06_41",
                                   "DO-153_GC02_63.5",
                                   "DO-154_GC04_170june",
                                   "DO-156_GC02_90.5",
                                   "DO-159_GC04_206",
                                   "DO-15_GC06_185june",
                                   "DO-161_GC04_171",
                                   "DO-162_GC04_191",
                                   "DO-17_GC06_186",
                                   "DO-19_GC06_176",
                                   "DO-200_MC01_23",
                                   "DO-202_MC01_18",
                                   "DO-203_MC01_12",
                                   "DO-206_MC01_8",
                                   "DO-209_MC01_2",
                                   "DO-20_GC06_101",
                                   "DO-211_MC01_4",
                                   "DO-21_GC06_51",
                                   "DO-23_GC06_91",
                                   "DO-2_GC06_21",
                                   "DO-30_GC02_114.5",
                                   "DO-31_GC04_201",
                                   "DO-35_GC02_33june",
                                   "DO-37_GC02_99.5",
                                   "DO-39_GC06_131",
                                   "DO-3_GC06_61",
                                   "DO-46_GC02_84.5",
                                   "DO-48_GC02_105.5",
                                   "DO-49_GC02_126.5",
                                   "DO-50_GC04_86",
                                   "DO-51_GC02_111.5",
                                   "DO-53_GC02_45.5",
                                   "DO-54_GC02_57.5",
                                   "DO-55_GC06_5june",
                                   "DO-58_GC02_72.5",
                                   "DO-5_GC06_31",
                                   "DO-60_GC02_81.5",
                                   "DO-61_GC02_18.5",
                                   "DO-63_GC02_51.5",
                                   "DO-6_GC06_181",
                                   "DO-70_GC02_6.5",
                                   "DO-72_GC02_39.5",
                                   "DO-73_GC04_186",
                                   "DO-74_GC04_166",
                                   "DO-77_GC04_196",
                                   "DO-7_GC06_16",
                                   "DO-81_GC02_21.5",
                                   "DO-83_GC06_161",
                                   "DO-86_GC04_146",
                                   "DO-88_GC04_71",
                                   "DO-8_GC06_141",
                                   "DO-98_GC04_110june",
                                   "DO-99_GC04_106",
                                   "DO-a_GC04_96",
                                   "DO-01_GC06_126",
                                   "DO-105_GC04_26",
                                   "DO-107_GC02_15.5",
                                   "DO-14_MC01_27",
                                   "DO-152_GC02_60.5",
                                   "DO-152_GC02_93.5",
                                   "DO-154_GC04_176",
                                   "DO-157_GC04_16",
                                   "DO-158_GC02_117.5",
                                   "DO-15_MC01_31",
                                   "DO-16_MC01_36",
                                   "DO-18_GC06_56",
                                   "DO-28_GC04_76",
                                   "DO-35_GC02_36.5",
                                   "DO-44_GC02_96.5",
                                   "DO-47_GC04_41",
                                   "DO-55_GC06_11",
                                   "DO-59_GC02_87.5",
                                   "DO-68_GC04_51",
                                   "DO-69_GC04_46",
                                   "DO-71_GC04_61",
                                   "DO-78_GC02_9.5",
                                   "DO-79_GC04_66",
                                   "DO-84_GC04_11",
                                   "DO-85_GC04_126",
                                   "DO-87_GC02_12.5",
                                   "DO-93_GC04_31",
                                   "DO-98_GC04_116",
                                   "DO-ne_GCnone_35",
                                   "DO-ne_GCnone_45",
                                   "DO-np_GCnone_35",
                                   "DO-np_GCnone_45"))

#Make replacements in count table
colnames (counts_tab) <- c("DO-101_GC04_136",
                           "DO-102_GC06_111",
                           "DO-10_GC06_81",
                           "DO-110_GC06_156",
                           "DO-113_GC04_156",
                           "DO-114_GC06_166",
                           "DO-115_GC06_171",
                           "DO-116_GC06_146",
                           "DO-11_GC06_71",
                           "DO-12_GC06_121",
                           "DO-13_GC06_41",
                           "DO-153_GC02_63.5",
                           "DO-154_GC04_170june",
                           "DO-156_GC02_90.5",
                           "DO-159_GC04_206",
                           "DO-15_GC06_185june",
                           "DO-161_GC04_171",
                           "DO-162_GC04_191",
                           "DO-17_GC06_186",
                           "DO-19_GC06_176",
                           "DO-200_MC01_23",
                           "DO-202_MC01_18",
                           "DO-203_MC01_12",
                           "DO-206_MC01_8",
                           "DO-209_MC01_2",
                           "DO-20_GC06_101",
                           "DO-211_MC01_4",
                           "DO-21_GC06_51",
                           "DO-23_GC06_91",
                           "DO-2_GC06_21",
                           "DO-30_GC02_114.5",
                           "DO-31_GC04_201",
                           "DO-35_GC02_33june",
                           "DO-37_GC02_99.5",
                           "DO-39_GC06_131",
                           "DO-3_GC06_61",
                           "DO-46_GC02_84.5",
                           "DO-48_GC02_105.5",
                           "DO-49_GC02_126.5",
                           "DO-50_GC04_86",
                           "DO-51_GC02_111.5",
                           "DO-53_GC02_45.5",
                           "DO-54_GC02_57.5",
                           "DO-55_GC06_5june",
                           "DO-58_GC02_72.5",
                           "DO-5_GC06_31",
                           "DO-60_GC02_81.5",
                           "DO-61_GC02_18.5",
                           "DO-63_GC02_51.5",
                           "DO-6_GC06_181",
                           "DO-70_GC02_6.5",
                           "DO-72_GC02_39.5",
                           "DO-73_GC04_186",
                           "DO-74_GC04_166",
                           "DO-77_GC04_196",
                           "DO-7_GC06_16",
                           "DO-81_GC02_21.5",
                           "DO-83_GC06_161",
                           "DO-86_GC04_146",
                           "DO-88_GC04_71",
                           "DO-8_GC06_141",
                           "DO-98_GC04_110june",
                           "DO-99_GC04_106",
                           "DO-a_GC04_96",
                           "DO-01_GC06_126",
                           "DO-105_GC04_26",
                           "DO-107_GC02_15.5",
                           "DO-14_MC01_27",
                           "DO-152_GC02_60.5",
                           "DO-152_GC02_93.5",
                           "DO-154_GC04_176",
                           "DO-157_GC04_16",
                           "DO-158_GC02_117.5",
                           "DO-15_MC01_31",
                           "DO-16_MC01_36",
                           "DO-18_GC06_56",
                           "DO-28_GC04_76",
                           "DO-35_GC02_36.5",
                           "DO-44_GC02_96.5",
                           "DO-47_GC04_41",
                           "DO-55_GC06_11",
                           "DO-59_GC02_87.5",
                           "DO-68_GC04_51",
                           "DO-69_GC04_46",
                           "DO-71_GC04_61",
                           "DO-78_GC02_9.5",
                           "DO-79_GC04_66",
                           "DO-84_GC04_11",
                           "DO-85_GC04_126",
                           "DO-87_GC02_12.5",
                           "DO-93_GC04_31",
                           "DO-98_GC04_116",
                           "DO-ne_GCnone_35",
                           "DO-ne_GCnone_45",
                           "DO-np_GCnone_35",
                           "DO-np_GCnone_45")

#Make replacements in count table
colnames (counts_tab) <- c("DO-101_GC04_136",
                           "DO-102_GC06_111",
                           "DO-10_GC06_81",
                           "DO-110_GC06_156",
                           "DO-113_GC04_156",
                           "DO-114_GC06_166",
                           "DO-115_GC06_171",
                           "DO-116_GC06_146",
                           "DO-11_GC06_71",
                           "DO-12_GC06_121",
                           "DO-13_GC06_41",
                           "DO-153_GC02_63.5",
                           "DO-154_GC04_170june",
                           "DO-156_GC02_90.5",
                           "DO-159_GC04_206",
                           "DO-15_GC06_185june",
                           "DO-161_GC04_171",
                           "DO-162_GC04_191",
                           "DO-17_GC06_186",
                           "DO-19_GC06_176",
                           "DO-200_MC01_23",
                           "DO-202_MC01_18",
                           "DO-203_MC01_12",
                           "DO-206_MC01_8",
                           "DO-209_MC01_2",
                           "DO-20_GC06_101",
                           "DO-211_MC01_4",
                           "DO-21_GC06_51",
                           "DO-23_GC06_91",
                           "DO-2_GC06_21",
                           "DO-30_GC02_114.5",
                           "DO-31_GC04_201",
                           "DO-35_GC02_33june",
                           "DO-37_GC02_99.5",
                           "DO-39_GC06_131",
                           "DO-3_GC06_61",
                           "DO-46_GC02_84.5",
                           "DO-48_GC02_105.5",
                           "DO-49_GC02_126.5",
                           "DO-50_GC04_86",
                           "DO-51_GC02_111.5",
                           "DO-53_GC02_45.5",
                           "DO-54_GC02_57.5",
                           "DO-55_GC06_5june",
                           "DO-58_GC02_72.5",
                           "DO-5_GC06_31",
                           "DO-60_GC02_81.5",
                           "DO-61_GC02_18.5",
                           "DO-63_GC02_51.5",
                           "DO-6_GC06_181",
                           "DO-70_GC02_6.5",
                           "DO-72_GC02_39.5",
                           "DO-73_GC04_186",
                           "DO-74_GC04_166",
                           "DO-77_GC04_196",
                           "DO-7_GC06_16",
                           "DO-81_GC02_21.5",
                           "DO-83_GC06_161",
                           "DO-86_GC04_146",
                           "DO-88_GC04_71",
                           "DO-8_GC06_141",
                           "DO-98_GC04_110june",
                           "DO-99_GC04_106",
                           "DO-a_GC04_96",
                           "DO-01_GC06_126",
                           "DO-105_GC04_26",
                           "DO-107_GC02_15.5",
                           "DO-14_MC01_27",
                           "DO-152_GC02_60.5",
                           "DO-152_GC02_93.5",
                           "DO-154_GC04_176",
                           "DO-157_GC04_16",
                           "DO-158_GC02_117.5",
                           "DO-15_MC01_31",
                           "DO-16_MC01_36",
                           "DO-18_GC06_56",
                           "DO-28_GC04_76",
                           "DO-35_GC02_36.5",
                           "DO-44_GC02_96.5",
                           "DO-47_GC04_41",
                           "DO-55_GC06_11",
                           "DO-59_GC02_87.5",
                           "DO-68_GC04_51",
                           "DO-69_GC04_46",
                           "DO-71_GC04_61",
                           "DO-78_GC02_9.5",
                           "DO-79_GC04_66",
                           "DO-84_GC04_11",
                           "DO-85_GC04_126",
                           "DO-87_GC02_12.5",
                           "DO-93_GC04_31",
                           "DO-98_GC04_116",
                           "DO-ne_GCnone_35",
                           "DO-ne_GCnone_45",
                           "DO-np_GCnone_35",
                           "DO-np_GCnone_45")

#Dropping the blank and june samples
counts_tab <- counts_tab [!grepl("np|ne", colnames(counts_tab))]
colnames(counts_tab)

#Make replacements in proteobacteria count table
colnames (proteobacteria_counts_tab) <- c("DO-101_GC04_136",
                                          "DO-102_GC06_111",
                                          "DO-10_GC06_81",
                                          "DO-110_GC06_156",
                                          "DO-113_GC04_156",
                                          "DO-114_GC06_166",
                                          "DO-115_GC06_171",
                                          "DO-116_GC06_146",
                                          "DO-11_GC06_71",
                                          "DO-12_GC06_121",
                                          "DO-13_GC06_41",
                                          "DO-153_GC02_63.5",
                                          "DO-154_GC04_170june",
                                          "DO-156_GC02_90.5",
                                          "DO-159_GC04_206",
                                          "DO-15_GC06_185june",
                                          "DO-161_GC04_171",
                                          "DO-162_GC04_191",
                                          "DO-17_GC06_186",
                                          "DO-19_GC06_176",
                                          "DO-200_MC01_23",
                                          "DO-202_MC01_18",
                                          "DO-203_MC01_12",
                                          "DO-206_MC01_8",
                                          "DO-209_MC01_2",
                                          "DO-20_GC06_101",
                                          "DO-211_MC01_4",
                                          "DO-21_GC06_51",
                                          "DO-23_GC06_91",
                                          "DO-2_GC06_21",
                                          "DO-30_GC02_114.5",
                                          "DO-31_GC04_201",
                                          "DO-35_GC02_33june",
                                          "DO-37_GC02_99.5",
                                          "DO-39_GC06_131",
                                          "DO-3_GC06_61",
                                          "DO-46_GC02_84.5",
                                          "DO-48_GC02_105.5",
                                          "DO-49_GC02_126.5",
                                          "DO-50_GC04_86",
                                          "DO-51_GC02_111.5",
                                          "DO-53_GC02_45.5",
                                          "DO-54_GC02_57.5",
                                          "DO-55_GC06_5june",
                                          "DO-58_GC02_72.5",
                                          "DO-5_GC06_31",
                                          "DO-60_GC02_81.5",
                                          "DO-61_GC02_18.5",
                                          "DO-63_GC02_51.5",
                                          "DO-6_GC06_181",
                                          "DO-70_GC02_6.5",
                                          "DO-72_GC02_39.5",
                                          "DO-73_GC04_186",
                                          "DO-74_GC04_166",
                                          "DO-77_GC04_196",
                                          "DO-7_GC06_16",
                                          "DO-81_GC02_21.5",
                                          "DO-83_GC06_161",
                                          "DO-86_GC04_146",
                                          "DO-88_GC04_71",
                                          "DO-8_GC06_141",
                                          "DO-98_GC04_110june",
                                          "DO-99_GC04_106",
                                          "DO-a_GC04_96",
                                          "DO-01_GC06_126",
                                          "DO-105_GC04_26",
                                          "DO-107_GC02_15.5",
                                          "DO-14_MC01_27",
                                          "DO-152_GC02_60.5",
                                          "DO-152_GC02_93.5",
                                          "DO-154_GC04_176",
                                          "DO-157_GC04_16",
                                          "DO-158_GC02_117.5",
                                          "DO-15_MC01_31",
                                          "DO-16_MC01_36",
                                          "DO-18_GC06_56",
                                          "DO-28_GC04_76",
                                          "DO-35_GC02_36.5",
                                          "DO-44_GC02_96.5",
                                          "DO-47_GC04_41",
                                          "DO-55_GC06_11",
                                          "DO-59_GC02_87.5",
                                          "DO-68_GC04_51",
                                          "DO-69_GC04_46",
                                          "DO-71_GC04_61",
                                          "DO-78_GC02_9.5",
                                          "DO-79_GC04_66",
                                          "DO-84_GC04_11",
                                          "DO-85_GC04_126",
                                          "DO-87_GC02_12.5",
                                          "DO-93_GC04_31",
                                          "DO-98_GC04_116",
                                          "DO-ne_GCnone_35",
                                          "DO-ne_GCnone_45",
                                          "DO-np_GCnone_35",
                                          "DO-np_GCnone_45")
#Dropping the blank and june samples
proteobacteria_counts_tab <- proteobacteria_counts_tab [!grepl("np|ne", colnames(proteobacteria_counts_tab))]
colnames(proteobacteria_counts_tab)


# Make replacements in unique_IDs
samples <- as.character (Replaces[match(samples, Replaces$find), "replace"])

# dropping the "blank" sample
samples <- samples [!grepl("np|ne", samples)]
samples


#Make a directory to save plots
# here we are making a directory to hold the output quality-filtered reads
dir.create("plots")

#Removing likely contaminants
predicted_controls <- grepl("np|ne", colnames(counts_tab))
predicted_controls_names <- setNames(predicted_controls, colnames(counts_tab))
conc <- c() # c(410,310,220,400,394,120)  # MAKE SURE IN SAME ORDER AS counts_tab COLUMNS

if (any(predicted_controls) | (length(conc) > 1)){
  if (length(conc) == length(colnames(counts_tab))) {
    print('Decontam using "concentration" method')
    contam_predict <- decontam::isContaminant(t(counts_tab), method='frequency', conc=conc) 
    contam_asvs <- rownames(contam_predict)[contam_predict$contaminant] 
    p <- melt(as.data.frame(cbind(t(relative_tab[contam_asvs,]), conc=conc)), id.vars = 'conc', 
              value.name = 'relabund', variable.name = 'contam_asv') %>% filter(relabund > 0) %>% 
      ggplot(aes(x=conc, y=relabund, color=contam_asv)) + scale_y_continuous(expand=c(0,0)) + scale_x_log10() + 
      geom_smooth(method=stats::lm, se = F, linetype='dashed', size=0.5) + geom_point() + theme_classic()
    print(p)
  } else {
    contam_predict <- decontam::isContaminant(t(counts_tab), method='prevalence', neg=predicted_controls, threshold=0.5) 
    contam_asvs <- rownames(contam_predict)[contam_predict$contaminant] 
  }
  print(apply(tax_tab[contam_asvs,], 1, function(x) paste0(x, collapse=";")))
  
  relative_tab %>% melt(varnames = c('ASV','Sample'), value.name = 'relabund') %>% 
    mutate(decontam=factor(setNames(contam_predict$contaminant, rownames(contam_predict))[ASV], levels=c('TRUE','FALSE'))) %>% 
    group_by(decontam, Sample) %>% summarise(relabund=sum(relabund)) %>% ungroup() %>%
    mutate(library_type=replace(rep("Real", length(Sample)), predicted_controls, "Control")) %>% 
    mutate(Sample = factor(Sample, levels=unique(Sample[order(decontam,relabund)]))) %>%
    ggplot(aes(x=Sample, y=relabund, fill=decontam)) + geom_col() + scale_fill_manual(values=c("TRUE"='#ca2a00', 'FALSE'='black')) +
    facet_grid(.~library_type, scales='free_x', space='free_x') + scale_y_continuous(expand=c(0,0)) +
    theme_classic() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), axis.ticks.x = element_blank())
  
  ## THIS WILL UPDATE THE COUNTS - running this block more than once will produce different results
  counts_tab <- counts_tab[!(rownames(counts_tab) %in% contam_asvs),]
  relative_tab <- apply(counts_tab, 2, function(x) 100*(x/sum(x)))
}


# Here we are making our samples “vector” into a column of a dataframe, 
#then splitting each sample ID into multiple columns with headers we are 
#providing, then pulling out the column we want:

#Pulling out core
core <- data.frame(samples) %>% 
  separate(samples, c("DO", "ID", "GC", "Depth_cm")) %>% 
  pull(GC)

core
#Number of samples per core
as.data.frame(cbind(vec = unique(core), n = tabulate(match(core, unique(core)))))

#Pulling out depth
depth <- data.frame(samples) %>% 
  separate(samples, c("DO", "ID", "GC", "Depth_cm")) %>% 
  pull(Depth_cm)

depth

#Now we are making a new dataframe that holds our sample IDs, cores, and depths, 
# and assigning it to a new variable:
sample_info_tab_1 <- data.frame("Sample" = samples, "Core" = core, "Depth" = depth)
sample_info_tab_1

#It will help with some plotting later if we add colors associated with 
# different characteristics now. So here we will do that:
# first let's see how many unique values there are for Core and Depth
sample_info_tab_1 %>% pull(Core) %>% unique %>% length
# 4 unique ones for core
sample_info_tab_1 %>% pull(Depth) %>% unique %>% length
# 72 unique values for depth
# let's see what the distribution looks like for depth
#hist(sample_info_tab %>% pull(Depth) %>% as.numeric)

# We’ll make 4 bins for depth: 0-20, 21-50, 51-100, and 101-220

# first let's get 4 color codes we can use
four_colors <-colorblind.pal(4, "Dark2")
four_colors


# reminding what the values are
sample_info_tab_1 %>% pull(Core) %>% unique()

# now we're going to loop through the values in our Core column
# and set the color value based on the core
# making an empty vector for the core colors
core_colors <- c()

for ( entry in sample_info_tab_1 %>% pull(Core)) {
  
  if ( entry == "MC01" ) {
    
    core_colors <- c(core_colors, four_colors[1])
    
  } else if ( entry == "GC02" ) {
    
    core_colors <- c(core_colors, four_colors[2])
    
  } else if ( entry == "GC04" ) {
    
    core_colors <- c(core_colors, four_colors[3])
    
  } else {
    
    core_colors <- c(core_colors, four_colors[4])
  }
  
}

core_colors

# now doing something similar for Depth, but also making a new column
# holding the group labels
depth_colors <- c()
depth_bins <- c()

for ( entry in sample_info_tab_1 %>% pull(Depth)) {
  
  if ( as.numeric(entry) %>% between(0,20)) {
    
    #depth_colors <- c(depth_colors, four_colors[1])
    depth_bins <- c(depth_bins, "0-20")
    
  } else if ( as.numeric(entry) %>% between(21,40)) {
    
    #depth_colors <- c(depth_colors, four_colors[2])
    depth_bins <- c(depth_bins, "21-40")
    
  } else if ( as.numeric(entry) %>% between(41,60)) {
    
    #depth_colors <- c(depth_colors, four_colors[2])
    depth_bins <- c(depth_bins, "41-60")
    
  } else if ( as.numeric(entry) %>% between(61,100)) {
    
    #depth_colors <- c(depth_colors, four_colors[3])
    depth_bins <- c(depth_bins, "61-100")
    
  } else if ( as.numeric(entry) %>% between(101,140)) {
    
    #depth_colors <- c(depth_colors, four_colors[3])
    depth_bins <- c(depth_bins, "101-140")
    
  } else {
    
    #depth_colors <- c(depth_colors, four_colors[4])
    depth_bins <- c(depth_bins, "141-220")
    
  }
  
}

#depth_colors
depth_bins

#And now we can add those color columns to our sample_info_tab:
sample_info_tab_1$core_colors <- core_colors
#sample_info_tab_1$depth_colors <- depth_colors
sample_info_tab_1$Depth_bin <- depth_bins

sample_info_tab_1

#Beta diversity

#We’re going to use Euclidean distances to generate some exploratory 
#visualizations of our samples. Since differences in sampling depths 
#between samples can influence distance/dissimilarity metrics, we first 
# need to do something to normalize our samples.
# One better method of normalizing across samples is to use a variance 
#stabilizing transformation – which fortunately we can do with the DESeq2
#package that we already have loaded 

#Normalize data
# first we need to make a deseq2 data object

#deseq_counts <- DESeqDataSetFromMatrix(counts_tab, colData = sample_info_tab, 
#                                    design = ~1)
deseq_counts <- DESeqDataSetFromMatrix(counts_tab, data.frame(samples=colnames(counts_tab)), 
                                       design=~samples)

deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")

# we have to include the "colData" and "design" arguments because they are 
# required, as they are needed for further downstream processing by DESeq2, 
# but for our purposes of simply transforming the data right now, they don't 
# matter

deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# and here is how we can pull out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)
# and this is a regular table now
head(vst_trans_count_tab)
vst_non_neg <- vst_trans_count_tab
vst_non_neg[vst_non_neg<0] <- 0
head(vst_non_neg)
#Export table for NMDS analysis
write.table(vst_non_neg, "ASVs-counts_norm.tsv", sep = "\t", quote = FALSE, col.names = NA)
#Hierarchical clustering
# We’re going to turn our normalized table into a distance matrix now, which 
#will be the foundation of what we use for hierarchical clustering:

# the dist() program works on rows, so we transpose our 
# table with the t() function
euc_dist <- dist(t(vst_trans_count_tab))

# hclust() is the function we use to cluster our distance matrix
#euc_clust <- hclust(euc_dist, method = "ward.D2")
euc_dend <- as.dendrogram(hclust(euc_dist, method="ward.D2"), hang=0.1)

# and hclust objects can be plotted with the generic plot() function
plot(euc_clust)

# but I like changing them to dendrograms as they are easier to
# color, and they can be rotated with rotate() if wanted

# making one for core and depth so we can color both ways
core_euc_dend <- as.dendrogram(euc_clust, hang = 0.1)
depth_euc_dend <- as.dendrogram(euc_clust, hang = 0.1)

# this line is getting a vector of colors in the order the sample names are in the dendrogram
core_dend_colors <- as.character(sample_info_tab$core_color[order.dendrogram(core_euc_dend)])

# now doing the same for depth
depth_dend_colors <- as.character(sample_info_tab$depth_color[order.dendrogram(depth_euc_dend)])

# Coloring the labels in the dendrogram objects for each
labels_colors(core_euc_dend) <- core_dend_colors
labels_colors(depth_euc_dend) <- depth_dend_colors

# now here's plotting the one colored by core
plot(core_euc_dend, ylab="VST Euc. dist.", title(main = "Colored by Core"))

# Adjusting margins so the labels are cut off from the bottom
par(mar = c(10, 4, 4, 2))
plot(core_euc_dend, ylab="VST Euc. dist.")
# Adding a legend
legend("topleft", legend = c("MC01", "GC02", "GC04", "GC06"), fill = four_colors,
       pt.cex = 0.5, title = "Core")

# and here's the one colored by depth
plot(depth_euc_dend, ylab="VST Euc. dist.", main = "Colored by Depth")
# adding a legend
legend("topleft", legend = c("0-20", "11-50", "51-100","101-210"), fill = four_colors,
       pt.cex = 0.5, title = "Depth")

## setting margins back to the default
# we can find that info in ?par
par(mar = c(5, 4, 4, 2))

#Multidimensional analysis (MDA)
library(reshape2)

makeMDSbyLevel <- function(counts=vst_trans_count_tab, tax=tax_tab, level='ASV', distance='bray', try=100, trymax=500) {
  counts <- cbind(tax[rownames(counts),], counts)
  counts <- dcast(melt(counts, variable.name='sample'), sample~get(level), value.var = 'value', fun.aggregate = sum)
  rownames(counts) <- counts$sample
  
  mds <- metaMDS(counts[,2:ncol(counts)], distance = distance, try = try, trymax = trymax)
  return(mds)
}

mds <- makeMDSbyLevel(counts=vst_non_neg, level='genus',try=100, trymax=200)

# if you have sample metadata, good place to merge that onto mds table to use for coloring points
mds_plot <- data.frame(mds$points)
mds_plot$sample <- rownames(mds_plot)

#Creating GC and Depth columns in the mds_plot df
mds_plot <- data.frame(mds_plot,do.call(rbind,str_split(mds_plot$sample,"_")))
mds_plot$X1 <- NULL
colnames(mds_plot)[colnames(mds_plot) == 'X2'] <- 'Core'
colnames(mds_plot)[colnames(mds_plot) == 'X3'] <- 'Depth'

# Make color bins for depth
# holding the group labels
depth_colors <- c()
depth_bins <- c()

for ( entry in mds_plot %>% pull(Depth)) {
  
  if ( as.numeric(entry) %>% between(0,20)) {
    
    depth_colors <- c(depth_colors, four_colors[1])
    depth_bins <- c(depth_bins, "0-20")
    
  } else if ( as.numeric(entry) %>% between(21,50)) {
    
    depth_colors <- c(depth_colors, four_colors[2])
    depth_bins <- c(depth_bins, "21-50")
    
  } else if ( as.numeric(entry) %>% between(51,100)) {
    
    depth_colors <- c(depth_colors, four_colors[3])
    depth_bins <- c(depth_bins, "51-100")
    
  } else {
    
    depth_colors <- c(depth_colors, four_colors[4])
    depth_bins <- c(depth_bins, "101-220")
    
  }
  
}

depth_colors
depth_bins

#And now we can add those color columns to our mds_plot dataframe:
mds_plot$Depth_bin <- depth_bins
mds_plot

#Make plot coloured by core
ggplot(mds_plot, aes(x=MDS1, y=MDS2, color = Core)) + 
  scale_colour_manual(values=cbbPalette) +
  geom_point(aes(shape=Core), size=2) + ggtitle(paste0('stress: ', round(mds$stress,3))) +
  #ggrepel::geom_text_repel(aes(label=sample)) +  # if a lot of samples, comment this off
  coord_equal() + theme_bw()

ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mds_ASV_counts_core.pdf", width=10, height=6)

#Make plot coloured by core and depth
ggplot(mds_plot, aes(x=MDS1, y=MDS2, color = as.numeric(Depth))) + 
  scale_colour_gradient(high = "#132B43", low = "#56B1F7") +
  geom_point(aes(shape=Core),size=2) + ggtitle(paste0('stress: ', round(mds$stress,3))) +
  #ggrepel::geom_text_repel(aes(label=sample)) +  # if a lot of samples, comment this off
  coord_equal() + theme_bw()

ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mds_ASV_counts_core_depth_genus.pdf", width=10, height=6)

#Ordination
#Generally speaking, ordinations provide visualizations of sample-relatedness 
#based on dimension reduction. The ‘dimensions’ can be based on whatever we 
#measured in each sample, in our case counts of ASVs. Principle coordinates 
#analysis (PCoA) is a type of ordination that operates on dissimilarities or distances.

library(ggrepel)

# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows = TRUE)
# phyloseq expects our sample ID to be the row names in the sample_info_tab we give it,
# so making that change here where we are specifying the input table
sample_info_tab_phy <- sample_data(sample_info_tab_1 %>% column_to_rownames("Sample"))

vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues 
# pulling these out helps us to scale the axes according to their 
# magnitude of separating apart the samples

# colored by core plot
plot_ordination(vst_physeq, vst_pcoa, color = "Core") + 
  geom_point(size = 1) + 
  geom_text_repel(label = sample_info_tab$Sample) + 
  coord_fixed(sqrt(eigen_vals[2] / eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values = unique(sample_info_tab$core_colors[order(sample_info_tab$Core)])) + 
  theme_bw()

# it's very busy with the labels, so we can drop them by removing 
# the geom_text_repel line
plot_ordination(vst_physeq, vst_pcoa, color = "Core") + 
  geom_point(size = 1) + 
  coord_fixed(sqrt(eigen_vals[2] / eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values = unique(sample_info_tab$core_colors[order(sample_info_tab$Core)])) + 
  theme_bw()

# now colored by depth
plot_ordination(vst_physeq, vst_pcoa, color = "Depth_bin") + 
  geom_point(size = 1) + 
  coord_fixed(sqrt(eigen_vals[2] / eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values = unique(sample_info_tab$depth_colors[order(sample_info_tab$Depth)])) + 
  theme_bw()

#Alpha diversity
# Rarefaction curves
#Useful overview of diversity of sequences (not organisms)
# colored by core

# colored by depth
rarecurve(t(counts_tab), step = 100, col = sample_info_tab$depth_color, lwd = 2, ylab = "ASVs", label = FALSE)
# adding a legend
legend("bottomright", legend = c("0-20", "21-50", "51-100","101-205"), fill = four_colors,
       pt.cex = 0.5, title = "Depth")

#The interpretation of these plots is that the lower the ASVs, the lower the diversity

#Richness and diversity estimates
# Here we’re going to plot Chao1 richness estimates and Shannon diversity values.
# first we need to create a phyloseq object using our not-transformed count table
count_tab_phy <- otu_table(counts_tab, taxa_are_rows = TRUE)
tax_tab_phy <- tax_table(as.matrix(tax_tab))
data_x <- as.data.frame((tax_tab_phy))
unique(data_x$phylum)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

# Define functions from the microViz package that are useful for filtering

# internal helper that get phyloseq sample_data as plain dataframe
# without changing invalid colnames (like microbiome::meta does)
# or losing rownames / sample_names (like data.frame() with defaults does)
samdatAsDataframe <- function(ps) {
  samdat <- phyloseq::sample_data(ps)
  df <- data.frame(samdat, check.names = FALSE)
  return(df)
}
# Function that filters phyloseq samples by sample_data variables
ps_filter <- function(ps,
                      ...,
                      .target = "sample_data",
                      .keep_all_taxa = FALSE) {
  if (!inherits(ps, "phyloseq")) {
    stop("ps must be a phyloseq object. It is of class: ", class(ps))
  }
  
  if (!identical(.target, "sample_data")) {
    stop("Only .target = 'sample_data', has been implemented so far.")
  }
  # TODO: see if it is useful to facilitate
  # filtering by variables in other phyloseq slots
  
  df <- samdatAsDataframe(ps)
  df <- dplyr::filter(df, ...)
  phyloseq::sample_data(ps) <- df
  
  # remove taxa that now have zero counts (or relative abundance)
  # across all remaining samples
  if (isFALSE(.keep_all_taxa)) ps <- tax_filter_zeros(ps)
  return(ps)
}

# helper function used here and in ps_join
# removes all taxa which sum to zero across all samples
# (phyloseq::taxa_sums(ps) == 0)
# provides helpful warning if otu_table contains negative values
tax_filter_zeros <- function(ps) {
  # remove taxa that now have zero counts (or relative abundance)
  # across all remaining samples
  if (any(phyloseq::otu_table(ps) < 0)) {
    warning(
      "Removing taxa whose abundance across filtered samples is equal to zero.",
      "\nThis may not result in the desired outcome, ",
      "as some values in the otu_table are negative.",
      "\nAvoid performing transformations, ",
      "e.g. clr, before using `ps_filter()`, or set .keep_all_taxa = TRUE "
    )
  }
  return(phyloseq::prune_taxa(taxa = phyloseq::taxa_sums(ps) != 0, x = ps))
}

#We'll filter the phyloseq object by core
ASV_physeq_mc01 <- ps_filter(ASV_physeq, Core == 'MC01')
ASV_physeq_gc02 <- ps_filter(ASV_physeq, Core == 'GC02')
ASV_physeq_gc04 <- ps_filter(ASV_physeq, Core == 'GC04')
ASV_physeq_gc06 <- ps_filter(ASV_physeq, Core == 'GC06')
#ASV_physeq_gcnone <- ps_filter(ASV_physeq, Core == 'GCnone')

# Colorblind palette with black:
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#009E73", "#0072B2", "#D55E00", "#CC79A7")

# and now we can call the plot_richness() function on our phyloseq object
# doing one based on Core
p<-plot_richness(ASV_physeq, color = "Core", measures = c("Chao1", "Shannon")) + 
  scale_color_manual(values = unique(sample_info_tab$core_colors[order(sample_info_tab$Sample)])) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p+scale_colour_manual(values=cbbPalette)
#Numbers on the left in the Chao1 plot (numbers of unique ASVs)


# and now one colored based on Depth
plot_richness(ASV_physeq, color = "Depth_bin", measures = c("Chao1", "Shannon")) + 
  scale_color_manual(values = unique(sample_info_tab$depth_colors[order(sample_info_tab$Sample)])) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#We can also easily group them with phyloseq, here’s by Core:
plot_richness(ASV_physeq, x = "Core", color = "Core", measures = c("Chao1", "Shannon")) + 
  scale_colour_manual(values=cbbPalette) +
  #scale_color_manual(values = unique(sample_info_tab$core_color[order(sample_info_tab$Sample)])) +
  geom_point(aes(shape=Core), size=2) +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/alpha_diversity_core.pdf")

# And by depth
plot_richness(ASV_physeq, x = "Depth_bin", color = "depth_colors", measures = c("Chao1", "Shannon")) + 
  scale_colour_manual(values=cbbPalette) +
  #scale_color_manual(values = unique(sample_info_tab$depth_color[order(sample_info_tab$Sample)])) +
  geom_point(aes(shape=Core), size=2) +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
pdf(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/alpha_diversity_depth.pdf")

#Or by core and depth
plot_richness(ASV_physeq, x = "Core", measures = c("Chao1", "Shannon")) + 
  scale_colour_gradient(high = "#132B43", low = "#56B1F7") +
  #scale_color_manual(values = unique(sample_info_tab$core_color[order(sample_info_tab$Sample)])) +
  geom_point(aes(shape=Core, colour = as.numeric(Depth)), size=2) +
  theme_bw() + theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/alpha_div_core_depth.pdf", width=6, height=8, dpi=300)

#Making ridges plot
#get diversity data
df_richness<-estimate_richness(ASV_physeq, split = TRUE, measures = NULL)
#Create sample column
df_richness$sample<-rownames(df_richness)
#Creating GC and Depth columns in the df_richness df
df_richness <- data.frame(df_richness,do.call(rbind,str_split(df_richness$sample,"_")))
df_richness$X1 <- NULL
colnames(df_richness)[colnames(df_richness) == 'X2'] <- 'Core'
colnames(df_richness)[colnames(df_richness) == 'X3'] <- 'Depth'
head(df_richness)

#Make ridges plot based on Shannon diversity
library(ggridges)
library("ggplot2")
library("ggpp")
ggplot(df_richness, aes(x = Shannon, y = Core)) +
  #geom_density_ridges_gradient(scale = 0.7, rel_min_height = 0.005) +
  geom_density_ridges(
    aes(point_color = as.numeric(Depth), point_fill = as.numeric(Depth)),
    jittered_points = TRUE, alpha=0, point_alpha=1, scale = 0.9) +
  scale_point_colour_gradient(high = "#132B43", low = "#56B1F7")+
  geom_boxplot(aes(x = as.numeric(Shannon)+0.25, y = Core),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  geom_text(aes(label = Depth), vjust = -1, angle=45) + 
  theme_ridges(font_size = 13, line_size = 0.5, grid = TRUE, center_axis_labels = TRUE) 

ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/Shannon_ridges_labels.pdf", width=12, height=5, dpi=300)

#Taxonomic summaries
#Taxonomic summary figures usually aren’t all that exciting, but they can still be useful, 
#and they’re certainly still a reasonable way to present a summary of our data (all caveats included – 
#like recovered amplicon sequences do not equal organisms, or genomes, or cells)

#Make a summary of all major taxa proportions across all samples
#To start, we need to parse the count matrix by taxonomy.
# using phyloseq to make a count table that has summed all ASVs
# that were in the same phylum
phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank = "phylum")) 
head(phyla_counts_tab)
write.table(asv_tab, "phyla_counts_tab.tsv", sep = "\t", quote = FALSE, col.names = NA)
#for order
order_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank = "order"))
head(order_counts_tab)
write.table(asv_tab, "order_counts_tab.tsv", sep = "\t", quote = FALSE, col.names = NA)
#for genus
genus_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank = "genus"))
write.table(asv_tab, "genus_counts_tab.tsv", sep = "\t", quote = FALSE, col.names = NA)

#For Proteobacteria
ASV_proteo_physeq <- subset_taxa(ASV_physeq, phylum=="Proteobacteria")
#Make a count table by class
ASV_proteo_physeq_counts_tab <- otu_table(tax_glom(ASV_proteo_physeq, taxrank = "class"))
ASV_proteo_physeq_counts_tab_order <- otu_table(tax_glom(ASV_proteo_physeq, taxrank = "order"))

#For Firmicutes
ASV_firm_physeq <- subset_taxa(ASV_physeq, phylum=="Firmicutes")
#Make a count table by genus
ASV_firm_physeq_counts_tab <- otu_table(tax_glom(ASV_firm_physeq, taxrank = "genus"))

#For other sulfur transforming microbes
ASV_sulfurmicrobes_physeq <-subset_taxa(ASV_physeq, phylum=="Bacteroidota"
                                        | phylum=="Chloroflexi"
                                        | phylum=="Dadabacteria"
                                        | phylum=="Desulfobacterota_1"
                                        | phylum=="Desulfobacterota_2"
                                        | phylum=="Nitrospirota"
                                        | phylum=="Halobacterota"
                                        | phylum=="Actinobacteriota"
                                        | phylum=="Planctomycetota"
                                        | phylum=="Crenarchaeota")
#Make a count table by phyla
ASV_sulfurmicrobes_physeq_counts_tab <- otu_table(tax_glom(ASV_sulfurmicrobes_physeq, taxrank = "phylum"))
#Fuse both count tables
class_proteobacteria_counts_tab <- merge_phyloseq(ASV_proteo_physeq_counts_tab,ASV_sulfurmicrobes_physeq_counts_tab)
#Fuse both ASV tables
ASV_proteobacteria_physeq <- merge_phyloseq(ASV_proteo_physeq, ASV_sulfurmicrobes_physeq)

#We'll filter each phyloseq object by core 
#for Proteobacteria
ASV_proteobacteria_physeq_mc01 <- ps_filter(ASV_proteo_physeq, Core == 'MC01')
ASV_proteobacteria_physeq_gc02 <- ps_filter(ASV_proteo_physeq, Core == 'GC02')
ASV_proteobacteria_physeq_gc04 <- ps_filter(ASV_proteo_physeq, Core == 'GC04')
ASV_proteobacteria_physeq_gc06 <- ps_filter(ASV_proteo_physeq, Core == 'GC06')
#Make count table by class
class_counts_tab_proteo_mc01 <- otu_table(tax_glom(ASV_proteobacteria_physeq_mc01, taxrank = "class")) 
class_counts_tab_proteo_gc02 <- otu_table(tax_glom(ASV_proteobacteria_physeq_gc02, taxrank = "class")) 
class_counts_tab_proteo_gc04 <- otu_table(tax_glom(ASV_proteobacteria_physeq_gc04, taxrank = "class")) 
class_counts_tab_proteo_gc06 <- otu_table(tax_glom(ASV_proteobacteria_physeq_gc06, taxrank = "class"))
#Make count table by order
order_counts_tab_proteo_mc01 <- otu_table(tax_glom(ASV_proteobacteria_physeq_mc01, taxrank = "order")) 
order_counts_tab_proteo_gc02 <- otu_table(tax_glom(ASV_proteobacteria_physeq_gc02, taxrank = "order")) 
order_counts_tab_proteo_gc04 <- otu_table(tax_glom(ASV_proteobacteria_physeq_gc04, taxrank = "order")) 
order_counts_tab_proteo_gc06 <- otu_table(tax_glom(ASV_proteobacteria_physeq_gc06, taxrank = "order"))


#For sulfur transforming microbes
ASV_sulfurmicrobes_physeq_mc01 <- ps_filter(ASV_sulfurmicrobes_physeq, Core == 'MC01')
ASV_sulfurmicrobes_physeq_gc02 <- ps_filter(ASV_sulfurmicrobes_physeq, Core == 'GC02')
ASV_sulfurmicrobes_physeq_gc04 <- ps_filter(ASV_sulfurmicrobes_physeq, Core == 'GC04')
ASV_sulfurmicrobes_physeq_gc06 <- ps_filter(ASV_sulfurmicrobes_physeq, Core == 'GC06')

#Make count table phylum
phylum_counts_tab_sulfurmicrobes_mc01 <- otu_table(tax_glom(ASV_sulfurmicrobes_physeq_mc01, taxrank = "phylum")) 
phylum_counts_tab_sulfurmicrobes_gc02 <- otu_table(tax_glom(ASV_sulfurmicrobes_physeq_gc02, taxrank = "phylum")) 
phylum_counts_tab_sulfurmicrobes_gc04 <- otu_table(tax_glom(ASV_sulfurmicrobes_physeq_gc04, taxrank = "phylum")) 
phylum_counts_tab_sulfurmicrobes_gc06 <- otu_table(tax_glom(ASV_sulfurmicrobes_physeq_gc06, taxrank = "phylum"))

#Fuse both count tables
class_counts_tab_proteobacteria_mc01 <- merge_phyloseq(class_counts_tab_proteo_mc01,phylum_counts_tab_sulfurmicrobes_mc01)
class_counts_tab_proteobacteria_gc02 <- merge_phyloseq(class_counts_tab_proteo_gc02,phylum_counts_tab_sulfurmicrobes_gc02)
class_counts_tab_proteobacteria_gc04 <- merge_phyloseq(class_counts_tab_proteo_gc04,phylum_counts_tab_sulfurmicrobes_gc04)
class_counts_tab_proteobacteria_gc06 <- merge_phyloseq(class_counts_tab_proteo_gc06,phylum_counts_tab_sulfurmicrobes_gc06)

#For Firmicutes
ASV_firm_physeq_mc01 <- ps_filter(ASV_firm_physeq, Core == 'MC01')
ASV_firm_physeq_gc02 <- ps_filter(ASV_firm_physeq, Core == 'GC02')
ASV_firm_physeq_gc04 <- ps_filter(ASV_firm_physeq, Core == 'GC04')
ASV_firm_physeq_gc06 <- ps_filter(ASV_firm_physeq, Core == 'GC06')
#Make count table by genus
genus_counts_tab_firm_mc01 <- otu_table(tax_glom(ASV_firm_physeq_mc01, taxrank = "genus")) 
genus_counts_tab_firm_gc02 <- otu_table(tax_glom(ASV_firm_physeq_gc02, taxrank = "genus")) 
genus_counts_tab_firm_gc04 <- otu_table(tax_glom(ASV_firm_physeq_gc04, taxrank = "genus")) 
genus_counts_tab_firm_gc06 <- otu_table(tax_glom(ASV_firm_physeq_gc06, taxrank = "genus"))


#Make a vector of names to set as rownames
#For classes of Proteobacteria 
class_tax_vec_proteo_mc01 <- as.vector(tax_table(tax_glom(ASV_proteobacteria_physeq_mc01, taxrank = "class"))[,"class"])
class_tax_vec_proteo_gc02 <- as.vector(tax_table(tax_glom(ASV_proteobacteria_physeq_gc02, taxrank = "class"))[,"class"])
class_tax_vec_proteo_gc04 <- as.vector(tax_table(tax_glom(ASV_proteobacteria_physeq_gc04, taxrank = "class"))[,"class"])
class_tax_vec_proteo_gc06 <- as.vector(tax_table(tax_glom(ASV_proteobacteria_physeq_gc06, taxrank = "class"))[,"class"])
#For orders of proteobacteria
order_tax_vec_proteo_mc01 <- as.vector(tax_table(tax_glom(ASV_proteobacteria_physeq_mc01, taxrank = "order"))[,"order"])
order_tax_vec_proteo_gc02 <- as.vector(tax_table(tax_glom(ASV_proteobacteria_physeq_gc02, taxrank = "order"))[,"order"])
order_tax_vec_proteo_gc04 <- as.vector(tax_table(tax_glom(ASV_proteobacteria_physeq_gc04, taxrank = "order"))[,"order"])
order_tax_vec_proteo_gc06 <- as.vector(tax_table(tax_glom(ASV_proteobacteria_physeq_gc06, taxrank = "order"))[,"order"])

# For phyla of sulfur transforming microbes
phylum_tax_vec_sulfurmicrobes_mc01 <- as.vector(tax_table(tax_glom(ASV_sulfurmicrobes_physeq_mc01, taxrank = "phylum"))[,"phylum"])
phylum_tax_vec_sulfurmicrobes_gc02 <- as.vector(tax_table(tax_glom(ASV_sulfurmicrobes_physeq_gc02, taxrank = "phylum"))[,"phylum"])
phylum_tax_vec_sulfurmicrobes_gc04 <- as.vector(tax_table(tax_glom(ASV_sulfurmicrobes_physeq_gc04, taxrank = "phylum"))[,"phylum"])
phylum_tax_vec_sulfurmicrobes_gc06 <- as.vector(tax_table(tax_glom(ASV_sulfurmicrobes_physeq_gc06, taxrank = "phylum"))[,"phylum"])
#Fuse both vectors
class_tax_vec_proteobacteria_mc01 <- c(class_tax_vec_proteo_mc01,phylum_tax_vec_sulfurmicrobes_mc01)
class_tax_vec_proteobacteria_gc02 <- c(class_tax_vec_proteo_gc02,phylum_tax_vec_sulfurmicrobes_gc02)
class_tax_vec_proteobacteria_gc04 <- c(class_tax_vec_proteo_gc04,phylum_tax_vec_sulfurmicrobes_gc04)
class_tax_vec_proteobacteria_gc06 <- c(class_tax_vec_proteo_gc06,phylum_tax_vec_sulfurmicrobes_gc06)

#For genuses of Firmicutes
genus_tax_vec_firm_mc01 <- as.vector(tax_table(tax_glom(ASV_firm_physeq_mc01, taxrank = "genus"))[,"genus"])
genus_tax_vec_firm_gc02 <- as.vector(tax_table(tax_glom(ASV_firm_physeq_gc02, taxrank = "genus"))[,"genus"])
genus_tax_vec_firm_gc04 <- as.vector(tax_table(tax_glom(ASV_firm_physeq_gc04, taxrank = "genus"))[,"genus"])
genus_tax_vec_firm_gc06 <- as.vector(tax_table(tax_glom(ASV_firm_physeq_gc06, taxrank = "genus"))[,"genus"])

#Replace rownames in the table of each core
#For proteobacteria and sulfur transforming microbes
rownames(class_counts_tab_proteobacteria_mc01) <- class_tax_vec_proteobacteria_mc01
rownames(class_counts_tab_proteobacteria_gc02) <- class_tax_vec_proteobacteria_gc02
rownames(class_counts_tab_proteobacteria_gc04) <- class_tax_vec_proteobacteria_gc04
rownames(class_counts_tab_proteobacteria_gc06) <- class_tax_vec_proteobacteria_gc06
#For orders of Proteobacteria
rownames(order_counts_tab_proteo_mc01) <- order_tax_vec_proteo_mc01
rownames(order_counts_tab_proteo_gc02) <- order_tax_vec_proteo_gc02
rownames(order_counts_tab_proteo_gc04) <- order_tax_vec_proteo_gc04
rownames(order_counts_tab_proteo_gc06) <- order_tax_vec_proteo_gc06
#For genuses of Firmicutes
rownames(genus_counts_tab_firm_mc01) <- genus_tax_vec_firm_mc01
rownames(genus_counts_tab_firm_gc02) <- genus_tax_vec_firm_gc02
rownames(genus_counts_tab_firm_gc04) <- genus_tax_vec_firm_gc04
rownames(genus_counts_tab_firm_gc06) <- genus_tax_vec_firm_gc06

#We'll filter the phyloseq object by core
ASV_physeq_mc01 <- ps_filter(ASV_physeq, Core == 'MC01')
ASV_physeq_gc02 <- ps_filter(ASV_physeq, Core == 'GC02')
ASV_physeq_gc04 <- ps_filter(ASV_physeq, Core == 'GC04')
ASV_physeq_gc06 <- ps_filter(ASV_physeq, Core == 'GC06')

#We'll create a table for each core
#For phylum
phyla_counts_tab_mc01 <- otu_table(tax_glom(ASV_physeq_mc01, taxrank = "phylum")) 
phyla_counts_tab_gc02 <- otu_table(tax_glom(ASV_physeq_gc02, taxrank = "phylum")) 
phyla_counts_tab_gc04 <- otu_table(tax_glom(ASV_physeq_gc04, taxrank = "phylum")) 
phyla_counts_tab_gc06 <- otu_table(tax_glom(ASV_physeq_gc06, taxrank = "phylum")) 
#phyla_counts_tab_gcnone <- otu_table(tax_glom(ASV_physeq_gcnone, taxrank = "phylum")) 
#For order
order_counts_tab_mc01 <- otu_table(tax_glom(ASV_physeq_mc01, taxrank = "order")) 
order_counts_tab_gc02 <- otu_table(tax_glom(ASV_physeq_gc02, taxrank = "order")) 
order_counts_tab_gc04 <- otu_table(tax_glom(ASV_physeq_gc04, taxrank = "order")) 
order_counts_tab_gc06 <- otu_table(tax_glom(ASV_physeq_gc06, taxrank = "order")) 
#family_counts_tab_gcnone <- otu_table(tax_glom(ASV_physeq_gcnone, taxrank = "family")) 
#For genus
genus_counts_tab_mc01 <- otu_table(tax_glom(ASV_physeq_mc01, taxrank = "genus")) 
genus_counts_tab_gc02 <- otu_table(tax_glom(ASV_physeq_gc02, taxrank = "genus")) 
genus_counts_tab_gc04 <- otu_table(tax_glom(ASV_physeq_gc04, taxrank = "genus")) 
genus_counts_tab_gc06 <- otu_table(tax_glom(ASV_physeq_gc06, taxrank = "genus"))
#genus_counts_tab_gcnone <- otu_table(tax_glom(ASV_physeq_gcnone, taxrank = "genus")) 

# That kept some representative sequence headers as the row names though
head(phyla_counts_tab[ ,1:5])
head(phyla_counts_tab_gc02[ ,1:5])
head(family_counts_tab[ ,1:5])
head(genus_counts_tab[ ,1:5])

# Make a vector of phyla names to set as row names next
#For phylum
phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank = "phylum"))[,"phylum"])
#For order
order_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank = "order"))[,"order"])
#For genus
genus_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank = "genus"))[,"genus"])
#Same for each core
#for phylum
phyla_tax_vec_mc01 <- as.vector(tax_table(tax_glom(ASV_physeq_mc01, taxrank = "phylum"))[,"phylum"])
phyla_tax_vec_gc02 <- as.vector(tax_table(tax_glom(ASV_physeq_gc02, taxrank = "phylum"))[,"phylum"])
phyla_tax_vec_gc04 <- as.vector(tax_table(tax_glom(ASV_physeq_gc04, taxrank = "phylum"))[,"phylum"])
phyla_tax_vec_gc06 <- as.vector(tax_table(tax_glom(ASV_physeq_gc06, taxrank = "phylum"))[,"phylum"])
#phyla_tax_vec_gcnone <- as.vector(tax_table(tax_glom(ASV_physeq_gcnone, taxrank = "phylum"))[,"phylum"])
#for order
order_tax_vec_mc01 <- as.vector(tax_table(tax_glom(ASV_physeq_mc01, taxrank = "order"))[,"order"])
order_tax_vec_gc02 <- as.vector(tax_table(tax_glom(ASV_physeq_gc02, taxrank = "order"))[,"order"])
order_tax_vec_gc04 <- as.vector(tax_table(tax_glom(ASV_physeq_gc04, taxrank = "order"))[,"order"])
order_tax_vec_gc06 <- as.vector(tax_table(tax_glom(ASV_physeq_gc06, taxrank = "order"))[,"order"])
#family_tax_vec_gcnone <- as.vector(tax_table(tax_glom(ASV_physeq_gcnone, taxrank = "family"))[,"family"])

#For genus
genus_tax_vec_mc01 <- as.vector(tax_table(tax_glom(ASV_physeq_mc01, taxrank = "genus"))[,"genus"])
genus_tax_vec_gc02 <- as.vector(tax_table(tax_glom(ASV_physeq_gc02, taxrank = "genus"))[,"genus"])
genus_tax_vec_gc04 <- as.vector(tax_table(tax_glom(ASV_physeq_gc04, taxrank = "genus"))[,"genus"])
genus_tax_vec_gc06 <- as.vector(tax_table(tax_glom(ASV_physeq_gc06, taxrank = "genus"))[,"genus"])
#genus_tax_vec_gcnone <- as.vector(tax_table(tax_glom(ASV_physeq_gcnone, taxrank = "genus"))[,"genus"])


#Replace row names in the general table
#for phyla
rownames(phyla_counts_tab) <- phyla_tax_vec
head(phyla_counts_tab[ ,1:5])
#for order
rownames(family_counts_tab) <- order_tax_vec
head(order_counts_tab[ ,1:5])
#for genus
rownames(genus_counts_tab) <- genus_tax_vec
head(genus_counts_tab[ ,1:5])

#Replace row names in the table of each core
#For phylum
rownames(phyla_counts_tab_mc01) <- phyla_tax_vec_mc01
rownames(phyla_counts_tab_gc02) <- phyla_tax_vec_gc02
rownames(phyla_counts_tab_gc04) <- phyla_tax_vec_gc04
rownames(phyla_counts_tab_gc06) <- phyla_tax_vec_gc06
#rownames(phyla_counts_tab_gcnone) <- phyla_tax_vec_gcnone
#For order
rownames(order_counts_tab_mc01) <- order_tax_vec_mc01
rownames(order_counts_tab_gc02) <- order_tax_vec_gc02
rownames(order_counts_tab_gc04) <- order_tax_vec_gc04
rownames(order_counts_tab_gc06) <- order_tax_vec_gc06
#rownames(family_counts_tab_gcnone) <- family_tax_vec_gcnone
#For genus
rownames(genus_counts_tab_mc01) <- genus_tax_vec_mc01
rownames(genus_counts_tab_gc02) <- genus_tax_vec_gc02
rownames(genus_counts_tab_gc04) <- genus_tax_vec_gc04
rownames(genus_counts_tab_gc06) <- genus_tax_vec_gc06
#rownames(genus_counts_tab_gcnone) <- genus_tax_vec_gcnone

#Rearrange columns in the tables in order of depth
#For phylum
phyla_counts_tab_gc02 <-phyla_counts_tab_gc02[,order(as.numeric(sapply(strsplit(colnames(phyla_counts_tab_gc02), "_"), "[[", 3)))]
phyla_counts_tab_gc04 <-phyla_counts_tab_gc04[,order(as.numeric(sapply(strsplit(colnames(phyla_counts_tab_gc04), "_"), "[[", 3)))]
phyla_counts_tab_gc06 <-phyla_counts_tab_gc06[,order(as.numeric(sapply(strsplit(colnames(phyla_counts_tab_gc06), "_"), "[[", 3)))]
phyla_counts_tab_mc01 <-phyla_counts_tab_mc01[,order(as.numeric(sapply(strsplit(colnames(phyla_counts_tab_mc01), "_"), "[[", 3)))]
#phyla_counts_tab_gcnone <-phyla_counts_tab_gcnone[,order(as.numeric(sapply(strsplit(colnames(phyla_counts_tab_gcnone), "_"), "[[", 3)))]

#For order
order_counts_tab_gc02 <-order_counts_tab_gc02[,order(as.numeric(sapply(strsplit(colnames(order_counts_tab_gc02), "_"), "[[", 3)))]
order_counts_tab_gc04 <-order_counts_tab_gc04[,order(as.numeric(sapply(strsplit(colnames(order_counts_tab_gc04), "_"), "[[", 3)))]
order_counts_tab_gc06 <-order_counts_tab_gc06[,order(as.numeric(sapply(strsplit(colnames(order_counts_tab_gc06), "_"), "[[", 3)))]
order_counts_tab_mc01 <-order_counts_tab_mc01[,order(as.numeric(sapply(strsplit(colnames(order_counts_tab_mc01), "_"), "[[", 3)))]
colnames(order_counts_tab_mc01)
#family_counts_tab_gcnone <-family_counts_tab_gcnone[,order(as.numeric(sapply(strsplit(colnames(family_counts_tab_gcnone), "_"), "[[", 3)))]

#For genus
genus_counts_tab_gc02 <-genus_counts_tab_gc02[,order(as.numeric(sapply(strsplit(colnames(genus_counts_tab_gc02), "_"), "[[", 3)))]
genus_counts_tab_gc04 <-genus_counts_tab_gc04[,order(as.numeric(sapply(strsplit(colnames(genus_counts_tab_gc04), "_"), "[[", 3)))]
genus_counts_tab_gc06 <-genus_counts_tab_gc06[,order(as.numeric(sapply(strsplit(colnames(genus_counts_tab_gc06), "_"), "[[", 3)))]
genus_counts_tab_mc01 <-genus_counts_tab_mc01[,order(as.numeric(sapply(strsplit(colnames(genus_counts_tab_mc01), "_"), "[[", 3)))]
#genus_counts_tab_gcnone <-genus_counts_tab_gcnone[,order(as.numeric(sapply(strsplit(colnames(genus_counts_tab_gcnone), "_"), "[[", 3)))]

#For classes of proteobacteria and phyla of sulfur transforming microbes
class_counts_tab_proteobacteria_mc01 <-class_counts_tab_proteobacteria_mc01[,order(as.numeric(sapply(strsplit(colnames(class_counts_tab_proteobacteria_mc01), "_"), "[[", 3)))]
class_counts_tab_proteobacteria_gc02 <-class_counts_tab_proteobacteria_gc02[,order(as.numeric(sapply(strsplit(colnames(class_counts_tab_proteobacteria_gc02), "_"), "[[", 3)))]
class_counts_tab_proteobacteria_gc04 <-class_counts_tab_proteobacteria_gc04[,order(as.numeric(sapply(strsplit(colnames(class_counts_tab_proteobacteria_gc04), "_"), "[[", 3)))]
class_counts_tab_proteobacteria_gc06 <-class_counts_tab_proteobacteria_gc06[,order(as.numeric(sapply(strsplit(colnames(class_counts_tab_proteobacteria_gc06), "_"), "[[", 3)))]

#For orders of proteobacteria 
order_counts_tab_proteo_mc01 <-order_counts_tab_proteo_mc01[,order(as.numeric(sapply(strsplit(colnames(order_counts_tab_proteo_mc01), "_"), "[[", 3)))]
order_counts_tab_proteo_gc02 <-order_counts_tab_proteo_gc02[,order(as.numeric(sapply(strsplit(colnames(order_counts_tab_proteo_gc02), "_"), "[[", 3)))]
order_counts_tab_proteo_gc04 <-order_counts_tab_proteo_gc04[,order(as.numeric(sapply(strsplit(colnames(order_counts_tab_proteo_gc04), "_"), "[[", 3)))]
order_counts_tab_proteo_gc06 <-order_counts_tab_proteo_gc06[,order(as.numeric(sapply(strsplit(colnames(order_counts_tab_proteo_gc06), "_"), "[[", 3)))]

#For genuses of Firmicutes
genus_counts_tab_firm_mc01 <-genus_counts_tab_firm_mc01[,order(as.numeric(sapply(strsplit(colnames(genus_counts_tab_firm_mc01), "_"), "[[", 3)))]
genus_counts_tab_firm_gc02 <-genus_counts_tab_firm_gc02[,order(as.numeric(sapply(strsplit(colnames(genus_counts_tab_firm_gc02), "_"), "[[", 3)))]
genus_counts_tab_firm_gc04 <-genus_counts_tab_firm_gc04[,order(as.numeric(sapply(strsplit(colnames(genus_counts_tab_firm_gc04), "_"), "[[", 3)))]
genus_counts_tab_firm_gc06 <-genus_counts_tab_firm_gc06[,order(as.numeric(sapply(strsplit(colnames(genus_counts_tab_firm_gc06), "_"), "[[", 3)))]

#Account for sequences that were not assigned any taxonomy at the phylum level
#For phylum
unclassified_tax_counts <- colSums(counts_tab) - colSums(phyla_counts_tab)
#For order
unclassified_tax_counts_o <- colSums(counts_tab) - colSums(order_counts_tab)
#For genus
unclassified_tax_counts_g <- colSums(counts_tab) - colSums(genus_counts_tab)

# Add the unclassifieds as a row to the phylum count table:
#For phylum
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified" = unclassified_tax_counts)
#For order
order_and_unidentified_counts_tab <- rbind(order_counts_tab, "Unclassified" = unclassified_tax_counts_o)
#For genus
genus_and_unidentified_counts_tab <- rbind(genus_counts_tab, "Unclassified" = unclassified_tax_counts_g)
#For each core
#Make count tables for each core
counts_tab_mc01 <- counts_tab[,grepl("MC01", names(counts_tab))]
counts_tab_gc02 <- counts_tab[,grepl("GC02", names(counts_tab))]
counts_tab_gc04 <- counts_tab[,grepl("GC04", names(counts_tab))]
counts_tab_gc06 <- counts_tab[,grepl("GC06", names(counts_tab))]

#Make count tables for each core for Proteobacteria
proteobacteria_counts_tab_mc01 <- proteobacteria_counts_tab[,grepl("MC01", names(proteobacteria_counts_tab))]
proteobacteria_counts_tab_gc02 <- proteobacteria_counts_tab[,grepl("GC02", names(proteobacteria_counts_tab))]
proteobacteria_counts_tab_gc04 <- proteobacteria_counts_tab[,grepl("GC04", names(proteobacteria_counts_tab))]
proteobacteria_counts_tab_gc06 <- proteobacteria_counts_tab[,grepl("GC06", names(proteobacteria_counts_tab))]

#Rearrange columns in the count tables in order of depth
counts_tab_gc02 <-counts_tab_gc02[,order(as.numeric(sapply(strsplit(colnames(counts_tab_gc02), "_"), "[[", 3)))]
counts_tab_gc04 <-counts_tab_gc04[,order(as.numeric(sapply(strsplit(colnames(counts_tab_gc04), "_"), "[[", 3)))]
counts_tab_gc06 <-counts_tab_gc06[,order(as.numeric(sapply(strsplit(colnames(counts_tab_gc06), "_"), "[[", 3)))]
counts_tab_mc01 <-counts_tab_mc01[,order(as.numeric(sapply(strsplit(colnames(counts_tab_mc01), "_"), "[[", 3)))]

#Rearrange columns in the count tables in order of depth for proteobacteria
proteobacteria_counts_tab_gc02 <-proteobacteria_counts_tab_gc02[,order(as.numeric(sapply(strsplit(colnames(proteobacteria_counts_tab_gc02), "_"), "[[", 3)))]
proteobacteria_counts_tab_gc04 <-proteobacteria_counts_tab_gc04[,order(as.numeric(sapply(strsplit(colnames(proteobacteria_counts_tab_gc04), "_"), "[[", 3)))]
proteobacteria_counts_tab_gc06 <-proteobacteria_counts_tab_gc06[,order(as.numeric(sapply(strsplit(colnames(proteobacteria_counts_tab_gc06), "_"), "[[", 3)))]
proteobacteria_counts_tab_mc01 <-proteobacteria_counts_tab_mc01[,order(as.numeric(sapply(strsplit(colnames(proteobacteria_counts_tab_mc01), "_"), "[[", 3)))]
head(proteobacteria_counts_tab_mc01)
#Account for sequences that were not assigned any taxonomy at the phylum level
unclassified_tax_counts_mc01 <- colSums(counts_tab_mc01) - colSums(phyla_counts_tab_mc01)
unclassified_tax_counts_gc02 <- colSums(counts_tab_gc02) - colSums(phyla_counts_tab_gc02)
unclassified_tax_counts_gc04 <- colSums(counts_tab_gc04) - colSums(phyla_counts_tab_gc04)
unclassified_tax_counts_gc06 <- colSums(counts_tab_gc06) - colSums(phyla_counts_tab_gc06)
#unclassified_tax_counts_gcnone <- colSums(counts_tab_gcnone) - colSums(phyla_counts_tab_gcnone)

#Account for sequences that were not assigned any taxonomy at the order level
unclassified_tax_counts_f_mc01 <- colSums(counts_tab_mc01) - colSums(order_counts_tab_mc01)
unclassified_tax_counts_f_gc02 <- colSums(counts_tab_gc02) - colSums(order_counts_tab_gc02)
unclassified_tax_counts_f_gc04 <- colSums(counts_tab_gc04) - colSums(order_counts_tab_gc04)
unclassified_tax_counts_f_gc06 <- colSums(counts_tab_gc06) - colSums(order_counts_tab_gc06)
#unclassified_tax_counts_f_gcnone <- colSums(counts_tab_gcnone) - colSums(family_counts_tab_gcnone)

#Account for sequences that were not assigned any taxonomy at the genus level
unclassified_tax_counts_g_mc01 <- colSums(counts_tab_mc01) - colSums(genus_counts_tab_mc01)
unclassified_tax_counts_g_gc02 <- colSums(counts_tab_gc02) - colSums(genus_counts_tab_gc02)
unclassified_tax_counts_g_gc04 <- colSums(counts_tab_gc04) - colSums(genus_counts_tab_gc04)
unclassified_tax_counts_g_gc06 <- colSums(counts_tab_gc06) - colSums(genus_counts_tab_gc06)
#unclassified_tax_counts_g_gcnone <- colSums(counts_tab_gcnone) - colSums(genus_counts_tab_gcnone)

#Account for sequences that were not assigned any taxonomy at the class level for proteobacteria and phyla for sulfur transforming microbes
proteobacteria_unclassified_tax_counts_g_mc01 <- colSums(proteobacteria_counts_tab_mc01) - colSums(class_counts_tab_proteobacteria_mc01)
proteobacteria_unclassified_tax_counts_g_gc02 <- colSums(proteobacteria_counts_tab_gc02) - colSums(class_counts_tab_proteobacteria_gc02)
proteobacteria_unclassified_tax_counts_g_gc04 <- colSums(proteobacteria_counts_tab_gc04) - colSums(class_counts_tab_proteobacteria_gc04)
proteobacteria_unclassified_tax_counts_g_gc06 <- colSums(proteobacteria_counts_tab_gc06) - colSums(class_counts_tab_proteobacteria_gc06)

#Account for sequences that were not assigned any taxonomy at the order level for proteobacteria 
proteobacteria_unclassified_tax_counts_o_mc01 <- colSums(proteobacteria_counts_tab_mc01) - colSums(order_counts_tab_proteo_mc01)
proteobacteria_unclassified_tax_counts_o_gc02 <- colSums(proteobacteria_counts_tab_gc02) - colSums(order_counts_tab_proteo_gc02)
proteobacteria_unclassified_tax_counts_o_gc04 <- colSums(proteobacteria_counts_tab_gc04) - colSums(order_counts_tab_proteo_gc04)
proteobacteria_unclassified_tax_counts_o_gc06 <- colSums(proteobacteria_counts_tab_gc06) - colSums(order_counts_tab_proteo_gc06)

#Account for sequences that were not assigned any taxonomy at the genus level for Firmicutes 
#Didn't do it. Not necessary.

# Add the unclassifieds as a row to the phylum count table:
phyla_and_unidentified_counts_tab_mc01 <- rbind(phyla_counts_tab_mc01, "Unclassified" = unclassified_tax_counts_mc01)
phyla_and_unidentified_counts_tab_gc02 <- rbind(phyla_counts_tab_gc02, "Unclassified" = unclassified_tax_counts_gc02)
phyla_and_unidentified_counts_tab_gc04 <- rbind(phyla_counts_tab_gc04, "Unclassified" = unclassified_tax_counts_gc04)
phyla_and_unidentified_counts_tab_gc06 <- rbind(phyla_counts_tab_gc06, "Unclassified" = unclassified_tax_counts_gc06)
#phyla_and_unidentified_counts_tab_gcnone <- rbind(phyla_counts_tab_gcnone, "Unclassified" = unclassified_tax_counts_gcnone)

# Add the unclassifieds as a row to the order count table:
order_and_unidentified_counts_tab_mc01 <- rbind(order_counts_tab_mc01, "Unclassified" = unclassified_tax_counts_f_mc01)
order_and_unidentified_counts_tab_gc02 <- rbind(order_counts_tab_gc02, "Unclassified" = unclassified_tax_counts_f_gc02)
order_and_unidentified_counts_tab_gc04 <- rbind(order_counts_tab_gc04, "Unclassified" = unclassified_tax_counts_f_gc04)
order_and_unidentified_counts_tab_gc06 <- rbind(order_counts_tab_gc06, "Unclassified" = unclassified_tax_counts_f_gc06)
#order_and_unidentified_counts_tab_gcnone <- rbind(order_counts_tab_gcnone, "Unclassified" = unclassified_tax_counts_f_gcnone)

# Add the unclassifieds as a row to the genus count table:
genus_and_unidentified_counts_tab_mc01 <- rbind(genus_counts_tab_mc01, "Unclassified" = unclassified_tax_counts_g_mc01)
genus_and_unidentified_counts_tab_gc02 <- rbind(genus_counts_tab_gc02, "Unclassified" = unclassified_tax_counts_g_gc02)
genus_and_unidentified_counts_tab_gc04 <- rbind(genus_counts_tab_gc04, "Unclassified" = unclassified_tax_counts_g_gc04)
genus_and_unidentified_counts_tab_gc06 <- rbind(genus_counts_tab_gc06, "Unclassified" = unclassified_tax_counts_g_gc06)
#genus_and_unidentified_counts_tab_gcnone <- rbind(genus_counts_tab_gcnone, "Unclassified" = unclassified_tax_counts_g_gcnone)

# Add the unclassifieds as a row to the class count table for classes of proteobacteria and sulfur transforming microbes
class_and_unidentified_counts_tab_proteobacteria_mc01 <- rbind(class_counts_tab_proteobacteria_mc01, "Unclassified" = proteobacteria_unclassified_tax_counts_g_mc01)
class_and_unidentified_counts_tab_proteobacteria_gc02 <- rbind(class_counts_tab_proteobacteria_gc02, "Unclassified" = proteobacteria_unclassified_tax_counts_g_gc02)
class_and_unidentified_counts_tab_proteobacteria_gc04 <- rbind(class_counts_tab_proteobacteria_gc04, "Unclassified" = proteobacteria_unclassified_tax_counts_g_gc04)
class_and_unidentified_counts_tab_proteobacteria_gc06 <- rbind(class_counts_tab_proteobacteria_gc06, "Unclassified" = proteobacteria_unclassified_tax_counts_g_gc06)
#genus_and_unidentified_counts_tab_gcnone <- rbind(genus_counts_tab_gcnone, "Unclassified" = unclassified_tax_counts_g_gcnone)

# Add the unclassifieds as a row to the class count table for orders of proteobacteria 
order_and_unidentified_counts_tab_proteobacteria_mc01 <- rbind(order_counts_tab_proteo_mc01, "Unclassified" = proteobacteria_unclassified_tax_counts_o_mc01)
order_and_unidentified_counts_tab_proteobacteria_gc02 <- rbind(order_counts_tab_proteo_gc02, "Unclassified" = proteobacteria_unclassified_tax_counts_o_gc02)
order_and_unidentified_counts_tab_proteobacteria_gc04 <- rbind(order_counts_tab_proteo_gc04, "Unclassified" = proteobacteria_unclassified_tax_counts_o_gc04)
order_and_unidentified_counts_tab_proteobacteria_gc06 <- rbind(order_counts_tab_proteo_gc06, "Unclassified" = proteobacteria_unclassified_tax_counts_o_gc06)
#genus_and_unidentified_counts_tab_gcnone <- rbind(genus_counts_tab_gcnone, "Unclassified" = unclassified_tax_counts_g_gcnone)

# To check we didn't miss any other sequences somehow, we can compare the column
# sums to see if they are the same
# if "TRUE", we know nothing fell through the cracks
#For phylum
identical(colSums(phyla_and_unidentified_counts_tab), colSums(counts_tab)) 
identical(colSums(phyla_and_unidentified_counts_tab_mc01), colSums(counts_tab_mc01)) 
identical(colSums(phyla_and_unidentified_counts_tab_gc02), colSums(counts_tab_gc02)) 
identical(colSums(phyla_and_unidentified_counts_tab_gc04), colSums(counts_tab_gc04)) 
identical(colSums(phyla_and_unidentified_counts_tab_gc06), colSums(counts_tab_gc06)) 
#For order
identical(colSums(order_and_unidentified_counts_tab), colSums(counts_tab)) 
identical(colSums(order_and_unidentified_counts_tab_mc01), colSums(counts_tab_mc01)) 
identical(colSums(order_and_unidentified_counts_tab_gc02), colSums(counts_tab_gc02)) 
identical(colSums(order_and_unidentified_counts_tab_gc04), colSums(counts_tab_gc04)) 
identical(colSums(order_and_unidentified_counts_tab_gc06), colSums(counts_tab_gc06)) 
#For genus
identical(colSums(genus_and_unidentified_counts_tab), colSums(counts_tab)) 
identical(colSums(genus_and_unidentified_counts_tab_mc01), colSums(counts_tab_mc01)) 
identical(colSums(genus_and_unidentified_counts_tab_gc02), colSums(counts_tab_gc02)) 
identical(colSums(genus_and_unidentified_counts_tab_gc04), colSums(counts_tab_gc04)) 
identical(colSums(genus_and_unidentified_counts_tab_gc06), colSums(counts_tab_gc06)) 

#For classes of proteobacteria and phyla of sulfur transforming microbes
identical(colSums(class_and_unidentified_counts_tab_proteobacteria_mc01), colSums(proteobacteria_counts_tab_mc01)) 
identical(colSums(class_and_unidentified_counts_tab_proteobacteria_gc02), colSums(proteobacteria_counts_tab_gc02)) 
identical(colSums(class_and_unidentified_counts_tab_proteobacteria_gc04), colSums(proteobacteria_counts_tab_gc04)) 
identical(colSums(class_and_unidentified_counts_tab_proteobacteria_gc06), colSums(proteobacteria_counts_tab_gc06)) 

#For orders of proteobacteria 
identical(colSums(order_and_unidentified_counts_tab_proteobacteria_mc01), colSums(proteobacteria_counts_tab_mc01)) 
identical(colSums(order_and_unidentified_counts_tab_proteobacteria_gc02), colSums(proteobacteria_counts_tab_gc02)) 
identical(colSums(order_and_unidentified_counts_tab_proteobacteria_gc04), colSums(proteobacteria_counts_tab_gc04)) 
identical(colSums(order_and_unidentified_counts_tab_proteobacteria_gc06), colSums(proteobacteria_counts_tab_gc06)) 

# now we'll generate a percent table for summarizing:
#For phylum
major_taxa_proportions_tab <- apply(phyla_and_unidentified_counts_tab, 2, function(x) x/sum(x)*100)
major_taxa_proportions_tab_mc01 <- apply(phyla_and_unidentified_counts_tab_mc01, 2, function(x) x/sum(x)*100)
major_taxa_proportions_tab_gc02 <- apply(phyla_and_unidentified_counts_tab_gc02, 2, function(x) x/sum(x)*100)
major_taxa_proportions_tab_gc04 <- apply(phyla_and_unidentified_counts_tab_gc04, 2, function(x) x/sum(x)*100)
major_taxa_proportions_tab_gc06 <- apply(phyla_and_unidentified_counts_tab_gc06, 2, function(x) x/sum(x)*100)
#major_taxa_proportions_tab_gcnone <- apply(phyla_and_unidentified_counts_tab_gcnone, 2, function(x) x/sum(x)*100)
head(major_taxa_proportions_tab_mc01)
major_taxa_proportions_tab_mc01[(rownames=='Firmicutes'),]
data_mod <- major_taxa_proportions_tab_mc01[rownames(major_taxa_proportions_tab_mc01) %in% "Firmicutes", ]
#For order
major_taxa_proportions_tab_o <- apply(order_and_unidentified_counts_tab, 2, function(x) x/sum(x)*100)
major_taxa_proportions_tab_o_mc01 <- apply(order_and_unidentified_counts_tab_mc01, 2, function(x) x/sum(x)*100)
major_taxa_proportions_tab_o_gc02 <- apply(order_and_unidentified_counts_tab_gc02, 2, function(x) x/sum(x)*100)
major_taxa_proportions_tab_o_gc04 <- apply(order_and_unidentified_counts_tab_gc04, 2, function(x) x/sum(x)*100)
major_taxa_proportions_tab_o_gc06 <- apply(order_and_unidentified_counts_tab_gc06, 2, function(x) x/sum(x)*100)
#major_taxa_proportions_tab_o_gcnone <- apply(order_and_unidentified_counts_tab_gcnone, 2, function(x) x/sum(x)*100)

#For genus
major_taxa_proportions_tab_g <- apply(genus_and_unidentified_counts_tab, 2, function(x) x/sum(x)*100)
major_taxa_proportions_tab_g_mc01 <- apply(genus_and_unidentified_counts_tab_mc01, 2, function(x) x/sum(x)*100)
major_taxa_proportions_tab_g_gc02 <- apply(genus_and_unidentified_counts_tab_gc02, 2, function(x) x/sum(x)*100)
major_taxa_proportions_tab_g_gc04 <- apply(genus_and_unidentified_counts_tab_gc04, 2, function(x) x/sum(x)*100)
major_taxa_proportions_tab_g_gc06 <- apply(genus_and_unidentified_counts_tab_gc06, 2, function(x) x/sum(x)*100)
#major_taxa_proportions_tab_g_gcnone <- apply(genus_and_unidentified_counts_tab_gcnone, 2, function(x) x/sum(x)*100)

#For classes of proteobacteria and phyla of sulfur transforming microbes
#major_taxa_proportions_tab_g <- apply(genus_and_unidentified_counts_tab, 2, function(x) x/sum(x)*100)
proteobacteria_major_taxa_proportions_tab_g_mc01 <- apply(class_and_unidentified_counts_tab_proteobacteria_mc01, 2, function(x) x/sum(x)*100)
proteobacteria_major_taxa_proportions_tab_g_gc02 <- apply(class_and_unidentified_counts_tab_proteobacteria_gc02, 2, function(x) x/sum(x)*100)
proteobacteria_major_taxa_proportions_tab_g_gc04 <- apply(class_and_unidentified_counts_tab_proteobacteria_gc04, 2, function(x) x/sum(x)*100)
proteobacteria_major_taxa_proportions_tab_g_gc06 <- apply(class_and_unidentified_counts_tab_proteobacteria_gc06, 2, function(x) x/sum(x)*100)
colSums(proteobacteria_major_taxa_proportions_tab_g_gc02)
colSums(proteobacteria_major_taxa_proportions_tab_g_gc04)
colSums(proteobacteria_major_taxa_proportions_tab_g_gc06)

#For orders of proteobacteria 
#major_taxa_proportions_tab_g <- apply(genus_and_unidentified_counts_tab, 2, function(x) x/sum(x)*100)
proteobacteria_major_taxa_proportions_tab_o_mc01 <- apply(order_and_unidentified_counts_tab_proteobacteria_mc01, 2, function(x) x/sum(x)*100)
proteobacteria_major_taxa_proportions_tab_o_gc02 <- apply(order_and_unidentified_counts_tab_proteobacteria_gc02, 2, function(x) x/sum(x)*100)
proteobacteria_major_taxa_proportions_tab_o_gc04 <- apply(order_and_unidentified_counts_tab_proteobacteria_gc04, 2, function(x) x/sum(x)*100)
proteobacteria_major_taxa_proportions_tab_o_gc06 <- apply(order_and_unidentified_counts_tab_proteobacteria_gc06, 2, function(x) x/sum(x)*100)
colSums(proteobacteria_major_taxa_proportions_tab_o_gc02)
colSums(proteobacteria_major_taxa_proportions_tab_o_gc04)
colSums(proteobacteria_major_taxa_proportions_tab_o_gc06)

#For genuses of Firmicutes
firmicutes_major_taxa_proportions_tab_g_mc01 <- apply(genus_counts_tab_firm_mc01, 2, function(x) x/sum(x)*100)
firmicutes_major_taxa_proportions_tab_g_gc02 <- apply(genus_counts_tab_firm_gc02, 2, function(x) x/sum(x)*100)
firmicutes_major_taxa_proportions_tab_g_gc04 <- apply(genus_counts_tab_firm_gc04, 2, function(x) x/sum(x)*100)
firmicutes_major_taxa_proportions_tab_g_gc06 <- apply(genus_counts_tab_firm_gc06, 2, function(x) x/sum(x)*100)
colSums(firmicutes_major_taxa_proportions_tab_g_gc02)


# If we check the dimensions of this table at this point
dim(major_taxa_proportions_tab)
dim(major_taxa_proportions_tab_f)
dim(major_taxa_proportions_tab_g)
head(major_taxa_proportions_tab_g[ ,1:5])
# We see there are currently 53 rows, which might be a little busy for a summary figure
#Dimensions of the core tables
#For phylum
dim(major_taxa_proportions_tab_mc01)
dim(major_taxa_proportions_tab_gc02)
dim(major_taxa_proportions_tab_gc04)
dim(major_taxa_proportions_tab_gc06)
#For family
dim(major_taxa_proportions_tab_f_mc01)
dim(major_taxa_proportions_tab_f_gc02)
dim(major_taxa_proportions_tab_f_gc04)
dim(major_taxa_proportions_tab_f_gc06)
#For genus
dim(major_taxa_proportions_tab_g_mc01) #67 rows
dim(major_taxa_proportions_tab_g_gc02) #84 rows
dim(major_taxa_proportions_tab_g_gc04) #79 rows
dim(major_taxa_proportions_tab_g_gc06) #92 rows

#Get the names of all genuses identified
row.names(major_taxa_proportions_tab_g_mc01)
row.names(major_taxa_proportions_tab_g_gc02)
row.names(major_taxa_proportions_tab_g_gc04)
row.names(major_taxa_proportions_tab_g_gc06)
#row.names(major_taxa_proportions_tab_g_gcnone)

#Get the genuses that are in mc01 and not in gc02, gc04, gc06
setdiff(row.names(major_taxa_proportions_tab_g_mc01),row.names(major_taxa_proportions_tab_g_gc02))
setdiff(row.names(major_taxa_proportions_tab_g_mc01),row.names(major_taxa_proportions_tab_g_gc04))
setdiff(row.names(major_taxa_proportions_tab_g_mc01),row.names(major_taxa_proportions_tab_g_gc06))
setdiff(row.names(major_taxa_proportions_tab_g_mc01),row.names(major_taxa_proportions_tab_g_gcnone))

#Get the genuses that are in gc02 and not in mc01, gc04, gc06
setdiff(row.names(major_taxa_proportions_tab_g_gc02),row.names(major_taxa_proportions_tab_g_mc01))
setdiff(row.names(major_taxa_proportions_tab_g_gc02),row.names(major_taxa_proportions_tab_g_gc04))
setdiff(row.names(major_taxa_proportions_tab_g_gc02),row.names(major_taxa_proportions_tab_g_gc06))
setdiff(row.names(major_taxa_proportions_tab_g_gc02),row.names(major_taxa_proportions_tab_g_gcnone))

#Get the genuses that are in gc04 and not in mc01, gc02, gc06
setdiff(row.names(major_taxa_proportions_tab_g_gc04),row.names(major_taxa_proportions_tab_g_mc01))
setdiff(row.names(major_taxa_proportions_tab_g_gc04),row.names(major_taxa_proportions_tab_g_gc02))
setdiff(row.names(major_taxa_proportions_tab_g_gc04),row.names(major_taxa_proportions_tab_g_gc06))
setdiff(row.names(major_taxa_proportions_tab_g_gc04),row.names(major_taxa_proportions_tab_g_gcnone))

#Get the genuses that are in gc06 and not in mc01, gc02, gc04
setdiff(row.names(major_taxa_proportions_tab_g_gc06),row.names(major_taxa_proportions_tab_g_mc01))
setdiff(row.names(major_taxa_proportions_tab_g_gc06),row.names(major_taxa_proportions_tab_g_gc02))
setdiff(row.names(major_taxa_proportions_tab_g_gc06),row.names(major_taxa_proportions_tab_g_gc04))
setdiff(row.names(major_taxa_proportions_tab_g_gc06),row.names(major_taxa_proportions_tab_g_gcnone))

#Get the genuses that are in the neg controls and not in mc01, gc02, gc04, gc06
setdiff(row.names(major_taxa_proportions_tab_g_gcnone),row.names(major_taxa_proportions_tab_g_mc01))
setdiff(row.names(major_taxa_proportions_tab_g_gcnone),row.names(major_taxa_proportions_tab_g_gc02))
setdiff(row.names(major_taxa_proportions_tab_g_gcnone),row.names(major_taxa_proportions_tab_g_gc04))
setdiff(row.names(major_taxa_proportions_tab_g_gcnone),row.names(major_taxa_proportions_tab_g_gc06))

#######

#Create dataframe with unique taxa for each core
comp<-list(list_mc01 = row.names(major_taxa_proportions_tab_g_mc01),
           list_gc02 = row.names(major_taxa_proportions_tab_g_gc02),
           list_gc04 = row.names(major_taxa_proportions_tab_g_gc04),
           list_gc06 = row.names(major_taxa_proportions_tab_g_gc06),
           #list_gcnone = row.names(major_taxa_proportions_tab_g_gcnone),
)
df_unique_taxa<-lapply(1:length(comp), function(n) setdiff(comp[[n]], unlist(comp[-n])))
df_unique_taxa

# Many of these taxa make up a very small percentage, so we're going to filter some out
# this is a completely arbitrary decision solely to ease visualization and
# is entirely up to your data and you
# Here, we'll only keep rows (taxa) that make up greater than 5% in any
# individual sample
#For phylum
temp_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 5, ])
temp_filt_major_taxa_proportions_tab_mc01 <- data.frame(major_taxa_proportions_tab_mc01[apply(major_taxa_proportions_tab_mc01, 1, max) > 5, ])
temp_filt_major_taxa_proportions_tab_gc02 <- data.frame(major_taxa_proportions_tab_gc02[apply(major_taxa_proportions_tab_gc02, 1, max) > 5, ])
temp_filt_major_taxa_proportions_tab_gc04 <- data.frame(major_taxa_proportions_tab_gc04[apply(major_taxa_proportions_tab_gc04, 1, max) > 5, ])
temp_filt_major_taxa_proportions_tab_gc06 <- data.frame(major_taxa_proportions_tab_gc06[apply(major_taxa_proportions_tab_gc06, 1, max) > 5, ])
#For family
temp_filt_major_taxa_proportions_tab_f <- data.frame(major_taxa_proportions_tab_f[apply(major_taxa_proportions_tab_f, 1, max) > 5, ])
temp_filt_major_taxa_proportions_tab_f_mc01 <- data.frame(major_taxa_proportions_tab_f_mc01[apply(major_taxa_proportions_tab_f_mc01, 1, max) > 5, ])
temp_filt_major_taxa_proportions_tab_f_gc02 <- data.frame(major_taxa_proportions_tab_f_gc02[apply(major_taxa_proportions_tab_f_gc02, 1, max) > 5, ])
temp_filt_major_taxa_proportions_tab_f_gc04 <- data.frame(major_taxa_proportions_tab_f_gc04[apply(major_taxa_proportions_tab_f_gc04, 1, max) > 5, ])
temp_filt_major_taxa_proportions_tab_f_gc06 <- data.frame(major_taxa_proportions_tab_f_gc06[apply(major_taxa_proportions_tab_f_gc06, 1, max) > 5, ])
#For genus
temp_filt_major_taxa_proportions_tab_g <- data.frame(major_taxa_proportions_tab_g[apply(major_taxa_proportions_tab_g, 1, max) > 5, ])
temp_filt_major_taxa_proportions_tab_g_mc01 <- data.frame(major_taxa_proportions_tab_g_mc01[apply(major_taxa_proportions_tab_g_mc01, 1, max) > 5, ])
temp_filt_major_taxa_proportions_tab_g_gc02 <- data.frame(major_taxa_proportions_tab_g_gc02[apply(major_taxa_proportions_tab_g_gc02, 1, max) > 5, ])
temp_filt_major_taxa_proportions_tab_g_gc04 <- data.frame(major_taxa_proportions_tab_g_gc04[apply(major_taxa_proportions_tab_g_gc04, 1, max) > 5, ])
temp_filt_major_taxa_proportions_tab_g_gc06 <- data.frame(major_taxa_proportions_tab_g_gc06[apply(major_taxa_proportions_tab_g_gc06, 1, max) > 5, ])
#And make dataframe for all the others for phylum
other_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) < 5, ])
other_filt_major_taxa_proportions_tab_mc01 <- data.frame(major_taxa_proportions_tab_mc01[apply(major_taxa_proportions_tab_mc01, 1, max) < 5, ])
other_filt_major_taxa_proportions_tab_gc02 <- data.frame(major_taxa_proportions_tab_gc02[apply(major_taxa_proportions_tab_gc02, 1, max) < 5, ])
other_filt_major_taxa_proportions_tab_gc04 <- data.frame(major_taxa_proportions_tab_gc04[apply(major_taxa_proportions_tab_gc04, 1, max) < 5, ])
other_filt_major_taxa_proportions_tab_gc06 <- data.frame(major_taxa_proportions_tab_gc06[apply(major_taxa_proportions_tab_gc06, 1, max) < 5, ])

options(max.print = 1000000)
other_filt_major_taxa_proportions_tab_g_mc01
# checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) # now we have 16, much more manageable for an overview figure
dim(temp_filt_major_taxa_proportions_tab_f) 
dim(temp_filt_major_taxa_proportions_tab_g_gc06) 

#Check for the cores
#For phylum
dim(temp_filt_major_taxa_proportions_tab_mc01) # 5 main rows (taxa)
dim(temp_filt_major_taxa_proportions_tab_gc02) # 11 main rows (taxa)
dim(temp_filt_major_taxa_proportions_tab_gc04) # 8 main rows (taxa)
dim(temp_filt_major_taxa_proportions_tab_gc06) # 14 main rows (taxa)

#For family
dim(temp_filt_major_taxa_proportions_tab_f_mc01) # 8 main rows (taxa)
dim(temp_filt_major_taxa_proportions_tab_f_gc02) # 16 main rows (taxa)
dim(temp_filt_major_taxa_proportions_tab_f_gc04) # 15 main rows (taxa)
dim(temp_filt_major_taxa_proportions_tab_f_gc06) # 20 main rows (taxa)

#For genus
dim(temp_filt_major_taxa_proportions_tab_g_mc01) # 7 main rows (taxa)
dim(temp_filt_major_taxa_proportions_tab_g_gc02) # 14 main rows (taxa)
dim(temp_filt_major_taxa_proportions_tab_g_gc04) # 12 main rows (taxa)
dim(temp_filt_major_taxa_proportions_tab_g_gc06) # 17 main rows (taxa)

#Get the names of the taxa that comprise the "Other" pool for each core
#For phyla
others_phyla <- data.frame(major_taxa_proportions_tab) %>% filter_all(all_vars(.< 5.0))
others_phyla_mc01 <- data.frame(major_taxa_proportions_tab_mc01) %>% filter_all(all_vars(.< 5.0))
others_phyla_gc02 <- data.frame(major_taxa_proportions_tab_gc02) %>% filter_all(all_vars(.< 5.0))
others_phyla_gc04 <- data.frame(major_taxa_proportions_tab_gc04) %>% filter_all(all_vars(.< 5.0))
others_phyla_gc06 <- data.frame(major_taxa_proportions_tab_gc06) %>% filter_all(all_vars(.< 5.0))

#For family
others_family_f <- data.frame(major_taxa_proportions_tab_f) %>% filter_all(all_vars(.< 5.0))
others_family_mc01 <- data.frame(major_taxa_proportions_tab_f_mc01) %>% filter_all(all_vars(.< 5.0))
others_family_gc02 <- data.frame(major_taxa_proportions_tab_f_gc02) %>% filter_all(all_vars(.< 5.0))
others_family_gc04 <- data.frame(major_taxa_proportions_tab_f_gc04) %>% filter_all(all_vars(.< 5.0))
others_family_gc06 <- data.frame(major_taxa_proportions_tab_f_gc06) %>% filter_all(all_vars(.< 5.0))

#For genus
others_genus_f <- data.frame(major_taxa_proportions_tab_g) %>% filter_all(all_vars(.< 5.0))
others_genus_mc01 <- data.frame(major_taxa_proportions_tab_g_mc01) %>% filter_all(all_vars(.< 5.0))
others_genus_gc02 <- data.frame(major_taxa_proportions_tab_g_gc02) %>% filter_all(all_vars(.< 5.0))
others_genus_gc04 <- data.frame(major_taxa_proportions_tab_g_gc04) %>% filter_all(all_vars(.< 5.0))
others_genus_gc06 <- data.frame(major_taxa_proportions_tab_g_gc06) %>% filter_all(all_vars(.< 5.0))

#Create dataframe with unique "Other" taxa for each core
comp_other<-list(list_mc01 = row.names(others_genus_mc01),
                 list_gc02 = row.names(others_genus_gc02),
                 list_gc04 = row.names(others_genus_gc04),
                 list_gc06 = row.names(others_genus_gc06)
)
df_unique_other_taxa<-lapply(1:length(comp_other), function(n) setdiff(comp_other[[n]], unlist(comp_other[-n])))
df_unique_other_taxa

# though each of the filtered taxa made up less than 5% alone, together they
# may add up, and they should still be included in the overall summary
# so we're going to add a row called "Other" that keeps track of how much we
# filtered out (which will also keep our totals at 100%)
#For phylum
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other" = filtered_proportions)
filt_major_taxa_proportions_tab
#For family
filtered_proportions_f <- colSums(major_taxa_proportions_tab_f) - colSums(temp_filt_major_taxa_proportions_tab_f)
filt_major_taxa_proportions_tab_f <- rbind(temp_filt_major_taxa_proportions_tab_f, "Other" = filtered_proportions_f)
filt_major_taxa_proportions_tab_f
#For genus
filtered_proportions_g <- colSums(major_taxa_proportions_tab_g) - colSums(temp_filt_major_taxa_proportions_tab_g)
filt_major_taxa_proportions_tab_g <- rbind(temp_filt_major_taxa_proportions_tab_g, "Other" = filtered_proportions_g)
filt_major_taxa_proportions_tab_g

#For the cores 
#By phylum
filtered_proportions_mc01 <- colSums(major_taxa_proportions_tab_mc01) - colSums(temp_filt_major_taxa_proportions_tab_mc01)
filtered_proportions_gc02 <- colSums(major_taxa_proportions_tab_gc02) - colSums(temp_filt_major_taxa_proportions_tab_gc02)
filtered_proportions_gc04 <- colSums(major_taxa_proportions_tab_gc04) - colSums(temp_filt_major_taxa_proportions_tab_gc04)
filtered_proportions_gc06 <- colSums(major_taxa_proportions_tab_gc06) - colSums(temp_filt_major_taxa_proportions_tab_gc06)

filt_major_taxa_proportions_tab_mc01 <- rbind(temp_filt_major_taxa_proportions_tab_mc01, "Other" = filtered_proportions_mc01)
filt_major_taxa_proportions_tab_gc02 <- rbind(temp_filt_major_taxa_proportions_tab_gc02, "Other" = filtered_proportions_gc02)
filt_major_taxa_proportions_tab_gc04 <- rbind(temp_filt_major_taxa_proportions_tab_gc04, "Other" = filtered_proportions_gc04)
filt_major_taxa_proportions_tab_gc06 <- rbind(temp_filt_major_taxa_proportions_tab_gc06, "Other" = filtered_proportions_gc06)

filt_major_taxa_proportions_tab_mc01
filt_major_taxa_proportions_tab_gc02
filt_major_taxa_proportions_tab_gc04
filt_major_taxa_proportions_tab_gc06

#By family
filtered_proportions_f_mc01 <- colSums(major_taxa_proportions_tab_f_mc01) - colSums(temp_filt_major_taxa_proportions_tab_f_mc01)
filtered_proportions_f_gc02 <- colSums(major_taxa_proportions_tab_f_gc02) - colSums(temp_filt_major_taxa_proportions_tab_f_gc02)
filtered_proportions_f_gc04 <- colSums(major_taxa_proportions_tab_f_gc04) - colSums(temp_filt_major_taxa_proportions_tab_f_gc04)
filtered_proportions_f_gc06 <- colSums(major_taxa_proportions_tab_f_gc06) - colSums(temp_filt_major_taxa_proportions_tab_f_gc06)

filt_major_taxa_proportions_tab_f_mc01 <- rbind(temp_filt_major_taxa_proportions_tab_f_mc01, "Other" = filtered_proportions_f_mc01)
filt_major_taxa_proportions_tab_f_gc02 <- rbind(temp_filt_major_taxa_proportions_tab_f_gc02, "Other" = filtered_proportions_f_gc02)
filt_major_taxa_proportions_tab_f_gc04 <- rbind(temp_filt_major_taxa_proportions_tab_f_gc04, "Other" = filtered_proportions_f_gc04)
filt_major_taxa_proportions_tab_f_gc06 <- rbind(temp_filt_major_taxa_proportions_tab_f_gc06, "Other" = filtered_proportions_f_gc06)

filt_major_taxa_proportions_tab_f_mc01
filt_major_taxa_proportions_tab_f_gc02
filt_major_taxa_proportions_tab_f_gc04
filt_major_taxa_proportions_tab_f_gc06

#By genus
filtered_proportions_g_mc01 <- colSums(major_taxa_proportions_tab_g_mc01) - colSums(temp_filt_major_taxa_proportions_tab_g_mc01)
filtered_proportions_g_gc02 <- colSums(major_taxa_proportions_tab_g_gc02) - colSums(temp_filt_major_taxa_proportions_tab_g_gc02)
filtered_proportions_g_gc04 <- colSums(major_taxa_proportions_tab_g_gc04) - colSums(temp_filt_major_taxa_proportions_tab_g_gc04)
filtered_proportions_g_gc06 <- colSums(major_taxa_proportions_tab_g_gc06) - colSums(temp_filt_major_taxa_proportions_tab_g_gc06)

filt_major_taxa_proportions_tab_g_mc01 <- rbind(temp_filt_major_taxa_proportions_tab_g_mc01, "Other" = filtered_proportions_g_mc01)
filt_major_taxa_proportions_tab_g_gc02 <- rbind(temp_filt_major_taxa_proportions_tab_g_gc02, "Other" = filtered_proportions_g_gc02)
filt_major_taxa_proportions_tab_g_gc04 <- rbind(temp_filt_major_taxa_proportions_tab_g_gc04, "Other" = filtered_proportions_g_gc04)
filt_major_taxa_proportions_tab_g_gc06 <- rbind(temp_filt_major_taxa_proportions_tab_g_gc06, "Other" = filtered_proportions_g_gc06)

head(filt_major_taxa_proportions_tab_g_mc01)
filt_major_taxa_proportions_tab_g_gc02
filt_major_taxa_proportions_tab_g_gc04
filt_major_taxa_proportions_tab_g_gc06

#Now we can make some plots
# first let's make a copy of our table that's safe for manipulating
#For phylum
filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab
#for each core
filt_major_taxa_proportions_tab_for_plot_mc01 <- filt_major_taxa_proportions_tab_mc01
filt_major_taxa_proportions_tab_for_plot_gc02 <- filt_major_taxa_proportions_tab_gc02
filt_major_taxa_proportions_tab_for_plot_gc04 <- filt_major_taxa_proportions_tab_gc04
filt_major_taxa_proportions_tab_for_plot_gc06 <- filt_major_taxa_proportions_tab_gc06
#For family
filt_major_taxa_proportions_tab_for_plot_f <- filt_major_taxa_proportions_tab_f
#for each core
filt_major_taxa_proportions_tab_for_plot_f_mc01 <- filt_major_taxa_proportions_tab_f_mc01
filt_major_taxa_proportions_tab_for_plot_f_gc02 <- filt_major_taxa_proportions_tab_f_gc02
filt_major_taxa_proportions_tab_for_plot_f_gc04 <- filt_major_taxa_proportions_tab_f_gc04
filt_major_taxa_proportions_tab_for_plot_f_gc06 <- filt_major_taxa_proportions_tab_f_gc06
#For genus
filt_major_taxa_proportions_tab_for_plot_g <- filt_major_taxa_proportions_tab_g
#for each core
filt_major_taxa_proportions_tab_for_plot_g_mc01 <- filt_major_taxa_proportions_tab_g_mc01
filt_major_taxa_proportions_tab_for_plot_g_gc02 <- filt_major_taxa_proportions_tab_g_gc02
filt_major_taxa_proportions_tab_for_plot_g_gc04 <- filt_major_taxa_proportions_tab_g_gc04
filt_major_taxa_proportions_tab_for_plot_g_gc06 <- filt_major_taxa_proportions_tab_g_gc06

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
#For phylum
filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab_for_plot %>% rownames_to_column("Phylum")

phyla_and_unidentified_counts_tab <- data.frame(phyla_and_unidentified_counts_tab)
phyla_and_unidentified_counts_tab_export <- phyla_and_unidentified_counts_tab %>% rownames_to_column("Phylum")
# Export the table outside of R
write.table(phyla_and_unidentified_counts_tab_export, "SR2113_phylum.tsv", sep = "\t", quote = FALSE, col.names = NA)

order_and_unidentified_counts_tab <- data.frame(order_and_unidentified_counts_tab)
order_and_unidentified_counts_tab_export <- order_and_unidentified_counts_tab %>% rownames_to_column("Order")
# Export the table outside of R
write.table(order_and_unidentified_counts_tab_export, "SR2113_order.tsv", sep = "\t", quote = FALSE, col.names = NA)

genus_and_unidentified_counts_tab <- data.frame(genus_and_unidentified_counts_tab)
genus_and_unidentified_counts_tab_export <- genus_and_unidentified_counts_tab %>% rownames_to_column("Genus")
# Export the table outside of R
write.table(genus_and_unidentified_counts_tab_export, "SR2113_genus.tsv", sep = "\t", quote = FALSE, col.names = NA)



#for the cores
filt_major_taxa_proportions_tab_for_plot_mc01 <- filt_major_taxa_proportions_tab_for_plot_mc01 %>% rownames_to_column("Phylum")
filt_major_taxa_proportions_tab_for_plot_gc02 <- filt_major_taxa_proportions_tab_for_plot_gc02 %>% rownames_to_column("Phylum")
filt_major_taxa_proportions_tab_for_plot_gc04 <- filt_major_taxa_proportions_tab_for_plot_gc04 %>% rownames_to_column("Phylum")
filt_major_taxa_proportions_tab_for_plot_gc06 <- filt_major_taxa_proportions_tab_for_plot_gc06 %>% rownames_to_column("Phylum")

#For family
filt_major_taxa_proportions_tab_for_plot_f <- filt_major_taxa_proportions_tab_for_plot_f %>% rownames_to_column("Family")
phyla_and_unidentified_counts_tab <- data.frame(phyla_and_unidentified_counts_tab)
phyla_and_unidentified_counts_tab_export <- phyla_and_unidentified_counts_tab %>% rownames_to_column("Phylum")
#phyla_and_unidentified_counts_tab_export <- phyla_and_unidentified_counts_tab_export %>%
#filter(Phylum != "Firmicutes")
# Export the table outside of R
write.table(phyla_and_unidentified_counts_tab_export, "SR2113_phylum.tsv", sep = "\t", quote = FALSE, col.names = NA)

#for the cores
filt_major_taxa_proportions_tab_for_plot_f_mc01 <- filt_major_taxa_proportions_tab_for_plot_f_mc01 %>% rownames_to_column("Family")
filt_major_taxa_proportions_tab_for_plot_f_gc02 <- filt_major_taxa_proportions_tab_for_plot_f_gc02 %>% rownames_to_column("Family")
filt_major_taxa_proportions_tab_for_plot_f_gc04 <- filt_major_taxa_proportions_tab_for_plot_f_gc04 %>% rownames_to_column("Family")
filt_major_taxa_proportions_tab_for_plot_f_gc06 <- filt_major_taxa_proportions_tab_for_plot_f_gc06 %>% rownames_to_column("Family")

#For genus
filt_major_taxa_proportions_tab_for_plot_g <- filt_major_taxa_proportions_tab_for_plot_g %>% rownames_to_column("Genus")
#for the cores
filt_major_taxa_proportions_tab_for_plot_g_mc01 <- filt_major_taxa_proportions_tab_for_plot_g_mc01 %>% rownames_to_column("Genus")
filt_major_taxa_proportions_tab_for_plot_g_gc02 <- filt_major_taxa_proportions_tab_for_plot_g_gc02 %>% rownames_to_column("Genus")
filt_major_taxa_proportions_tab_for_plot_g_gc04 <- filt_major_taxa_proportions_tab_for_plot_g_gc04 %>% rownames_to_column("Genus")
filt_major_taxa_proportions_tab_for_plot_g_gc06 <- filt_major_taxa_proportions_tab_for_plot_g_gc06 %>% rownames_to_column("Genus")

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
#For phylum
filt_major_taxa_proportions_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot %>% pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)
#for the cores
filt_major_taxa_proportions_tab_for_plot_mc01.g <- filt_major_taxa_proportions_tab_for_plot_mc01 %>% pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)
filt_major_taxa_proportions_tab_for_plot_gc02.g <- filt_major_taxa_proportions_tab_for_plot_gc02 %>% pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)
filt_major_taxa_proportions_tab_for_plot_gc04.g <- filt_major_taxa_proportions_tab_for_plot_gc04 %>% pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)
filt_major_taxa_proportions_tab_for_plot_gc06.g <- filt_major_taxa_proportions_tab_for_plot_gc06 %>% pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)

#For family
filt_major_taxa_proportions_tab_for_plot_f.g <- filt_major_taxa_proportions_tab_for_plot_f %>% pivot_longer(!Family, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)
#for the cores
filt_major_taxa_proportions_tab_for_plot_f_mc01.g <- filt_major_taxa_proportions_tab_for_plot_f_mc01 %>% pivot_longer(!Family, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)
filt_major_taxa_proportions_tab_for_plot_f_gc02.g <- filt_major_taxa_proportions_tab_for_plot_f_gc02 %>% pivot_longer(!Family, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)
filt_major_taxa_proportions_tab_for_plot_f_gc04.g <- filt_major_taxa_proportions_tab_for_plot_f_gc04 %>% pivot_longer(!Family, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)
filt_major_taxa_proportions_tab_for_plot_f_gc06.g <- filt_major_taxa_proportions_tab_for_plot_f_gc06 %>% pivot_longer(!Family, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)

#For genus
filt_major_taxa_proportions_tab_for_plot_g.g <- filt_major_taxa_proportions_tab_for_plot_g %>% pivot_longer(!Genus, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)
#for the cores
filt_major_taxa_proportions_tab_for_plot_g_mc01.g <- filt_major_taxa_proportions_tab_for_plot_g_mc01 %>% pivot_longer(!Genus, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)
filt_major_taxa_proportions_tab_for_plot_g_gc02.g <- filt_major_taxa_proportions_tab_for_plot_g_gc02 %>% pivot_longer(!Genus, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)
filt_major_taxa_proportions_tab_for_plot_g_gc04.g <- filt_major_taxa_proportions_tab_for_plot_g_gc04 %>% pivot_longer(!Genus, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)
filt_major_taxa_proportions_tab_for_plot_g_gc06.g <- filt_major_taxa_proportions_tab_for_plot_g_gc06 %>% pivot_longer(!Genus, names_to = "Sample", values_to = "Proportion") %>% data.frame(check.names = FALSE)

# take a look at the new table and compare it with the old one
#For phylum
head(filt_major_taxa_proportions_tab_for_plot.g)
head(filt_major_taxa_proportions_tab_for_plot)
#for the cores
head(filt_major_taxa_proportions_tab_for_plot_mc01.g)
head(filt_major_taxa_proportions_tab_for_plot_mc01)
#for family
head(filt_major_taxa_proportions_tab_for_plot_f.g)
head(filt_major_taxa_proportions_tab_for_plot_f)
#for the cores
head(filt_major_taxa_proportions_tab_for_plot_f_mc01.g)
head(filt_major_taxa_proportions_tab_for_plot_f_mc01)
#for genus
head(filt_major_taxa_proportions_tab_for_plot_g.g)
head(filt_major_taxa_proportions_tab_for_plot_g)
#for the cores
head(filt_major_taxa_proportions_tab_for_plot_g_mc01.g)
head(filt_major_taxa_proportions_tab_for_plot_g_mc01)

# manipulating tables like this is something you may need to do frequently in R

# it will also help us to have some of our sample information included in our
# plotting table, we can take advantage of the following function to merge the two
# tables for use based on the Sample ID
# first, we do need to modify our sample names a bit to match how they were
# changed by phyloseq (really the data.frame() function)
mod_sample_info_tab <- sample_info_tab_1
#Make sample info tables for each core
mod_sample_info_tab_mc01 <- filter(mod_sample_info_tab, Core == "MC01")
mod_sample_info_tab_gc02 <- filter(mod_sample_info_tab, Core == "GC02")
mod_sample_info_tab_gc04 <- filter(mod_sample_info_tab, Core == "GC04")
mod_sample_info_tab_gc06 <- filter(mod_sample_info_tab, Core == "GC06")

#Change the depth column to integers in the core tables
mod_sample_info_tab_mc01$Depth <- as.integer(mod_sample_info_tab_mc01$Depth)
mod_sample_info_tab_gc02$Depth <- as.integer(mod_sample_info_tab_gc02$Depth)
mod_sample_info_tab_gc04$Depth <- as.integer(mod_sample_info_tab_gc04$Depth)
mod_sample_info_tab_gc06$Depth <- as.integer(mod_sample_info_tab_gc06$Depth)

#Rearrange the tables of each core by depth
mod_sample_info_tab_mc01 <- mod_sample_info_tab_mc01[order(mod_sample_info_tab_mc01$Depth),]
mod_sample_info_tab_gc02 <- mod_sample_info_tab_gc02[order(mod_sample_info_tab_gc02$Depth),]
mod_sample_info_tab_gc04 <- mod_sample_info_tab_gc04[order(mod_sample_info_tab_gc04$Depth),]
mod_sample_info_tab_gc06 <- mod_sample_info_tab_gc06[order(mod_sample_info_tab_gc06$Depth),]

#Change dashes to periods
mod_sample_info_tab$Sample <- gsub("-", ".", mod_sample_info_tab$Sample)
#For phylum
filt_major_taxa_proportions_tab_for_plot.g2 <- filt_major_taxa_proportions_tab_for_plot.g %>% left_join(mod_sample_info_tab)
#For family
filt_major_taxa_proportions_tab_for_plot_f.g2 <- filt_major_taxa_proportions_tab_for_plot_f.g %>% left_join(mod_sample_info_tab)
#For genus
filt_major_taxa_proportions_tab_for_plot_g.g2 <- filt_major_taxa_proportions_tab_for_plot_g.g %>% left_join(mod_sample_info_tab)

#For the cores
#MC01
mod_sample_info_tab_mc01$Sample <- gsub("-", ".", mod_sample_info_tab_mc01$Sample)
#For phylum
filt_major_taxa_proportions_tab_for_plot_mc01.g2 <- filt_major_taxa_proportions_tab_for_plot_mc01.g %>% left_join(mod_sample_info_tab_mc01)
#For family
filt_major_taxa_proportions_tab_for_plot_f_mc01.g2 <- filt_major_taxa_proportions_tab_for_plot_f_mc01.g %>% left_join(mod_sample_info_tab_mc01)
#For genus
filt_major_taxa_proportions_tab_for_plot_g_mc01.g2 <- filt_major_taxa_proportions_tab_for_plot_g_mc01.g %>% left_join(mod_sample_info_tab_mc01)


#GC02
mod_sample_info_tab_gc02$Sample <- gsub("-", ".", mod_sample_info_tab_gc02$Sample)
#For phylum
filt_major_taxa_proportions_tab_for_plot_gc02.g2 <- filt_major_taxa_proportions_tab_for_plot_gc02.g %>% left_join(mod_sample_info_tab_gc02)
#For family
filt_major_taxa_proportions_tab_for_plot_f_gc02.g2 <- filt_major_taxa_proportions_tab_for_plot_f_gc02.g %>% left_join(mod_sample_info_tab_gc02)
#For genus
filt_major_taxa_proportions_tab_for_plot_g_gc02.g2 <- filt_major_taxa_proportions_tab_for_plot_g_gc02.g %>% left_join(mod_sample_info_tab_gc02)

#GC04
mod_sample_info_tab_gc04$Sample <- gsub("-", ".", mod_sample_info_tab_gc04$Sample)
#For phylum
filt_major_taxa_proportions_tab_for_plot_gc04.g2 <- filt_major_taxa_proportions_tab_for_plot_gc04.g %>% left_join(mod_sample_info_tab_gc04)
#For family
filt_major_taxa_proportions_tab_for_plot_f_gc04.g2 <- filt_major_taxa_proportions_tab_for_plot_f_gc04.g %>% left_join(mod_sample_info_tab_gc04)
#For genus
filt_major_taxa_proportions_tab_for_plot_g_gc04.g2 <- filt_major_taxa_proportions_tab_for_plot_g_gc04.g %>% left_join(mod_sample_info_tab_gc04)

#GCO6
mod_sample_info_tab_gc06$Sample <- gsub("-", ".", mod_sample_info_tab_gc06$Sample)
#For phylum
filt_major_taxa_proportions_tab_for_plot_gc06.g2 <- filt_major_taxa_proportions_tab_for_plot_gc06.g %>% left_join(mod_sample_info_tab_gc06)
#For family
filt_major_taxa_proportions_tab_for_plot_f_gc06.g2 <- filt_major_taxa_proportions_tab_for_plot_f_gc06.g %>% left_join(mod_sample_info_tab_gc06)
#For genus
filt_major_taxa_proportions_tab_for_plot_g_gc06.g2 <- filt_major_taxa_proportions_tab_for_plot_g_gc06.g %>% left_join(mod_sample_info_tab_gc06)

head(filt_major_taxa_proportions_tab_for_plot_gc06.g2)
#Create colorblind friendly palette
# The palette with black:
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills
#scale_fill_manual(values=cbbPalette)
# To use for line and point colors
#scale_colour_manual(values=cbbPalette)

#Now we have a table that’s easy to plug into ggplot2. One common way to visualize this is with stacked bar charts:
#Define color palette
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E69F00","#004949","#009292","#ff6db6","#924900",
                "#490092","#006ddb","#b66dff","#ffb6db", "#6db6ff",
                "#000000","#ffff6d",'#808080',"tomato", '#999999', '#f0e442', '#FFFFFF',
                "#920000","#924900","#24ff24","#ffff6d", '#999999', '#f0e442', '#FFFFFF')

stock_colors <- c("#FF0000","#F6A300","#0068CC","#6600AA","#AC0088","#AA33FF","#00FFFF","#00CC00","#006611","#00AC99",
                  "#AC6844","#FFFF00","#991100","#ACAC11","#a0f0aa","#FF00FF","#FF8611","#B9F6F6","#001166","#AC9A00",
                  "#994141","#ff1169","#0AF622","#119924","#Ac3311","#004A9A","#AcAc99","turquoise","tomato","sienna1",
                  "rosybrown","peachpuff","olivedrab3","mistyrose1","mediumorchid","indianred2","#114914","#660011",
                  "ivory3","deeppink","#331111", "#FF0000","#F6A300","#0068CC","#6600AA","#AC0088","#AA33FF","#00FFFF","#00CC00","#006611","#00AC99",
                  "#AC6844","#FFFF00","#991100","#ACAC11","#a0f0aa","#FF00FF","#FF8611","#B9F6F6","#001166","#AC9A00",
                  "#994141","#ff1169","#0AF622","#119924","#Ac3311","#004A9A","#AcAc99","turquoise","tomato","sienna1",
                  "rosybrown","peachpuff","olivedrab3","mistyrose1","mediumorchid","indianred2","#114914","#660011",
                  "ivory3","deeppink","#331111")
bold_mc01 <- c('#7F3C8D','#11A579','#3969AC','#F2B701','#E73F74',
               '#80BA5A','#E68310','#CF1C90','#008695','#f97b72',"#006ddb","#924900", "#FFE800", "#999999", "#E69F00", "#56B4E9", "#009E73")

#General plot
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(x = Sample, y = Proportion, fill = Phylum)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1), legend.title = element_blank()) +
  labs(x = "Sample", y = "% of 16S rRNA gene copies recovered", title = "All samples")

#Plots for S-cycling microbes By Phyla
#Create a copy of the tables for plots
major_taxa_proportions_tab_mc01_for_plot <- data.frame(major_taxa_proportions_tab_mc01)
major_taxa_proportions_tab_gc02_for_plot <- data.frame(major_taxa_proportions_tab_gc02)
major_taxa_proportions_tab_gc04_for_plot <- data.frame(major_taxa_proportions_tab_gc04)
major_taxa_proportions_tab_gc06_for_plot <- data.frame(major_taxa_proportions_tab_gc06)
head(major_taxa_proportions_tab_gc04)

#Simplify sample names
colnames(major_taxa_proportions_tab_mc01_for_plot) <- gsub("^[^_]*_", "", colnames(major_taxa_proportions_tab_mc01_for_plot))
colnames(major_taxa_proportions_tab_gc02_for_plot) <- gsub("^[^_]*_", "", colnames(major_taxa_proportions_tab_gc02_for_plot))
colnames(major_taxa_proportions_tab_gc04_for_plot) <- gsub("^[^_]*_", "", colnames(major_taxa_proportions_tab_gc04_for_plot))
colnames(major_taxa_proportions_tab_gc06_for_plot) <- gsub("^[^_]*_", "", colnames(major_taxa_proportions_tab_gc06_for_plot))

#Filter out samples with less than 600 total AVSs
major_taxa_proportions_tab_mc01_for_plot<-
  major_taxa_proportions_tab_mc01_for_plot[,colnames(major_taxa_proportions_tab_mc01_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
major_taxa_proportions_tab_gc02_for_plot<-
  major_taxa_proportions_tab_gc02_for_plot[,colnames(major_taxa_proportions_tab_gc02_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
major_taxa_proportions_tab_gc04_for_plot<-
  major_taxa_proportions_tab_gc04_for_plot[,colnames(major_taxa_proportions_tab_gc04_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
major_taxa_proportions_tab_gc06_for_plot<-
  major_taxa_proportions_tab_gc06_for_plot[,colnames(major_taxa_proportions_tab_gc06_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
head(major_taxa_proportions_tab_gc04_for_plot)

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
major_taxa_proportions_tab_mc01_for_plot <- 
  major_taxa_proportions_tab_mc01_for_plot %>% rownames_to_column("Taxa")
major_taxa_proportions_tab_gc02_for_plot <- 
  major_taxa_proportions_tab_gc02_for_plot %>% rownames_to_column("Taxa")
major_taxa_proportions_tab_gc04_for_plot <- 
  major_taxa_proportions_tab_gc04_for_plot %>% rownames_to_column("Taxa")
major_taxa_proportions_tab_gc06_for_plot <- 
  major_taxa_proportions_tab_gc06_for_plot %>% rownames_to_column("Taxa")

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
major_taxa_proportions_tab_mc01_for_plot <- 
  major_taxa_proportions_tab_mc01_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
major_taxa_proportions_tab_gc02_for_plot <- 
  major_taxa_proportions_tab_gc02_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
major_taxa_proportions_tab_gc04_for_plot <- 
  major_taxa_proportions_tab_gc04_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
major_taxa_proportions_tab_gc06_for_plot <- 
  major_taxa_proportions_tab_gc06_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)

#Add depth and core to tables
major_taxa_proportions_tab_mc01_for_plot <- 
  major_taxa_proportions_tab_mc01_for_plot %>% left_join(mod_sample_info_tab_mc01)
major_taxa_proportions_tab_gc02_for_plot <- 
  major_taxa_proportions_tab_gc02_for_plot %>% left_join(mod_sample_info_tab_gc02)
major_taxa_proportions_tab_gc04_for_plot <- 
  major_taxa_proportions_tab_gc04_for_plot %>% left_join(mod_sample_info_tab_gc04)
major_taxa_proportions_tab_gc06_for_plot <- 
  major_taxa_proportions_tab_gc06_for_plot %>% left_join(mod_sample_info_tab_gc06)

#SCB 01
#Organize things in plot

major_taxa_proportions_tab_mc01_for_plot$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_mc01_for_plot$Sample)
major_taxa_proportions_tab_mc01_for_plot$Sample <- factor(major_taxa_proportions_tab_mc01_for_plot$Sample, levels=unique(major_taxa_proportions_tab_mc01_for_plot$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_mc01_for_plot$Sample))))])
#major_taxa_proportions_tab_mc01_for_plot$Genus <- with(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2,factor(Class,Genus = rev(sort(unique(Genus)))))

#Make plot
#Include only S cycling taxa
major_taxa_proportions_tab_mc01_for_plot.1 <- major_taxa_proportions_tab_mc01_for_plot %>%
  filter(Taxa %in%  c("Bacteroidota",
                      "Spirochaetota",
                      "Marinimicrobia (SAR406 clade)",
                      "Planctomycetota",
                      "Actinobacteriota",
                      "Crenarchaeota",
                      "Thermoplasmatota",
                      "Verrucomicrobiota",
                      "Desulfobacterota_1",
                      "Desulfobacterota_2",
                      "Dadabacteria",
                      "Chloroflexi",
                      "Zixibacteria",
                      "Halobacterota",
                      "Nitrospirota",
                      "Acidobacteriota"))
major_taxa_proportions_tab_mc01_for_plot.2 <- major_taxa_proportions_tab_mc01_for_plot.1 %>% 
  mutate(Taxa = if_else(Taxa %ni% c("Desulfobacterota_1",
                                    "Desulfobacterota_2") , Taxa, "Desulfobacterota"))
#For phyla
major_taxa_proportions_tab_mc01_for_plot.2 %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "SCB 01")
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_sulfur_phyla_alldepths.pdf", width=12, height=5, dpi=300)


#Get averages for depth_bin
major_taxa_proportions_tab_mc01_for_plot.3<- major_taxa_proportions_tab_mc01_for_plot.2 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
major_taxa_proportions_tab_mc01_for_plot.3 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "SCB 01")+
  scale_x_discrete(limits = c("0-20", "21-40")) +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_sulfur_phyla.pdf", width=12, height=5, dpi=300)

#GC02
#Organize things in plot
major_taxa_proportions_tab_gc02_for_plot$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_gc02_for_plot$Sample)
major_taxa_proportions_tab_gc02_for_plot$Sample <- factor(major_taxa_proportions_tab_gc02_for_plot$Sample, levels=unique(major_taxa_proportions_tab_gc02_for_plot$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_gc02_for_plot$Sample))))])
#major_taxa_proportions_tab_gc02_for_plot$Genus <- with(proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2,factor(Class,Genus = rev(sort(unique(Genus)))))

#Make plot
#Include only S cycling taxa
major_taxa_proportions_tab_gc02_for_plot.1 <- major_taxa_proportions_tab_gc02_for_plot %>%
  filter(Taxa %in%  c("Bacteroidota",
                      "Spirochaetota",
                      "Marinimicrobia (SAR406 clade)",
                      "Planctomycetota",
                      "Actinobacteriota",
                      "Crenarchaeota",
                      "Thermoplasmatota",
                      "Verrucomicrobiota",
                      "Desulfobacterota_1",
                      "Desulfobacterota_2",
                      "Dadabacteria",
                      "Chloroflexi",
                      "Zixibacteria",
                      "Halobacterota",
                      "Nitrospirota",
                      "Acidobacteriota"))
major_taxa_proportions_tab_gc02_for_plot.2 <- major_taxa_proportions_tab_gc02_for_plot.1 %>% 
  mutate(Taxa = if_else(Taxa %ni% c("Desulfobacterota_1",
                                    "Desulfobacterota_2") , Taxa, "Desulfobacterota"))
#For phyla
major_taxa_proportions_tab_gc02_for_plot.2 %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR 02")
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_sulfur_phyla_alldepths.pdf", width=12, height=5, dpi=300)


#Get averages for depth_bin
major_taxa_proportions_tab_gc02_for_plot.3<- major_taxa_proportions_tab_gc02_for_plot.2 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
major_taxa_proportions_tab_gc02_for_plot.3 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR 02")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140')) +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_sulfur_phyla.pdf", width=12, height=5, dpi=300)


#GC04
#Organize things in plot
major_taxa_proportions_tab_gc04_for_plot$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_gc04_for_plot$Sample)
major_taxa_proportions_tab_gc04_for_plot$Sample <- factor(major_taxa_proportions_tab_gc04_for_plot$Sample, levels=unique(major_taxa_proportions_tab_gc04_for_plot$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_gc04_for_plot$Sample))))])
#major_taxa_proportions_tab_gc04_for_plot$Genus <- with(proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2,factor(Class,Genus = rev(sort(unique(Genus)))))

#Make plot
#Include only S cycling taxa
major_taxa_proportions_tab_gc04_for_plot.1 <- major_taxa_proportions_tab_gc04_for_plot %>%
  filter(Taxa %in%  c("Bacteroidota",
                      "Spirochaetota",
                      "Marinimicrobia (SAR406 clade)",
                      "Planctomycetota",
                      "Actinobacteriota",
                      "Crenarchaeota",
                      "Thermoplasmatota",
                      "Verrucomicrobiota",
                      "Desulfobacterota_1",
                      "Desulfobacterota_2",
                      "Dadabacteria",
                      "Chloroflexi",
                      "Zixibacteria",
                      "Halobacterota",
                      "Nitrospirota",
                      "Acidobacteriota"))
major_taxa_proportions_tab_gc04_for_plot.2 <- major_taxa_proportions_tab_gc04_for_plot.1 %>% 
  mutate(Taxa = if_else(Taxa %ni% c("Desulfobacterota_1",
                                    "Desulfobacterota_2") , Taxa, "Desulfobacterota"))
#For phyla
major_taxa_proportions_tab_gc04_for_plot.2 %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR 03")
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_sulfur_phyla_alldepths.pdf", width=12, height=5, dpi=300)


#Get averages for depth_bin
major_taxa_proportions_tab_gc04_for_plot.3<- major_taxa_proportions_tab_gc04_for_plot.2 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
major_taxa_proportions_tab_gc04_for_plot.3 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR 03")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140','141-220')) +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_sulfur_phyla.pdf", width=12, height=5, dpi=300)

#GC06
#Organize things in plot
major_taxa_proportions_tab_gc06_for_plot$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_gc06_for_plot$Sample)
major_taxa_proportions_tab_gc06_for_plot$Sample <- factor(major_taxa_proportions_tab_gc06_for_plot$Sample, levels=unique(major_taxa_proportions_tab_gc06_for_plot$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_gc06_for_plot$Sample))))])
#major_taxa_proportions_tab_gc06_for_plot$Genus <- with(proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2,factor(Class,Genus = rev(sort(unique(Genus)))))

#Make plot
#Include only S cycling taxa
major_taxa_proportions_tab_gc06_for_plot.1 <- major_taxa_proportions_tab_gc06_for_plot %>%
  filter(Taxa %in%  c("Bacteroidota",
                      "Spirochaetota",
                      "Marinimicrobia (SAR406 clade)",
                      "Planctomycetota",
                      "Actinobacteriota",
                      "Crenarchaeota",
                      "Thermoplasmatota",
                      "Verrucomicrobiota",
                      "Desulfobacterota_1",
                      "Desulfobacterota_2",
                      "Dadabacteria",
                      "Chloroflexi",
                      "Zixibacteria",
                      "Halobacterota",
                      "Nitrospirota",
                      "Acidobacteriota"))
major_taxa_proportions_tab_gc06_for_plot.2 <- major_taxa_proportions_tab_gc06_for_plot.1 %>% 
  mutate(Taxa = if_else(Taxa %ni% c("Desulfobacterota_1",
                                    "Desulfobacterota_2") , Taxa, "Desulfobacterota"))
#For phyla
major_taxa_proportions_tab_gc06_for_plot.2 %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR 04")
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_sulfur_phyla_alldepths.pdf", width=12, height=5, dpi=300)


#Get averages for depth_bin
major_taxa_proportions_tab_gc06_for_plot.3<- major_taxa_proportions_tab_gc06_for_plot.2 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
major_taxa_proportions_tab_gc06_for_plot.3 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR 04")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140','141-220')) +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_sulfur_phyla.pdf", width=12, height=5, dpi=300)

#Plots for genuses of Bacillus that may be sulfur oxidizers
#Create a copy of the tables for plots
major_taxa_proportions_tab_g_mc01_for_plot <- data.frame(major_taxa_proportions_tab_g_mc01)
major_taxa_proportions_tab_g_gc02_for_plot <- data.frame(major_taxa_proportions_tab_g_gc02)
major_taxa_proportions_tab_g_gc04_for_plot <- data.frame(major_taxa_proportions_tab_g_gc04)
major_taxa_proportions_tab_g_gc06_for_plot <- data.frame(major_taxa_proportions_tab_g_gc06)
head(major_taxa_proportions_tab_g_mc01)

#Simplify sample names
colnames(major_taxa_proportions_tab_g_mc01_for_plot) <- gsub("^[^_]*_", "", colnames(major_taxa_proportions_tab_g_mc01_for_plot))
colnames(major_taxa_proportions_tab_g_gc02_for_plot) <- gsub("^[^_]*_", "", colnames(major_taxa_proportions_tab_g_gc02_for_plot))
colnames(major_taxa_proportions_tab_g_gc04_for_plot) <- gsub("^[^_]*_", "", colnames(major_taxa_proportions_tab_g_gc04_for_plot))
colnames(major_taxa_proportions_tab_g_gc06_for_plot) <- gsub("^[^_]*_", "", colnames(major_taxa_proportions_tab_g_gc06_for_plot))

#Filter out samples with less than 600 total AVSs
major_taxa_proportions_tab_g_mc01_for_plot<-
  major_taxa_proportions_tab_g_mc01_for_plot[,colnames(major_taxa_proportions_tab_g_mc01_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
major_taxa_proportions_tab_g_gc02_for_plot<-
  major_taxa_proportions_tab_g_gc02_for_plot[,colnames(major_taxa_proportions_tab_g_gc02_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
major_taxa_proportions_tab_g_gc04_for_plot<-
  major_taxa_proportions_tab_g_gc04_for_plot[,colnames(major_taxa_proportions_tab_g_gc04_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
major_taxa_proportions_tab_g_gc06_for_plot<-
  major_taxa_proportions_tab_g_gc06_for_plot[,colnames(major_taxa_proportions_tab_g_gc06_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
head(major_taxa_proportions_tab_g_gc04_for_plot)

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
major_taxa_proportions_tab_g_mc01_for_plot <- 
  major_taxa_proportions_tab_g_mc01_for_plot %>% rownames_to_column("Taxa")
major_taxa_proportions_tab_g_gc02_for_plot <- 
  major_taxa_proportions_tab_g_gc02_for_plot %>% rownames_to_column("Taxa")
major_taxa_proportions_tab_g_gc04_for_plot <- 
  major_taxa_proportions_tab_g_gc04_for_plot %>% rownames_to_column("Taxa")
major_taxa_proportions_tab_g_gc06_for_plot <- 
  major_taxa_proportions_tab_g_gc06_for_plot %>% rownames_to_column("Taxa")

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
major_taxa_proportions_tab_g_mc01_for_plot <- 
  major_taxa_proportions_tab_g_mc01_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
major_taxa_proportions_tab_g_gc02_for_plot <- 
  major_taxa_proportions_tab_g_gc02_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
major_taxa_proportions_tab_g_gc04_for_plot <- 
  major_taxa_proportions_tab_g_gc04_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
major_taxa_proportions_tab_g_gc06_for_plot <- 
  major_taxa_proportions_tab_g_gc06_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)

#Add depth and core to tables
major_taxa_proportions_tab_g_mc01_for_plot <- 
  major_taxa_proportions_tab_g_mc01_for_plot %>% left_join(mod_sample_info_tab_mc01)
major_taxa_proportions_tab_g_gc02_for_plot <- 
  major_taxa_proportions_tab_g_gc02_for_plot %>% left_join(mod_sample_info_tab_gc02)
major_taxa_proportions_tab_g_gc04_for_plot <- 
  major_taxa_proportions_tab_g_gc04_for_plot %>% left_join(mod_sample_info_tab_gc04)
major_taxa_proportions_tab_g_gc06_for_plot <- 
  major_taxa_proportions_tab_g_gc06_for_plot %>% left_join(mod_sample_info_tab_gc06)

#SCB 01
#Organize things in plot
major_taxa_proportions_tab_g_mc01_for_plot$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_g_mc01_for_plot$Sample)
major_taxa_proportions_tab_g_mc01_for_plot$Sample <- factor(major_taxa_proportions_tab_g_mc01_for_plot$Sample, levels=unique(major_taxa_proportions_tab_g_mc01_for_plot$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_g_mc01_for_plot$Sample))))])

#Make plot
#Include only S cycling taxa
major_taxa_proportions_tab_g_mc01_for_plot.1 <- major_taxa_proportions_tab_g_mc01_for_plot %>%
  filter(Taxa %in%  c("unclassified_Hyphomicrobiaceae",
                      "unclassified_Bacillaceaeae_16",
                      "unclassified_Bacilli",
                      "Omnitrophales",
                      "Bacillus_13",
                      "Bacillus_1",
                      "Candidatus.Omnitrophus"))

# Plot for all depths
major_taxa_proportions_tab_g_mc01_for_plot %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "SCB 01")
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_bacillus_alldepths.pdf", width=12, height=5, dpi=300)
unique(major_taxa_proportions_tab_g_mc01_for_plot$Taxa)

#Get averages for depth_bin
major_taxa_proportions_tab_g_mc01_for_plot.2<- major_taxa_proportions_tab_g_mc01_for_plot.1 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
major_taxa_proportions_tab_g_mc01_for_plot.2 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "SCB 01")+
  scale_x_discrete(limits = c("0-20", "21-40")) +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_bacillus_genus.pdf", width=12, height=5, dpi=300)


#GC 02
#Organize things in plot
major_taxa_proportions_tab_g_gc02_for_plot$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_g_gc02_for_plot$Sample)
major_taxa_proportions_tab_g_gc02_for_plot$Sample <- factor(major_taxa_proportions_tab_g_gc02_for_plot$Sample, levels=unique(major_taxa_proportions_tab_g_gc02_for_plot$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_g_gc02_for_plot$Sample))))])

#Make plot
#Include only S cycling taxa
major_taxa_proportions_tab_g_gc02_for_plot.1 <- major_taxa_proportions_tab_g_gc02_for_plot %>%
  filter(Taxa %in%  c("unclassified_Hyphomicrobiaceae",
                      "unclassified_Bacillaceaeae_16",
                      "unclassified_Bacilli",
                      "Omnitrophales",
                      "Bacillus_13",
                      "Bacillus_1",
                      "Candidatus.Omnitrophus"))

#Get averages for depth_bin
major_taxa_proportions_tab_g_gc02_for_plot.2<- major_taxa_proportions_tab_g_gc02_for_plot.1 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
head(major_taxa_proportions_tab_g_gc02_for_plot.2)
#plot
major_taxa_proportions_tab_g_gc02_for_plot.2 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR 02")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140')) +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_bacillus_genus.pdf", width=12, height=5, dpi=300)

#GC 04
#Organize things in plot
major_taxa_proportions_tab_g_gc04_for_plot$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_g_gc04_for_plot$Sample)
major_taxa_proportions_tab_g_gc04_for_plot$Sample <- factor(major_taxa_proportions_tab_g_gc04_for_plot$Sample, levels=unique(major_taxa_proportions_tab_g_gc04_for_plot$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_g_gc04_for_plot$Sample))))])

#Make plot
#Include only S cycling taxa
major_taxa_proportions_tab_g_gc04_for_plot.1 <- major_taxa_proportions_tab_g_gc04_for_plot %>%
  filter(Taxa %in%  c("unclassified_Hyphomicrobiaceae",
                      "unclassified_Bacillaceaeae_16",
                      "unclassified_Bacilli",
                      "Omnitrophales",
                      "Bacillus_13",
                      "Bacillus_1",
                      "Candidatus.Omnitrophus"))

#Get averages for depth_bin
major_taxa_proportions_tab_g_gc04_for_plot.2<- major_taxa_proportions_tab_g_gc04_for_plot.1 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
major_taxa_proportions_tab_g_gc04_for_plot.2 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR 03")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140', '141-220')) +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_bacillus_genus.pdf", width=12, height=5, dpi=300)

#GC 06
#Organize things in plot
major_taxa_proportions_tab_g_gc06_for_plot$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_g_gc06_for_plot$Sample)
major_taxa_proportions_tab_g_gc06_for_plot$Sample <- factor(major_taxa_proportions_tab_g_gc06_for_plot$Sample, levels=unique(major_taxa_proportions_tab_g_gc06_for_plot$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_g_gc06_for_plot$Sample))))])

#Make plot
#Include only S cycling taxa
major_taxa_proportions_tab_g_gc06_for_plot.1 <- major_taxa_proportions_tab_g_gc06_for_plot %>%
  filter(Taxa %in%  c("unclassified_Hyphomicrobiaceae",
                      "unclassified_Bacillaceaeae_16",
                      "unclassified_Bacilli",
                      "Omnitrophales",
                      "Bacillus_13",
                      "Bacillus_1",
                      "Candidatus.Omnitrophus"))

#Get averages for depth_bin
major_taxa_proportions_tab_g_gc06_for_plot.2<- major_taxa_proportions_tab_g_gc06_for_plot.1 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
major_taxa_proportions_tab_g_gc06_for_plot.2 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR 04")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140', '141-220')) +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_bacillus_genus.pdf", width=12, height=5, dpi=300)


#Plots for S-cycling microbes By Order
#Create a copy of the tables for plots
major_taxa_proportions_tab_o_mc01_for_plot <- data.frame(major_taxa_proportions_tab_o_mc01)
major_taxa_proportions_tab_o_gc02_for_plot <- data.frame(major_taxa_proportions_tab_o_gc02)
major_taxa_proportions_tab_o_gc04_for_plot <- data.frame(major_taxa_proportions_tab_o_gc04)
major_taxa_proportions_tab_o_gc06_for_plot <- data.frame(major_taxa_proportions_tab_o_gc06)
head(major_taxa_proportions_tab_o_gc04)

#Simplify sample names
colnames(major_taxa_proportions_tab_o_mc01_for_plot) <- gsub("^[^_]*_", "", colnames(major_taxa_proportions_tab_o_mc01_for_plot))
colnames(major_taxa_proportions_tab_o_gc02_for_plot) <- gsub("^[^_]*_", "", colnames(major_taxa_proportions_tab_o_gc02_for_plot))
colnames(major_taxa_proportions_tab_o_gc04_for_plot) <- gsub("^[^_]*_", "", colnames(major_taxa_proportions_tab_o_gc04_for_plot))
colnames(major_taxa_proportions_tab_o_gc06_for_plot) <- gsub("^[^_]*_", "", colnames(major_taxa_proportions_tab_o_gc06_for_plot))

#Filter out samples with less than 600 total AVSs
major_taxa_proportions_tab_mc01_for_plot<-
  major_taxa_proportions_tab_mc01_for_plot[,colnames(major_taxa_proportions_tab_mc01_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
major_taxa_proportions_tab_gc02_for_plot<-
  major_taxa_proportions_tab_gc02_for_plot[,colnames(major_taxa_proportions_tab_gc02_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
major_taxa_proportions_tab_gc04_for_plot<-
  major_taxa_proportions_tab_gc04_for_plot[,colnames(major_taxa_proportions_tab_gc04_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
major_taxa_proportions_tab_gc06_for_plot<-
  major_taxa_proportions_tab_gc06_for_plot[,colnames(major_taxa_proportions_tab_gc06_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
head(major_taxa_proportions_tab_gc04_for_plot)

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
major_taxa_proportions_tab_o_mc01_for_plot <- 
  major_taxa_proportions_tab_o_mc01_for_plot %>% rownames_to_column("Taxa")
major_taxa_proportions_tab_o_gc02_for_plot <- 
  major_taxa_proportions_tab_o_gc02_for_plot %>% rownames_to_column("Taxa")
major_taxa_proportions_tab_o_gc04_for_plot <- 
  major_taxa_proportions_tab_o_gc04_for_plot %>% rownames_to_column("Taxa")
major_taxa_proportions_tab_o_gc06_for_plot <- 
  major_taxa_proportions_tab_o_gc06_for_plot %>% rownames_to_column("Taxa")

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
major_taxa_proportions_tab_o_mc01_for_plot <- 
  major_taxa_proportions_tab_o_mc01_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
major_taxa_proportions_tab_o_gc02_for_plot <- 
  major_taxa_proportions_tab_o_gc02_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
major_taxa_proportions_tab_o_gc04_for_plot <- 
  major_taxa_proportions_tab_o_gc04_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
major_taxa_proportions_tab_o_gc06_for_plot <- 
  major_taxa_proportions_tab_o_gc06_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)

#Add depth and core to tables
major_taxa_proportions_tab_o_mc01_for_plot <- 
  major_taxa_proportions_tab_o_mc01_for_plot %>% left_join(mod_sample_info_tab_mc01)
major_taxa_proportions_tab_o_gc02_for_plot <- 
  major_taxa_proportions_tab_o_gc02_for_plot %>% left_join(mod_sample_info_tab_gc02)
major_taxa_proportions_tab_o_gc04_for_plot <- 
  major_taxa_proportions_tab_o_gc04_for_plot %>% left_join(mod_sample_info_tab_gc04)
major_taxa_proportions_tab_o_gc06_for_plot <- 
  major_taxa_proportions_tab_o_gc06_for_plot %>% left_join(mod_sample_info_tab_gc06)

#SCB 01
#Organize things in plot
major_taxa_proportions_tab_o_mc01_for_plot$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_o_mc01_for_plot$Sample)
major_taxa_proportions_tab_o_mc01_for_plot$Sample <- factor(major_taxa_proportions_tab_o_mc01_for_plot$Sample, levels=unique(major_taxa_proportions_tab_o_mc01_for_plot$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_o_mc01_for_plot$Sample))))])
#major_taxa_proportions_tab_mc01_for_plot$Genus <- with(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2,factor(Class,Genus = rev(sort(unique(Genus)))))

#Make plot
#Include only S cycling taxa
major_taxa_proportions_tab_o_mc01_for_plot.1 <- major_taxa_proportions_tab_o_mc01_for_plot %>%
  filter(Taxa %in%  c("Rhodobacterales",
                      "Rhodospirillales",
                      "Sphingomonadales",
                      "unclassified_Alphaproteobacteria",
                      "unclassified_Gammaproteobacteria",
                      "unclassified_Proteobacteria",
                      "Oceanospirillales_12",
                      "Oceanospirillales_8",
                      "Acetobacterales",
                      "Burkholderiales",
                      "Thiomicrospirales",
                      "Pseudomonadales"))
major_taxa_proportions_tab_o_mc01_for_plot.2 <- major_taxa_proportions_tab_o_mc01_for_plot.1 %>% 
  mutate(Taxa = if_else(Taxa %ni% c("Oceanospirillales_12",
                                    "Oceanospirillales_8") , Taxa, "Oceanospirillales"))
#For order
major_taxa_proportions_tab_o_mc01_for_plot.1 %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "SCB 01")
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_sulfur_proteo_orders_alldepths.pdf", width=12, height=5, dpi=300)
a<-unique(major_taxa_proportions_tab_o_mc01_for_plot$Taxa)
a<-relist(sort(unlist(a)), a)

#Get averages for depth_bin
major_taxa_proportions_tab_o_mc01_for_plot.3<- major_taxa_proportions_tab_o_mc01_for_plot.2 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
major_taxa_proportions_tab_o_mc01_for_plot.3 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "SCB 01")+
  scale_x_discrete(limits = c("0-20", "21-40")) +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_sulfur_proteo_orders.pdf", width=12, height=5, dpi=300)

#Plot for clades inside Phylum, all depths
#Filter by Actinobacteria orders
major_taxa_proportions_tab_o_mc01_actinobacteria <- major_taxa_proportions_tab_o_mc01_for_plot %>%
  filter(Taxa %in%  c("unclassified_Actinobacteria",
                      "unclassified_Actinobacteriota",
                      "unclassified_Coriobacteriia",
                      "Actinomycetales",
                      "Actinomarinales",
                      "Corynebacteriales",
                      "Frankiales_1",
                      "Gaiellales",
                      "Propionibacteriales",
                      "Pseudonocardiales",
                      "Rubrobacterales",
                      "Solirubrobacterales",
                      "FS118.23B.02",
                      "FS117.23B.02",
                      "CG2.30.50.142",
                      "Streptosporangiales",
                      "WCHB1.81",
                      "Acidimicrobiales"))
#Make plot
major_taxa_proportions_tab_o_mc01_actinobacteria %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Actinobacteria orders", title = "SCB 01") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_actinobacteria_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Planctomycetota orders
major_taxa_proportions_tab_o_mc01_planctomycetota <- major_taxa_proportions_tab_o_mc01_for_plot %>%
  filter(Taxa %in%  c("Gemmatales",
                      "Planctomycetales",
                      "Phycisphaerales",
                      "Pla1.lineage",
                      "Pla4.lineage",
                      "Pirellulales",
                      "MSBL9",
                      "SPG12.343.353.B69"))
#Make plot
major_taxa_proportions_tab_o_mc01_planctomycetota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Planctomycetota orders", title = "SCB 01") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_planctomycetota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Chloroflexi orders
major_taxa_proportions_tab_o_mc01_chloroflexi <- major_taxa_proportions_tab_o_mc01_for_plot %>%
  filter(Taxa %in%  c("unclassified_Chloroflexi",
                      "unclassified_Dehalococcoidia",
                      "unclassified_Anaerolineae",
                      "Anaerolineales",
                      "Ardenticatenales",
                      "Caldilineales",
                      "Chloroflexales",
                      "SBR1031",
                      "RBG.13.54.9",
                      "Napoli.4B.65",
                      "S085",
                      "Sh765B.AG.111",
                      "DscP2",
                      "Thermoflexales",
                      "vadinHA49"))
#Make plot
major_taxa_proportions_tab_o_mc01_chloroflexi %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Chloroflexi orders", title = "SCB 01") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_chloroflexi_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Crenarchaeota orders
major_taxa_proportions_tab_o_mc01_crenarchaeota <- major_taxa_proportions_tab_o_mc01_for_plot %>%
  filter(Taxa %in%  c("Nitrososphaerales",
                      "Bathyarchaeia",
                      "Caldiarchaeales",
                      "Thermoproteales", 
                      "Sulfolobales",
                      "Desulfurococcales",
                      "Thermococcales",
                      "Thermoplasmatales"))
#Make plot
major_taxa_proportions_tab_o_mc01_crenarchaeota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Crenarchaeota orders", title = "SCB 01") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_crenarchaeota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Acidobacteriota orders
major_taxa_proportions_tab_o_mc01_acidobacteriota <- major_taxa_proportions_tab_o_mc01_for_plot %>%
  filter(Taxa %in%  c("Pyrinomonadales",
                      "Subgroup.9",
                      "Thermoanaerobaculales",
                      "Acidobacteriales",
                      "Acidoferrales",
                      "Bryobacterales",
                      "Blastocatellales",
                      "Vicinamibacteriales"))
#Make plot
major_taxa_proportions_tab_o_mc01_acidobacteriota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Acidobacteriota orders", title = "SCB 01") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_acidobacteriota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Bacteroidota orders
major_taxa_proportions_tab_o_mc01_bacteroidota <- major_taxa_proportions_tab_o_mc01_for_plot %>%
  filter(Taxa %in%  c("Bacteroidales_1",
                      "Cytophagales_1",
                      "Flavobacteriales",
                      "Ignavibacteriales",
                      "Sphingobacteriales",
                      "Chitinophagales",
                      "Melioribacterales",
                      "Rhodotermales",
                      "Chlorobiales"))
#Make plot
major_taxa_proportions_tab_o_mc01_bacteroidota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Bacteroidota orders", title = "SCB 01") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_bacteroidota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Firmicutes orders
major_taxa_proportions_tab_o_mc01_firmicutes <- major_taxa_proportions_tab_o_mc01_for_plot %>%
  filter(Taxa %in%  c("Clostridiales",
                      "Bacillales_13",
                      "Bacillales_14",
                      "Erysipelotrichales_2",
                      "Lachnospirales",
                      "Lactobacillales",
                      "H3.93"))
#Make plot
major_taxa_proportions_tab_o_mc01_firmicutes %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Firmicutes orders", title = "SCB 01") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_firmicutes_orders_alldepths.pdf", width=12, height=5, dpi=300)

#GC 02
#Organize things in plot
major_taxa_proportions_tab_o_gc02_for_plot$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_o_gc02_for_plot$Sample)
major_taxa_proportions_tab_o_gc02_for_plot$Sample <- factor(major_taxa_proportions_tab_o_gc02_for_plot$Sample, levels=unique(major_taxa_proportions_tab_o_gc02_for_plot$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_o_gc02_for_plot$Sample))))])
#major_taxa_proportions_tab_gc02_for_plot$Genus <- with(proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2,factor(Class,Genus = rev(sort(unique(Genus)))))

#Make plot
#Include only S cycling taxa
major_taxa_proportions_tab_o_gc02_for_plot.1 <- major_taxa_proportions_tab_o_gc02_for_plot %>%
  filter(Taxa %in%  c("Rhodobacterales",
                      "Rhodospirillales",
                      "Sphingomonadales",
                      "unclassified_Alphaproteobacteria",
                      "unclassified_Gammaproteobacteria",
                      "unclassified_Proteobacteria",
                      "Oceanospirillales_12",
                      "Oceanospirillales_8",
                      "Acetobacterales",
                      "Burkholderiales",
                      "Thiomicrospirales",
                      "Pseudomonadales"))
major_taxa_proportions_tab_o_gc02_for_plot.2 <- major_taxa_proportions_tab_o_gc02_for_plot.1 %>% 
  mutate(Taxa = if_else(Taxa %ni% c("Oceanospirillales_12",
                                    "Oceanospirillales_8") , Taxa, "Oceanospirillales"))
#For order
major_taxa_proportions_tab_o_gc02_for_plot.1 %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR 02")
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_sulfur_proteo_orders_alldepths.pdf", width=12, height=5, dpi=300)


#Get averages for depth_bin
major_taxa_proportions_tab_o_gc02_for_plot.3<- major_taxa_proportions_tab_o_gc02_for_plot.2 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
major_taxa_proportions_tab_o_gc02_for_plot.3 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_gc02) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "GC 02")+
  scale_x_discrete(limits = c("0-20", "21-40")) +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_sulfur_proteo_orders.pdf", width=12, height=5, dpi=300)

#Plot for clades inside Phylum, all depths
orders_gc02<-unique(major_taxa_proportions_tab_o_gc02_for_plot$Taxa)
orders_gc02<-relist(sort(unlist(orders_gc02)), orders_gc02)
#Filter by Actinobacteria orders
major_taxa_proportions_tab_o_gc02_actinobacteria <- major_taxa_proportions_tab_o_gc02_for_plot %>%
  filter(Taxa %in%  c("unclassified_Actinobacteria",
                      "unclassified_Actinobacteriota",
                      "unclassified_Coriobacteriia",
                      "Actinomycetales",
                      "Actinomarinales",
                      "Corynebacteriales",
                      "Frankiales_1",
                      "Frankiales_2",
                      "Gaiellales",
                      "Propionibacteriales",
                      "Pseudonocardiales",
                      "Rubrobacterales",
                      "Solirubrobacterales",
                      "FS118.23B.02",
                      "FS117.23B.02",
                      "CG2.30.50.142",
                      "Streptosporangiales",
                      "WCHB1.81",
                      "Acidimicrobiales"))
#Make plot
major_taxa_proportions_tab_o_gc02_actinobacteria %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Actinobacteria orders", title = "CR 02") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_actinobacteria_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Planctomycetota orders
major_taxa_proportions_tab_o_gc02_planctomycetota <- major_taxa_proportions_tab_o_gc02_for_plot %>%
  filter(Taxa %in%  c("Gemmatales",
                      "Planctomycetales",
                      "Phycisphaerales",
                      "Pla1.lineage",
                      "Pla3.lineage",
                      "Pla4.lineage",
                      "Pirellulales",
                      "MSBL9",
                      "SPG12.343.353.B69",
                      "Brocadiales",
                      "unclassified_Planctomycetes",
                      "unclassified_Phycisphaerae",
                      "unclassified_Planctomycetota"))
#Make plot
major_taxa_proportions_tab_o_gc02_planctomycetota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Planctomycetota orders", title = "CR 02") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_planctomycetota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Chloroflexi orders
major_taxa_proportions_tab_o_gc02_chloroflexi <- major_taxa_proportions_tab_o_gc02_for_plot %>%
  filter(Taxa %in%  c("unclassified_Chloroflexi",
                      "unclassified_Dehalococcoidia",
                      "unclassified_Anaerolineae",
                      "Anaerolineales",
                      "Ardenticatenales",
                      "Caldilineales",
                      "Chloroflexales",
                      "SBR1031",
                      "RBG.13.54.9",
                      "Napoli.4B.65",
                      "S085",
                      "Sh765B.AG.111",
                      "DscP2",
                      "Thermoflexales",
                      "vadinHA49",
                      "ADurb.Bin180",
                      "Dehalococcoidia"))
#Make plot
major_taxa_proportions_tab_o_gc02_chloroflexi %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Chloroflexi orders", title = "CR 02") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_chloroflexi_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Crenarchaeota orders
major_taxa_proportions_tab_o_gc02_crenarchaeota <- major_taxa_proportions_tab_o_gc02_for_plot %>%
  filter(Taxa %in%  c("Nitrososphaerales",
                      "Bathyarchaeia",
                      "Caldiarchaeales",
                      "Thermoproteales", 
                      "Sulfolobales",
                      "Desulfurococcales",
                      "Thermococcales",
                      "Thermoplasmatales"))
#Make plot
major_taxa_proportions_tab_o_gc02_crenarchaeota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Crenarchaeota orders", title = "CR 02") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_crenarchaeota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Acidobacteriota orders
major_taxa_proportions_tab_o_gc02_acidobacteriota <- major_taxa_proportions_tab_o_gc02_for_plot %>%
  filter(Taxa %in%  c("Pyrinomonadales",
                      "Subgroup.9",
                      "Subgroup.15",
                      "Subgroup.17",
                      "Subgroup.2" ,
                      "Subgroup.22",
                      "Subgroup.7",
                      "Subgroup.26",
                      "Subgroup.21",
                      "Thermoanaerobaculales",
                      "Acidobacteriales",
                      "Acidoferrales",
                      "Bryobacterales",
                      "Blastocatellales",
                      "Vicinamibacteriales",
                      "unclassified_Vicinamibacteria",
                      "unclassified_Acidobacteriota"))
#Make plot
major_taxa_proportions_tab_o_gc02_acidobacteriota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Acidobacteriota orders", title = "CR 02") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_acidobacteriota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Bacteroidota orders
major_taxa_proportions_tab_o_gc02_bacteroidota <- major_taxa_proportions_tab_o_gc02_for_plot %>%
  filter(Taxa %in%  c("Bacteroidales_1",
                      "Bacteroidales_2",
                      "Cytophagales_1",
                      "Flavobacteriales",
                      "Ignavibacteriales",
                      "Sphingobacteriales",
                      "Chitinophagales",
                      "Melioribacterales",
                      "Rhodotermales",
                      "Chlorobiales"))
#Make plot
major_taxa_proportions_tab_o_gc02_bacteroidota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Bacteroidota orders", title = "CR 02") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_bacteroidota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Firmicutes orders
major_taxa_proportions_tab_o_gc02_firmicutes <- major_taxa_proportions_tab_o_gc02_for_plot %>%
  filter(Taxa %in%  c("Clostridiales",
                      "Bacillales_13",
                      "Bacillales_14",
                      "Bacillales_15",
                      "Erysipelotrichales_2",
                      "Lachnospirales",
                      "Lactobacillales",
                      "H3.93",
                      "C86",
                      "unclassified_Firmicutes"))
#Make plot
major_taxa_proportions_tab_o_gc02_firmicutes %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Firmicutes orders", title = "CR 02") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_firmicutes_orders_alldepths.pdf", width=12, height=5, dpi=300)

#GC 04
#Organize things in plot
major_taxa_proportions_tab_o_gc04_for_plot$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_o_gc04_for_plot$Sample)
major_taxa_proportions_tab_o_gc04_for_plot$Sample <- factor(major_taxa_proportions_tab_o_gc04_for_plot$Sample, levels=unique(major_taxa_proportions_tab_o_gc04_for_plot$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_o_gc04_for_plot$Sample))))])
#major_taxa_proportions_tab_gc04_for_plot$Genus <- with(proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2,factor(Class,Genus = rev(sort(unique(Genus)))))

#Make plot
#Include only S cycling taxa
major_taxa_proportions_tab_o_gc04_for_plot.1 <- major_taxa_proportions_tab_o_gc04_for_plot %>%
  filter(Taxa %in%  c("Rhodobacterales",
                      "Rhodospirillales",
                      "Sphingomonadales",
                      "unclassified_Alphaproteobacteria",
                      "unclassified_Gammaproteobacteria",
                      "unclassified_Proteobacteria",
                      "Oceanospirillales_12",
                      "Oceanospirillales_8",
                      "Acetobacterales",
                      "Burkholderiales",
                      "Thiomicrospirales",
                      "Pseudomonadales"))
major_taxa_proportions_tab_o_gc04_for_plot.2 <- major_taxa_proportions_tab_o_gc04_for_plot.1 %>% 
  mutate(Taxa = if_else(Taxa %ni% c("Oceanospirillales_12",
                                    "Oceanospirillales_8") , Taxa, "Oceanospirillales"))
#For order
major_taxa_proportions_tab_o_gc04_for_plot.2 %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR 03")
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_sulfur_proteo_orders_alldepths.pdf", width=12, height=5, dpi=300)


#Get averages for depth_bin
major_taxa_proportions_tab_o_gc04_for_plot.3<- major_taxa_proportions_tab_o_gc04_for_plot.2 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
major_taxa_proportions_tab_o_gc04_for_plot.3 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_gc04) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "GC 04")+
  scale_x_discrete(limits = c("0-20", "21-40")) +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_sulfur_proteo_orders.pdf", width=12, height=5, dpi=300)

#Plot for clades inside Phylum, all depths
orders_gc04<-unique(major_taxa_proportions_tab_o_gc04_for_plot$Taxa)
orders_gc04<-relist(sort(unlist(orders_gc04)), orders_gc04)
#Filter by Actinobacteria orders
major_taxa_proportions_tab_o_gc04_actinobacteria <- major_taxa_proportions_tab_o_gc04_for_plot %>%
  filter(Taxa %in%  c("unclassified_Actinobacteria",
                      "unclassified_Actinobacteriota",
                      "unclassified_Coriobacteriia",
                      "Actinomycetales",
                      "Actinomarinales",
                      "Corynebacteriales",
                      "Frankiales_1",
                      "Frankiales_2",
                      "Gaiellales",
                      "Propionibacteriales",
                      "Pseudonocardiales",
                      "Rubrobacterales",
                      "Solirubrobacterales",
                      "FS118.23B.02",
                      "FS117.23B.02",
                      "CG2.30.50.142",
                      "Streptosporangiales",
                      "WCHB1.81",
                      "Acidimicrobiales"))
#Make plot
major_taxa_proportions_tab_o_gc04_actinobacteria %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Actinobacteria orders", title = "CR 03") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_actinobacteria_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Planctomycetota orders
major_taxa_proportions_tab_o_gc04_planctomycetota <- major_taxa_proportions_tab_o_gc04_for_plot %>%
  filter(Taxa %in%  c("Gemmatales",
                      "Planctomycetales",
                      "Phycisphaerales",
                      "Pla1.lineage",
                      "Pla3.lineage",
                      "Pla4.lineage",
                      "Pirellulales",
                      "MSBL9",
                      "SPG12.343.353.B69",
                      "Brocadiales",
                      "unclassified_Planctomycetes",
                      "unclassified_Phycisphaerae",
                      "unclassified_Planctomycetota"))
#Make plot
major_taxa_proportions_tab_o_gc04_planctomycetota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Planctomycetota orders", title = "CR 03") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_planctomycetota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Chloroflexi orders
major_taxa_proportions_tab_o_gc04_chloroflexi <- major_taxa_proportions_tab_o_gc04_for_plot %>%
  filter(Taxa %in%  c("unclassified_Chloroflexi",
                      "unclassified_Dehalococcoidia",
                      "unclassified_Anaerolineae",
                      "Anaerolineales",
                      "Ardenticatenales",
                      "Caldilineales",
                      "Chloroflexales",
                      "SBR1031",
                      "RBG.13.54.9",
                      "Napoli.4B.65",
                      "S085",
                      "Sh765B.AG.111",
                      "DscP2",
                      "Thermoflexales",
                      "vadinHA49",
                      "ADurb.Bin180",
                      "Dehalococcoidia"))
#Make plot
major_taxa_proportions_tab_o_gc04_chloroflexi %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Chloroflexi orders", title = "CR 03") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_chloroflexi_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Crenarchaeota orders
major_taxa_proportions_tab_o_gc04_crenarchaeota <- major_taxa_proportions_tab_o_gc04_for_plot %>%
  filter(Taxa %in%  c("Nitrososphaerales",
                      "Bathyarchaeia",
                      "Caldiarchaeales",
                      "Thermoproteales", 
                      "Sulfolobales",
                      "Desulfurococcales",
                      "Thermococcales",
                      "Thermoplasmatales"))
#Make plot
major_taxa_proportions_tab_o_gc04_crenarchaeota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Crenarchaeota orders", title = "CR 03") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_crenarchaeota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Acidobacteriota orders
major_taxa_proportions_tab_o_gc04_acidobacteriota <- major_taxa_proportions_tab_o_gc04_for_plot %>%
  filter(Taxa %in%  c("Pyrinomonadales",
                      "Subgroup.9",
                      "Subgroup.15",
                      "Subgroup.17",
                      "Subgroup.2" ,
                      "Subgroup.22",
                      "Subgroup.7",
                      "Subgroup.26",
                      "Subgroup.21",
                      "Subgroup.13",
                      "Thermoanaerobaculales",
                      "Acidobacteriales",
                      "Acidoferrales",
                      "Bryobacterales",
                      "Blastocatellales",
                      "Vicinamibacteriales",
                      "unclassified_Vicinamibacteria",
                      "unclassified_Acidobacteriota"))
#Make plot
major_taxa_proportions_tab_o_gc04_acidobacteriota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Acidobacteriota orders", title = "CR 03") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_acidobacteriota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Bacteroidota orders
major_taxa_proportions_tab_o_gc04_bacteroidota <- major_taxa_proportions_tab_o_gc04_for_plot %>%
  filter(Taxa %in%  c("Bacteroidales_1",
                      "Bacteroidales_2",
                      "Cytophagales_1",
                      "Flavobacteriales",
                      "Ignavibacteriales",
                      "Sphingobacteriales",
                      "Chitinophagales",
                      "Melioribacterales",
                      "Rhodotermales",
                      "Chlorobiales"))
#Make plot
major_taxa_proportions_tab_o_gc04_bacteroidota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Bacteroidota orders", title = "CR 03") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_bacteroidota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Firmicutes orders
major_taxa_proportions_tab_o_gc04_firmicutes <- major_taxa_proportions_tab_o_gc04_for_plot %>%
  filter(Taxa %in%  c("Clostridiales",
                      "Bacillales_13",
                      "Bacillales_14",
                      "Bacillales_15",
                      "Bacillales_11",
                      "Bacillales_6",
                      "Erysipelotrichales_2",
                      "Lachnospirales",
                      "Lactobacillales",
                      "H3.93",
                      "C86",
                      "unclassified_Firmicutes",
                      "Oscillospirales"))
#Make plot
major_taxa_proportions_tab_o_gc04_firmicutes %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Firmicutes orders", title = "CR 03") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_firmicutes_orders_alldepths.pdf", width=12, height=5, dpi=300)

#GC 06
#Organize things in plot
major_taxa_proportions_tab_o_gc06_for_plot$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_o_gc06_for_plot$Sample)
major_taxa_proportions_tab_o_gc06_for_plot$Sample <- factor(major_taxa_proportions_tab_o_gc06_for_plot$Sample, levels=unique(major_taxa_proportions_tab_o_gc06_for_plot$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_o_gc06_for_plot$Sample))))])
#major_taxa_proportions_tab_gc06_for_plot$Genus <- with(proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2,factor(Class,Genus = rev(sort(unique(Genus)))))

#Make plot
#Include only S cycling taxa
major_taxa_proportions_tab_o_gc06_for_plot.1 <- major_taxa_proportions_tab_o_gc06_for_plot %>%
  filter(Taxa %in%  c("Rhodobacterales",
                      "Rhodospirillales",
                      "Sphingomonadales",
                      "unclassified_Alphaproteobacteria",
                      "unclassified_Gammaproteobacteria",
                      "unclassified_Proteobacteria",
                      "Oceanospirillales_12",
                      "Oceanospirillales_8",
                      "Acetobacterales",
                      "Burkholderiales",
                      "Thiomicrospirales",
                      "Pseudomonadales"))
major_taxa_proportions_tab_o_gc06_for_plot.2 <- major_taxa_proportions_tab_o_gc06_for_plot.1 %>% 
  mutate(Taxa = if_else(Taxa %ni% c("Oceanospirillales_12",
                                    "Oceanospirillales_8") , Taxa, "Oceanospirillales"))
#For order
major_taxa_proportions_tab_o_gc06_for_plot.2 %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR 04")
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_sulfur_proteo_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Get averages for depth_bin
major_taxa_proportions_tab_o_gc06_for_plot.3<- major_taxa_proportions_tab_o_gc06_for_plot.2 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
major_taxa_proportions_tab_o_gc06_for_plot.3 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_gc06) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "GC 06")+
  scale_x_discrete(limits = c("0-20", "21-40")) +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_sulfur_proteo_orders.pdf", width=12, height=5, dpi=300)

#Plot for clades inside Phylum, all depths
orders_gc06<-unique(major_taxa_proportions_tab_o_gc06_for_plot$Taxa)
orders_gc06<-relist(sort(unlist(orders_gc06)), orders_gc06)
#Filter by Actinobacteria orders
major_taxa_proportions_tab_o_gc06_actinobacteria <- major_taxa_proportions_tab_o_gc06_for_plot %>%
  filter(Taxa %in%  c("unclassified_Actinobacteria",
                      "unclassified_Actinobacteriota",
                      "unclassified_Coriobacteriia",
                      "Actinomycetales",
                      "Actinomarinales",
                      "Corynebacteriales",
                      "Frankiales_1",
                      "Frankiales_2",
                      "Gaiellales",
                      "Propionibacteriales",
                      "Pseudonocardiales",
                      "Rubrobacterales",
                      "Solirubrobacterales",
                      "FS118.23B.02",
                      "FS117.23B.02",
                      "CG2.30.50.142",
                      "Streptosporangiales",
                      "WCHB1.81",
                      "Acidimicrobiales"))
#Make plot
major_taxa_proportions_tab_o_gc06_actinobacteria %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Actinobacteria orders", title = "CR 04") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_actinobacteria_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Planctomycetota orders
major_taxa_proportions_tab_o_gc06_planctomycetota <- major_taxa_proportions_tab_o_gc06_for_plot %>%
  filter(Taxa %in%  c("Gemmatales",
                      "Planctomycetales",
                      "Phycisphaerales",
                      "Pla1.lineage",
                      "Pla3.lineage",
                      "Pla4.lineage",
                      "Pirellulales",
                      "MSBL9",
                      "SPG12.343.353.B69",
                      "Brocadiales",
                      "unclassified_Planctomycetes",
                      "unclassified_Phycisphaerae",
                      "unclassified_Planctomycetota"))
#Make plot
major_taxa_proportions_tab_o_gc06_planctomycetota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Planctomycetota orders", title = "CR 04") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_planctomycetota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Chloroflexi orders
major_taxa_proportions_tab_o_gc06_chloroflexi <- major_taxa_proportions_tab_o_gc06_for_plot %>%
  filter(Taxa %in%  c("unclassified_Chloroflexi",
                      "unclassified_Dehalococcoidia",
                      "unclassified_Anaerolineae",
                      "Anaerolineales",
                      "Ardenticatenales",
                      "Caldilineales",
                      "Chloroflexales",
                      "SBR1031",
                      "RBG.13.54.9",
                      "Napoli.4B.65",
                      "S085",
                      "Sh765B.AG.111",
                      "DscP2",
                      "Thermoflexales",
                      "vadinHA49",
                      "ADurb.Bin180",
                      "Dehalococcoidia"))
#Make plot
major_taxa_proportions_tab_o_gc06_chloroflexi %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Chloroflexi orders", title = "CR 04") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_chloroflexi_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Crenarchaeota orders
major_taxa_proportions_tab_o_gc06_crenarchaeota <- major_taxa_proportions_tab_o_gc06_for_plot %>%
  filter(Taxa %in%  c("Nitrososphaerales",
                      "Bathyarchaeia",
                      "Caldiarchaeales",
                      "Thermoproteales", 
                      "Sulfolobales",
                      "Desulfurococcales",
                      "Thermococcales",
                      "Thermoplasmatales"))
#Make plot
major_taxa_proportions_tab_o_gc06_crenarchaeota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Crenarchaeota orders", title = "CR 04") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_crenarchaeota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Acidobacteriota orders
major_taxa_proportions_tab_o_gc06_acidobacteriota <- major_taxa_proportions_tab_o_gc06_for_plot %>%
  filter(Taxa %in%  c("Pyrinomonadales",
                      "Subgroup.9",
                      "Subgroup.15",
                      "Subgroup.17",
                      "Subgroup.2" ,
                      "Subgroup.22",
                      "Subgroup.7",
                      "Subgroup.26",
                      "Subgroup.21",
                      "Subgroup.13",
                      "Thermoanaerobaculales",
                      "Acidobacteriales",
                      "Acidoferrales",
                      "Bryobacterales",
                      "Blastocatellales",
                      "Vicinamibacteriales",
                      "unclassified_Vicinamibacteria",
                      "unclassified_Acidobacteriota"))
#Make plot
major_taxa_proportions_tab_o_gc06_acidobacteriota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Acidobacteriota orders", title = "CR 04") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_acidobacteriota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Bacteroidota orders
major_taxa_proportions_tab_o_gc06_bacteroidota <- major_taxa_proportions_tab_o_gc06_for_plot %>%
  filter(Taxa %in%  c("Bacteroidales_1",
                      "Bacteroidales_2",
                      "Cytophagales_1",
                      "Flavobacteriales",
                      "Ignavibacteriales",
                      "Sphingobacteriales",
                      "Chitinophagales",
                      "Melioribacterales",
                      "Rhodotermales",
                      "Chlorobiales"))
#Make plot
major_taxa_proportions_tab_o_gc06_bacteroidota %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Bacteroidota orders", title = "CR 04") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_bacteroidota_orders_alldepths.pdf", width=12, height=5, dpi=300)

#Filter by Firmicutes orders
major_taxa_proportions_tab_o_gc06_firmicutes <- major_taxa_proportions_tab_o_gc06_for_plot %>%
  filter(Taxa %in%  c("Clostridiales",
                      "Bacillales_13",
                      "Bacillales_14",
                      "Bacillales_15",
                      "Bacillales_11",
                      "Bacillales_6",
                      "Bacillales_1",
                      "Erysipelotrichales_2",
                      "Lachnospirales",
                      "Lactobacillales",
                      "H3.93",
                      "C86",
                      "unclassified_Firmicutes",
                      "Oscillospirales",
                      "Peptostreptococcales.Tissierellales"))
#Make plot
major_taxa_proportions_tab_o_gc06_firmicutes %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from Firmicutes orders", title = "CR 04") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_firmicutes_orders_alldepths.pdf", width=12, height=5, dpi=300)







#For family
#Make plot
filt_major_taxa_proportions_tab_for_plot_f_mc01.g2 %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
  ggplot(aes(x = Sample, y = Proportion, fill = Family)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = stock_colors) +
  theme_bw() +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=14), legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% of 16S rRNA gene copies recovered", title = "San Clemente Basin multicore") +
  scale_x_discrete(labels = c('ext_35c', 'ext_45c', 'pcr_35c', 'pcr_45c')) 
#For genus
#Group all the unclassifieds and uncultured together
`%ni%` <- Negate(`%in%`)
mc_01_genus_unc <- filt_major_taxa_proportions_tab_for_plot_g_mc01.g2 %>% 
  mutate(Genus = if_else(Genus %ni% c("uncultured.1",
                                      "uncultured.2",
                                      "unclassified_Root") , Genus, "Unclassified"))
mc_01_genus_unc_g1 <- mc_01_genus_unc %>% 
  mutate(Genus = if_else(Genus %ni% c("AT.s2.59",
                                      "Filomicrobium_1",
                                      "Hydrothermarchaeales",
                                      "S085",
                                      "SG8.4"
  ) , Genus, "Other"))
mc_01_genus_unc_g2 <- mc_01_genus_unc_g1 %>% 
  mutate(Genus = if_else(Genus %ni% c("Aminicenantales",
                                      "JS1") , Genus, "Aminicenantales+JS1"))
mc_01_genus_unc_g3 <- mc_01_genus_unc_g2 %>% 
  mutate(Genus = if_else(Genus %ni% c("Pir4.lineage",
                                      "unclassified_Pirellulaceae"), Genus, "Pir4 lineage"))

mc_01_genus_unc_g4 <- mc_01_genus_unc_g3 %>% 
  mutate(Genus = if_else(Genus %ni% c("uncultured"), Genus, "Pedomicrobium"))

mc_01_genus_unc_g5 <- mc_01_genus_unc_g4 %>% 
  mutate(Genus = if_else(Genus %ni% c('Nitrosomonas',
                                      'Nitrosopumilaceae',
                                      'Candidatus.Nitrosopumilus'), Genus, "Nitrosomonas + Nitrosopumilus"))

mc_01_genus_unc_g6 <- mc_01_genus_unc_g5 %>% 
  mutate(Genus = if_else(Genus %ni% c('Subgroup.10',
                                      'Subgroup.22',
                                      'Subgroup.23',
                                      'Subgroup.9'), Genus, "Acidobacteria"))

mc_01_genus_unc_g7 <- mc_01_genus_unc_g6 %>% 
  mutate(Genus = if_else(Genus %ni% c('Bathyarchaeia',
                                      'Lokiarchaeia'), Genus, "Bathyarchaeia + Lokiarchaeia"))
#Make plot
#Filter unclassifieds out
mc_01_genus_unc_g7 <- mc_01_genus_unc_g7 %>%
  #filter(Genus != "Unclassified")
  filt_major_taxa_proportions_tab_for_plot_g_mc01 %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
  ggplot(aes(x = Sample, y = Proportion, fill = Genus)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = bold_mc01) +
  theme_bw() +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=14), legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% of 16S rRNA gene copies recovered", title = "San Clemente Basin multicore") +
  scale_x_discrete(labels = c(2, 4, 8, 12, 18, 23)) 

filt_major_taxa_proportions_tab_for_plot_g_mc01.g2 %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
  ggplot(aes(x = Sample, y = Proportion, fill = Genus)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = stock_colors) +
  theme_bw() +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=14), legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% of 16S rRNA gene copies recovered", title = "San Clemente Basin multicore") +
  scale_x_discrete(labels = c('ext_35c', 'ext_45c', 'pcr_35c', 'pcr_45c')) 

#For GC02
#Define color palette
gc_02_Palette <- c("#ff6db6", "#ffff6d", "#b66dff", "#004949", "#b6dbff", "#6db6ff",
                   "#006ddb", "#009292","#db6d00","#ffb6db", "#490092", "#924900")

#Make plot
filt_major_taxa_proportions_tab_for_plot_gc02.g2 %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
  ggplot(aes(x = Sample, y = Proportion, fill = Phylum)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = stock_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% of 16S rRNA gene copies recovered", title = "GC02_CR02") +
  geom_bar(stat = "identity", width = 0.5) +
  scale_x_discrete(labels = c(3, 15, 18, 33, 36, 42, 48, 54, 60, 69, 78, 81, 87, 96, 102, 108, 111, 123)) 
#For family
filt_major_taxa_proportions_tab_for_plot_f_gc02.g2 %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
  ggplot(aes(x = Sample, y = Proportion, fill = Family)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = stock_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% of 16S rRNA gene copies recovered", title = "GC02_CR02") +
  geom_bar(stat = "identity", width = 0.5) +
  scale_x_discrete(labels = c(3, 15, 18, 33, 36, 42, 48, 54, 60, 69, 78, 81, 87, 96, 102, 108, 111, 123)) 
#For genus
#Group all the unclassifieds and uncultured together
`%ni%` <- Negate(`%in%`)
gc_02_genus_unc <- filt_major_taxa_proportions_tab_for_plot_g_gc02.g2 %>% 
  mutate(Genus = if_else(Genus %ni% c("uncultured.1",
                                      "uncultured.2",
                                      "uncultured.3",
                                      "uncultured",
                                      "unclassified_Root",
                                      "unclassified_Bacteria",
                                      "unclassified_Halobacterales",
                                      "Unclassified") , Genus, "Unclassified"))
gc_02_genus_unc_g1 <- gc_02_genus_unc %>% 
  mutate(Genus = if_else(Genus %ni% c("DscP2",
                                      "SG8.4",
                                      "WCHB1.81") , Genus, "Other"))
gc_02_genus_unc_g2 <- gc_02_genus_unc_g1 %>% 
  mutate(Genus = if_else(Genus %ni% c("Aminicenantales",
                                      "JS1",
                                      "Aerophobales",
                                      "BHI80.139") , Genus, "Aminicenantales group"))
gc_02_genus_unc_g3 <- gc_02_genus_unc_g2 %>% 
  mutate(Genus = if_else(Genus %ni% c('Other'), Genus, "Uz_Other"))
#Filter unclassifieds out
gc_02_genus_unc_g3 <- gc_02_genus_unc_g3 %>%
  filter(Genus != "Unclassified")


bold_gc02 <- c('#11A579','#3969AC','#f97b72','#D55E00', "#6db6ff",
               '#80BA5A','#E68310','#CF1C90','#008695','#E73F74',
               "#000000", '#7F3C8D', "#006ddb","#924900", "#FFE800", 
               "#009E73", "#999999", "#E69F00", "#56B4E9")

#Filter unclassifieds out
gc_02_genus_unc_g7 <- gc_02_genus_unc_g7 %>%
  filter(Genus != "Unclassified")
filt_major_taxa_proportions_tab_for_plot_g_gc02.g2 %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
  ggplot(aes(x = Sample, y = Proportion, fill = Genus)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = stock_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% of 16S rRNA gene copies recovered", title = "GC02_CR02") +
  geom_bar(stat = "identity", width = 0.5) +
  scale_x_discrete(labels = c(3, 15, 18, 33, 36, 42, 48, 54, 60, 69, 78, 81, 87, 96, 102, 108, 111, 123)) 
getwd()
#For GC04
#Define color palette
gc_04_Palette <- c("#ff6db6", "#ffff6d", "#920000","#6db6ff",
                   "#006ddb", "#009292","#db6d00","#ffb6db", "#490092")
#Make plot
filt_major_taxa_proportions_tab_for_plot_gc04.g2 %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
  ggplot(aes(x = Sample, y = Proportion, fill = Phylum)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% of 16S rRNA gene copies recovered", title = "GC04_CR03") +
  scale_x_discrete(labels = c(65, 80, 90, 100, 110, 130, 140, 150, 160, 165, 170, 180, 185, 190, 195, 200)) 

#For genus
#Group all the unclassifieds and uncultured together
`%ni%` <- Negate(`%in%`)
gc_04_genus_unc <- filt_major_taxa_proportions_tab_for_plot_g_gc04.g2 %>% 
  mutate(Genus = if_else(Genus %ni% c("uncultured.1",
                                      "uncultured.2",
                                      "uncultured",
                                      "unclassified_Root",
                                      "unclassified_Bacteria") , Genus, "Unclassified"))

gc_04_genus_unc_g1 <- gc_04_genus_unc %>% 
  mutate(Genus = if_else(Genus %ni% c("TA06",
                                      "Azotobacter_5") , Genus, "Other"))

gc_04_genus_unc_g2 <- gc_04_genus_unc_g1 %>% 
  mutate(Genus = if_else(Genus %ni% c("Aminicenantales",
                                      "JS1",
                                      "Aerophobales",
                                      "BHI80.139") , Genus, "Aminicenantales group"))

gc_04_genus_unc_g3 <- gc_04_genus_unc_g2 %>% 
  mutate(Genus = if_else(Genus %ni% c("Other"), Genus, "Uz_Other"))

#Filter unclassifieds out
gc_04_genus_unc_g3 <- gc_04_genus_unc_g3 %>%
  filter(Genus != "Unclassified")


bold_gc04 <- c('#F2B701','#11A579','#3969AC',"#6db6ff",'#E68310',
               '#E73F74', "#000000", '#7F3C8D')


#Make plot
filt_major_taxa_proportions_tab_for_plot_g_gc04.g2 %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
  ggplot(aes(x = Sample, y = Proportion, fill = Genus)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = stock_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% of 16S rRNA gene copies recovered", title = "GC04_CR03") +
  scale_x_discrete(labels = c(65, 80, 90, 100, 110, 130, 140, 150, 160, 165, 170, 180, 185, 190, 195, 200)) 

#For GC06
#Define color palette
gc_06_Palette <- c("#000000", "#ff6db6", "#ffff6d", "#b66dff", "#004949", "#b6dbff", 
                   "#6db6ff", "#24ff24", "#006ddb", "#999999", "#009292","#db6d00",
                   "#ffb6db", "#490092", "#E69F00")

#Make plot
filt_major_taxa_proportions_tab_for_plot_gc06.g2 %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
  ggplot(aes(x = Sample, y = Proportion, fill = Phylum)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = gc_06_Palette) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% of 16S rRNA gene copies recovered", title = "GC06_CR04.5") +
  scale_x_discrete(labels = c(5, 10, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 
                              140, 150, 155, 160, 165, 170, 175, 180, 185)) 

#For genus
#Filter unclassifieds out
#Group all unclassifieds together
gc_06_genus_unc <- filt_major_taxa_proportions_tab_for_plot_g_gc06.g2 %>% 
  mutate(Genus = if_else(Genus %ni% c("uncultured.1",
                                      "uncultured.2",
                                      "uncultured.3",
                                      "uncultured.4",
                                      "unclassified_Root",
                                      "unclassified_Bacteria",
                                      "unclassified_Halobacterales",
                                      "unclassified_Burkholderiales",
                                      "unclassified_Anaerolineae",
                                      "unclassified_Pseudonocardiales",
                                      "unclassified_Xanthomonadaceae_1",
                                      'unclassified_Streptosporangiales',
                                      'Incertae.Sedis') , Genus, "Unclassified"))
gc_06_genus_unc_g1 <- gc_06_genus_unc %>% 
  mutate(Genus = if_else(Genus %ni% c("Pla1.lineage",
                                      "Chloroplast",
                                      'Pla3.lineage',
                                      "DscP2",
                                      "SG8.4",
                                      "NKB15",
                                      'Candidatus.Woesebacteria',
                                      'Candidatus.Chisholmbacteria',
                                      'Zixibacteria',
                                      'A0839',
                                      'ANME.2a.2b',
                                      'ANME.2b',
                                      'Berkelbacteria',
                                      'C0119',
                                      'X4.29',
                                      'Vermiphilaceae',
                                      'Candidatus.Chloroploca',
                                      'Chloroflexus',
                                      'Pseudoxanthomonas',
                                      'Phormidesmis.ANT.L52.6',
                                      'Paracoccus_1',
                                      'DG.56',
                                      'Geitlerinema.PCC.7105',
                                      'H3.93',
                                      'Jeotgalibacillus',
                                      'JG30.KF.CM45',
                                      'Leptococcus.JA.3.3Ab',
                                      'MBMPE27',
                                      'ODP1230B23.02',
                                      'Sumerlaea',
                                      'SEEP.SRB1',
                                      'SAR324.clade.Marine.group.B.'), Genus, "Other"))
gc_06_genus_unc_g2 <- gc_06_genus_unc_g1 %>% 
  mutate(Genus = if_else(Genus %ni% c("Aminicenantales",
                                      "JS1",
                                      "Aerophobales",
                                      "BHI80.139",
                                      "TA06"), Genus, "Aminicenantales+JS1"))
gc_06_genus_unc_g3 <- gc_06_genus_unc_g2 %>% 
  mutate(Genus = if_else(Genus %ni% c("SCGC.AB.539.J10",
                                      "unclassified_Dehalococcoidia",
                                      'MSB.5B2'), Genus, "Dehalococcoides"))

gc_06_genus_unc_g4 <- gc_06_genus_unc_g3 %>% 
  mutate(Genus = if_else(Genus %ni% c("uncultured"), Genus, "Pedomicrobium"))

gc_06_genus_unc_g5 <- gc_06_genus_unc_g4 %>% 
  mutate(Genus = if_else(Genus %ni% c('Novosphingobium_2',
                                      'Prauserella_1'), Genus, "Hydrocarbon degraders"))

gc_06_genus_unc_g6 <- gc_06_genus_unc_g5 %>% 
  mutate(Genus = if_else(Genus %ni% c('Bathyarchaeia',
                                      'Lokiarchaeia'), Genus, "Bathyarchaeia + Lokiarchaeia"))

gc_06_genus_unc_g7 <- gc_06_genus_unc_g6 %>% 
  mutate(Genus = if_else(Genus %ni% c('Candidatus.Nitrosopumilus',
                                      'Nitrospira'), Genus, "Nitrosopumilus + Nitrospira"))

gc_06_genus_unc_g8 <- gc_06_genus_unc_g7 %>% 
  mutate(Genus = if_else(Genus %ni% c('Subgroup.26',
                                      'Subgroup.22'), Genus, "Acidobacteria"))

gc_06_genus_unc_g9 <- gc_06_genus_unc_g8 %>% 
  mutate(Genus = if_else(Genus %ni% c('AKAU3564.sediment.group',
                                      'unclassified_Pirellulaceae'), Genus, "Unclassified Pirellulaceae"))

gc_06_genus_unc_g10 <- gc_06_genus_unc_g9 %>% 
  mutate(Genus = if_else(Genus %ni% c('unclassified_Actinobacteria',
                                      'unclassified_Actinobacteriota',
                                      'WCHB1.81',
                                      'WS1'), Genus, "Unclassified Actinobacteria"))

bold_gc06 <- c('#7F3C8D','#F2B701','#11A579','#3969AC','#f97b72',
               "#999999",'#D55E00', "#6db6ff", '#80BA5A','#E68310',
               '#CF1C90','#008695','#E73F74',"#000000", "#006ddb",
               '#7F3C8D', "#924900", "#FFE800", "#009E73","#56B4E9")

#Make plot
gc_06_genus_unc_g10 <- gc_06_genus_unc_g10 %>%
  filter(Genus != "Unclassified")
filt_major_taxa_proportions_tab_for_plot_g_gc06.g2 %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
  ggplot(aes(x = Sample, y = Proportion, fill = Genus)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = stock_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% of 16S rRNA gene copies recovered", title = "GC06_CR04.5") +
  scale_x_discrete(labels = c(5, 10, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 
                              140, 150, 155, 160, 165, 170, 175, 180, 185)) 

# Another way to visualize taxonomy could be using boxplots where each box is a major taxon, with each point being colored/shaped based on core:
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(Phylum, Proportion)) +
  geom_jitter(aes(color = Core, shape = Core), size = 2, width = 0.15, height = 0) +
  scale_color_manual(values = unique(filt_major_taxa_proportions_tab_for_plot.g2$core_colors[order(filt_major_taxa_proportions_tab_for_plot.g2$Core)])) +
  geom_boxplot(fill = NA, outlier.color = NA) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x="Phylum", y = "% of 16S rRNA gene copies recovered", title = "All samples")

# Or making the point color/shape based on our Depth bins:
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(Phylum, Proportion)) +
  geom_jitter(aes(color = Depth_bin, shape = Depth_bin), size = 2, width = 0.15, height = 0) +
  scale_color_manual(values = unique(filt_major_taxa_proportions_tab_for_plot.g2$depth_colors[order(filt_major_taxa_proportions_tab_for_plot.g2$Depth_bin)])) +
  geom_boxplot(fill = NA, outlier.color = NA) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x="Phylum", y = "% of 16S rRNA gene copies recovered", title = "All samples")

# And the last bit of code we’ll include here for example purposes is making a pie chart for each sample:
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(x = " ", y = Proportion, group = Phylum, fill = Phylum)) +
  geom_bar(width = 1, stat = "identity") + scale_fill_brewer(palette = "Paired") +
  coord_polar("y", start = 0) + facet_wrap(~ Sample) + theme_void()

#Make a summary of all major taxa proportions by Family

# Betadisperser and permutational ANOVA
# One way we can integrate this info is with a permutational ANOVA test to see 
# if any of the available information correlates with overall community structure (based on our ASVs, not taxonomy).
# Here we are going to test if there is a statistically significant difference between our dates or our depth bins.
#Calculate significance difference between intragroup variants
anova(betadisper(euc_dist, sample_info_tab$core))
#Checking by Date, we get a significant result (2.2 e-16) from the betadisper test. 
#This tells us that there is a difference between group dispersions, which means 
#that we can’t trust the results of an adonis (permutational anova) test on this, 
#because the assumption of homogenous within-group dispersions is not met.

#Let’s check based on Depth bin:
anova(betadisper(euc_dist, sample_info_tab$Depth_bin))
# Going with a standard 0.05 cutoff, this is significant too (1.34 e-6), and tells us that 
# within-group dispersion of our Depth bins isn’t low enough for us to trust the 
# results of running the adonis permutational anova.

#Differential abundance analysis
#We are going to take advantage of another phyloseq convenience, and use the 
#phyloseq_to_deseq2() function to make our DESeq2 object:
ASV_deseq <- phyloseq_to_deseq2(ASV_physeq, ~ Core + Depth_bin)

# and running deseq standard analysis:
ASV_deseq <- DESeq(ASV_deseq)

#We can now access the results.
#First, we’ll look at Dates.We can get a look at the names DESeq has for each contrast with the following:
resultsNames(ASV_deseq)

# And we can pull out the results for one of those contrasts like so:
deseq_res_Core_GC06_vs_GC02 <- results(ASV_deseq, alpha = 0.01, name = "Core_GC06_vs_GC02")

#And we can use the summary() function to get an overview of that results object:
summary(deseq_res_Core_GC06_vs_GC02)

#Let’s subset our results table to only include those that pass our 0.01 adjusted p-value cutoff:
sigtab_deseq_res_Core_GC06_vs_GC02 <- deseq_res_Core_GC06_vs_GC02[which(deseq_res_Core_GC06_vs_GC02$padj < 0.01), ]

# Now we can see this table only contains those we consider significantly differentially abundant:
summary(sigtab_deseq_res_Core_GC06_vs_GC02) 

# Next let’s stitch that table together with these ASVs’ taxonomic classifications for a quick look at both together:
sigtab_deseq_res_Core_GC06_vs_GC02_with_tax <- cbind(as(sigtab_deseq_res_Core_GC06_vs_GC02, "data.frame"), as(tax_table(ASV_physeq)[row.names(sigtab_deseq_res_Core_GC06_vs_GC02), ], "matrix"))

# and now let's sort that table by the baseMean column
sigtab_deseq_res_Core_GC06_vs_GC02_with_tax <- sigtab_deseq_res_Core_GC06_vs_GC02_with_tax[order(sigtab_deseq_res_Core_GC06_vs_GC02$baseMean, decreasing=T), ]

head(sigtab_deseq_res_Core_GC06_vs_GC02_with_tax)

# Create relative abundance dataframes for genus for each core
relative_tab_genus_mc01 <- apply(counts_tab_mc01, 2, function(x) 100*(x/sum(x)))
relative_tab_genus_gc02 <- apply(counts_tab_gc02, 2, function(x) 100*(x/sum(x)))
relative_tab_genus_gc04 <- apply(counts_tab_gc04, 2, function(x) 100*(x/sum(x)))
relative_tab_genus_gc06 <- apply(counts_tab_gc06, 2, function(x) 100*(x/sum(x)))

#### Relative abundance stackbar
Make a stackbar, as one does:
  ```{r eval=F}
relative_long <- melt(relative_tab, varnames = c("ASV","sample"), value.name = "relabund")
relative_long <- cbind(relative_long, tax_tab[as.character(relative_long$ASV), 
                                              grep("ASV", colnames(tax_tab),invert=T, value=T)])
#For genus for each core
relative_long_genus_mc01 <- melt(relative_tab_genus_mc01, varnames = c("ASV","sample"), value.name = "relabund")
relative_long_genus_mc01 <- cbind(relative_long_genus_mc01, tax_tab[as.character(relative_long_genus_mc01$ASV), 
                                                                    grep("ASV", colnames(tax_tab),invert=T, value=T)])
relative_long_genus_gc02 <- melt(relative_tab_genus_gc02, varnames = c("ASV","sample"), value.name = "relabund")
relative_long_genus_gc02 <- cbind(relative_long_genus_gc02, tax_tab[as.character(relative_long_genus_gc02$ASV), 
                                                                    grep("ASV", colnames(tax_tab),invert=T, value=T)])
relative_long_genus_gc04 <- melt(relative_tab_genus_gc04, varnames = c("ASV","sample"), value.name = "relabund")
relative_long_genus_gc04 <- cbind(relative_long_genus_gc04, tax_tab[as.character(relative_long_genus_gc04$ASV), 
                                                                    grep("ASV", colnames(tax_tab),invert=T, value=T)])
relative_long_genus_gc06 <- melt(relative_tab_genus_gc06, varnames = c("ASV","sample"), value.name = "relabund")
relative_long_genus_gc06 <- cbind(relative_long_genus_gc06, tax_tab[as.character(relative_long_genus_gc06$ASV), 
                                                                    grep("ASV", colnames(tax_tab),invert=T, value=T)])


```
Can subset (or not)
```{r eval=F}
sub_relative_long <- relative_long[grep("", relative_long$phylum),] %>% group_by(sample, phylum, class) %>% summarise(relabund = sum(value))
#For genus for each core
sub_relative_long_genus_mc01 <- relative_long_genus_mc01[grep("", relative_long_genus_mc01$genus),] %>% group_by(sample, phylum, class) %>% summarise(relabund = sum(value))
sub_relative_long_genus_gc02 <- relative_long_genus_gc02[grep("", relative_long_genus_gc02$phylum),] %>% group_by(sample, phylum, class) %>% summarise(relabund = sum(value))
sub_relative_long_genus_gc04 <- relative_long_genus_gc04[grep("", relative_long_genus_gc04$phylum),] %>% group_by(sample, phylum, class) %>% summarise(relabund = sum(value))
sub_relative_long_genus_gc06 <- relative_long_genus_gc06[grep("", relative_long_genus_gc06$phylum),] %>% group_by(sample, phylum, class) %>% summarise(relabund = sum(value))


#Make heatmaps

#Plot for MC01
#Filters
`%ni%` <- Negate(`%in%`)
mc_01_genus_unc <- filt_major_taxa_proportions_tab_for_plot_g_mc01.g2 %>% 
  mutate(Genus = if_else(Genus %ni% c("uncultured.1",
                                      "uncultured.2",
                                      "uncultured",
                                      "unclassified_Root",
                                      "unclassified_Bacilli") , Genus, "Unclassified"))

mc_01_genus_unc_g1 <- mc_01_genus_unc %>% 
  mutate(Genus = if_else(Genus %ni% c('Other'), Genus, "Uz_Other"))

#Organize things in plot
mc_01_genus_unc_g1$Sample <-gsub("^.*_", "", mc_01_genus_unc_g1$Sample)
mc_01_genus_unc_g1$Sample <- factor(mc_01_genus_unc_g1$Sample, levels=unique(mc_01_genus_unc_g1$Sample)[order(as.numeric(gsub("^.*_","", unique(mc_01_genus_unc_g1$Sample))))])
mc_01_genus_unc_g1$Genus <- with(mc_01_genus_unc_g1,factor(Genus,levels = rev(sort(unique(Genus)))))

#Make plot
ggplot(mc_01_genus_unc_g1, aes(x=Sample, y = Genus, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3","#EA5835"),
    values=rescale(c(0,5,95)),
    limits=c(0,95)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/MC01_heatmap.pdf", width=12, height=5, dpi=300)

#Phylum
#Filters
`%ni%` <- Negate(`%in%`)
mc01_phylum_unc <- filt_major_taxa_proportions_tab_for_plot_mc01.g2 %>% 
  mutate(Phylum = if_else(Phylum %ni% c('Other'), Phylum, "z_Other"))

#Organize things in plot
mc01_phylum_unc$Sample <-gsub("^.*_", "", mc01_phylum_unc$Sample)
mc01_phylum_unc$Sample <- factor(mc01_phylum_unc$Sample, levels=unique(mc01_phylum_unc$Sample)[order(as.numeric(gsub("^.*_","", unique(mc01_phylum_unc$Sample))))])
mc01_phylum_unc$Phylum <- with(mc01_phylum_unc,factor(Phylum,levels = rev(sort(unique(Phylum)))))

#Make plot
ggplot(mc01_phylum_unc, aes(x=Sample, y = Phylum, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3","#EA5835"),
    values=rescale(c(0,5,95)),
    limits=c(0,95)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/MC01_heatmap_phylum.pdf", width=12, height=5, dpi=300)

#Plot for gc02

#Genus
#Filters
`%ni%` <- Negate(`%in%`)
gc_02_genus_unc <- filt_major_taxa_proportions_tab_for_plot_g_gc02.g2 %>% 
  mutate(Genus = if_else(Genus %ni% c("uncultured.1",
                                      "uncultured.2",
                                      "uncultured.3",
                                      "uncultured",
                                      "unclassified_Root",
                                      "unclassified_Bacteria",
                                      "unclassified_Halobacterales",
                                      "Unclassified") , Genus, "Unclassified"))
gc_02_genus_unc_g1 <- gc_02_genus_unc %>% 
  mutate(Genus = if_else(Genus %ni% c("DscP2",
                                      "SG8.4",
                                      "WCHB1.81") , Genus, "Other"))
gc_02_genus_unc_g2 <- gc_02_genus_unc_g1 %>% 
  mutate(Genus = if_else(Genus %ni% c("Aminicenantales",
                                      "JS1",
                                      "Aerophobales",
                                      "BHI80.139") , Genus, "Aminicenantales group"))
gc_02_genus_unc_g3 <- gc_02_genus_unc_g2 %>% 
  mutate(Genus = if_else(Genus %ni% c('Other'), Genus, "Uz_Other"))
#Filter unclassifieds out
gc_02_genus_unc_g3 <- gc_02_genus_unc_g3 %>%
  filter(Genus != "Unclassified")

#Organize things in plot
gc_02_genus_unc_g3$Sample <-gsub("^.*_", "", gc_02_genus_unc_g3$Sample)
gc_02_genus_unc_g3$Sample <- factor(gc_02_genus_unc_g3$Sample, levels=unique(gc_02_genus_unc_g3$Sample)[order(as.numeric(gsub("^.*_","", unique(gc_02_genus_unc_g3$Sample))))])
gc_02_genus_unc_g3$Genus <- with(gc_02_genus_unc_g3,factor(Genus,levels = rev(sort(unique(Genus)))))

filt_major_taxa_proportions_tab_for_plot_g_gc02.g2$Sample <-gsub("^.*_", "", filt_major_taxa_proportions_tab_for_plot_g_gc02.g2$Sample)
filt_major_taxa_proportions_tab_for_plot_g_gc02.g2$Sample <- factor(filt_major_taxa_proportions_tab_for_plot_g_gc02.g2$Sample, levels=unique(filt_major_taxa_proportions_tab_for_plot_g_gc02.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(filt_major_taxa_proportions_tab_for_plot_g_gc02.g2$Sample))))])
filt_major_taxa_proportions_tab_for_plot_g_gc02.g2$Genus <- with(filt_major_taxa_proportions_tab_for_plot_g_gc02.g2,factor(Genus,levels = rev(sort(unique(Genus)))))

#Make plot
ggplot(filt_major_taxa_proportions_tab_for_plot_g_gc02.g2, aes(x=Sample, y = Genus, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3","#EA5835"),
    values=rescale(c(0,5,95)),
    limits=c(0,95)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/GC02_heatmap.pdf", width=12, height=5, dpi=300)

#Phylum
#Filters
`%ni%` <- Negate(`%in%`)
gc_02_phylum_unc <- filt_major_taxa_proportions_tab_for_plot_gc02.g2 %>% 
  mutate(Phylum = if_else(Phylum %ni% c('Other'), Phylum, "z_Other"))

#Organize things in plot
gc_02_phylum_unc$Sample <-gsub("^.*_", "", gc_02_phylum_unc$Sample)
gc_02_phylum_unc$Sample <- factor(gc_02_phylum_unc$Sample, levels=unique(gc_02_phylum_unc$Sample)[order(as.numeric(gsub("^.*_","", unique(gc_02_phylum_unc$Sample))))])
gc_02_phylum_unc$Phylum <- with(gc_02_phylum_unc,factor(Phylum,levels = rev(sort(unique(Phylum)))))

#Make plot
ggplot(gc_02_phylum_unc, aes(x=Sample, y = Phylum, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3","#EA5835"),
    values=rescale(c(0,5,95)),
    limits=c(0,95)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/GC02_heatmap_phylum.pdf", width=12, height=5, dpi=300)

#Plot for gc04
#Genus
#Filters
`%ni%` <- Negate(`%in%`)
gc_04_genus_unc <- filt_major_taxa_proportions_tab_for_plot_g_gc04.g2 %>% 
  mutate(Genus = if_else(Genus %ni% c("uncultured.1",
                                      "uncultured.2",
                                      "uncultured",
                                      "unclassified_Root",
                                      "unclassified_Bacteria") , Genus, "Unclassified"))

gc_04_genus_unc_g1 <- gc_04_genus_unc %>% 
  mutate(Genus = if_else(Genus %ni% c("TA06",
                                      "Azotobacter_5") , Genus, "Other"))

gc_04_genus_unc_g2 <- gc_04_genus_unc_g1 %>% 
  mutate(Genus = if_else(Genus %ni% c("Aminicenantales",
                                      "JS1",
                                      "Aerophobales",
                                      "BHI80.139") , Genus, "Aminicenantales group"))

gc_04_genus_unc_g3 <- gc_04_genus_unc_g2 %>% 
  mutate(Genus = if_else(Genus %ni% c("Other"), Genus, "Uz_Other"))

#Filter unclassifieds out
gc_04_genus_unc_g3 <- gc_04_genus_unc_g3 %>%
  filter(Genus != "Unclassified")

#Organize things in plot
gc_04_genus_unc_g3$Sample <-gsub("^.*_", "", gc_04_genus_unc_g3$Sample)
gc_04_genus_unc_g3$Sample <- factor(gc_04_genus_unc_g3$Sample, levels=unique(gc_04_genus_unc_g3$Sample)[order(as.numeric(gsub("^.*_","", unique(gc_04_genus_unc_g3$Sample))))])
gc_04_genus_unc_g3$Genus <- with(gc_04_genus_unc_g3,factor(Genus,levels = rev(sort(unique(Genus)))))

#Make plot
ggplot(gc_04_genus_unc_g3, aes(x=Sample, y = Genus, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.5) +
  scale_fill_gradientn(
    colors=c("white","#2699D3","#EA5835"),
    values=rescale(c(0,5,95)),
    limits=c(0,95)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/GC04_heatmap.pdf", width=12, height=5, dpi=300)

#Phylum
#Filters
`%ni%` <- Negate(`%in%`)
gc_04_phylum_unc <- filt_major_taxa_proportions_tab_for_plot_gc04.g2 %>% 
  mutate(Phylum = if_else(Phylum %ni% c('Other'), Phylum, "z_Other"))

#Organize things in plot
gc_04_phylum_unc$Sample <-gsub("^.*_", "", gc_04_phylum_unc$Sample)
gc_04_phylum_unc$Sample <- factor(gc_04_phylum_unc$Sample, levels=unique(gc_04_phylum_unc$Sample)[order(as.numeric(gsub("^.*_","", unique(gc_04_phylum_unc$Sample))))])
gc_04_phylum_unc$Phylum <- with(gc_04_phylum_unc,factor(Phylum,levels = rev(sort(unique(Phylum)))))

#Make plot
ggplot(gc_04_phylum_unc, aes(x=Sample, y = Phylum, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.5) +
  scale_fill_gradientn(
    colors=c("white","#2699D3","#EA5835"),
    values=rescale(c(0,5,95)),
    limits=c(0,95)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/GC04_heatmap_phylum.pdf", width=12, height=5, dpi=300)


#Plot for gc06
#Genus
#Filters
#Group all unclassifieds together
gc_06_genus_unc <- filt_major_taxa_proportions_tab_for_plot_g_gc06.g2 %>% 
  mutate(Genus = if_else(Genus %ni% c("uncultured.1",
                                      "uncultured.2",
                                      "uncultured.3",
                                      "uncultured.4",
                                      "uncultured",
                                      "unclassified_Bacteria",
                                      "unclassified_Root") , Genus, "Unclassified"))

gc_06_genus_unc_g1 <- gc_06_genus_unc %>% 
  mutate(Genus = if_else(Genus %ni% c('Other'), Genus, "Uz_Other"))

#Filter unclassifieds out
gc_06_genus_unc_g1 <- gc_06_genus_unc_g1 %>%
  filter(Genus != "Unclassified")

#Organize things in plot
gc_06_genus_unc_g1$Sample <-gsub("^.*_", "", gc_06_genus_unc_g1$Sample)
gc_06_genus_unc_g1$Sample <- factor(gc_06_genus_unc_g1$Sample, levels=unique(gc_06_genus_unc_g1$Sample)[order(as.numeric(gsub("^.*_","", unique(gc_06_genus_unc_g1$Sample))))])
gc_06_genus_unc_g1$Genus <- with(gc_06_genus_unc_g1,factor(Genus,levels = rev(sort(unique(Genus)))))

#Make plot
ggplot(gc_06_genus_unc_g1, aes(x=Sample, y = Genus, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3","#EA5835"),
    values=rescale(c(0,5,95)),
    limits=c(0,95)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/GC06_heatmap.pdf", width=12, height=5, dpi=300)

#Phylum
#Filters
`%ni%` <- Negate(`%in%`)
gc_06_phylum_unc <- filt_major_taxa_proportions_tab_for_plot_gc06.g2 %>% 
  mutate(Phylum = if_else(Phylum %ni% c('Other'), Phylum, "z_Other"))

#Organize things in plot
gc_06_phylum_unc$Sample <-gsub("^.*_", "", gc_06_phylum_unc$Sample)
gc_06_phylum_unc$Sample <- factor(gc_06_phylum_unc$Sample, levels=unique(gc_06_phylum_unc$Sample)[order(as.numeric(gsub("^.*_","", unique(gc_06_phylum_unc$Sample))))])
gc_06_phylum_unc$Phylum <- with(gc_06_phylum_unc,factor(Phylum,levels = rev(sort(unique(Phylum)))))

#Make plot
ggplot(gc_06_phylum_unc, aes(x=Sample, y = Phylum, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.5) +
  scale_fill_gradientn(
    colors=c("white","#2699D3","#EA5835"),
    values=rescale(c(0,5,95)),
    limits=c(0,95)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/GC06_heatmap_phylum.pdf", width=12, height=5, dpi=300)

#Making heatplots without filtering phyla with less than 5% representation out
#Create a copy of the tables for plots
major_taxa_proportions_tab_mc01_for_plot <- data.frame(major_taxa_proportions_tab_mc01)
major_taxa_proportions_tab_gc02_for_plot <- data.frame(major_taxa_proportions_tab_gc02)
major_taxa_proportions_tab_gc04_for_plot <- data.frame(major_taxa_proportions_tab_gc04)
major_taxa_proportions_tab_gc06_for_plot <- data.frame(major_taxa_proportions_tab_gc06)

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
#For phylum
major_taxa_proportions_tab_mc01_for_plot <- 
  major_taxa_proportions_tab_mc01_for_plot %>% rownames_to_column("Phylum")
major_taxa_proportions_tab_gc02_for_plot <- 
  major_taxa_proportions_tab_gc02_for_plot %>% rownames_to_column("Phylum")
major_taxa_proportions_tab_gc04_for_plot <- 
  major_taxa_proportions_tab_gc04_for_plot %>% rownames_to_column("Phylum")
major_taxa_proportions_tab_gc06_for_plot <- 
  major_taxa_proportions_tab_gc06_for_plot %>% rownames_to_column("Phylum")

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
#For phylum
major_taxa_proportions_tab_mc01_for_plot.g <- 
  major_taxa_proportions_tab_mc01_for_plot %>% 
  pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
major_taxa_proportions_tab_gc02_for_plot.g <- 
  major_taxa_proportions_tab_gc02_for_plot %>% 
  pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
major_taxa_proportions_tab_gc04_for_plot.g <- 
  major_taxa_proportions_tab_gc04_for_plot %>% 
  pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
major_taxa_proportions_tab_gc06_for_plot.g <- 
  major_taxa_proportions_tab_gc06_for_plot %>% 
  pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)

#Add depth and core to tables
major_taxa_proportions_tab_mc01_for_plot.g2 <- 
  major_taxa_proportions_tab_mc01_for_plot.g %>% left_join(mod_sample_info_tab_mc01)
major_taxa_proportions_tab_gc02_for_plot.g2 <- 
  major_taxa_proportions_tab_gc02_for_plot.g %>% left_join(mod_sample_info_tab_gc02)
major_taxa_proportions_tab_gc04_for_plot.g2 <- 
  major_taxa_proportions_tab_gc04_for_plot.g %>% left_join(mod_sample_info_tab_gc04)
major_taxa_proportions_tab_gc06_for_plot.g2 <- 
  major_taxa_proportions_tab_gc06_for_plot.g %>% left_join(mod_sample_info_tab_gc06)

#MC01
#Combine Desulfobacterota
major_taxa_proportions_tab_mc01_for_plot.g2 <- major_taxa_proportions_tab_mc01_for_plot.g2 %>% 
  mutate(Phylum = if_else(Phylum  %ni% c('Desulfobacterota_1','Desulfobacterota_2'), Phylum , "Desulfobacterota"))
#Organize things in plot
major_taxa_proportions_tab_mc01_for_plot.g2$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_mc01_for_plot.g2$Sample)
major_taxa_proportions_tab_mc01_for_plot.g2$Sample <- factor(major_taxa_proportions_tab_mc01_for_plot.g2$Sample, levels=unique(major_taxa_proportions_tab_mc01_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_mc01_for_plot.g2$Sample))))])
major_taxa_proportions_tab_mc01_for_plot.g2$Phylum <- with(major_taxa_proportions_tab_mc01_for_plot.g2,factor(Phylum,levels = rev(sort(unique(Phylum)))))

#Make plot
ggplot(major_taxa_proportions_tab_mc01_for_plot.g2, aes(x=Sample, y = Phylum, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3","#EA5835"),
    values=rescale(c(0,5,100)),
    limits=c(0,100)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/MC01_heatmap_phylum_all.pdf", width=12, height=5, dpi=300)

#Filter out

#GC02
#Combine Desulfobacterota
major_taxa_proportions_tab_gc02_for_plot.g2 <- major_taxa_proportions_tab_gc02_for_plot.g2 %>% 
  mutate(Phylum = if_else(Phylum  %ni% c('Desulfobacterota_1','Desulfobacterota_2'), Phylum , "Desulfobacterota"))
#Combine Bdellovibrionota
major_taxa_proportions_tab_gc02_for_plot.g2 <- major_taxa_proportions_tab_gc02_for_plot.g2 %>% 
  mutate(Phylum = if_else(Phylum  %ni% c('Bdellovibrionota_1','Bdellovibrionota_2'), Phylum , "Bdellovibrionota"))
#Organize things in plot
major_taxa_proportions_tab_gc02_for_plot.g2$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_gc02_for_plot.g2$Sample)
major_taxa_proportions_tab_gc02_for_plot.g2$Sample <- factor(major_taxa_proportions_tab_gc02_for_plot.g2$Sample, levels=unique(major_taxa_proportions_tab_gc02_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_gc02_for_plot.g2$Sample))))])
major_taxa_proportions_tab_gc02_for_plot.g2$Phylum <- with(major_taxa_proportions_tab_gc02_for_plot.g2,factor(Phylum,levels = rev(sort(unique(Phylum)))))

#Make plot
ggplot(major_taxa_proportions_tab_gc02_for_plot.g2, aes(x=Sample, y = Phylum, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3","#EA5835"),
    values=rescale(c(0,5,100)),
    limits=c(0,100)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/GC02_heatmap_phylum_all.pdf", width=12, height=5, dpi=300)

#GC04
#Combine Desulfobacterota
major_taxa_proportions_tab_gc04_for_plot.g2 <- major_taxa_proportions_tab_gc04_for_plot.g2 %>% 
  mutate(Phylum = if_else(Phylum  %ni% c('Desulfobacterota_1','Desulfobacterota_2'), Phylum , "Desulfobacterota"))
#Combine Bdellovibrionota
major_taxa_proportions_tab_gc04_for_plot.g2 <- major_taxa_proportions_tab_gc04_for_plot.g2 %>% 
  mutate(Phylum = if_else(Phylum  %ni% c('Bdellovibrionota_1','Bdellovibrionota_2'), Phylum , "Bdellovibrionota"))
#Organize things in plot
major_taxa_proportions_tab_gc04_for_plot.g2$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_gc04_for_plot.g2$Sample)
major_taxa_proportions_tab_gc04_for_plot.g2$Sample <- factor(major_taxa_proportions_tab_gc04_for_plot.g2$Sample, levels=unique(major_taxa_proportions_tab_gc04_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_gc04_for_plot.g2$Sample))))])
major_taxa_proportions_tab_gc04_for_plot.g2$Phylum <- with(major_taxa_proportions_tab_gc04_for_plot.g2,factor(Phylum,levels = rev(sort(unique(Phylum)))))
#Make plot
ggplot(major_taxa_proportions_tab_gc04_for_plot.g2, aes(x=Sample, y = Phylum, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3","#EA5835"),
    values=rescale(c(0,5,100)),
    limits=c(0,100)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/GC04_heatmap_phylum_all.pdf", width=12, height=5, dpi=300)

#GC06
#Combine Desulfobacterota
major_taxa_proportions_tab_gc06_for_plot.g2 <- major_taxa_proportions_tab_gc06_for_plot.g2 %>% 
  mutate(Phylum = if_else(Phylum  %ni% c('Desulfobacterota_1','Desulfobacterota_2'), Phylum , "Desulfobacterota"))
#Combine Bdellovibrionota
major_taxa_proportions_tab_gc06_for_plot.g2 <- major_taxa_proportions_tab_gc06_for_plot.g2 %>% 
  mutate(Phylum = if_else(Phylum  %ni% c('Bdellovibrionota_1','Bdellovibrionota_2'), Phylum , "Bdellovibrionota"))
#Organize things in plot
major_taxa_proportions_tab_gc06_for_plot.g2$Sample <-gsub("^.*_", "", major_taxa_proportions_tab_gc06_for_plot.g2$Sample)
major_taxa_proportions_tab_gc06_for_plot.g2$Sample <- factor(major_taxa_proportions_tab_gc06_for_plot.g2$Sample, levels=unique(major_taxa_proportions_tab_gc06_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(major_taxa_proportions_tab_gc06_for_plot.g2$Sample))))])
major_taxa_proportions_tab_gc06_for_plot.g2$Phylum <- with(major_taxa_proportions_tab_gc06_for_plot.g2,factor(Phylum,levels = rev(sort(unique(Phylum)))))
#Make plot
ggplot(major_taxa_proportions_tab_gc06_for_plot.g2, aes(x=Sample, y = Phylum, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3","#EA5835"),
    values=rescale(c(0,5,100)),
    limits=c(0,100)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/GC06_heatmap_phylum_all.pdf", width=12, height=5, dpi=300)

#Making heatplots for Others

#Create a copy of the tables for plots
other_filt_major_taxa_proportions_tab_mc01_for_plot <- other_filt_major_taxa_proportions_tab_mc01
other_filt_major_taxa_proportions_tab_gc02_for_plot <- other_filt_major_taxa_proportions_tab_gc02
other_filt_major_taxa_proportions_tab_gc04_for_plot <- other_filt_major_taxa_proportions_tab_gc04
other_filt_major_taxa_proportions_tab_gc06_for_plot <- other_filt_major_taxa_proportions_tab_gc06

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
#For phylum
other_filt_major_taxa_proportions_tab_mc01_for_plot <- 
  other_filt_major_taxa_proportions_tab_mc01_for_plot %>% rownames_to_column("Phylum")
other_filt_major_taxa_proportions_tab_gc02_for_plot <- 
  other_filt_major_taxa_proportions_tab_gc02_for_plot %>% rownames_to_column("Phylum")
other_filt_major_taxa_proportions_tab_gc04_for_plot <- 
  other_filt_major_taxa_proportions_tab_gc04_for_plot %>% rownames_to_column("Phylum")
other_filt_major_taxa_proportions_tab_gc06_for_plot <- 
  other_filt_major_taxa_proportions_tab_gc06_for_plot %>% rownames_to_column("Phylum")

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
#For phylum
other_filt_major_taxa_proportions_tab_mc01_for_plot.g <- 
  other_filt_major_taxa_proportions_tab_mc01_for_plot %>% 
  pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
other_filt_major_taxa_proportions_tab_gc02_for_plot.g <- 
  other_filt_major_taxa_proportions_tab_gc02_for_plot %>% 
  pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
other_filt_major_taxa_proportions_tab_gc04_for_plot.g <- 
  other_filt_major_taxa_proportions_tab_gc04_for_plot %>% 
  pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
other_filt_major_taxa_proportions_tab_gc06_for_plot.g <- 
  other_filt_major_taxa_proportions_tab_gc06_for_plot %>% 
  pivot_longer(!Phylum, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)

#Add depth and core to tables
other_filt_major_taxa_proportions_tab_mc01_for_plot.g2 <- 
  other_filt_major_taxa_proportions_tab_mc01_for_plot.g %>% left_join(mod_sample_info_tab_mc01)
other_filt_major_taxa_proportions_tab_gc02_for_plot.g2 <- 
  other_filt_major_taxa_proportions_tab_gc02_for_plot.g %>% left_join(mod_sample_info_tab_gc02)
other_filt_major_taxa_proportions_tab_gc04_for_plot.g2 <- 
  other_filt_major_taxa_proportions_tab_gc04_for_plot.g %>% left_join(mod_sample_info_tab_gc04)
other_filt_major_taxa_proportions_tab_gc06_for_plot.g2 <- 
  other_filt_major_taxa_proportions_tab_gc06_for_plot.g %>% left_join(mod_sample_info_tab_gc06)

#MC01
#Organize things in plot
other_filt_major_taxa_proportions_tab_mc01_for_plot.g2$Sample <-gsub("^.*_", "", other_filt_major_taxa_proportions_tab_mc01_for_plot.g2$Sample)
other_filt_major_taxa_proportions_tab_mc01_for_plot.g2$Sample <- factor(other_filt_major_taxa_proportions_tab_mc01_for_plot.g2$Sample, levels=unique(other_filt_major_taxa_proportions_tab_mc01_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(other_filt_major_taxa_proportions_tab_mc01_for_plot.g2$Sample))))])
other_filt_major_taxa_proportions_tab_mc01_for_plot.g2$Phylum <- with(other_filt_major_taxa_proportions_tab_mc01_for_plot.g2,factor(Phylum,levels = rev(sort(unique(Phylum)))))
#Make plot
ggplot(other_filt_major_taxa_proportions_tab_mc01_for_plot.g2, aes(x=Sample, y = Phylum, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3"),
    values=rescale(c(0,5)),
    limits=c(0,5)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/others_MC01_heatmap_phylum.pdf", width=12, height=5, dpi=300)

#GC02
#Organize things in plot
other_filt_major_taxa_proportions_tab_gc02_for_plot.g2$Sample <-gsub("^.*_", "", other_filt_major_taxa_proportions_tab_gc02_for_plot.g2$Sample)
other_filt_major_taxa_proportions_tab_gc02_for_plot.g2$Sample <- factor(other_filt_major_taxa_proportions_tab_gc02_for_plot.g2$Sample, levels=unique(other_filt_major_taxa_proportions_tab_gc02_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(other_filt_major_taxa_proportions_tab_gc02_for_plot.g2$Sample))))])
other_filt_major_taxa_proportions_tab_gc02_for_plot.g2$Phylum <- with(other_filt_major_taxa_proportions_tab_gc02_for_plot.g2,factor(Phylum,levels = rev(sort(unique(Phylum)))))

#Make plot
ggplot(other_filt_major_taxa_proportions_tab_gc02_for_plot.g2, aes(x=Sample, y = Phylum, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3"),
    values=rescale(c(0,5)),
    limits=c(0,5)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/others_GC02_heatmap_phylum.pdf", width=12, height=5, dpi=300)

#GC04
#Organize things in plot
other_filt_major_taxa_proportions_tab_gc04_for_plot.g2$Sample <-gsub("^.*_", "", other_filt_major_taxa_proportions_tab_gc04_for_plot.g2$Sample)
other_filt_major_taxa_proportions_tab_gc04_for_plot.g2$Sample <- factor(other_filt_major_taxa_proportions_tab_gc04_for_plot.g2$Sample, levels=unique(other_filt_major_taxa_proportions_tab_gc04_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(other_filt_major_taxa_proportions_tab_gc04_for_plot.g2$Sample))))])
other_filt_major_taxa_proportions_tab_gc04_for_plot.g2$Phylum <- with(other_filt_major_taxa_proportions_tab_gc04_for_plot.g2,factor(Phylum,levels = rev(sort(unique(Phylum)))))
#Make plot
ggplot(other_filt_major_taxa_proportions_tab_gc04_for_plot.g2, aes(x=Sample, y = Phylum, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3"),
    values=rescale(c(0,5)),
    limits=c(0,5)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/others_GC04_heatmap_phylum.pdf", width=12, height=5, dpi=300)

#GC06
#Organize things in plot
other_filt_major_taxa_proportions_tab_gc06_for_plot.g2$Sample <-gsub("^.*_", "", other_filt_major_taxa_proportions_tab_gc06_for_plot.g2$Sample)
other_filt_major_taxa_proportions_tab_gc06_for_plot.g2$Sample <- factor(other_filt_major_taxa_proportions_tab_gc06_for_plot.g2$Sample, levels=unique(other_filt_major_taxa_proportions_tab_gc06_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(other_filt_major_taxa_proportions_tab_gc06_for_plot.g2$Sample))))])
other_filt_major_taxa_proportions_tab_gc06_for_plot.g2$Phylum <- with(other_filt_major_taxa_proportions_tab_gc06_for_plot.g2,factor(Phylum,levels = rev(sort(unique(Phylum)))))
#Make plot
ggplot(other_filt_major_taxa_proportions_tab_gc06_for_plot.g2, aes(x=Sample, y = Phylum, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3"),
    values=rescale(c(0,5)),
    limits=c(0,5)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/others_GC06_heatmap_phylum.pdf", width=12, height=5, dpi=300)

#Making pie charts for Proteobacteria and sulfur transforming microbes


#Create a copy of the tables for plots
proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot <- data.frame(proteobacteria_major_taxa_proportions_tab_g_mc01)
proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot <- data.frame(proteobacteria_major_taxa_proportions_tab_g_gc02)
proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot <- data.frame(proteobacteria_major_taxa_proportions_tab_g_gc04)
proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot <- data.frame(proteobacteria_major_taxa_proportions_tab_g_gc06)
head(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot)

#Simplify sample names
colnames(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot) <- gsub("^[^_]*_", "", colnames(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot))
colnames(proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot) <- gsub("^[^_]*_", "", colnames(proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot))
colnames(proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot) <- gsub("^[^_]*_", "", colnames(proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot))
colnames(proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot) <- gsub("^[^_]*_", "", colnames(proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot))

#Filter out samples with less than 600 total AVSs
proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot<-
  proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot[,colnames(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot<-
  proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot[,colnames(proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot<-
  proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot[,colnames(proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot<-
  proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot[,colnames(proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot <- 
  proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot %>% rownames_to_column("Taxa")
proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot <- 
  proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot %>% rownames_to_column("Taxa")
proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot <- 
  proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot %>% rownames_to_column("Taxa")
proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot <- 
  proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot %>% rownames_to_column("Taxa")

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g <- 
  proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g <- 
  proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g <- 
  proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g <- 
  proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)

#Simplify sample names in sample info tabs
mod_sample_info_tab_mc01$Sample <- gsub("^[^_]*_", "", mod_sample_info_tab_mc01$Sample)
mod_sample_info_tab_gc02$Sample <- gsub("^[^_]*_", "", mod_sample_info_tab_gc02$Sample)
mod_sample_info_tab_gc04$Sample <- gsub("^[^_]*_", "", mod_sample_info_tab_gc04$Sample)
mod_sample_info_tab_gc06$Sample <- gsub("^[^_]*_", "", mod_sample_info_tab_gc06$Sample)

#Add depth and core to tables
proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2 <- 
  proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g %>% left_join(mod_sample_info_tab_mc01)
proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2 <- 
  proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g %>% left_join(mod_sample_info_tab_gc02)
proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2 <- 
  proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g %>% left_join(mod_sample_info_tab_gc04)
proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2 <- 
  proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g %>% left_join(mod_sample_info_tab_gc06)

#Organize things in plot
proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2$Sample <-gsub("^.*_", "", proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2$Sample)
proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2$Sample <- factor(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2$Sample, levels=unique(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2$Sample))))])
proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2$Genus <- with(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2,factor(Class,Genus = rev(sort(unique(Genus)))))

#Make plot

head(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2)
ggplot(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2, aes(x=Sample, y = Class, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3", "#EA5835"),
    values=rescale(c(0,100)),
    limits=c(0,100)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/others_MC01_heatmap_phylum.pdf", width=12, height=5, dpi=300)

# Colorblind palettes:
Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#44BB99','#000000','#EE8866', '#FFAABB','#BBCC33', '#AAAA00', '#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#DDDDDD')

#Make pie chart

#Filter NA out
mc01_prot_plot <- na.omit(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2)
#Filter unclassifieds out
mc01_prot_plot <- mc01_prot_plot %>%
  filter(Taxa != "Unclassified")
#Group Desulfobacterota together
mc01_prot_plot <- mc01_prot_plot %>% 
  mutate(Taxa = if_else(Taxa %ni% c('Desulfobacterota_1','Desulfobacterota_2'), Taxa, "Desulfobacterota"))
#Get averages for depth_bin
mc01_prot_plot_depth_bin<- mc01_prot_plot %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion)) 
#Make bar plot
mc01_prot_plot_depth_bin %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "MC01")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140','141-220'))

#Make pie chart
ggplot(mc01_prot_plot_depth_bin, aes(x = " ", y =ave, group = Taxa, fill = Taxa)) +
  geom_bar(width = 1, stat = "identity") + scale_fill_manual(values=Tol_muted) +
  coord_polar("y", start = 0) + facet_wrap(~ Depth_bin) + theme_void()
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/MC01_sulfurmicrobes.pdf", width=12, height=5, dpi=300)


#GC02
#Organize things in plot
proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2$Sample <-gsub("^.*_", "", proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2$Sample)
proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2$Sample <- factor(proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2$Sample, levels=unique(proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2$Sample))))])
proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2$Class <- with(proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2,factor(Class,levels = rev(sort(unique(Class)))))

#Make plot
ggplot(proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2, aes(x=Sample, y = Class, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3", "#EA5835"),
    values=rescale(c(0,100)),
    limits=c(0,100)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/others_GC02_heatmap_phylum.pdf", width=12, height=5, dpi=300)

#Make pie chart

#Filter NA out
gc02_prot_plot <- na.omit(proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2)
gc02_prot_plot <- proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2
#Filter unclassifieds out
gc02_prot_plot <- gc02_prot_plot %>%
  filter(Taxa != "Unclassified")
#Group Desulfobacterota together
gc02_prot_plot <- gc02_prot_plot %>% 
  mutate(Taxa = if_else(Taxa %ni% c('Desulfobacterota_1','Desulfobacterota_2'), Taxa, "Desulfobacterota"))
head(gc02_prot_plot)
#Get averages for depth_bin
gc02_prot_plot_depth_bin<- gc02_prot_plot %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion)) 
agg_df <- aggregate(gc02_prot_plot_depth_bin$ave, by=list(gc02_prot_plot_depth_bin$Depth_bin), FUN=sum)
tail(gc02_prot_plot_depth_bin)
#Make bar plot
gc02_prot_plot_depth_bin %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR02")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140','141-220'))

#Make pie chart 
ggplot(gc02_prot_plot_depth_bin, aes(x = " ", y =ave, group = Taxa, fill = Taxa)) +
  geom_bar(width = 1, stat = "identity") + scale_fill_manual(values=Tol_muted) +
  coord_polar("y", start = 0) + facet_wrap(~ Depth_bin) + theme_void()
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/GC02_sulfurmicrobes.pdf", width=12, height=5, dpi=300)

#GC04
#Organize things in plot
proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2$Sample <-gsub("^.*_", "", proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2$Sample)
proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2$Sample <- factor(proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2$Sample, levels=unique(proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2$Sample))))])
proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2$Taxa <- with(proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2,factor(Taxa,levels = rev(sort(unique(Taxa)))))

#Filter NA out
gc04_prot_plot <- na.omit(proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2)
gc04_prot_plot <- proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2
#Filter unclassifieds out
gc04_prot_plot <- gc04_prot_plot %>%
  filter(Taxa != "Unclassified")
#Group Desulfobacterota together
gc04_prot_plot <- gc04_prot_plot %>% 
  mutate(Taxa = if_else(Taxa %ni% c('Desulfobacterota_1','Desulfobacterota_2'), Taxa, "Desulfobacterota"))
head(gc04_prot_plot)
#Normalize for bar plot
gc04_prot_plot_depth_bin<- gc04_prot_plot %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion)) 
agg_df <- aggregate(gc04_prot_plot_depth_bin$ave, by=list(gc04_prot_plot_depth_bin$Depth_bin), FUN=sum)

#Make bar plot
#Make plot
gc04_prot_plot_depth_bin %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR03")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140','141-220'))
#Make pie chart
ggplot(gc04_prot_plot, aes(x = " ", y =Proportion, group = Taxa, fill = Taxa)) +
  geom_bar(width = 1, stat = "identity") + scale_fill_manual(values=Tol_muted) +
  coord_polar("y", start = 0) + facet_wrap(~ Depth_bin) + theme_void()
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_sulfurmicrobes.pdf", width=12, height=5, dpi=300)

#GC06
#Organize things in plot
proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2$Sample <-gsub("^.*_", "", proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2$Sample)
proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2$Sample <- factor(proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2$Sample, levels=unique(proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2$Sample))))])
proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2$Class <- with(proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2,factor(Class,levels = rev(sort(unique(Class)))))
#Make heatmap
ggplot(proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2, aes(x=Sample, y = Class, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3", "#EA5835"),
    values=rescale(c(0,100)),
    limits=c(0,100)
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))

#Make pie chart
#Filter NA out
gc06_prot_plot <- na.omit(proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2)
gc06_prot_plot <- proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2
#Filter unclassifieds out
gc06_prot_plot <- gc06_prot_plot %>%
  filter(Taxa != "Unclassified")
#Group Desulfobacterota together
gc06_prot_plot <- gc06_prot_plot %>% 
  mutate(Taxa = if_else(Taxa %ni% c('Desulfobacterota_1','Desulfobacterota_2'), Taxa, "Desulfobacterota"))
type(gc06_prot_plot)
#Get averages for depth_bin
gc06_prot_plot_depth_bin<- gc06_prot_plot %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion)) 
agg_df <- aggregate(gc06_prot_plot_depth_bin$ave, by=list(gc06_prot_plot_depth_bin$Depth_bin), FUN=sum)
#Make bar plot
gc06_prot_plot_depth_bin %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR04")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140','141-220'))

#Make plot
ggplot(gc06_prot_plot_depth_bin, aes(x = " ", y =ave, group = Taxa, fill = Taxa)) +
  geom_bar(width = 1, stat = "identity") + scale_fill_manual(values=Tol_muted) +
  coord_polar("y", start = 0) + facet_wrap(~ Depth_bin) + theme_void()
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_sulfurmicrobes.pdf", width=12, height=5, dpi=300)

#Making bar charts for orders of Proteobacteria


#Create a copy of the tables for plots
proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot <- data.frame(proteobacteria_major_taxa_proportions_tab_o_mc01)
proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot <- data.frame(proteobacteria_major_taxa_proportions_tab_o_gc02)
proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot <- data.frame(proteobacteria_major_taxa_proportions_tab_o_gc04)
proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot <- data.frame(proteobacteria_major_taxa_proportions_tab_o_gc06)
colSums(proteobacteria_major_taxa_proportions_tab_o_mc01)

#Simplify sample names
colnames(proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot) <- gsub("^[^_]*_", "", colnames(proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot))
colnames(proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot) <- gsub("^[^_]*_", "", colnames(proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot))
colnames(proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot) <- gsub("^[^_]*_", "", colnames(proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot))
colnames(proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot) <- gsub("^[^_]*_", "", colnames(proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot))

#Filter out samples with less than 600 total AVSs
proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot<-
  proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot[,colnames(proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot<-
  proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot[,colnames(proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot<-
  proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot[,colnames(proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot<-
  proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot[,colnames(proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot <- 
  proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot %>% rownames_to_column("Taxa")
proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot <- 
  proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot %>% rownames_to_column("Taxa")
proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot <- 
  proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot %>% rownames_to_column("Taxa")
proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot <- 
  proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot %>% rownames_to_column("Taxa")

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot.g <- 
  proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot.g <- 
  proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot.g <- 
  proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g <- 
  proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)

#Simplify sample names in sample info tabs
mod_sample_info_tab_mc01$Sample <- gsub("^[^_]*_", "", mod_sample_info_tab_mc01$Sample)
mod_sample_info_tab_gc02$Sample <- gsub("^[^_]*_", "", mod_sample_info_tab_gc02$Sample)
mod_sample_info_tab_gc04$Sample <- gsub("^[^_]*_", "", mod_sample_info_tab_gc04$Sample)
mod_sample_info_tab_gc06$Sample <- gsub("^[^_]*_", "", mod_sample_info_tab_gc06$Sample)

#Add depth and core to tables
proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot.g2 <- 
  proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot.g %>% left_join(mod_sample_info_tab_mc01)
proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot.g2 <- 
  proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot.g %>% left_join(mod_sample_info_tab_gc02)
proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot.g2 <- 
  proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot.g %>% left_join(mod_sample_info_tab_gc04)
proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g2 <- 
  proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g %>% left_join(mod_sample_info_tab_gc06)

#Organize things in plot
proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot.g2$Sample <-gsub("^.*_", "", proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot.g2$Sample)
proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot.g2$Sample <- factor(proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot.g2$Sample, levels=unique(proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(proteobacteria_major_taxa_proportions_tab_g_mc01_for_plot.g2$Sample))))])
proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot.g2$Taxa <- with(proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot.g2,factor(Class,Genus = rev(sort(unique(Taxa)))))

#Make bar plot
#For MC01
#Get averages for depth_bin
proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot.g2<- proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot.g2 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
proteobacteria_major_taxa_proportions_tab_o_mc01_for_plot.g2 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "MC01")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140','141-220'))

#For CR02
#Organize things in plot
proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot.g2$Sample <-gsub("^.*_", "", proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot.g2$Sample)
proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot.g2$Sample <- factor(proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot.g2$Sample, levels=unique(proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(proteobacteria_major_taxa_proportions_tab_g_gc02_for_plot.g2$Sample))))])
proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot.g2$Taxa <- with(proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot.g2,factor(Class,Genus = rev(sort(unique(Taxa)))))

#Get averages for depth_bin
proteobacteria_order_gc02<- proteobacteria_major_taxa_proportions_tab_o_gc02_for_plot.g2 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
#Include only S cycling taxa
proteobacteria_order_gc02_plot <- proteobacteria_order_gc02 %>%
  filter(Taxa %in%  c("Rhodobacterales",
                      "Rhodospirillales",
                      "Sphingomonadales",
                      "unclassified_Alphaproteobacteria",
                      "unclassified_Gammaproteobacteria",
                      "unclassified_Proteobacteria",
                      "Oceanospirilalles_12",
                      "Acetobacterales",
                      "Burkholderiales"))
proteobacteria_order_gc02_plot %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "GC02")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140','141-220'))

#For CR04
#Organize things in plot
proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot.g2$Sample <-gsub("^.*_", "", proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot.g2$Sample)
proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot.g2$Sample <- factor(proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot.g2$Sample, levels=unique(proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(proteobacteria_major_taxa_proportions_tab_g_gc04_for_plot.g2$Sample))))])
proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot.g2$Taxa <- with(proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot.g2,factor(Class,Genus = rev(sort(unique(Taxa)))))

#Get averages for depth_bin
proteobacteria_order_gc04<- proteobacteria_major_taxa_proportions_tab_o_gc04_for_plot.g2 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
#Include only S cycling taxa
proteobacteria_order_gc04_plot <- proteobacteria_order_gc04 %>%
  filter(Taxa %in% c("Rhodobacterales",
                     "Rhodospirillales",
                     'Sphingomonadales',
                     'unclassified_Alphaproteobacteria',
                     'unclassified_Gammaproteobacteria',
                     'unclassified_Proteobacteria',
                     'Oceanospirillales_8',
                     'Oceanospirillales_12',
                     'Acetobacterales',
                     'Burkholderiales',
                     'Thiomicrospirales'))
proteobacteria_order_gc04_plot %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "GC04")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140','141-220'))

#For CR06
#Organize things in plot
proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g2$Sample <-gsub("^.*_", "", proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g2$Sample)
proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g2$Sample <- factor(proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g2$Sample, levels=unique(proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(proteobacteria_major_taxa_proportions_tab_g_gc06_for_plot.g2$Sample))))])
proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g2$Taxa <- with(proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g2,factor(Class,Genus = rev(sort(unique(Taxa)))))

#Get averages for depth_bin
proteobacteria_order_gc06<- proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g2 %>% 
  group_by(Depth_bin, Taxa) %>% summarise(ave = mean(Proportion))
#plot
proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g3 <- proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g2 %>%
  filter(Taxa %in% c("Rhodobacterales",
                     "Rhodospirillales",
                     'Sphingomonadales',
                     'unclassified_Alphaproteobacteria',
                     'unclassified_Gammaproteobacteria',
                     'unclassified_Proteobacteria',
                     'Oceanospirillales_8',
                     'Oceanospirillales_12',
                     'Acetobacterales',
                     'Burkholderiales',
                     'Thiomicrospirales',
                     'Pseudomonadales',
                     'uncultured'))
#Include only S cycling taxa
proteobacteria_order_gc06_plot <- proteobacteria_order_gc06 %>%
  filter(Taxa %in% c("Rhodobacterales",
                     "Rhodospirillales",
                     'Sphingomonadales',
                     'unclassified_Alphaproteobacteria',
                     'unclassified_Gammaproteobacteria',
                     'unclassified_Proteobacteria',
                     'Oceanospirillales_8',
                     'Oceanospirillales_12',
                     'Acetobacterales',
                     'Burkholderiales',
                     'Thiomicrospirales',
                     'Pseudomonadales',
                     'uncultured'))
# By depth bins
proteobacteria_order_gc06 %>% 
  mutate(Depth_bin = factor(Depth_bin, levels = unique(Depth_bin))) %>% 
  ggplot(aes(x = Depth_bin, y = ave, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "GC06")+
  scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140','141-220'))
# Not by depth bins
proteobacteria_major_taxa_proportions_tab_o_gc06_for_plot.g3 %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "GC06")
#scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140','141-220'))

#Making bar charts for genuses of Firmicutes
#Create a copy of the tables for plots
firmicutes_major_taxa_proportions_tab_g_mc01_for_plot <- data.frame(firmicutes_major_taxa_proportions_tab_g_mc01)
firmicutes_major_taxa_proportions_tab_g_gc02_for_plot <- data.frame(firmicutes_major_taxa_proportions_tab_g_gc02)
firmicutes_major_taxa_proportions_tab_g_gc04_for_plot <- data.frame(firmicutes_major_taxa_proportions_tab_g_gc04)
firmicutes_major_taxa_proportions_tab_g_gc06_for_plot <- data.frame(firmicutes_major_taxa_proportions_tab_g_gc06)

#Simplify sample names
colnames(firmicutes_major_taxa_proportions_tab_g_mc01_for_plot) <- gsub("^[^_]*_", "", colnames(firmicutes_major_taxa_proportions_tab_g_mc01_for_plot))
colnames(firmicutes_major_taxa_proportions_tab_g_gc02_for_plot) <- gsub("^[^_]*_", "", colnames(firmicutes_major_taxa_proportions_tab_g_gc02_for_plot))
colnames(firmicutes_major_taxa_proportions_tab_g_gc04_for_plot) <- gsub("^[^_]*_", "", colnames(firmicutes_major_taxa_proportions_tab_g_gc04_for_plot))
colnames(firmicutes_major_taxa_proportions_tab_g_gc06_for_plot) <- gsub("^[^_]*_", "", colnames(firmicutes_major_taxa_proportions_tab_g_gc06_for_plot))

#Filter out samples with less than 600 total AVSs
firmicutes_major_taxa_proportions_tab_g_mc01_for_plot<-
  firmicutes_major_taxa_proportions_tab_g_mc01_for_plot[,colnames(firmicutes_major_taxa_proportions_tab_g_mc01_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
firmicutes_major_taxa_proportions_tab_g_gc02_for_plot<-
  firmicutes_major_taxa_proportions_tab_g_gc02_for_plot[,colnames(firmicutes_major_taxa_proportions_tab_g_gc02_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
firmicutes_major_taxa_proportions_tab_g_gc04_for_plot<-
  firmicutes_major_taxa_proportions_tab_g_gc04_for_plot[,colnames(firmicutes_major_taxa_proportions_tab_g_gc04_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]
firmicutes_major_taxa_proportions_tab_g_gc06_for_plot<-
  firmicutes_major_taxa_proportions_tab_g_gc06_for_plot[,colnames(firmicutes_major_taxa_proportions_tab_g_gc06_for_plot) %in% colnames(filtered_absolute_asvs_for_manip)]

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
firmicutes_major_taxa_proportions_tab_g_mc01_for_plot <- 
  firmicutes_major_taxa_proportions_tab_g_mc01_for_plot %>% rownames_to_column("Taxa")
firmicutes_major_taxa_proportions_tab_g_gc02_for_plot <- 
  firmicutes_major_taxa_proportions_tab_g_gc02_for_plot %>% rownames_to_column("Taxa")
firmicutes_major_taxa_proportions_tab_g_gc04_for_plot <- 
  firmicutes_major_taxa_proportions_tab_g_gc04_for_plot %>% rownames_to_column("Taxa")
firmicutes_major_taxa_proportions_tab_g_gc06_for_plot <- 
  firmicutes_major_taxa_proportions_tab_g_gc06_for_plot %>% rownames_to_column("Taxa")

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
firmicutes_major_taxa_proportions_tab_g_mc01_for_plot.g <- 
  firmicutes_major_taxa_proportions_tab_g_mc01_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
firmicutes_major_taxa_proportions_tab_g_gc02_for_plot.g <- 
  firmicutes_major_taxa_proportions_tab_g_gc02_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
firmicutes_major_taxa_proportions_tab_g_gc04_for_plot.g <- 
  firmicutes_major_taxa_proportions_tab_g_gc04_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
firmicutes_major_taxa_proportions_tab_g_gc06_for_plot.g <- 
  firmicutes_major_taxa_proportions_tab_g_gc06_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)

#Simplify sample names in sample info tabs
mod_sample_info_tab_mc01$Sample <- gsub("^[^_]*_", "", mod_sample_info_tab_mc01$Sample)
mod_sample_info_tab_gc02$Sample <- gsub("^[^_]*_", "", mod_sample_info_tab_gc02$Sample)
mod_sample_info_tab_gc04$Sample <- gsub("^[^_]*_", "", mod_sample_info_tab_gc04$Sample)
mod_sample_info_tab_gc06$Sample <- gsub("^[^_]*_", "", mod_sample_info_tab_gc06$Sample)

#Add depth and core to tables
firmicutes_major_taxa_proportions_tab_g_mc01_for_plot.g2 <- 
  firmicutes_major_taxa_proportions_tab_g_mc01_for_plot.g %>% left_join(mod_sample_info_tab_mc01)
firmicutes_major_taxa_proportions_tab_g_gc02_for_plot.g2 <- 
  firmicutes_major_taxa_proportions_tab_g_gc02_for_plot.g %>% left_join(mod_sample_info_tab_gc02)
firmicutes_major_taxa_proportions_tab_g_gc04_for_plot.g2 <- 
  firmicutes_major_taxa_proportions_tab_g_gc04_for_plot.g %>% left_join(mod_sample_info_tab_gc04)
firmicutes_major_taxa_proportions_tab_g_gc06_for_plot.g2 <- 
  firmicutes_major_taxa_proportions_tab_g_gc06_for_plot.g %>% left_join(mod_sample_info_tab_gc06)

#Make bar plot

#For MC01
#Organize things in plot
firmicutes_major_taxa_proportions_tab_g_mc01_for_plot.g2$Sample <-gsub("^.*_", "", firmicutes_major_taxa_proportions_tab_g_mc01_for_plot.g2$Sample)
firmicutes_major_taxa_proportions_tab_g_mc01_for_plot.g2$Sample <- factor(firmicutes_major_taxa_proportions_tab_g_mc01_for_plot.g2$Sample, levels=unique(firmicutes_major_taxa_proportions_tab_g_mc01_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(firmicutes_major_taxa_proportions_tab_g_mc01_for_plot.g2$Sample))))])
firmicutes_major_taxa_proportions_tab_g_mc01_for_plot.g2$Taxa <- with(firmicutes_major_taxa_proportions_tab_g_mc01_for_plot.g2,factor(Class,Genus = rev(sort(unique(Taxa)))))
head(firmicutes_major_taxa_proportions_tab_g_mc01_for_plot.g2)
# Not by depth bins
firmicutes_major_taxa_proportions_tab_g_mc01_for_plot.g2 %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "SCB 01")
#scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140','141-220'))

#For CR 02
#Organize things in plot
firmicutes_major_taxa_proportions_tab_g_gc06_for_plot.g2$Sample <-gsub("^.*_", "", firmicutes_major_taxa_proportions_tab_g_gc06_for_plot.g2$Sample)
firmicutes_major_taxa_proportions_tab_g_gc06_for_plot.g2$Sample <- factor(firmicutes_major_taxa_proportions_tab_g_gc06_for_plot.g2$Sample, levels=unique(firmicutes_major_taxa_proportions_tab_g_gc06_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(firmicutes_major_taxa_proportions_tab_g_gc06_for_plot.g2$Sample))))])
firmicutes_major_taxa_proportions_tab_g_gc06_for_plot.g2$Taxa <- with(firmicutes_major_taxa_proportions_tab_g_gc06_for_plot.g2,factor(Class,Genus = rev(sort(unique(Taxa)))))
head(firmicutes_major_taxa_proportions_tab_g_gc06_for_plot.g2)
# Not by depth bins
firmicutes_major_taxa_proportions_tab_g_gc06_for_plot.g2 %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Proportion, fill = Taxa)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_manual(values = Tol_muted) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "% ASVs from S cycling taxa", title = "CR 02")
#scale_x_discrete(limits = c("0-20", "21-40", "41-60", "61-100",'101-140','141-220'))




#For absolute abundances

#Create a copy of the tables for plots
phyla_and_unidentified_counts_tab_mc01_for_plot <- data.frame(phyla_and_unidentified_counts_tab_mc01)
phyla_and_unidentified_counts_tab_gc02_for_plot <- data.frame(phyla_and_unidentified_counts_tab_gc02)
phyla_and_unidentified_counts_tab_gc04_for_plot <- data.frame(phyla_and_unidentified_counts_tab_gc04)
phyla_and_unidentified_counts_tab_gc06_for_plot <- data.frame(phyla_and_unidentified_counts_tab_gc06)

colnames(phyla_and_unidentified_counts_tab_mc01_for_plot) <- gsub("^[^_]*_", "", colnames(phyla_and_unidentified_counts_tab_mc01_for_plot))
colnames(phyla_and_unidentified_counts_tab_gc02_for_plot) <- gsub("^[^_]*_", "", colnames(phyla_and_unidentified_counts_tab_gc02_for_plot))
colnames(phyla_and_unidentified_counts_tab_gc04_for_plot) <- gsub("^[^_]*_", "", colnames(phyla_and_unidentified_counts_tab_gc04_for_plot))
colnames(phyla_and_unidentified_counts_tab_gc06_for_plot) <- gsub("^[^_]*_", "", colnames(phyla_and_unidentified_counts_tab_gc06_for_plot))


# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
phyla_and_unidentified_counts_tab_mc01_for_plot <- 
  phyla_and_unidentified_counts_tab_mc01_for_plot %>% rownames_to_column("Taxa")
phyla_and_unidentified_counts_tab_gc02_for_plot <- 
  phyla_and_unidentified_counts_tab_gc02_for_plot %>% rownames_to_column("Taxa")
phyla_and_unidentified_counts_tab_gc04_for_plot <- 
  phyla_and_unidentified_counts_tab_gc04_for_plot %>% rownames_to_column("Taxa")
phyla_and_unidentified_counts_tab_gc06_for_plot <- 
  phyla_and_unidentified_counts_tab_gc06_for_plot %>% rownames_to_column("Taxa")

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
phyla_and_unidentified_counts_tab_mc01_for_plot.g <- 
  phyla_and_unidentified_counts_tab_mc01_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
phyla_and_unidentified_counts_tab_gc02_for_plot.g <- 
  phyla_and_unidentified_counts_tab_gc02_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
phyla_and_unidentified_counts_tab_gc04_for_plot.g <- 
  phyla_and_unidentified_counts_tab_gc04_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)
phyla_and_unidentified_counts_tab_gc06_for_plot.g <- 
  phyla_and_unidentified_counts_tab_gc06_for_plot %>% 
  pivot_longer(!Taxa, names_to = "Sample", values_to = "Proportion") %>% 
  data.frame(check.names = FALSE)

#Add depth and core to tables
phyla_and_unidentified_counts_tab_mc01_for_plot.g2 <- 
  phyla_and_unidentified_counts_tab_mc01_for_plot.g %>% left_join(mod_sample_info_tab_mc01)
phyla_and_unidentified_counts_tab_gc02_for_plot.g2 <- 
  phyla_and_unidentified_counts_tab_gc02_for_plot.g %>% left_join(mod_sample_info_tab_gc02)
phyla_and_unidentified_counts_tab_gc04_for_plot.g2 <- 
  phyla_and_unidentified_counts_tab_gc04_for_plot.g %>% left_join(mod_sample_info_tab_gc04)
phyla_and_unidentified_counts_tab_gc06_for_plot.g2 <- 
  phyla_and_unidentified_counts_tab_gc06_for_plot.g %>% left_join(mod_sample_info_tab_gc06)

#MC01
#Organize things in plot
phyla_and_unidentified_counts_tab_mc01_for_plot.g2$Sample <-gsub("^.*_", "", phyla_and_unidentified_counts_tab_mc01_for_plot.g2$Sample)
phyla_and_unidentified_counts_tab_mc01_for_plot.g2$Sample <- factor(phyla_and_unidentified_counts_tab_mc01_for_plot.g2$Sample, levels=unique(phyla_and_unidentified_counts_tab_mc01_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(phyla_and_unidentified_counts_tab_mc01_for_plot.g2$Sample))))])
phyla_and_unidentified_counts_tab_mc01_for_plot.g2$Taxa <- with(phyla_and_unidentified_counts_tab_mc01_for_plot.g2,factor(Taxa, levels = rev(sort(unique(Taxa)))))

#Get total Phylum ASVs for depth
df_mc01_sum <-
  phyla_and_unidentified_counts_tab_mc01_for_plot.g2 %>%
  group_by(Depth) %>%
  summarise(Freq = sum(Proportion))
#Make plot
df_mc01_sum %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Freq)) +
  geom_bar(width = 0.6, stat = "identity") 
theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "Total ASVs") 
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/mc01_total_AVSs_phylum.pdf", width=12, height=5, dpi=300)

#Make plot of ASVs by Phylum by depth
#Group Desulfobacterota together
phyla_and_unidentified_counts_tab_mc01_for_plot.g2 <- phyla_and_unidentified_counts_tab_mc01_for_plot.g2 %>% 
  mutate(Taxa = if_else(Taxa %ni% c('Desulfobacterota_1','Desulfobacterota_2'), Taxa, "Desulfobacterota"))
type(phyla_and_unidentified_counts_tab_mc01_for_plot.g2)
ggplot(phyla_and_unidentified_counts_tab_mc01_for_plot.g2, aes(x=Sample, y = Taxa, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3", "#EA5835"),
    values=rescale(c(0,max(phyla_and_unidentified_counts_tab_mc01_for_plot.g2$Proportion))),
    limits=c(0,max(phyla_and_unidentified_counts_tab_mc01_for_plot.g2$Proportion))
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/MC01_heatmap_phylum_filtered.pdf", width=12, height=5, dpi=300)

#GC02
#Organize things in plot
phyla_and_unidentified_counts_tab_gc02_for_plot.g2$Sample <-gsub("^.*_", "", phyla_and_unidentified_counts_tab_gc02_for_plot.g2$Sample)
phyla_and_unidentified_counts_tab_gc02_for_plot.g2$Sample <- factor(phyla_and_unidentified_counts_tab_gc02_for_plot.g2$Sample, levels=unique(phyla_and_unidentified_counts_tab_gc02_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(phyla_and_unidentified_counts_tab_gc02_for_plot.g2$Sample))))])
phyla_and_unidentified_counts_tab_gc02_for_plot.g2$Taxa <- with(phyla_and_unidentified_counts_tab_gc02_for_plot.g2,factor(Taxa, levels = rev(sort(unique(Taxa)))))

#Get total Phylum ASVs for depth
df_gc02_sum <-
  phyla_and_unidentified_counts_tab_gc02_for_plot.g2 %>%
  group_by(Depth) %>%
  summarise(Freq = sum(Proportion))
#Make plot
df_gc02_sum %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Freq)) +
  geom_bar(width = 0.6, stat = "identity") 
theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "Total ASVs") 
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_total_AVSs_phylum.pdf", width=12, height=5, dpi=300)


#Make plot of ASVs by Phylum by depth
#Group Desulfobacterota together
phyla_and_unidentified_counts_tab_gc02_for_plot.g2 <- phyla_and_unidentified_counts_tab_gc02_for_plot.g2 %>% 
  mutate(Taxa = if_else(Taxa %ni% c('Desulfobacterota_1','Desulfobacterota_2'), Taxa, "Desulfobacterota"))
type(phyla_and_unidentified_counts_tab_gc02_for_plot.g2)
ggplot(phyla_and_unidentified_counts_tab_gc02_for_plot.g2, aes(x=Sample, y = Taxa, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3", "#EA5835"),
    values=rescale(c(0,max(phyla_and_unidentified_counts_tab_gc02_for_plot.g2$Proportion))),
    limits=c(0,max(phyla_and_unidentified_counts_tab_gc02_for_plot.g2$Proportion))
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/others_gc02_heatmap_phylum.pdf", width=12, height=5, dpi=300)

#Filter out depths with less than 750 total ASVs
depths_to_keep_gc02 <- subset(df_gc02_sum, as.numeric(Freq) >= 750)
filtered_absolute_asvs_gc02 <- filter(phyla_and_unidentified_counts_tab_gc02_for_plot.g2,
                                      phyla_and_unidentified_counts_tab_gc02_for_plot.g2$Depth %in% depths_to_keep_gc02$Depth)
#Get relative abundance table for the filtered depths
filtered_relative_asvs_gc02 <- 
  filtered_absolute_asvs_gc02 %>%
  group_by(Depth) %>%
  mutate(freq = Proportion*100 / sum(Proportion))
head(filtered_relative_asvs_gc02)


#Make plot of ASVs by Phylum by depth
#Group Desulfobacterota together
filtered_relative_asvs_gc02 <- filtered_relative_asvs_gc02 %>% 
  mutate(Taxa = if_else(Taxa %ni% c('Desulfobacterota_1','Desulfobacterota_2'), Taxa, "Desulfobacterota"))
ggplot(filtered_relative_asvs_gc02, aes(x=Sample, y = Taxa, fill=freq)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3", "#EA5835"),
    values=rescale(c(0,max(filtered_relative_asvs_gc02$freq))),
    limits=c(0,max(filtered_relative_asvs_gc02$freq))
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc02_heatmap_phylum_filtered.pdf", width=12, height=5, dpi=300)


#GC04
#Organize things in plot
phyla_and_unidentified_counts_tab_gc04_for_plot.g2$Sample <-gsub("^.*_", "", phyla_and_unidentified_counts_tab_gc04_for_plot.g2$Sample)
phyla_and_unidentified_counts_tab_gc04_for_plot.g2$Sample <- factor(phyla_and_unidentified_counts_tab_gc04_for_plot.g2$Sample, levels=unique(phyla_and_unidentified_counts_tab_gc04_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(phyla_and_unidentified_counts_tab_gc04_for_plot.g2$Sample))))])
phyla_and_unidentified_counts_tab_gc04_for_plot.g2$Taxa <- with(phyla_and_unidentified_counts_tab_gc04_for_plot.g2,factor(Taxa, levels = rev(sort(unique(Taxa)))))

#Get total Phylum ASVs for depth
df_gc04_sum <-
  phyla_and_unidentified_counts_tab_gc04_for_plot.g2 %>%
  group_by(Depth) %>%
  summarise(Freq = sum(Proportion))
#Make plot
df_gc04_sum %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Freq)) +
  geom_bar(width = 0.6, stat = "identity") 
theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "Total ASVs") 
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_total_AVSs_phylum.pdf", width=12, height=5, dpi=300)


#Make plot of ASVs by Phylum by depth
#Group Desulfobacterota together
phyla_and_unidentified_counts_tab_gc04_for_plot.g2 <- phyla_and_unidentified_counts_tab_gc04_for_plot.g2 %>% 
  mutate(Taxa = if_else(Taxa %ni% c('Desulfobacterota_1','Desulfobacterota_2'), Taxa, "Desulfobacterota"))
type(phyla_and_unidentified_counts_tab_gc04_for_plot.g2)
ggplot(phyla_and_unidentified_counts_tab_gc04_for_plot.g2, aes(x=Sample, y = Taxa, fill=Proportion)) + 
  geom_tile(width=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3", "#EA5835"),
    values=rescale(c(0,max(phyla_and_unidentified_counts_tab_gc04_for_plot.g2$Proportion))),
    limits=c(0,max(phyla_and_unidentified_counts_tab_gc04_for_plot.g2$Proportion))
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_heatmap_phylum_absolute.pdf", width=12, height=5, dpi=300)

column <- "Taxa"
gc04_abs_Desulfobacterota <- 
  #phyla_and_unidentified_counts_tab_gc04_for_plot.g2[phyla_and_unidentified_counts_tab_gc04_for_plot.g2$Taxa == c('Desulfobacterota_1','Desulfobacterota_2')]
  phyla_and_unidentified_counts_tab_gc04_for_plot.g2 %>% filter_(paste(column, "==", c('Desulfobacterota_1','Desulfobacterota_2')))
#Filter out depths with less than 500 total ASVs
depths_to_keep_gc04 <- subset(df_gc04_sum, as.numeric(Freq) >= 750)
filtered_absolute_asvs_gc04 <- filter(phyla_and_unidentified_counts_tab_gc04_for_plot.g2,
                                      phyla_and_unidentified_counts_tab_gc04_for_plot.g2$Depth %in% depths_to_keep_gc04$Depth)
#Get relative abundance table for the filtered depths
filtered_relative_asvs_gc04 <- 
  filtered_absolute_asvs_gc04 %>%
  group_by(Depth) %>%
  mutate(freq = Proportion*100 / sum(Proportion))
head(filtered_relative_asvs_gc04)


#Make plot of ASVs by Phylum by depth
#Group Desulfobacterota together
filtered_relative_asvs_gc04 <- filtered_relative_asvs_gc04 %>% 
  mutate(Taxa = if_else(Taxa %ni% c('Desulfobacterota_1','Desulfobacterota_2'), Taxa, "Desulfobacterota"))
ggplot(filtered_relative_asvs_gc04, aes(x=Sample, y = Taxa, fill=freq)) + 
  geom_tile(w
            idth=0.8, height=0.8) +
  scale_fill_gradientn(
    colors=c("white","#2699D3", "#EA5835"),
    values=rescale(c(0,max(filtered_relative_asvs_gc04$freq))),
    limits=c(0,max(filtered_relative_asvs_gc04$freq))
  ) +
  guides(fill=guide_colorbar(ticks.colour = NA, barheight = 15)) +
  theme_classic() + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=0.5))
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc04_heatmap_phylum_filtered.pdf", width=12, height=5, dpi=300)

#GC06
#Organize things in plot
phyla_and_unidentified_counts_tab_gc06_for_plot.g2$Sample <-gsub("^.*_", "", phyla_and_unidentified_counts_tab_gc06_for_plot.g2$Sample)
phyla_and_unidentified_counts_tab_gc06_for_plot.g2$Sample <- factor(phyla_and_unidentified_counts_tab_gc06_for_plot.g2$Sample, levels=unique(phyla_and_unidentified_counts_tab_gc06_for_plot.g2$Sample)[order(as.numeric(gsub("^.*_","", unique(phyla_and_unidentified_counts_tab_gc06_for_plot.g2$Sample))))])
phyla_and_unidentified_counts_tab_gc06_for_plot.g2$Taxa <- with(phyla_and_unidentified_counts_tab_gc06_for_plot.g2,factor(Taxa, levels = rev(sort(unique(Taxa)))))

#Get total Phylum ASVs for depth
df_gc06_sum <-
  phyla_and_unidentified_counts_tab_gc06_for_plot.g2 %>%
  group_by(Depth) %>%
  summarise(Freq = sum(Proportion))
#Make plot
df_gc06_sum %>% 
  mutate(Depth = factor(Depth, levels = unique(Depth))) %>% 
  ggplot(aes(x = Depth, y = Freq)) +
  geom_bar(width = 0.6, stat = "identity") 
theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x = "Depth (cm)", y = "Total ASVs") 
ggsave(file="/export/data1/projects/daniela/SR2113/June_Nov_2022_16S/plots/gc06_total_AVSs_phylum.pdf", width=12, height=5, dpi=300)

