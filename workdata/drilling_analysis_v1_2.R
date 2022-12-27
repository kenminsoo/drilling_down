##--##
#Drilling down analysis v1.1
##Most basic querying feature has been completed on python side
###Analysis of basic query begins here
##Expected features and goal with v1
##--##
##FEATURES##

#Read in data from python query
#Clean data from low expression samples 
#Create heatmaps based upon categories with unmodified current data
#Create heatmaps with ratio normalization method
#Create boxplots comparing phenotypes
#Run T-Test denoting significant differences
#Create KM survival plots
#Conversions from RPKM to TPM & subseqent analysis with new, introduced datasets that will go thru pipe
#Create a file with figures named easily
#Create a file with an overall summary of statistical tests

##FEATURES END##
#import packages 
library(readr)
library(magrittr)
library(dplyr)
library(janitor)
library(ggplot2)
library(hrbrthemes)
library(tidyr)
library(viridis)
library(tidyverse)
library(survival)
library(gplots)
library(patchwork)

##Define commonly used functions##

#export list of plots as pdfs
export_figures <- function(list_of_plots, plot_type = 'none'){
  n = 1
  dir.create(paste("figures/",plot_type,sep=''))
  for (plotted in list_of_plots){
    ggsave(paste('figures/',plot_type,'/plot', n, ".pdf", sep = ""), plot = plotted)
    n = n + 1
  }
}

export_figures_ngg <- function(list_of_plots, plot_type = 'none'){
  n = 1
  dir.create(paste("figures/",plot_type,sep=''))
  for (plotted in list_of_plots){
    ggsave(paste('figures/',plot_type,'/plot', n, ".pdf", sep = ""), plot = plotted)
    n = n + 1
  }
}

##Define commonly used functions##

#import data
dataset_names <- scan("exported_data.csv", what ="", sep=',') #generated from python
dataset_files <- list()
for (dataset in dataset_names){
  temp <- paste(dataset, '.csv',sep ='')
  dataset_files <- append(dataset_files, temp)
  rm(temp)
  rm(dataset)
}

n = 1
for (dataset in dataset_files){
  assign(dataset_names[n], read.csv(paste(dataset, sep="")))
  n = n + 1
}

##close out import sequence
rm(dataset)
rm(dataset_names)
rm(n)                     
rm(dataset_files)

#next step is to produce basic histograms regarding the data
#import a config file from python, where to cut off values
rm_low_exp <- "Y"
num_genes <- 17
universal_units <- "log2(tpm + .001)"
gois_hks <- colnames(genomic_data[2:17])
gois <- colnames(genomic_data[2:8])
combined_data$mean_expression <- rowMeans(genomic_data[2:num_genes+1])
genomic_data$mean_expression <- rowMeans(genomic_data[2:num_genes+1])

#plot histogram
histograms <- list()
histograms[['Inital_Hist']] <- ggplot(data = combined_data, aes(x = mean_expression, group = X_study, fill = X_study)) +
  geom_density(adjust=1.5, alpha=0.4) + ggtitle("Average Expression for Selected Genes") +
  xlab(universal_units) + labs(fill = 'Study Name')
#filter dataset for mean expression <1
if (rm_low_exp == "Y"){
  combined_data <- combined_data[combined_data$mean_expression > 1,]
  genomic_data <- genomic_data[genomic_data$mean_expression > 1, ]
} else {
  print("Will not filter")
}
  
#remove duplicate columns
dup_col <- duplicated(as.list(phenotypic_data))
phenotypic_data <- phenotypic_data[!dup_col]
dup_col <- duplicated(as.list(combined_data))
combined_data <- combined_data[!dup_col]
rm(dup_col)

histograms[['Filtered_Hist']] <- ggplot(data = combined_data, aes(x = mean_expression, group = X_study, fill = X_study)) +
    geom_density(adjust=1.5, alpha=0.4) + ggtitle("Average Expression > 1 for Selected Genes") +
    xlab(universal_units) + labs(fill = 'Study Name')

export_figures(histograms, 'hist2')
#close off histograms
rm(histograms)

##analysis, two types of each, log2 & ratios
#heatmaps
#boxplots
#barcharts
#SD/SE^^
#t-tests
#survival analysis
##analysis##

#summary of data for plotting

mean_dfs <- list()
mean_dfs_asc <- list()
phenotypes <- colnames(phenotypic_data[2:6])
#create loop to sort datasets
for (phenotype in phenotypes){
  mean_dfs[[phenotype]] <- combined_data %>% aggregate(by = list(combined_data[[phenotype]]), FUN = mean)
  mean_dfs[[phenotype]] <- mean_dfs[[phenotype]] %>% arrange(., desc(PTRH2))
  mean_dfs_asc[[phenotype]] <- mean_dfs[[phenotype]] %>% arrange(., PTRH2)
}

#select the things that we want, genes
new_mean_dfs <- lapply(mean_dfs, function(x) x%>% select(all_of(gois_hks)))
goi_mean_dfs <- lapply(new_mean_dfs, function(x) x%>% select(all_of(gois)))
goi_mean_dfs_asc <- lapply(mean_dfs_asc, function(x) x%>% select(all_of(gois)))

#create openings for the incoming matracies
new_matrix <- list()
goi_matrix <- list()
asc_matrix <- list()

#transform dataframes into matracies
for (phenotype in phenotypes){
  new_matrix[[phenotype]] <- as.matrix(new_mean_dfs[[phenotype]])
  goi_matrix[[phenotype]] <- as.matrix(goi_mean_dfs[[phenotype]])
  asc_matrix[[phenotype]] <- as.matrix(goi_mean_dfs_asc[[phenotype]])
}

#give the rows names
#probably make this into a function at a later date since this is repeated
new_matrix_named <- list()

for (phenotype in phenotypes){
  name_bank <- mean_dfs[[phenotype]]
  temp_matrix <- new_matrix[[phenotype]]
  rnames <- as.list(name_bank[1])
  rownames(temp_matrix) <- c(rnames$Group.1)
  new_matrix_named[[phenotype]] <- temp_matrix
}

goi_matrix_names <- list()

for (phenotype in phenotypes){
  name_bank <- mean_dfs[[phenotype]]
  temp_matrix <- goi_matrix[[phenotype]]
  rnames <- as.list(name_bank[1])
  rownames(temp_matrix) <- c(rnames$Group.1)
  goi_matrix_names[[phenotype]] <- temp_matrix
}

goi_matrix_names_asc <- list()

for (phenotype in phenotypes){
  name_bank <- mean_dfs_asc[[phenotype]]
  temp_matrix <- asc_matrix[[phenotype]]
  rnames <- as.list(name_bank[1])
  rownames(temp_matrix) <- c(rnames$Group.1)
  goi_matrix_names_asc[[phenotype]] <- temp_matrix
}

#create plots with the matracies, heatplots
for (phenotype in phenotypes){
  pdf(file = paste("figures/Heatmaps/",phenotype,'.pdf',sep=''), pointsize = 12, width = 10, height = 15)
  
  heatmap.2(goi_matrix_names[[phenotype]], margins = c(6,12), Colv = FALSE, Rowv = FALSE,
            col = redgreen(75),
            key.xlab("log2(norm+1"), key = TRUE,
            main = paste("Avg Gene Expression for\n", phenotype, 'Descending'))
  
  dev.off()
}

for (phenotype in phenotypes){
  pdf(file = paste("figures/Heatmaps/",phenotype,'asc','.pdf',sep=''), pointsize = 12, width = 10, height = 15)
  
  heatmap.2(goi_matrix_names_asc[[phenotype]], margins = c(6,12), Colv = FALSE, Rowv = FALSE,
            col = redgreen(75),
            key.xlab("log2(norm+1"), key = TRUE,
            main = paste("Avg Gene Expression for\n", phenotype, 'Ascending'))
  
  dev.off()
}

for (phenotype in phenotypes){
  pdf(file = paste("figures/Heatmaps/",phenotype,'bp','.pdf',sep=''), pointsize = 12, width = 10, height = 15)
  
  heatmap.2(new_matrix_named[[phenotype]], margins = c(6,12), Colv = FALSE, Rowv = FALSE,
            col = redgreen(75),
            key.xlab("log2(norm+1"), key = TRUE,
            main = paste("Avg Gene Expression for\n", phenotype, 'Big Picture'))
  
  dev.off()
}

##Let's just do informative plots, KM & Box from this point on. 

plot_gene_goi <- function(gene){
  print(ggplot(combined_data, aes_(combined_data$X_sample_type, as.name(gene))) + 
          geom_boxplot() + coord_flip() + scale_fill_grey() + 
          labs(x = "Sample Type", y = "Expression (log2(norm_count + 1))", title = paste("Average Expression of Experimental", as.character(gene), "\nFor Different Sample Types"))) +
    ylim(0,18)
}

detailed_interest_l <- list(2,3,7,14,15,17,27,33,36,39,42,51,53,56,63)
for (num in detailed_interest_l){
  print(detailed_unique[num])
  detailed_interest[[detailed_unique[num]]] <- detailed_unique[num]
}


interest_dataset <- combined_data %>% 
  filter(combined_data$detailed_category %in% detailed_interest)
###

#create cohorts
#MEL 15, 39, 51
#NSCLC 14 17 33
#NB 42 56
#SAR 13 21 62 67
detailed_unique <- unique(combined_data$detailed_category)
mel_interest_l <- list(15, 39, 51)
mel_interest <- list()
for (num in mel_interest_l){
  print(detailed_unique[num])
  mel_interest[[detailed_unique[num]]] <- detailed_unique[num]
}

nsclc_interest_l <- list(14, 17, 33)
nsclc_interest <- list()
for (num in nsclc_interest_l){
  print(detailed_unique[num])
  nsclc_interest[[detailed_unique[num]]] <- detailed_unique[num]
}

nb_interest_l <- list(42, 56)
nb_interest <- list()
for (num in nb_interest_l){
  print(detailed_unique[num])
  nb_interest[[detailed_unique[num]]] <- detailed_unique[num]
}

sar_interest_l <- list(13, 21, 62, 67)
sar_interest <- list()
for (num in sar_interest_l){
  print(detailed_unique[num])
  sar_interest[[detailed_unique[num]]] <- detailed_unique[num]
}

rm(nb_interest_l)
rm(nsclc_interest_l)
rm(sar_interest_l)
rm(mel_interest_l)


nb_dataset <- combined_data %>% 
  filter(combined_data$detailed_category %in% nb_interest)

nsclc_dataset <- combined_data %>% 
  filter(combined_data$detailed_category %in% nsclc_interest)

sar_dataset <- combined_data %>% 
  filter(combined_data$detailed_category %in% sar_interest)

mel_dataset <- combined_data %>% 
  filter(combined_data$detailed_category %in% mel_interest)
###
plot_detailed <- function(dataset_i,gene){
  (ggplot(dataset_i, aes_(dataset_i$detailed_category, as.name(gene))) + 
          geom_boxplot() + coord_flip() + scale_fill_grey() + 
          labs(x = "Sample Type", y = "Expression (log2(tpm + 0.001))", title = paste(as.character(gene)))) +
    ylim(0,18)
}

ptrh2_nb <- plot_detailed(nb_dataset,'PTRH2')
tle5_nb <- plot_detailed(nb_dataset,'TLE5')
tle1_nb <- plot_detailed(nb_dataset,'TLE1')
tcf7_nb <- plot_detailed(nb_dataset,'TCF7')
rela_nb <- plot_detailed(nb_dataset,'RELA')


comb_nb <- ptrh2_nb / (tle5_nb|tle1_nb) / (tcf7_nb|rela_nb) + plot_annotation(title = 'Expression Profile of Apoptotic Factors PTRH2, TLE5, and TLE1 in Neuroblastoma Related Tissues',
                                                            caption = 'NB Data: TARGET \nHealthy Data: GTEX')
#############
ptrh2_nsclc <- plot_detailed(nsclc_dataset,'PTRH2')
tle5_nsclc <- plot_detailed(nsclc_dataset,'TLE5')
tle1_nsclc <- plot_detailed(nsclc_dataset,'TLE1')
tcf7_nsclc <- plot_detailed(nsclc_dataset,'TCF7')
rela_nsclc <- plot_detailed(nsclc_dataset,'RELA')

comb_nsclc <- ptrh2_nsclc / (tle5_nsclc|tle1_nsclc) / (tcf7_nsclc|rela_nsclc) + plot_annotation(title = 'Expression Profile of Apoptotic Factors PTRH2, TLE5, and TLE1 in NSCLC & Related Tissues',
                                                                      caption = 'NSCLC Data: TCGA \nHealthy Data: GTEX')
##############

ptrh2_sar <- plot_detailed(sar_dataset,'PTRH2')
tle5_sar <- plot_detailed(sar_dataset,'TLE5')
tle1_sar <- plot_detailed(sar_dataset,'TLE1')
tcf7_sar <- plot_detailed(sar_dataset,'TCF7')
rela_sar <- plot_detailed(sar_dataset,'RELA')

comb_sar <- ptrh2_sar / (tle5_sar|tle1_sar) / (tcf7_sar|rela_sar) + plot_annotation(title = 'Expression Profile of Apoptotic Factors PTRH2, TLE5, and TLE1 in Sarcoma & Related Tissues',
                                                              caption = 'Sarcoma Data: TCGA \nHealthy Data: GTEX')
##################

ptrh2_mel <- plot_detailed(mel_dataset,'PTRH2')
tle5_mel <- plot_detailed(mel_dataset,'TLE5')
tle1_mel <- plot_detailed(mel_dataset,'TLE1')
tcf7_mel <- plot_detailed(mel_dataset,'TCF7')
rela_mel <- plot_detailed(mel_dataset,'RELA')

comb_mel <- ptrh2_mel / (tle5_mel|tle1_mel) / (tcf7_mel|rela_mel) + plot_annotation(title = 'Expression Profile of Apoptotic Factors PTRH2, TLE5, and TLE1 in Melanomas & Related Tissues',
                                                              caption = 'MEL Data: TCGA \nHealthy Data: GTEX')

ggsave('exp_mel.pdf', plot = comb_mel, device = 'pdf', width = 12, height = 9)
ggsave('exp_nsclc.pdf', plot = comb_nsclc, device = 'pdf', width = 12, height = 9)
ggsave('exp_sar.pdf', plot = comb_sar, device = 'pdf', width = 12, height = 9)
ggsave('exp_nb.pdf', plot = comb_nb, device = 'pdf', width = 12, height = 9)

## Okay now let's prep for KM Plots



KM <- function(category_num, name_of_plot){
  temp_dataset <- combined_data %>% 
    filter(combined_data$detailed_category %in% detailed_unique[category_num])
  
  temp_dataset$score <- (temp_dataset$RELA)
  temp_dataset$assignment <- NA
  temp_dataset$assignment[temp_dataset$score >= (median(temp_dataset$score))] <- 1 #red
  temp_dataset$assignment[temp_dataset$score < (median(temp_dataset$score))] <- 0 #blue
  
  #view(temp_dataset)
  
  km_model_temp <- survfit(Surv(temp_dataset$OS.time,
                                temp_dataset$OS) ~ 
                             temp_dataset$assignment)
  
  temp_plot <- plot(km_model_temp, conf.int = F, xlab = "Time (Days)", ylab = "% Alive", main = paste("KM -", name_of_plot),
       col = c("Blue","Red"), las = 1, mark.time =TRUE)
  
  print(survdiff(Surv(temp_dataset$OS.time,
                temp_dataset$OS) ~ 
             temp_dataset$assignment))
  
  return(temp_plot)
}

KM_LK <- KM(17, 'Lung Adenocarcinoma')
KM_SAR <- KM(21, 'Sarcoma')
KM_SQ <- KM(33, 'Squamous Cell LK')
KM_MEL <- KM(39, 'Melanoma')
KM_MEL2 <- KM(51, 'Melanoma')
KM_SARUT <- KM(67, 'UtSAR')


