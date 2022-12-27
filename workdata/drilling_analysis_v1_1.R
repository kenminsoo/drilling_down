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
rm_low_exp <- readline(prompt = "Remove Low Expression? Y/N: ")
num_genes <- 24
universal_units <- "log2(norm_count + 1)"
gois_hks <- colnames(genomic_data[2:14])
gois <- colnames(genomic_data[2:6])
combined_data$mean_expression <- rowMeans(genomic_data[2:num_genes+1])
genomic_data$mean_expression <- rowMeans(genomic_data[2:num_genes+1])

#plot histogram
histograms <- list()
histograms[['Inital_Hist']] <- ggplot(data = combined_data, aes(x = mean_expression, group = X_study, fill = X_study)) +
  geom_density(adjust=1.5, alpha=0.4) + ggtitle("Average Expression for Selected Genes") +
  xlab(universal_units) + labs(fill = 'Study Name')
  theme_ipsum()
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
  theme_ipsum()

export_figures(histograms, 'hist')
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
new_matrix <- list()
goi_matrix <- list()
asc_matrix <- list()

for (phenotype in phenotypes){
  new_matrix[[phenotype]] <- as.matrix(new_mean_dfs[[phenotype]])
  goi_matrix[[phenotype]] <- as.matrix(goi_mean_dfs[[phenotype]])
  asc_matrix[[phenotype]] <- as.matrix(goi_mean_dfs_asc[[phenotype]])
}

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

##Finished with heatmaps for phenotype vs. gene expression
##Next is phenotype vs. phenotype


