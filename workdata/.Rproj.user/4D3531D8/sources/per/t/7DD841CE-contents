#import a list of the datasets we may need
datasets <- scan("data/exported_files.txt", what="", sep=',')
n <- length(datasets)
datasets_fixed <- datasets[-n]
rm(n)
rm(datasets)
datasets <- datasets_fixed
rm(datasets_fixed)

for (dataset in datasets){
  assign(dataset, read.csv(paste(dataset, sep="")))
}

phenotype <- `data/decoded_phenotype.csv`
rm(`data/decoded_phenotype.csv`)
genomic <- `data/goi_hk_prev_genomic_data.csv`
rm(`data/goi_hk_prev_genomic_data.csv`)
hk <- `data/house_keepers.csv`
rm(`data/house_keepers.csv`)
goi <- `data/goi.csv`
rm(`data/goi.csv`)
prev <- `data/prev_est.csv`
rm(`data/prev_est.csv`)
#this is stupid, fix it later

#import needed packages
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

#define functions

#--Data Cleaning--#

#summary information
sample_types <- aggregate(data.frame(count = phenotype$X_sample_type), list(value = phenotype$X_sample_type), length)
primary_sites <- aggregate(data.frame(count = phenotype$X_primary_site), list(value = phenotype$X_primary_site), length)
detailed_categories <- aggregate(data.frame(count = phenotype$detailed_category), list(value = phenotype$detailed_category), length)
studies <- aggregate(data.frame(count = phenotype$X_study), list(value = phenotype$X_study), length)

genomeXphenotype <- merge(genomic, phenotype, by ="sample")

cleaned_gXp <- genomeXphenotype[!(genomeXphenotype$X_gender == ""), ] #not completely necessary
cleaned_gXp <- cleaned_gXp[!(cleaned_gXp$sampleID == 'Sarcoma'), ]
cleaned_gXp <- cleaned_gXp[(cleaned_gXp$X_study == 'GTEX' | genomeXphenotype$X_study == 'TARGET' | genomeXphenotype$X_study == 'TCGA'),]
cleaned_gXp <- na.omit(cleaned_gXp)

#turn into cleaned
sample_types <- aggregate(data.frame(count = cleaned_gXp$X_sample_type), list(value = cleaned_gXp$X_sample_type), length)
primary_sites <- aggregate(data.frame(count = cleaned_gXp$X_primary_site), list(value = cleaned_gXp$X_primary_site), length)
detailed_categories <- aggregate(data.frame(count = cleaned_gXp$detailed_category), list(value = cleaned_gXp$detailed_category), length)
studies <- aggregate(data.frame(count = cleaned_gXp$X_study), list(value = cleaned_gXp$X_study), length)

#create mean expression levels <- automate later
headers <- colnames(genomeXphenotype)
cleaned_gXp$expression_mean <- rowMeans(cleaned_gXp[,c(headers[2:25])])

#prefilterhist
ggplot(data = cleaned_gXp, aes(x=expression_mean, group=X_study, fill=X_study)) +
  geom_density(adjust=1.5, alpha=.4) + ggtitle("Histogram of Average Expression for Selected Genes \nPrefilter") +
  xlab ("log2(norm_count + 1)") +
  theme_ipsum()

#remove values <4
cleaned_gXp <- cleaned_gXp[cleaned_gXp$expression_mean > 4, ]

#postfilter hist by study
ggplot(data = cleaned_gXp, aes(x=expression_mean, group=X_study, fill=X_study)) +
  geom_density(adjust=1.5, alpha=.4) + ggtitle("Histogram of Average Expression for Selected Genes \nBy Study, Filtered: Avg Exp > 4") +
  xlab ("log2(norm_count + 1)") +
  theme_ipsum()

#-----# 
#ANALYSIS#

#heatmaps

averages_by_detailed_cat <- (aggregate(cleaned_gXp, list(cat = cleaned_gXp$detailed_category), mean))
row.names(averages_by_detailed_cat) <-averages_by_detailed_cat$cat
averages_by_detailed_cat <- averages_by_detailed_cat[,3:26]
averages_by_detailed_cat_sort <- averages_by_detailed_cat[order(averages_by_detailed_cat$PTRH2),]
detailed_cat_matrix <- as.matrix(averages_by_detailed_cat_sort)
detailed_cat_matrix_goi <- detailed_cat_matrix[,1:5]

goi_det_cat_heatmap <- heatmap(detailed_cat_matrix, Colv=NA, Rowv =NA,
                               col = cm.colors(256), scale="column", margins=c(5,10),
                               main = "Relative Expression of All Genes")

goi_det_cat_heatmap <- heatmap(detailed_cat_matrix_goi, Colv=NA, Rowv =NA,
                               col = cm.colors(256), scale="column", margins=c(5,10),
                               main = "Relative Expression of GOI")

##Boxplot processsing
averages_by_sample_types <- aggregate(cbind(MYCN, RET, KIT, PROM1, PTRH2) ~ X_sample_type,
                                      data = cleaned_gXp, FUN = mean, na.rm = TRUE)

filt_cln_gXp <- cleaned_gXp %>%
  group_by(X_sample_type) %>% filter(n() > 15)

GOI_names <- (colnames(averages_by_sample_types)[2:6])

con_names <- (colnames(filt_cln_gXp[7:14]))

prev_names <- (colnames(filt_cln_gXp[15:25]))

goi_boxplots <- list()

con_boxplots <- list()

prev_boxplots <- list()

#----#
#plot boxplots for gene
plot_gene_goi <- function(gene){
  print(ggplot(filt_cln_gXp, aes_(filt_cln_gXp$X_sample_type, as.name(gene))) + 
          geom_boxplot() + coord_flip() + scale_fill_grey() + 
          labs(x = "Sample Type", y = "Expression (log2(norm_count + 1))", title = paste("Average Expression of Experimental", as.character(gene), "\nFor Different Sample Types"))) +
    ylim(0,18)
}

plot_gene_con <- function(gene){
  print(ggplot(filt_cln_gXp, aes_(filt_cln_gXp$X_sample_type, as.name(gene))) + 
          geom_boxplot() + coord_flip() + scale_fill_grey() + 
          labs(x = "Sample Type", y = "Expression (log2(norm_count + 1))", title = paste("Average Expression of House Keeper", as.character(gene), "\nFor Different Sample Types"))) +
    ylim(0,18)
}

plot_gene_prev <- function(gene){
  print(ggplot(filt_cln_gXp, aes_(filt_cln_gXp$X_sample_type, as.name(gene))) + 
          geom_boxplot() + coord_flip() + scale_fill_grey() + 
          labs(x = "Sample Type", y = "Expression (log2(norm_count + 1))", title = paste("Average Expression of Prev.Est.", as.character(gene), "\nFor Different Sample Types"))) +
    ylim(0,18)
}
#----#

for (gene_goi in GOI_names){
  goi_boxplots[[gene_goi]] <- plot_gene_goi(gene = gene_goi)
}

for (gene_con in con_names){
  con_boxplots[[gene_con]] <- plot_gene_con(gene = gene_con)
}

for (gene_prev in prev_names){
  prev_boxplots[[gene_prev]] <- plot_gene_prev(gene = gene_prev)
}

#----#
#testing change to above fx
nb_exp <- subset(cleaned_gXp, detailed_category == detailed_categories[64, 1])
nb_con <- subset(cleaned_gXp , detailed_category == detailed_categories[5,1])
nsclc_exp <- subset(cleaned_gXp, (detailed_category == detailed_categories[57, 1] | detailed_category == detailed_categories[58, 1]))
nsclc_con <- subset(cleaned_gXp, detailed_category == detailed_categories[56, 1])
mel <- subset(cleaned_gXp, (detailed_category == detailed_categories[76, 1] | detailed_category == detailed_categories[89, 1]))
mel_con <- subset(cleaned_gXp, (detailed_category == detailed_categories[74, 1] | detailed_category == detailed_categories[75, 1]))

plot_gene_pheno <- function(gene, phenotype){
  print(ggplot(samples_of_interest, aes_(as.name(phenotype), as.name(gene))) + 
          geom_boxplot() + coord_flip() + scale_fill_grey() + 
          labs(x = "Sample Type", y = "Expression (log2(norm_count + 1))", title = paste("Average Expression of", as.character(gene), "\nFor Detailed Categories"))) +
    ylim(0,18)
}
#----#
#boxplots for specific samples of interest
samples_of_interest <- rbind(mel, mel_con, nb_con, nb_exp, nsclc_con, nsclc_exp)

target_sample_box_goi <- list()

for (gene_goi in GOI_names){
  target_sample_box_goi[[gene_goi]] <- plot_gene_pheno(gene = gene_goi, phenotype = 'detailed_category')
}

#KM survival plots
all_survival_data <- bind_rows(TARGET_NBL_survival, TCGA_survival_data[1:3])
all_survival_data <- subset(all_survival_data, select = c(sample, OS, OS.time))
survival_x_genomic <- merge(cleaned_gXp, all_survival_data, by = "sample" )
survival_x_genomic <- survival_x_genomic %>% filter(OS == 0 | OS == 1)
cats_of_interest <- as.data.frame(unique(samples_of_interest$detailed_category))
survival_x_genomic_filt <- survival_x_genomic %>% filter(detailed_category %in% cats_of_interest$`unique(samples_of_interest$detailed_category)`)
#nb data <- DRILL DOWN

##boxplot fx##
plot_gene_pheno_data <- function(gene, phenotype, dataset = samples_of_interest){
  print(ggplot(dataset, aes_(as.name(phenotype), as.name(gene))) + 
          geom_boxplot() + coord_flip() + scale_fill_grey() + 
          labs(x = "Sample Type", y = "Expression (log2(norm_count + 1))", title = paste("Average Expression of", as.character(gene), "\nFor Detailed Categories"))) +
    ylim(0,18)
}
#nb
nb_exp$TARGET_USI <- substr(nb_exp$sample, 1, 16)
nb_exp_more_data <- merge(nb_exp, TARGET_NBL_ClinicalData, by.x =c('TARGET_USI'),
                          by.y = c('TARGET USI'))  
nb_exp_more_data$`Vital Status`[nb_exp_more_data$`Vital Status` == "Dead"] <- 1
nb_exp_more_data$`Vital Status`[nb_exp_more_data$`Vital Status` == "Alive"] <- 0
nb_exp_more_data$BIT1_Higher_Med <- 12
nb_exp_more_data$BIT1_Higher_Med[nb_exp_more_data$PTRH2 >= median(nb_exp_more_data$PTRH2)] <- 1
nb_exp_more_data$BIT1_Higher_Med[nb_exp_more_data$PTRH2 < median(nb_exp_more_data$PTRH2)] <- 0

nb_by_stage <- list()

for (gene_goi in GOI_names){
  nb_by_stage[[gene_goi]] <- plot_gene_pheno_data(gene = gene_goi, phenotype = 'INSS Stage', dataset = nb_exp_more_data)
}

for (gene_goi in GOI_names){
  nb_by_stage[[gene_goi]] <- plot_gene_pheno_data(gene = gene_goi, phenotype = 'Vital Status', dataset = nb_exp_more_data)
}

for (gene_goi in GOI_names){
  nb_by_stage[[gene_goi]] <- plot_gene_pheno_data(gene = gene_goi, phenotype = 'Ploidy', dataset = nb_exp_more_data)
}
#Survival Analysis NB
km_model_nb_overall <- survfit(Surv(nb_exp_more_data$`Overall Survival Time in Days`,
                                    nb_exp_more_data$`Vital Status`) ~
                                 nb_exp_more_data$BIT1_Higher_Med)

plot(km_model_nb_overall, conf.int = F, xlab = "Time (Days)", ylab = "% Alive", main = "KM - Neuroblastoma",
     col = c("Blue","Red"), las = 1)

legend(3200, .97, legend = c("BIT1 < Median (Low)", "BIT1 > Median (High)"), lty = 1, 
       col = c("Blue", "Red"))

summary(km_model_nb_overall)
km_model_nb_overall

km_model_nb_overall_p <- survdiff(Surv(nb_exp_more_data$`Overall Survival Time in Days`,
                                       nb_exp_more_data$`Vital Status`) ~
                                    nb_exp_more_data$BIT1_Higher_Med)

text(500, .2, "p = 0.004")
text(500, .15, "n = 162")


#Survival Analysis Rest of Cancers

survival_x_genomic_filt$OS <- as.numeric(survival_x_genomic_filt$OS)
survival_x_genomic_filt$BIT1_Higher_Med <- NA

survival_x_genomic_filt$BIT1_Higher_Med[survival_x_genomic_filt$PTRH2 >= median(survival_x_genomic_filt$PTRH2)] <- 1
survival_x_genomic_filt$BIT1_Higher_Med[survival_x_genomic_filt$PTRH2 < median(survival_x_genomic_filt$PTRH2)] <- 0

survival_x_genomic_filt_mel <- survival_x_genomic_filt %>% filter(detailed_category == "Uveal Melanoma" | detailed_category == "Skin Cutaneous Melanoma")
survival_x_genomic_filt_mel$BIT1_Higher_Med[survival_x_genomic_filt_mel$PTRH2 >= median(survival_x_genomic_filt_mel$PTRH2)] <- 1
survival_x_genomic_filt_mel$BIT1_Higher_Med[survival_x_genomic_filt_mel$PTRH2 < median(survival_x_genomic_filt_mel$PTRH2)] <- 0


km_model_mel <- survfit(Surv(survival_x_genomic_filt_mel$OS.time,
                             survival_x_genomic_filt_mel$OS) ~
                          survival_x_genomic_filt_mel$BIT1_Higher_Med)

plot(km_model_mel, conf.int = F, xlab = "Time (Days)", ylab = "% Alive", main = "KM - Cutaneous & Uveal Melanoma",
     col = c("Blue","Red"), las = 1, mark.time =TRUE)

km_model_mel_p <- survdiff(Surv(survival_x_genomic_filt_mel$OS.time,
                                survival_x_genomic_filt_mel$OS) ~
                             survival_x_genomic_filt_mel$BIT1_Higher_Med)

legend(6000, .97, legend = c("BIT1 < Median (Low)", "BIT1 > Median (High)"), lty = 1, 
       col = c("Blue", "Red"))

text(500, .2, "p = 0.02")
text(500, .15, "n = 524")
##
survival_x_genomic_filt_nsclc <- survival_x_genomic_filt %>% filter(detailed_category == "Lung Adenocarcinoma" | detailed_category == "Lung Squamous Cell Carcinoma")
survival_x_genomic_filt_nsclc$BIT1_Higher_Med[survival_x_genomic_filt_nsclc$PTRH2 >= median(survival_x_genomic_filt_nsclc$PTRH2)] <- 1
survival_x_genomic_filt_nsclc$BIT1_Higher_Med[survival_x_genomic_filt_nsclc$PTRH2 < median(survival_x_genomic_filt_nsclc$PTRH2)] <- 0


km_model_nsclc <- survfit(Surv(survival_x_genomic_filt_nsclc$OS.time,
                               survival_x_genomic_filt_nsclc$OS) ~
                            survival_x_genomic_filt_nsclc$BIT1_Higher_Med)

plot(km_model_nsclc, conf.int = F, xlab = "Time (Days)", ylab = "% Alive", main = "KM - Lung Squamous Cell Carcinoma & Adenocarcinoma",
     col = c("Blue","Red"), las = 1, mark.time =TRUE)

km_model_nsclc_p <- survdiff(Surv(survival_x_genomic_filt_nsclc$OS.time,
                                  survival_x_genomic_filt_nsclc$OS) ~
                               survival_x_genomic_filt_nsclc$BIT1_Higher_Med)

legend(4000, .97, legend = c("BIT1 < Median (Low)", "BIT1 > Median (High)"), lty = 1, 
       col = c("Blue", "Red"))

text(500, .2, "p = 0.6")
text(500, .15, "n = 1120")
##NOTE: Both NSCLC and MEL survival datasets were recconmended for use by study: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6066282/

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="/Users/kenminsoo/Desktop/Projects/JABSOM_2022/Dr_Jijiwa/ncc_project_files")

#sclc analysis
sclc_rpkm <- data_mrna_seq_rpkm
rm(data_mrna_seq_rpkm)
sclc_rpkm_goi <- sclc_rpkm %>% filter(sclc_rpkm$Hugo_Symbol %in% GOI_names)
sclc_rpkm_con <- sclc_rpkm %>% filter(sclc_rpkm$Hugo_Symbol %in% con_names)

t_sclc_rpkm_goi <- t(sclc_rpkm_goi)
colnames(t_sclc_rpkm_goi) <- t_sclc_rpkm_goi[1,]
t_sclc_rpkm_goi <- t_sclc_rpkm_goi[!(row.names(t_sclc_rpkm_goi) %in% "Hugo_Symbol"), ]

t_sclc_rpkm_con <- t(sclc_rpkm_con)
colnames(t_sclc_rpkm_con) <- t_sclc_rpkm_con[1,]
t_sclc_rpkm_con <- t_sclc_rpkm_con[!(row.names(t_sclc_rpkm_con) %in% "Hugo_Symbol"), ]

orig_t_sclc_rpkm_con <- t(sclc_rpkm_con)
colnames(orig_t_sclc_rpkm_con) <- orig_t_sclc_rpkm_con[1,]
orig_t_sclc_rpkm_con <- orig_t_sclc_rpkm_con[!(row.names(orig_t_sclc_rpkm_con) %in% "Hugo_Symbol"), ]

t_sclc_rpkm_con <- matrix(as.numeric(t_sclc_rpkm_con), ncol = ncol(t_sclc_rpkm_con))
t_sclc_rpkm_con <- as.data.frame((t_sclc_rpkm_con))
colnames(t_sclc_rpkm_con) <- colnames(orig_t_sclc_rpkm_con)
rownames(t_sclc_rpkm_con) <- rownames(orig_t_sclc_rpkm_con)
t_sclc_rpkm_con$rpkm_sum <- colSums(sclc_rpkm[2:82])

orig_t_sclc_rpkm_goi <- t(sclc_rpkm_goi)
colnames(orig_t_sclc_rpkm_goi) <- orig_t_sclc_rpkm_goi[1,]
orig_t_sclc_rpkm_goi <- orig_t_sclc_rpkm_goi[!(row.names(orig_t_sclc_rpkm_goi) %in% "Hugo_Symbol"), ]

t_sclc_rpkm_goi <- matrix(as.numeric(t_sclc_rpkm_goi), ncol = ncol(t_sclc_rpkm_goi))
t_sclc_rpkm_goi <- as.data.frame((t_sclc_rpkm_goi))
colnames(t_sclc_rpkm_goi) <- colnames(orig_t_sclc_rpkm_goi)
rownames(t_sclc_rpkm_goi) <- rownames(orig_t_sclc_rpkm_goi)
t_sclc_rpkm_goi$rpkm_sum <- colSums(sclc_rpkm[2:82])

#convert to tpm

t_sclc_tpm_con <- (t_sclc_rpkm_con/t_sclc_rpkm_con$rpkm_sum)*10^6
t_sclc_tpm_goi <- (t_sclc_rpkm_goi/t_sclc_rpkm_goi$rpkm_sum)*10^6

t_sclc_log_con <- t_sclc_tpm_con + 1
t_sclc_log_con <- log(t_sclc_log_con, base = exp(2))
t_sclc_log_goi <- t_sclc_tpm_goi + 1
t_sclc_log_goi <- log(t_sclc_log_goi, base = exp(2))

######---log2(norm+1) ---> norm 
##formula: (2^x)-1 = norm
for (gene_goi in GOI_names){
  norm_cleaned_gXp$as.name(gene_goi) <- 2^norm_cleaned_gXp$as.name(gene_goi)
  norm_cleaned_gXp$as.name(gene_goi) <- norm_cleaned_gXp$as.name(gene_goi) - 1
}

norm_cleaned_gXp <- cleaned_gXp
norm_cleaned_gXp$PTRH2 <- 2^norm_cleaned_gXp$PTRH2
norm_cleaned_gXp$PTRH2 <- norm_cleaned_gXp$PTRH2 - 1 

norm_cleaned_gXp[2:25] <- 2^norm_cleaned_gXp[2:24]
norm_cleaned_gXp[2:25] <- norm_cleaned_gXp[2:24] - 1

##find categorized means
averages_by_detailed_cat_norm <- (aggregate(norm_cleaned_gXp, list(cat = norm_cleaned_gXp$detailed_category), mean))
averages_by_detailed_cat_norm_log <- log(averages_by_detailed_cat_norm[3:25], 2)
row.names(averages_by_detailed_cat_norm_log) <- row.names(averages_by_detailed_cat_norm)
averages_by_detailed_cat_norm_log$type <- averages_by_detailed_cat_norm$cat