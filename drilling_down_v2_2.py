#import packages
from cbio_py import cbio_mod as cb
from select import select
from matplotlib.style import available
import numpy as np
import pandas as pd
import xenaPython as xena
import json
import csv
#---#
##version 2.2: 
#New Features to Implement
#########

#Gene alias search
#Translate genes to universal location code
#Addition of survival data
#Combining multiple pieces of data
#User input for queries
#Easier to work with data classest that are less fragile
#Testing edge cases where datasets may not be picked

#########
#---#
#this is an semi-artifical program meant to eventually act as a back end for an easy to use GUI interface
#this uses UCSC's XENA API... Thank you Banana Slugs. (will add official citation later...)

#------------------------------#

#DEFINE FUNCTIONS
##___##
#user to select data theyd like to download
def select_data_xena(para1, para2):
    #query storage for datasets
    potential_dataset = []
    datasets_to_use = []
    num_potential_datasets = 0
    query_para2 = para2 #input data type (in this case normalized) ***********
    
    #add relevant datasets
    for dataset in toil: 
        if para1.lower() in dataset['name'].lower():
            if query_para2.lower() in dataset['name'].lower():
                potential_dataset.append(dataset)
                num_potential_datasets += 1
        else:
            continue

    print("You have " + str(num_potential_datasets) + " potential genomic datasets")

    if num_potential_datasets == 0:
        return print("no datasets found")

    #have the user pick dataset they would like to use
    for dataset in potential_dataset: 
        print("Would you like to include the following dataset?")
        print(dataset['name'])
        condition = input("Y/N: ")
    if condition == 'Y' or 'y':
        datasets_to_use.append(dataset['name'])
    elif condition == 'N' or 'n':
        print("Skipping")
    else: 
        print("WARNING: Invalid input, skipping dataset")
    return datasets_to_use

#extract genomic data from the API
#data: [sample: gene1exp, gene2exp, gene3exp, ...]
def extract_xena_data_g(sample_dict, targets_gene):#i.e. samp or phenotypic
    extracted_data = {}
    extracted_data['sample'] = targets_gene
    for dataset in sample_dict: #extracts data as: [x1, x2, x3,...], [y1, y2, y3,...], ...
        temp_data_storage = xena.dataset_fetch(hub, dataset, sample_dict[dataset], targets_gene)
    for dataset in sample_dict:
        iter_num = 0
        for sample in sample_dict[dataset]:
            extracted_data[sample] = []
            for data in temp_data_storage:
                extracted_data[sample] += [data[iter_num]]
            iter_num += 1
    return extracted_data

#extract phenotypic data from the API, I don't think this scales to multiple clinical matracies implement
#data: [sample: pheno1, pheno2, pheno3, ...]
#necessary to have multiple because target phenotypes may differ between datasets while genes stay the same
def extract_xena_data_p(sample_dict, targets_pheno):#i.e. samp or phenotypic
    extracted_data = {}
    for dataset in sample_dict: #extracts data as: [x1, x2, x3,...], [y1, y2, y3,...], ...
        temp_data_storage = xena.dataset_fetch(hub, dataset, sample_dict[dataset], targets_pheno[dataset])
        extracted_data['sample'] = targets_pheno[dataset] #produce "headers"
    for dataset in sample_dict:
        iter_num = 0
        for sample in sample_dict[dataset]:
            extracted_data[sample] = []
            for data in temp_data_storage:
                extracted_data[sample] += [data[iter_num]]
            iter_num += 1    
    return extracted_data

def extract_xena_data_s(sample_dict):
    extracted_data = {}
    for dataset in survival_data:
        temp_data_storage = xena. 

#gives the user a choice to keep the three category of genes separate
#combines datasets based upon sample tag
def combine_dicts(_of_dict_to_combine): 
    combined_data = {}
    combined_data['sample'] = []
    #first, give headers
    for average_dict in _of_dict_to_combine: 
        combined_data['sample'] += average_dict['sample']
    #create empty list for each sample in the dictionary
    for sample in all_samples:
        combined_data[sample] = []
    #now, add in data from each dict based on key
    for average_dict in _of_dict_to_combine:
        average_dict_keys = list(average_dict.keys())
        for sample in all_samples:
            if sample in average_dict_keys:
                combined_data[sample] += average_dict[sample]
            else:
                print("WARNING: Skipping sample: {sample_txt}. Not in dataset.".format(sample_txt = sample))
    return combined_data

#tests if a certain value is a number
def is_num(n):
    if isinstance(n, int):
        return True
    if isinstance(n, float):
        return True
    return False

#decodes phenotypes, will eventually implement greater support for more clinical matracies
def xena_phenotype_decode(dictionary):
    decoded_dictionary = dictionary.copy()
    for dataset in phenotypic_fields:
        for sample in dictionary:
            n = 0
            for field in phenotypic_fields[dataset]:
                key = dictionary[sample][n]
                if is_num(key):
                    if key + 1 > len(phenotypic_codes[dataset][field]):
                        print("ERROR: Skipping. Info:")
                        print(key)
                        print(field)
                        print(sample)
                        input()
                        n += 1
                    else: 
                        decoded_dictionary[sample][n] = phenotypic_codes[dataset][field][key]
                        n += 1
                else:
                    n += 1
                    continue
    return decoded_dictionary

#turn dictionary to csv

exported_csv_dicts = []
def dict_to_csv(dictionary, name):
    with open('data/{file_name}.csv'.format(file_name = name), 'w') as f:
        n = 0
        for key in dictionary.keys():
            if n != 0:
                f.write("\n%s"%(key))
            if n == 0:
                f.write("%s"%(key))
                n += 1
            for item in dictionary[key]:
                f.write(",%s"%(item))
    exported_csv_dicts.append('data/{file_name}.csv'.format(file_name = name))

#create txt file with csv values for R to read
def finished_dict_export():
    with open('data/exported_files.txt', 'w') as f:
        for file_name in exported_csv_dicts:
            f.write(file_name)
            f.write(",")
#-----------------------------------------#

#pick datatype from xenahub
hub = "https://toil.xenahubs.net" #use toil hub with normalized rna-seq data, change later to apply to more data hubs
toil = np.array(xena.all_datasets(hub)) #create array with all datasets
cohort = '' #how can I apply cohorts?

#implement other API & eventual normalization
cbioportal = cb.getAllStudies(return_type = 'dict') #pause for now

#gene storage, input
genes_of_interest = ['MYCN', 'RET', 'KIT', "PROM1", 'PTRH2']
house_keepers = ['GAPDH', 'RPL30', 'PSMC5', 'WDTC1', 'RRP1', 'HNRNPL', 'PCBP1', 'RER1', 'A2ML1']
prev_est_GOI = ['ZFHX3', 'CACNA2D3', 'CDK1', 'CCNB1', 'EGFR', 'ROS1', 'ALK', 'ERBB2', 'FLG', 'DSG3']
all_genes = genes_of_interest + house_keepers + prev_est_GOI

#convert genes to their equivalent


#check genes exist/valid
gene_dataset = 'TcgaTargetGtex_RSEM_Hugo_norm_count'
is_genes = xena.dataset_field(hub, gene_dataset)
for list_1 in all_genes:
    for gene in list_1:
        if gene in is_genes:
            continue
        else:
            print("ERROR: Gene not found. Searching for aliases...") #implement at a later date
            print("Alias not found, omitting.") 
            list_1.remove(gene)

#convert genes to their ID
#implement at a later date
gene_probe_map = 'probeMap/gencode.v23.annotation.gene.probemap'
transcript_probe_map = 'probeMap/gencode.v23.annotation.transcript.probemap'

#query for genomic datasets
genomic_data = select_data_xena("N/A", 'N/A') #so for example, files with TcgaTargetGtex and norm

#query storage for phenotypic datasets
phenotypic_data = select_data_xena("N/A", "N/A")

#query storage for survival data (if applicable)
survival_data = select_data_xena("survival", "TCGA")

#extract phenotype fields into dictionary
phenotypic_fields = {}
try:
    for dataset in phenotypic_data:
        phenotypic_fields[dataset] = xena.dataset_field(hub, dataset)
except:
    print("No phenotypic fields found")

#extract phenotype codes into dictionary
phenotypic_codes_raw = {}
try:
    for dataset in phenotypic_data:
        phenotypic_codes_raw[dataset] = xena.field_codes(hub, dataset, phenotypic_fields[dataset])
except:
    print("No phenotypic codes")
phenotypic_codes = {}

#implement - generalize further to apply to other clinical matracies, i.e. if clinical matracies not str, search for decoder... etc.
try:
    for dataset in phenotypic_data:
        ph_n = 0
        phenotypic_codes[dataset] = {}
        for field in phenotypic_fields[dataset]:
            name = phenotypic_codes_raw[dataset][ph_n]['name']
            phenotypic_codes[dataset][name] = phenotypic_codes_raw[dataset][ph_n]['code'].split('\t')
            ph_n += 1
except:
    print("No Phenotypic data found")
#store samples for data, implement a way to add samples from multiple datasets (I think it can now)
genomic_samples = {}
phenotypic_samples = {}
all_samples = []
try:
    for dataset in genomic_data:
        genomic_samples[dataset] = xena.dataset_samples(hub, dataset, None)
        for sample in genomic_samples[dataset]:
            if sample not in all_samples:
                all_samples.append(sample)
            else:
                continue
except:
    print("No genomic data")

try:
    for dataset in phenotypic_data: #usually genomic data has a clinical matrix attached to it. automate finding it?
        phenotypic_samples[dataset] = xena.dataset_samples(hub, dataset, None)
        for sample in phenotypic_samples[dataset]:
            if sample not in all_samples:
                all_samples.append(sample)
            else:
                continue    
except:
    print("No phenotypic data")

#----#

phenotype = extract_xena_data_p(phenotypic_samples, phenotypic_fields) #phenotype data

#collected data from the above code
house_keepers = extract_xena_data_g(genomic_samples, house_keepers) #house keeper genes, separate data
goi = extract_xena_data_g(genomic_samples, genes_of_interest) #genes of interests, separate data
prev_est = extract_xena_data_g(genomic_samples, prev_est_GOI) #previously established DEGs to "test"

exp_gene_data = [goi, house_keepers,prev_est] #set up for fx

goi_hk_prev_genomic_data = combine_dicts(exp_gene_data) #combines all genomic data, allows one to eventually convert data from other APIs
decoded_phenotype = xena_phenotype_decode(phenotype) #decoded phenotypic data