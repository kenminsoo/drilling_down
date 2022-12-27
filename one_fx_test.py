from matplotlib.style import available
import numpy as np
import xenaPython as xena
from cbio_py import cbio_mod as cb

hub = "https://toil.xenahubs.net" #use toil hub with normalized rna-seq data, change later to apply to more data hubs
toil = np.array(xena.all_datasets(hub)) #create array with all datasets

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

    print("You have " + str(num_potential_datasets) + " potential datasets")

    #have the user pick dataset they would like to use
    for dataset in potential_dataset: 
        print("Would you like to include the following dataset?")
        print(dataset['name'])
        condition = ""
        condition = input("Y/N: ")
        print(condition)
        if condition == 'Y':
            datasets_to_use.append(dataset['name'])
        elif condition == 'N':
            print("Skipping")
        else: 
            print("WARNING: Invalid input, skipping dataset")
    return datasets_to_use

def extract_xena_data_s(sample_dict):
    extracted_data = {}
    for dataset in survival_data:
        fields = xena.dataset_field(hub, dataset)
        extracted_data['samples'] = fields #produce "headers"
        temp_data_storage = xena.dataset_fetch(hub, dataset, sample_dict[dataset], fields)
    for dataset in sample_dict:
        iter_num = 0
        for sample in sample_dict[dataset]:
            extracted_data[sample] = []
            for data in temp_data_storage:
                extracted_data[sample] += [data[iter_num]]
            iter_num += 1    
    return extracted_data

survival_data = select_data_xena("survival", "TCGA")

genomic_samples = {}
phenotypic_samples = {}
survival_samples = {}
all_samples = []

