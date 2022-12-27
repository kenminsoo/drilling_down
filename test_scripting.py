#import packages
from matplotlib.style import available
import numpy as np
import xenaPython as xena
import json
from cbio_py import cbio_mod as cb

#this is an semi-artifical program meant to eventually act as a back end for an easy to use GUI interface

hub = "https://toil.xenahubs.net" #use toil hub with normalized rna-seq data, change later to apply to more data hubs
toil = np.array(xena.all_datasets(hub)) #create array with all datasets

cbioportal = cb.getAllStudies(return_type = 'dict') #pause for now
datatype = 'tpm' #input datat type

#query storage
potential_genomic_dataset = []
genomic_datasets_inuse = []

#genes
genes_of_interest = ['MYCN', 'RET', 'KIT', "PROM1", 'PTRH2']
house_keepers = ['GAPDH', 'RPL30', 'PSMC5', 'WDTC1', 'RRP1', 'HNRNPL', 'PCBP1', 'RER1']
prev_est_GOI = ['ZFHX3', 'CACNA2D3', 'CDK1', 'CCNB1', 'EGFR', 'ROS1', 'ALK', 'ERBB2', 'FLG', 'DSG3']
all_genes = [genes_of_interest, house_keepers, prev_est_GOI]


gene_probe_map = 'gencode.v23.annotation.gene.probemap'
transcript_probe_map = 'probeMap/gencode.v23.annotation.transcript.probemap'

test_sample = 'GTEX-S4Q7-0003-SM-3NM8M'


hub = "https://toil.xenahubs.net"
dataset = 'TcgaTargetGtex_RSEM_Hugo_norm_count'
pheno_dataset = 'TcgaTargetGTEX_phenotype.txt'

genomic_samples = xena.dataset_samples(hub, 'TcgaTargetGtex_RSEM_Hugo_norm_count', 2)

test_data = xena.dataset_probe_values(hub, dataset, genomic_samples, ['MYCN', 'RET'])

pheno_fields = xena.dataset_field(hub, pheno_dataset)
field_codes = xena.field_codes(hub, pheno_dataset, pheno_fields[0:2])

test_data_2 = [genomic_samples] + xena.dataset_fetch(hub, dataset, genomic_samples, ['MYCN', 'RET']) + xena.dataset_fetch(hub, pheno_dataset, genomic_samples, pheno_fields[0:2])

field_names = []
n = 0

dictionary_storage_2 = {}

for field in field_codes:
    field_names += [field['name']]

for field in field_names:
    dictionary_storage_2[field] = field_codes[n]['code'].split('\t')
    n += 1

iter_num = 0

dictionary_storage = {}

for sample in test_data_2[0]:
    sample = []
    for data in test_data_2:
        sample += [data[iter_num]]
    iter_num += 1

iter_num = 0

for sample in test_data_2[0]:
    dictionary_storage[sample] = []
    for data in test_data_2:
        dictionary_storage[sample] += [data[iter_num]]
    iter_num += 1

dictionary_storage_copy = dictionary_storage.copy()

print(test_data_2)
print(dictionary_storage)
print(dictionary_storage_2)
print(pheno_fields)
for sample in dictionary_storage:
    num_gene = len(['MYCN', 'RET']) + 1
    for field in pheno_fields[0:2]:
        key = dictionary_storage[sample][num_gene]
        dictionary_storage_copy[sample][num_gene] = dictionary_storage_2[field][key]
        print(field)
        print(key)
        print(dictionary_storage_2[field][key])
        num_gene += 1


print(dictionary_storage_copy)
    