import numpy as np
import pandas as pd

## import DE_probes file into a dataframe
de_probes = pd.read_csv('./H5N1_VN1203_DE_Probes.txt', delimiter = '\t')
## create a list of DE gene, eliminate duplicate genes
de_genes = list(set(de_probes['gene']))

## import UNIVERSE_probes file into a dataframe
all_probes = pd.read_csv('./H5N1_VN1203_UNIVERSE_Probes.txt', delimiter = '\t')
## create a list of all genes, eliminate duplicate genes
all_genes = list(set(all_probes['gene']))

## import KEGG_Pathway_Genes file into a 1-column dataframe
pathway = pd.read_csv('./KEGG_Pathway_Genes.txt', squeeze=True)
## split the dataframe into 3 columns: ID, Title, Members (Members = string of all pathway genes)
pathway = pathway.str.split('\t', 2, expand=True)
pathway.columns = ['ID','Title', 'Members']
## replace tabs between pathway genes with space
pathway['Members'] = pathway['Members'].str.replace('\t',' ')

## add space before and after the string of pathway genes - useful for searching through the list of genes
pathway['Members'] = ' ' + pathway['Members'] + ' '
## create a string of genes in the pathway that are tested (in the list of all genes) 
pathway['Tested Members'] = ''
for gene in all_genes:
    pathway['Tested Members'] += np.where(pathway['Members'].str.contains(' ' + gene + ' '), gene + ' ', '')

## count number of genes that are tested in each pathway
pathway['Tested Members'] = pathway['Tested Members'].str.strip()
pathway['Mem count'] = pathway['Tested Members'].str.split().apply(len)

## count number of genes that are in differentially expressed (in the DE list)
pathway['PW DE count'] = 0
for gene in de_genes:
    pathway['PW DE count'] += np.where(pathway['Members'].str.contains(' ' + gene + ' '), 1, 0)

## calculate the rest of the gene counts to get odds ratios later
non_de_count = len(all_genes) - len(de_genes)
pathway['PW NDE count'] = pathway['Mem count'] - pathway['PW DE count']
pathway['NPW DE count'] = len(de_genes) - pathway['PW DE count']
pathway['NPW NDE count'] = non_de_count - pathway['PW NDE count']

## calculate odds ratios
pathway['Odds Ratio'] = (pathway['PW DE count']*(pathway['NPW NDE count']))/(pathway['PW NDE count']*(pathway['NPW DE count']))
## extract ID, Title, Odds Ratios only
final_pw = pathway.loc[:,['ID','Title','Odds Ratio', 'Tested Members']]

## save dataframe to txt file
final_pw.to_csv(path_or_buf = './pathway.txt', sep='\t', index = False, line_terminator='\n')