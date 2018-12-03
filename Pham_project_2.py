import numpy as np
import pandas as pd
from Bio import Entrez
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import ClustalwCommandline
from scipy import stats

## import list of DE genes
de_probes = pd.read_csv('./H5N1_VN1203_DE_Probes.txt', delimiter = '\t')
de_genes = list(set(de_probes['gene']))

## import pathway file into dataframe
pathway = pd.read_csv('./pathway.txt', delimiter = '\t')

## extract list of genes in the pathway of Small cell lung cancer and put it in dataframe
gene_list = pathway.loc[pathway.Title == 'Small cell lung cancer', 'Tested Members'].str.split(pat=' ').tolist()[0]
seqdf = pd.DataFrame({'ID':gene_list})

## email for Entrez
email = "phamnh@ohsu.edu"
Entrez.email = email

## list of organisms to search for sequence
organisms = ['Human', 'Mouse', 'Dog']

## search through the nuccore databse to get list of IDs of sequences
for animal in organisms:
    seq_id = []
    for gene in gene_list:
        handle = Entrez.esearch(db="nuccore", term = gene +" AND "+ animal + "[Organism] AND mRNA[Filter] AND Refseq[Filter]")
        record = Entrez.read(handle)
        handle.close()
        ## if the species has a sequence (ID length > 0), add the first ID to the dataframe
        if len(record["IdList"]) > 0:
            seq_id.append(record['IdList'][0])
        ## if the species does not have a sequence, add 0 to the data frame
        else:
            seq_id.append(0)
    seqdf[animal] = seq_id

## remove genes that do not have sequences for all 3 species
seqdf = seqdf.loc[(seqdf!=0).all(axis=1)]

clustalw = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"

## go through each gene, extract the sequences for all 3 species to a SeqFile and run clustalw to align the sequence
consensus_seq = []
consensus_rate = []
for index, row in seqdf.iterrows():
    handle = Entrez.efetch(db="nuccore", id=[row['Human'], row['Mouse'], row['Dog']], rettype="fasta", retmode="text")
    fasta_records = handle.read().strip()
    with open ('./SeqFile', 'w') as fh:
        fh.writelines(fasta_records)
    command = ClustalwCommandline(clustalw, infile='./SeqFile')
    command()
    ## read the aligned sequence into align object
    align = AlignIO.read(".aln", "clustal")
    ## get a gap_consensus sequence, only keeping the nucleotide if all 3 sequences agree (threshold = 1)
    summary = AlignInfo.SummaryInfo(align)
    cseq = summary.gap_consensus(threshold=1, ambiguous='X')
    
    ## Add consensus sequence to list
    consensus_seq.append(str(cseq))
    ## calculate rate of consensus (similarity)
    consensus_rate.append((cseq.count('A') + cseq.count('C') + cseq.count('G') + cseq.count('T'))/len(cseq))
    open('./SeqFile', 'w').close()

## update dataframe
seqdf['Consensus Sequence'] = consensus_seq
seqdf['Similarity'] = consensus_rate

## create flag for DE genes
de_flag = []
## create lists for similarity ratios between DE vs. Non-DE genes
sim_de = []
sim_nde = []

## create DE flag column and list of DE vs. Non-DE genes
for index, row in seqdf.iterrows():
    if row['ID'] in de_genes:
        de_flag.append(1)
        sim_de.append(row['Similarity'])
    else:
        de_flag.append(0)
        sim_nde.append(row['Similarity'])

seqdf['DE Flag'] = de_flag

## draw boxplot to compare similarity ratios between DE genes vs. Non-DE genes
bp = seqdf.boxplot(column='Similarity', by='DE Flag')
bp.set_xticklabels(['Non DE genes','DE genes']);

## Do Mann-Whitney U test
x1 = np.array(sim_de)
x2 = np.array(sim_nde)

u, p = stats.mannwhitneyu(x1, x2)
print("DE genes =", x1)
print("Non DE genes =", x2)
print("U =", u)
print("P-value =", p)