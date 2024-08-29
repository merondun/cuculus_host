cat *fa > All_Except_064.fa

# Align and keep the coordinates of the chr3 autosome 
mafft --thread 10 --auto --addfragments All_Except_064.fa --keeplength --reorder 064_CC_GRW_BGR_M__SRR11531726.fasta > NDUFAF4_chr3anchor.fa

# Align freely
mafft --thread 10 --auto All.fa > NDUFAF4_freealign.fa

# With python, split the MSA
from Bio import SeqIO

# Path to your MSA file in FASTA format
msa_file = 'NDUFAF4_chr3anchor.fa'

# Read the MSA file and create individual FASTA files
with open(msa_file, 'r') as msa:
    for record in SeqIO.parse(msa, 'fasta'):
        filename = f"{record.id}.fa"
        with open(filename, 'w') as output_file:
            SeqIO.write(record, output_file, 'fasta')

# And then identify gaps 
seqkit stats -a *

# Identify the samples with < 50% gaps, and create a tree:
for i in $(cat Keep.list); do cat ${i} >> NDUFAF4_Gapless50.fa; done

# Tree
iqtree --redo -keep-ident -T 20 -s NDUFAF4_Gapless50.fa --seqtype DNA -m "GTR" -B 1000