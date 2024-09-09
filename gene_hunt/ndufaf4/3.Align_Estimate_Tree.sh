# Align freely
mafft --thread 10 --auto Paralogs_G75.fa > NDUFAF4_freealign.fa

# Trim
trimal -automated1 -in NDUFAF4_freealign.fa -out NDUFAF4_freealign_trimalAuto.fa

# Tree
iqtree --redo -keep-ident -T 20 -s NDUFAF4_freealign_trimalAuto.fa --seqtype DNA -m "MFP" -B 1000

# And then identify gaps 
seqkit stats -a *

# Identify the samples with < 70% gaps, and create a tree:
for i in $(cat Keep.list); do cat ${i} >> NDUFAF4_freealign_trimalAuto.fa; done

# Tree
iqtree --redo -keep-ident -T 20 -s NDUFAF4_freealign_trimalAuto.fa --seqtype DNA -m "GTR" -B 1000
