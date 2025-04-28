# Relatedness among C. canorus

Relatedness was assessed with the KING coefficient among C. canorus individuals from all autosomal chromosomes. We then investigated egg phenotype matching among relatives. 

## Remove and Analyze Relatives

Calculate relatedness on the full set 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

#merge autosomes
bcftools concat --file-list Autosomes.list --threads 10 -Oz -o merged_full/autos.vcf.gz
bcftools index --threads 10 merged_full/autos.vcf.gz

#LD prune
~/modules/plink2 --threads 20 --vcf merged_full/autos.vcf.gz --keep Canorus_Samples_n300.list --allow-extra-chr --set-missing-var-ids @:# \
        --rm-dup --indep-pairwise 50 5 0.1 --maf 0.05 --hwe 1e-10 --max-alleles 2 --min-alleles 2 --out merged_full/autos_canorus_LD
#extract, also a vcf
~/modules/plink2 --threads 20 --vcf merged_full/autos.vcf.gz --keep Canorus_Samples_n300.list --allow-extra-chr --set-missing-var-ids @:# \
        --extract merged_full/autos_canorus_LD.prune.in \
        --make-bed --recode vcf bgz --pca --out merged_full/autos_canorus_LD
bcftools index --threads 5 merged_full/autos_canorus_LD.vcf.gz
sed -i 's/chr_//g' merged_full/autos_canorus_LD.bim

#run admixture 
for K in {2..10}; do
echo "Runnign admixture at K: ${K}"
admixture -j7 --cv=5 ../merged_full/autos_canorus_LD.bed ${K} > autos_canorus_LD.log${K}.out
done

#calculated relatedness 
zcat merged_full/autos_canorus_LD.vcf.gz | sed 's/VCFv4.3/VCFv4.2/g' | vcftools --maf 0.05 --vcf - --relatedness2 --out merged_full/autos_canorus_LD

#calculate missingness
zcat merged_full/autos_canorus_LD.vcf.gz | sed 's/VCFv4.3/VCFv4.2/g' | vcftools --maf 0.05 --vcf - --missing-indv --out merged_full/autos_canorus_LD

#also calculate IBS0
plink --vcf merged_full/autos_canorus_LD.vcf.gz --genome full --const-fid --allow-extra-chr
awk '{print $2, $4, $15}' plink.genome > merged_full/autos_canorus_LD.ibs
```

Calculate the number of SNPs between individuals to design maternal / paternal relatives. This will output double the amount of SNPs because it considers both genotypes if it is encoded as diploid.

Strict mapping and missingness, since we want to be incredibly accurate:

```bash
bcftools view chr_MT_All.SNP.DP3.vcf.gz --samples-file Samples.list | bcftools view --min-ac 1 -e 'F_MISSING > 0.05 || MQ < 40' -Oz -o chr_MT.DP3.MM05.MQ40.vcf.gz
```

Count snps: `python 02_Count_mtDNA_Mutations.py --vcf chr_MT.DP3.MM05.MQ40.vcf.gz --out chrMT_Pairwise_Differences.txt` 

```python
import allel
import itertools
import numpy as np
import gzip
import argparse

# Read VCF > genotypeArray
def read_vcf(filename):
    callset = allel.read_vcf(filename)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    return genotypes, callset['samples']

# Compute pairwise comparisons
def pairwise_comparisons(genotypes, samples):
    results = []
    ploidy = genotypes.ploidy  #ploidy 

    for id1, id2 in itertools.combinations(range(len(samples)), 2):
        g1 = genotypes[:, id1]
        g2 = genotypes[:, id2]

        # Filter out missing data
        non_missing = ~(g1.is_missing() | g2.is_missing())
        g1 = g1.compress(non_missing, axis=0)
        g2 = g2.compress(non_missing, axis=0)

        # Count the number of valid sites
        number_sites = non_missing.sum()

        #Count the number of SNPs correctly accounting for ploidy
        number_snps = np.count_nonzero(g1.to_n_alt() != g2.to_n_alt())

        results.append((samples[id1], samples[id2], number_sites, number_snps))
    return results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute pairwise SNP differences from a VCF file.')
    parser.add_argument('--vcf', required=True, help='Input VCF file (gzipped)')
    parser.add_argument('--out', required=True, help='Output file for pairwise comparisons')
    args = parser.parse_args()

    # Read VCF
    genotypes, samples = read_vcf(args.vcf)

    # Perform pairwise comparisons
    comparisons = pairwise_comparisons(genotypes, samples)

    # Output results
    with open(args.out, 'w') as f:
        for comparison in comparisons:
            f.write(f'{comparison[0]} {comparison[1]} {comparison[2]} {comparison[3]}\n')

    print(f'Pairwise comparisons saved to {args.out}')

```

### Plot relatedness pies

Identify relatives, shared phenotypes among them:

```bash
#### Determine relatedness  
setwd('~/EvoBioWolf/CUCKOO_gentes/relatedness/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(scales)

# Analyze LD-pruned relatedness from vcftools --relatedness2, KING coefficient 
rel <- read.table('autos_canorus_LD.relatedness2',header=T) %>% select(1,2,7) %>% as_tibble
names(rel) <- c('ID_A','ID_B','PHI')
rel = rel %>% mutate(PHI = pmax(0, pmin(0.5, PHI)))
rel = rel %>% filter(ID_A != ID_B)

# Values from here for designating relationships https://www.kingrelatedness.com/manual.shtml
fam = rel %>% mutate(Relationship = ifelse(PHI > 0.354, 'First Degree',
                                           ifelse(PHI > 0.177, 'First Degree',
                                                  ifelse(PHI > 0.0884,'Second Degree',
                                                         ifelse(PHI > 0.0442, 'Third Degree',
                                                                ifelse(PHI <= 0.0442,'Unrelated','Unassigned'))))))

fam %>% filter(PHI > 0.354) #usually > 0.354 is MZ twin / duplicate, but since our highest value is 0.393 and most around 0.37, seems more likely they are just first degree

# Add mtDNA differences, calculated externally which identifies pairwise SNPs between samples for mtDNA
mt = read.table('chrMT_Pairwise_Differences.txt')

names(mt) = c('ID_A','ID_B','Sites','SNPs')
# The script counted diploid genotypes as 2 SNPs, so divide by 2
#mt <- mt %>% mutate(SNPs = SNPs / 2)
fam2 = left_join(fam,mt) %>% drop_na(SNPs)

# Remove redundant comparisons, e.g. ID_A vs ID_B or pairwise redundancies
famrm = fam2 %>% 
  rowwise() %>% 
  mutate(pair = sort(c(ID_A,ID_B)) %>% paste(collapse = ",")) %>%
  group_by(pair,Sites,SNPs) %>%
  distinct(pair, .keep_all = T) %>% 
  separate(pair,into=c('ID_A','ID_B'),remove=F,sep=',') %>% ungroup() %>% select(-pair)

# Merge with metadata
md = read_tsv('Full_Metadata.txt')v
md <- read_tsv('../Metadata_Host.txt')
md = md %>% select(c(ID,Egg,GeographicGroup,Sampling_Year,Haplogroup,Sex,Age))
mda = md

# Add an '_A' and '_B' to each metadata field 
names(mda) = paste0(names(mda),'_A')
mdb = md
names(mdb) = paste0(names(mdb),'_B')
fam3 = left_join(famrm,mda) %>% 
  left_join(.,mdb)

# Inspect, look to see which Haplogroups don't align, etc 
fam3 %>% filter(Relationship == 'First Degree' & Haplogroup_A != Haplogroup_B) %>% data.frame
fam3 %>% filter(Relationship == 'First Degree' & GeographicGroup_A != GeographicGroup_B) %>% data.frame
fam3 %>% filter(Relationship == 'First Degree' & Egg_A != Egg_B) %>% data.frame

# Only retain those with known egg morph, which are NESTLINGS, and which are caught within 2 years ( no parent / sibling)
famd = fam3 %>% drop_na(Egg_A) %>% drop_na(Egg_B)
famd = famd %>% filter(GeographicGroup_A == GeographicGroup_B) # Only compare within the same distance clade
famd = famd %>% filter(Age_A == 'Young' & Age_B == 'Young') # Only compare within YOUNG NESTLINGS!
famd %>% count(Age_A,Age_B)
famd = famd %>% filter(abs(Sampling_Year_A - Sampling_Year_B) <= 2) # Only compare within the same sampling year (e.g. SIBLINGS)
famd %>% count(Age_A,Age_B) #n = 2806

# Within first degree relatives, what's the distribution of the # of mtDNA snps? 
mtdna_mutations <- famd %>% filter(Relationship == 'First Degree') %>% 
  ggplot(aes(x=SNPs))+geom_histogram(show.legend = F)+theme_bw()+
  scale_x_continuous(breaks=pretty_breaks(n=14))

png('20250207_mtDNAMutations.png',units='in',res=300,height=3,width=5)
mtdna_mutations
dev.off()

# Unhash this for making the plots with difference mtDNA difference thresholds
#for (snp_count in seq(0,3,1)){

# snp_count = 2 shows a similar number of unrelated paternal / maternal, and first degree more maternal (n=81 compared to 34 paternal), which we expect for cuckoos 
snp_count <- 0

# If mtDNA Haplogrouplotypes are identical, assign as maternal. 
famd = famd %>% mutate(Line = ifelse(SNPs <= snp_count,'Maternal','Paternal'))
famd %>% filter(Relationship == 'First Degree') %>%  count(Line)
# Line         n
#   1 Maternal    56
# 2 Paternal    59

famd %>% filter(Relationship != 'Unrelated') %>% count(Line)
# Line         n
#   1 Maternal   112
# 2 Paternal   163

#simply assign unrelated as paternal, move pie chart in final plot 
famd = famd %>% mutate(Line = ifelse(Relationship == 'Unrelated','Paternal',Line))

# and remove first degree relatives, OR NOT! 
#famd = famd %>% filter(Relationship != 'First Degree')

famd %>% filter(Relationship != 'Unrelated') %>% count(Relationship,Line)
famd %>% filter(Relationship != 'Unrelated') %>% nrow
# Relationship  Line         n
#   1 First Degree  Maternal    56
# 2 First Degree  Paternal    59
# 3 Second Degree Maternal    29
# 4 Second Degree Paternal    38
# 5 Third Degree  Maternal    27
# 6 Third Degree  Paternal    66

nrow(famd)
#[1] 2806 with first degree, or [1] 2948 with no age and year restrictions
#how many unique individuals?
length(unique(c(famd$ID_A,famd$ID_B)))

#function to get matched data
get_matched_data <- function(df) {
  df %>% 
    mutate(
      Year = ifelse(abs(Sampling_Year_A - Sampling_Year_B) <= 2, 'Matched', 'Unmatched'),
      Egg = ifelse(Egg_A == Egg_B, 'Matched', 'Unmatched'),
      Distance = ifelse(GeographicGroup_A == GeographicGroup_B, 'Matched', 'Unmatched'),
      Haplogrouplogroup = ifelse(Haplogroup_A == Haplogroup_B, 'Matched', 'Unmatched')
    ) %>% 
    gather(key = "Variable", value = "Matched", Year, Egg, Distance,Haplogrouplogroup) %>% 
    group_by(Relationship,Line,Variable) %>% 
    count(Matched)
}

#Divide data into frames and store in a list
data_frames <- list(
  fem = famd %>% filter(Sex_A == 'F' & Sex_B == 'F'),
  mal = famd %>% filter(Sex_A == 'M' & Sex_B == 'M'),
  fm = famd %>% filter(Sex_A != Sex_B)
)

#function to each data frame and store results in a list
results_list <- lapply(names(data_frames), function(frame_name) {
  result_df <- get_matched_data(data_frames[[frame_name]])
  result_df$frame <- frame_name
  result_df
})

#combine
mm <- bind_rows(results_list)
mm$Relationship <- factor(mm$Relationship,levels=c('First Degree','Second Degree','Third Degree','Unrelated'))
mm <- mm %>%
  mutate(frame = case_when(
    frame == "fem" ~ "Female-Female",
    frame == "mal" ~ "Male-Male",
    frame == "fm" ~ "Intersexual",
    TRUE ~ frame)) %>% dplyr::rename(Count = n)

#proportions for pie charts
piece = mm %>% select(-frame) %>% group_by(Relationship,Line,Variable,Matched) %>% na.omit %>% summarize(Count = sum(Count)) %>% ungroup %>% group_by(Relationship,Line,Variable) %>% 
  mutate(Total = sum(Count),
         Proportion = Count/Total,
         Percent = paste0(round(Count/Total,3)*100,'% (',Total,')'))
piece = piece %>% filter(!grepl('Year|Distance|Haplogroup',Variable))

#plot
relpie = piece %>% 
  ggplot(aes(x="",y=Proportion,fill=Matched))+
  geom_bar(stat='identity')+
  coord_polar("y", start=0)+xlab('')+ylab('')+
  facet_grid(Variable+Line~Relationship)+
  scale_fill_manual(values=c('forestgreen','grey95','black'))+
  theme_bw(base_size=7)+labs(fill='Phenotype Matching')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
#add labels showing the % matched, and then the total sample size 
labs = piece %>% 
  filter(Matched == 'Matched')
piesR = relpie + 
  geom_label(data = labs, 
             aes(y=Inf,x=-Inf,label=Percent),fill='white',
             size=2,vjust=.3,hjust=.5,alpha=0.8)
piesR

pdf(paste0('20250207-Relatives_Phenotype_Matching-',snp_count,'Mismatch-FIRSTDEGREE.pdf'),height=4.5,width=6)
print(piesR)
dev.off()

# png(paste0('20250207-Relatives_Phenotype_Matching-',snp_count,'Mismatch-FIRSTDEGREE.png'),units='in',res=300,height=4.5,width=6)
# print(piesR)
# dev.off()

#} # Unhash this for doing the loop based on mtDNA differences 

famd %>% filter(Relationship != 'Unrelated' & Line == 'Paternal') %>%
  count(Relationship,Egg_A,Egg_B)


famd %>% filter(Relationship != 'Unrelated' & Line == 'Maternal') %>% 
  count(Relationship,Egg_A,Egg_B)


famd %>% count(Relationship,Egg_A,Egg_B,Line)

## Aggregate
df_agg <- famd %>% 
  count(Relationship,Egg_A,Egg_B,Line) %>% 
  mutate(
    Egg = pmin(Egg_A, Egg_B),  # ensures consistent grouping
    Match = ifelse(Egg_A == Egg_B, n, 0),
    Unmatch = ifelse(Egg_A != Egg_B, n, 0)
  ) %>%
  group_by(Relationship, Egg, Line) %>%
  summarise(
    Match = sum(Match),
    Unmatch = sum(Unmatch),
    .groups = "drop"
  )

# Supplementary table 
write.table(df_agg,file='20250207_EggMismatch_Counts.txt',quote=F,sep='\t',row.names=F)

# Plot matched 
counts_mismatch <- df_agg %>% 
  filter(Relationship != 'Unrelated') %>% 
  pivot_longer(c(Match,Unmatch)) %>% 
  ggplot(aes(y=Egg,x=value,fill=name))+
  xlab('Count')+
  scale_fill_manual('Matched',values=c('forestgreen','grey95'))+
  geom_bar(col='black',stat='identity',position=position_dodge(width=0.5),lwd=0.25)+
  facet_grid(Relationship~Line,scales='free',space='free_y')+
  theme_bw(base_size=8)+
  scale_x_continuous(breaks = pretty_breaks(n = 3))
counts_mismatch

pdf('~/symlinks/host/figures/20250404-Relatives_EggMismatch.pdf',height=4.5,width=3)
counts_mismatch
dev.off()

# Perform Fisher's test for each Relationship category
results <- df_agg %>%
  filter(Relationship != 'Unrelated') %>% 
  group_by(Relationship) %>%
  summarise(
    p_value = list(
      fisher.test(matrix(c(sum(Match[Line == "Maternal"]),
                           sum(Unmatch[Line == "Maternal"]),
                           sum(Match[Line == "Paternal"]),
                           sum(Unmatch[Line == "Paternal"])),
                         nrow = 2))$p.value
    )
  ) %>%
  unnest(p_value) %>% 
  mutate(padj = p.adjust(p_value,method='bonferroni'))
options(scipen = 9999)
print(results)
# Relationship  p_value   padj
#   1 First Degree  1       1     
# 2 Second Degree 0.502   1     
# 3 Third Degree  0.00460 0.0138


#also add missingness, we will remove individual with more missing data if both are females or both are males 
miss <- read.table('autos_canorus_LD.imiss',head=T) %>% select(c(1,5))
names(miss) = c('ID_A','Missing_A')
miss$Missing_A = as.numeric(miss$Missing_A)
miss = miss %>% na.omit
fam4 = fam3 %>% left_join(.,miss) %>% 
  merge(.,miss %>% dplyr::rename(ID_B=ID_A,Missing_B=Missing_A)) 

#Remove unrelateds
reli = fam4 %>% filter(Relationship != 'Unrelated')
relatives = fam4

#Unique samples to test 
samples <- unique(c(reli$ID_A,reli$ID_B))
rm <- NULL 
choice <- reli
#Re-run this multiple times, until it stops removing any individuals. 
for (run in seq(1,5,1)) { 
  
  cat('Running ',run,'\n')
  #grab one sample at a time 
  for (samp in samples) {
    #grab all the records with this sample 
    sf.a <- choice[grepl(samp,choice$ID_A),] 
    sf.b <- choice[grepl(samp,choice$ID_B),] 
    sf <- rbind(sf.a,sf.b) %>% arrange(desc(PHI))
    if(nrow(sf) < 1) next
    #if one if male remove, otherwise, if sampled clade is higher, remove, otherwise if missing data is higher, remove, otherwise random (sample by ID# higher#)
    sf1 <- sf %>% mutate(Remove = ifelse(Sex_A == 'F' & Sex_B == 'M', ID_B,
                                         ifelse(Sex_A == 'M' & Sex_B == 'F', ID_A,
                                                ifelse(Missing_A > Missing_B, ID_A,
                                                       ifelse(Missing_A < Missing_B, ID_B,
                                                              'Problem')))))
    sf1 %>% select(contains(c('ID')))
    #take the ID from the sample to remove
    bad <- head(na.omit(sf1$Remove),1) 
    #unless there are no samples to remove.. 
    if(length(bad) < 1) next
    #remove that bad sample from the pool
    cat('Removing bad sample: ',bad,'\n')
    choice <- choice[!grepl(bad,choice$ID_A),]
    choice <- choice[!grepl(bad,choice$ID_B),]
    #and then restart the whole process 
    rm <- rbind(rm,bad)
    rm(bad)
  }
}

#these are the bad individuals we will remove
rms <- as.data.frame(rm)
names(rms) <- 'Remove'
row.names(rms) <- NULL

#remove those bad IDs from the full dataset, make sure to filter form both ID_A and ID_B columns 
keep <- relatives %>% filter(!ID_A %in% rms$Remove) %>% filter(!ID_B %in% rms$Remove)
#see how many individuals were in the full dataset (sanity check), and then how many are retained after filtering 
length(unique(c(relatives$ID_A,relatives$ID_B)))
length(unique(c(keep$ID_A,keep$ID_B)))
length(unique(c(relatives$ID_A,relatives$ID_B))) - length(unique(rms$Remove))

write.table(unique(c(keep$ID_A,keep$ID_B)),'Unrelated_2023OCT27.list',quote=F,row.names=F,sep='/t',col.names=F)
samps = read_tsv('Unrelated_2023OCT27.list',col_names = F)
```

