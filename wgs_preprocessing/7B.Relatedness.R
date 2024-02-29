#### Determine relatedness with PHI statistic 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/merged_full')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(adegenet)
library(vcfR)
library(ggpubr)

#for brevity analyze these chroms
rel <- read.table('autos_canorus_LD.relatedness2',header=T) %>% select(1,2,7) %>% as_tibble
names(rel) <- c('ID_A','ID_B','PHI')
rel = rel %>% mutate(PHI = pmax(0, pmin(0.5, PHI)))
rel = rel %>% filter(ID_A != ID_B)

#compare with plink IBS0
pk = read.table('autos_canorus_LD.ibs',sep=' ',header=TRUE); names(pk) = c('ID_A','ID_B','IBS0')
relpk = left_join(pk,rel)

#values from here https://www.kingrelatedness.com/manual.shtml
fam = relpk %>% mutate(Relationship = ifelse(PHI > 0.354, 'First Degree',
                                             ifelse(PHI > 0.177, 'First Degree',
                                                    ifelse(PHI > 0.0884,'Second Degree',
                                                           ifelse(PHI > 0.0442, 'Third Degree',
                                                                  ifelse(PHI <= 0.0442,'Unrelated','Unassigned'))))))
fam %>% filter(PHI > 0.354) #usually > 0.354 is MZ twin / duplicate, but since our highest value is 0.393 and most around 0.37, seems more likely they are just first degree
phi_ibs = fam %>% ggplot(aes(x=PHI,y=IBS0,col=Relationship))+
  geom_point()+
  scale_color_viridis(discrete=TRUE)+
  theme_bw()
pdf('../figures/Relatedness__PhivIBS_2023OCT27.pdf',height=5,width=6)
phi_ibs
dev.off()

#add mtDNA differences
mt = read.table('chrMT_Pairwise_Differences.txt')
names(mt) = c('ID_A','ID_B','Sites','SNPs')
fam2 = left_join(fam,mt) %>% drop_na(SNPs)

#remove redundant comparisons
famrm = fam2 %>% 
  select(-IBS0) %>% 
  rowwise() %>% 
  mutate(pair = sort(c(ID_A,ID_B)) %>% paste(collapse = ",")) %>%
  group_by(pair,Sites,SNPs) %>%
  distinct(pair, .keep_all = T) %>% 
  separate(pair,into=c('ID_A','ID_B'),remove=F,sep=',') %>% ungroup() %>% select(-pair)

#remove unrelated individuals ... first merge with metadata
md = read_tsv('../Cuckoo_Full_Metadata_2023OCT3.txt')
md = md %>% select(c(ID,HostParentShort,Habitat,BNB,Sampling_Year,KDist,Hap,Sex))
mda = md
names(mda) = paste0(names(mda),'_A')
mdb = md
names(mdb) = paste0(names(mdb),'_B')
fam3 = left_join(famrm,mda) %>% 
  left_join(.,mdb)

#inspect 
fam3 %>% filter(Relationship == 'First Degree' & Hap_A != Hap_B) %>% data.frame
fam3 %>% filter(Relationship == 'First Degree' & KDist_A != KDist_B) %>% data.frame
fam3 %>% filter(Relationship == 'First Degree' & BNB_A != BNB_B) %>% data.frame

famd = fam3 %>% filter(KDist_A == KDist_B) #only compare within the same distance clade
famd = famd %>% filter(abs(Sampling_Year_A - Sampling_Year_B) <= 2) #only compare within the same sampling year (e.g. SIBLINGS)
famd = famd %>% drop_na(Habitat_A) %>% drop_na(Habitat_B)
famd %>% filter(Relationship == 'First Degree') %>% ggplot(aes(x=SNPs,fill=ID_A))+geom_histogram(show.legend = F)+theme_bw()
famd = famd %>% mutate(Line = ifelse(SNPs <= 3,'Maternal','Paternal'))
famd %>% filter(Relationship == 'First Degree') %>%  count(Line)
# famd = famd %>% mutate(Relationship = ifelse(Relationship == 'First Degree',paste0(Relationship,' ',Line),
#                                              ifelse(Relationship == 'Second Degree',paste0(Relationship,' ',Line),Relationship)))
famd = famd %>% mutate(Line = ifelse(Relationship == 'Unrelated','Paternal',Line))
famd %>% filter(Relationship == 'Second Degree') %>% filter(Line == 'Maternal') %>% filter(HostParentShort_A != HostParentShort_B) %>% data.frame

#function to get matched data
get_matched_data <- function(df) {
  df %>% 
    mutate(
      Host = ifelse(HostParentShort_A == HostParentShort_B, 'Matched', 'Unmatched'),
      Habitat = ifelse(Habitat_A == Habitat_B, 'Matched', 'Unmatched'),
      Year = ifelse(abs(Sampling_Year_A - Sampling_Year_B) <= 2, 'Matched', 'Unmatched'),
      BNB = ifelse(BNB_A == BNB_B, 'Matched', 'Unmatched'),
      Distance = ifelse(KDist_A == KDist_B, 'Matched', 'Unmatched'),
      Haplogroup = ifelse(Hap_A == Hap_B, 'Matched', 'Unmatched')
    ) %>% 
    gather(key = "Variable", value = "Matched", Host, Habitat, Year, BNB, Distance,Haplogroup) %>% 
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

#proportions for pies
piece = mm %>% group_by(Relationship,Variable,frame) %>% na.omit %>% mutate(Total = sum(Count),
                                                                            Proportion = Count/Total,
                                                                            Percent = paste0(round(Count/Total,3)*100,'% (',Total,')'))

#or if you want to merge them 
piece = mm %>% select(-frame) %>% group_by(Relationship,Line,Variable,Matched) %>% na.omit %>% summarize(Count = sum(Count)) %>% ungroup %>% group_by(Relationship,Line,Variable) %>% 
  mutate(Total = sum(Count),
         Proportion = Count/Total,
         Percent = paste0(round(Count/Total,3)*100,'% (',Total,')'))
piece = piece %>% mutate(Variable = gsub('Host','Host Parent',Variable),
                         Variable = gsub('BNB','Is Blue',Variable)) %>% 
  filter(!grepl('Year|Blue|Distance|Hap',Variable))

#plot
relpie = piece %>% 
  #filter(grepl('Distance|Microhabitat|Year|Haplogroup',Variable)) %>% ##supplement 
  #filter(grepl('Inversion',Variable)) %>% ##supplement 
  #filter(frame == 'Female-Female') %>% 
  ggplot(aes(x="",y=Proportion,fill=Matched))+
  geom_bar(stat='identity')+
  coord_polar("y", start=0)+xlab('')+ylab('')+
  facet_grid(Relationship~Variable+Line)+
  scale_fill_manual(values=c('forestgreen','grey95','black'))+
  theme_bw()+labs(fill='Phenotype Matching')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
labs = piece %>% 
  #filter(grepl('Distance|Microhabitat|Year|Haplogroup',Variable)) %>%  #supplement
  #filter(grepl('Inversion',Variable)) %>%  #supplement
  #filter(grepl('Egg|Haplogroup',Variable)) %>% 
  filter(Matched == 'Matched')
piesR = relpie + 
  geom_label(data = labs, 
             aes(y=Inf,x=-Inf,label=Percent),fill='white',
             size=3,vjust=.3,hjust=.5,alpha=0.8)
piesR

pdf('../figures/Relatives_Phenotype_Matching-Sex-Supplement__2023NOV2.pdf',height=8,width=14)
pdf('../figures/Relatives_Phenotype_Matching-PAT-MAT__2023NOV7.pdf',height=7,width=6)
piesR
dev.off()

#also add missingness
miss <- read.table('autos_canorus_LD.imiss',head=T) %>% select(c(1,5))
names(miss) = c('ID_A','Missing_A')
miss$Missing_A = as.numeric(miss$Missing_A)
miss = miss %>% na.omit
fam4 = fam3 %>% left_join(.,miss) %>% 
  merge(.,miss %>% dplyr::rename(ID_B=ID_A,Missing_B=Missing_A)) 

reli = fam4 %>% filter(Relationship != 'Unrelated')
relatives = fam4
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

#calculate pairwise geographic distance
md = read_tsv('../Cuckoo_Full_Metadata_2023OCT3.txt')
ll = md %>% select(ID,Latitude,Longitude)
ll1 = geosphere::distm(ll %>% select(Longitude,Latitude)) %>% as.data.frame
names(ll1) = ll$ID
ll1$ID = ll$ID
ll2 = ll1 %>% pivot_longer(!(ID),names_to='ID_B',values_to = 'GDistance') %>% dplyr::rename(ID_A=ID)
ll2 = ll2 %>% mutate(GDistance = GDistance/1000)

#Add genetic distance
dp = left_join(fam4,ll2)