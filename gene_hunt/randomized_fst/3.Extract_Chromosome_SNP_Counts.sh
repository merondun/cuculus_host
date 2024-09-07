fstdat=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/randomizing_fst/fst/out

mkdir -p counts

for CHR in $(cat Chromosomes.list); do 

    cat $fstdat/${CHR}_rf10.fst.txt | wc -l | awk -v chr=${CHR} '{print $0, chr}' > counts/${CHR}_rf10.fst.count.txt

done
