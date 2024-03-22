#chr_W, run 2 iterations of this with a different SNP sample
bcftools view ../vcfs/chr_W_HostOG.SNP.DP3-AC2-MQ40.vcf.gz -Ov -o chr_W_HostOG.SNP.DP3-AC2-MQ40.vcf
ruby snapp_prep.rb -o chr_W -v chr_W_HostOG.SNP.DP3-AC2-MQ40.vcf -t chr_W_samples.txt -c constraints.txt -m 1000 -l 100000
mv snapp.xml chr_W.xml

#also re-run with 5K SNPs 
ruby snapp_prep.rb --xml chr_W-BC.xml -o chr_W-BC -v chr_W_HostOG.SNP.DP3-AC2-MQ40.vcf -t chr_W_samples.txt -c constraints_beast.txt -m 1000 -l 500000

#and second replicate
ruby snapp_prep.rb --xml chr_W-BC2.xml -o chr_W-BC2 -v chr_W_HostOG.SNP.DP3-AC2-MQ40.vcf -t chr_W_samples.txt -c constraints_beast.txt -m 1000 -l 500000

#and with double SNPs 
ruby snapp_prep.rb --xml chr_W-2KBC2.xml -o chr_W-2KBC2 -v chr_W_HostOG.SNP.DP3-AC2-MQ40.vcf -t chr_W_samples.txt -c constraints_beast.txt -m 2000 -l 500000