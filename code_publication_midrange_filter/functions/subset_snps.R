
#extracts SNPs randomly from file to keep for basis of imputation
#INPUT 
#input.data: filepath string to a vcf.gz file where SNPs should be deleted from (after simulation)
#exec.bcftools: filepath string to bcftools executable
#fraction to keep: numeric between 0 and 1 to indicate the percentage of how many SNPs to keep.
#     SNPs will be selected randomly
#output.data: filepath string to desired vcf file and location with reduced SNPs
#run.signifier: string that indicates current run, needs to be unique for each run
#seed.no: integer for seed setting
#OUTPUT
#return output.data
extract_random_snps<-function(input.data, exec.bcftools, fraction.to.keep, output.data, run.signifier, seed.no){
  
  #set directories for extracted SNP lists
  dir_list_complete_snps<- file.path(dirname(input.data),paste(run.signifier,"_complete_list_SNPs.txt",sep=""))
  dir_list_kept_snps<- file.path(dirname(input.data),paste(run.signifier, "_kept_list_SNPs.txt",sep=""))
  
  #set seed
  set.seed(seed.no)

  
  #extract list with all SNPs
  system2(exec.bcftools, c("query -f '%CHROM\t%POS\n'", input.data, ">", dir_list_complete_snps))
  
  #read in list and extract length=number of SNPs
  SNPs_complete<-read.delim(dir_list_complete_snps)
  number_snps<-dim(SNPs_complete)[1]
  
  #draw randomly the amout of SNPs to keep as intended with fraction.to.keep and sort them
  #(draws row numbers to keep)
  SNPs_random<-sort(sample(1:number_snps, round(number_snps*fraction.to.keep),replace = FALSE))
  
  #subset the SNP list according to random pattern
  SNPs_random<-SNPs_complete[SNPs_random, ]
  
  #save table with SNPs to keep
  write.table(SNPs_random, file=dir_list_kept_snps, row.names = FALSE, col.names = FALSE, sep="\t")
  
  #create new VCF-file with only kept SNPs
  system2(exec.bcftools, c("view -R", dir_list_kept_snps,
                               "-O z -o", output.data, input.data))
  
  #remove obsolete SNP list files
  file.remove(dir_list_complete_snps)
  file.remove(dir_list_kept_snps)
  
  #return path to output data
  return(output.data)
}




#extracts SNPs accordning to List to keep for basis of imputation(e.g. to simulate SNP array)
#INPUT 
#input.data: filepath string to a vcf.gz file where SNPs should be deleted from (after simulation)
#exec.bcftools: filepath string to bcftools executable
#input.snp.list: filepath string to simple txt file (one row per SNP), which contains SNPs to keep
#output.data: filepath string to desired vcf file and location with reduced SNPs
#OUTPUT
#return output.data
extract_listed_snps<-function(input.data, input.snp.list, exec.bcftools, output.data){
  
  #create new VCF-file with only kept SNPs
  system2(exec.bcftools, c("view -R", input.snp.list,
                           "-O z -o", output.data, input.data))
  
  #return path to output data
  return(output.data)
}




