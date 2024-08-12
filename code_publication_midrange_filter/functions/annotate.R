#INPUT
#input.data: string file path to data in vcf.gz format
#exec.bcftools: string file path to executable for bcftools
#output.data: string filepath to desired output data in vcf.gz format and location
#OUTPUT
#returns output.data
annotate_data <- function(input.data, exec.bcftools, output.data) {
  
  #set paths for intermediate files
  filtered <- file.path(dir_input, "data_filtered.vcf.gz")
  drop <- file.path(dir_input, "data_drop.vcf.gz")
  
  #filter out SNPs that does not "PASS" the quality criteria and have less than 10% missing entries
  system2(exec.bcftools, c("filter -i \"FILTER='PASS' && TYPE='snp' && F_PASS(GT='./.')<0.1\"", "-o", filtered, "-O z", input.data))
  # drop anything but biallelic SNPs
  system2(exec.bcftools, c("view -m2 -M2 -v snps", "-o", drop, "-O z", filtered))
  #annotate SNPsnames
  system2(exec.bcftools, c("annotate --set-id '%CHROM\\:%POS\\_%REF\\_%FIRST_ALT'", "-o", output.data, "-O z", drop))
  
  #remove obsolete files
  file.remove(filtered)
  file.remove(drop)
  
  #index final file
  system2(exec.bcftools, c("index -f", output.data))
  
  return(output.data)
}


#set registry names

#INPUT
#run.signifier: string that indicates current run, needs to be unique for each run
#OUTPUT
#returns vector of registries and their location marked with run signifier
get_registries<-function(run.signifier){
  
  return(c(paste(dir_reg_simulation, "_", run.signifier, sep=""),
           paste(dir_reg_reference, "_", run.signifier, sep=""),
           paste(dir_reg_imputation, "_", run.signifier, sep=""),
           paste(dir_reg_association, "_", run.signifier, sep=""),
           paste(dir_reg_char, "_", run.signifier, sep=""),
           file.path(dir_reg, paste0("set_sig_imp_", run.signifier))
           
           #paste(dir_reg_LD,"_", run.signifier, sep=""),
           #paste(dir_reg_split, "_", run.signifier, sep="")))  
))
}


