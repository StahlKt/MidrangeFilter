#problem_wrap
#deletion, QC, phasing


#flag.cases.problem=1 for cases, 0 for controls
#flag.random.deletion=1 means snps are randomly deleted, else provide list.deleted.snps
#data is the result from the batchmap simulation
#wrapperfunction that includes deletion, QC and Phasing before imputation
#INPUT
#data: filepath string to simulated testdata set in vcf.gz format
#exec.bcftools filepath string to bcftools executable
#exec.plink filepath string to plink bcftools executable
#exec.beagle filepath string to beagle executable
#map.plink filepath string to genetic map file in plink format
#flag.random.deletion: boolean or 1/0 to indicase if the deletion is random(1) 
#                      or according to list(0)
#fraction.to.keep: numeric between 0 and 1, indication the faction of simulated data
#                   used as a basis for imputation. 1-fraction.to.keep is deleted
#                   not needed if flag.random.deletion is 0
#list.snps.to.keep list of SNP ids in the data set that are to be kept for the imputation
#                   not needed if flag.random.deletion is 1
#run.signifier: string that indicates current run, needs to be unique for each run
#flag.cases.problem: boolean or 0/1 to indicate if the input data set contains cases(1)
#                     or controls(0)
#OUTPUT
#returns vector of fileüath strings. first one is the test data set after cq and with deleted SNPs
#for imputation, second is the test data set after qc with all SNPs
problem_wrap<-function(data,
                       job,
                       exec.bcftools,
                       exec.plink,
                       exec.beagle,
                       map.plink,
                       flag.random.deletion,
                       fraction.to.keep=NULL,
                       list.snps.to.keep=NULL,
                       run.signifier
                       ){
  
  draw_seed<-round(runif(1,1,1000000))
  set.seed(draw_seed)
  


  # check case/controls
  if(str_detect(data, "cases")){
    case_string<-"cases"
  }
  else{
    case_string<-"controls"
  }
  
 
  #delete SNPs
  
  dir_test_data_deleted<-paste(substr(data,1,nchar(data)-7),"_deleted.vcf.gz",sep="")
 
  if(flag.random.deletion==1){
    extract_random_snps(data, exec.bcftools, fraction.to.keep, dir_test_data_deleted, run.signifier, draw_seed)
  }
  else{
    extract_listed_snps(data, list.snps.to.keep, exec.bcftools,dir_test_data_deleted)
  }
 

  dir_vcf_deleted_phased<-beagle_phase(dir_test_data_deleted,
                                       exec.beagle,
                                       map.plink,
                                       file.path(dir_out_phasing,paste(
                                         run.signifier,
                                         "_", case_string,
                                         ".vcf.gz", sep="")),
                                       draw_seed)

  
  return(dir_vcf_deleted_phased)
}
