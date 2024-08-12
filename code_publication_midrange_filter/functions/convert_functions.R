#Convert between formats


#converts from vcf to hap/legend format
#INPUT
#input.data: filepath string to input file in vcf.gz format
#exec.bcftools: string file path to executable for bcftools
#output.prefix: filepath string to desired file and location without format suffix
#               both .hap and .legend are saved under that name
#OUTPUT
#returns output.prefix
vcf_to_hapleg<-function(input.data, exec.bcftools, output.prefix){
system2(exec.bcftools, c("convert -h", output.prefix , input.data))
return(output.prefix)
  }


#convert from vcf to bref2 format
#INPUT
#input.data: filepath string to input file in vcf.gz format
#exec.java filepath string to java executable
#exec.bref filepath string to bref executable
#OUTPUT
#filepath string to output data in .bref format
vcf_to_bref3 <- function(input.data, exec.java, exec.bref, output.data) {
  system2(exec.java, c("-jar", exec.bref, input.data, ">", output.data))
  return(output.data)
}


#convert from vcf to bref2 format
#INPUT
#input.prefix: filepath string to the hap input file without the suffixes
#exec.bcftools filepath string to bcftools executable
#output.data: string filepath to desired output data in vcf.gz format and location
#OUTPUT
#returns output.data
hapleg_to_vcf <- function(input.prefix, exec.bcftools, output.data) {
  
  hap_dir<-paste(input.prefix,".new.hap",sep="")
  
 
    leg_dir<-gsub(".cases","",input.prefix)
    leg_dir<-gsub(".controls","",leg_dir)
   leg_dir<-paste(leg_dir,".legend",sep="")
   
   sample_dir<-paste(input.prefix,".new.sample",sep="")
  
  
  #convert to vcf with bcftools
  system2(exec.bcftools, c("convert -o", output.data,
                           "-Oz -H", paste(hap_dir,leg_dir,sample_dir,sep=",")))
  
  #return converted file path
  return(output.data)
}



#convert from vcf to plink
#INPUT
#input.data: filepath string to the input file in vcf.gz
#exec.plink filepath string to plink executable
#output.prefix: filepath string to desired files in the same location without format suffixes
#OUTPUT
#returns output.prefix
vcf_to_plink<-function(input.data, exec.plink, output.prefix){
  
  system2(exec.plink, c("--vcf", input.data, "--make-bed --out", output.prefix))
  return(output.prefix)
}


#convert from plink bed/bim/fam to vcf 
#INPUT
#input.prefix: filepath string to the plink bed/bim/fam files without their suffixes
#exec.plink filepath string to plink executable
#output.data: string filepath to desired output data in vcf.gz format and location
#OUTPUT
#returns output.data
plink_to_vcf<-function(input.prefix, exec.plink, output.data){
  out_prefix<-gsub(".vcf.gz","", output.data)
  system2(exec.plink, c("--bed", paste(input.prefix,".bed", sep=""),
                        "--bim", paste(input.prefix,".bim", sep=""),
                        "--fam", paste(input.prefix,".fam", sep=""),
                        "--recode vcf --out", out_prefix))
  system2("gzip", paste(out_prefix, ".vcf", sep=""))
  
  unlink(paste(out_prefix, ".vcf", sep=""))
  return(output.data)
}


