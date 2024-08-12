#functions for using the beagle executable for phasing and imputation


#imputation with Beagle5.1
#INPUT
#data, job and instance as described in batchtools, 
#instance is vector of filepath strings, first entry to deleted and phased dataset
#                                 second entry to the complete simulated dataset
#flag.cases.algoritm: Boolean or 1,0 indicating case data(1) or control data (0)
#exec.beagle: filepath string to beagle executable
#map.plink: filepath string to genetic map data in plink format
#reference.bref3: filepath string to genetic reference data in bref3 format
#chunk.length: integer for how long each imputed chunk is going to be
#region.lower.bound: integer on which position imputation is going to start
##exec.bcftools: filepath string to bcftools executable
#run.signifier: string that indicates current run, needs to be unique for each run
#OUTPUT
#vector of 2 filepath strings, first one to imputed genetic data in vcf.gz
#                               second one echoes instance[2] to the complete simulated dataset
beagle_impute <- function(data, 
                          job,
                          instance,
                          flag.cases.algorithm,
                          exec.beagle,
                          map.plink,
                          reference.bref3,
                          chunk.length,
                          region.lower.bound,
                          exec.bcftools,
                          run.signifier
                          ) {
  
  #check run.signifier
  if(!grepl(run.signifier,instance)){
    print("Error: run.signifier does not match previous step!")
    return()
  }
  
  
  draw_seed<-round(runif(1,1,1000000))
  #define chunks
  set.seed(draw_seed)
  inte <- region.lower.bound + chunk.length - 1
  intstring1 <- toString(format(region.lower.bound, scientific = F))
  intstring2 <- toString(format(inte, scientific = F))
  
  
  #define output_name 
  if(str_detect(data, "cases")){
    case_string<-"cases"
  }
  else{
    case_string<-"controls"
  }
  output_prefix<-file.path(dir_out_imputation,paste(run.signifier, "_",
                                                  case_string, "_",
                                                  "imputed", "_",
                                                  region.lower.bound,
                                                  sep=""))
  
  output_data<-paste(output_prefix,".vcf.gz", sep="")
  
  #use Beagle5 for imputation
  system2(dir_exec_java,
          c("-jar",
            exec.beagle,
            paste("ref=", reference.bref3, sep = ""),
            paste("gt=",instance, sep = "" ),
            paste("out=", output_prefix, sep = ""),
            paste("map=", map.plink, sep = ""),
            paste("chrom=", 19,":", intstring1, "-", intstring2, sep = ""),
            paste("seed=", draw_seed, sep = ""),
            "nthreads=1 gp=true"))
  
  #return path to output file and simulated file for comparison
  system2(exec.bcftools, c("index -f", paste(output_prefix, ".vcf.gz", sep="")))
  

  #index annotated file
  system2(exec.bcftools, c("index -f", output_data))
  return(c(output_data))
}




beagle_impute_no_chunks <- function(data, 
                          job,
                          instance,
                          flag.cases.algorithm,
                          exec.beagle,
                          map.plink,
                          reference.bref3,
                          exec.bcftools,
                          run.signifier
) {
  
  #check run.signifier
  if(!grepl(run.signifier,instance)){
    print("Error: run.signifier does not match previous step!")
    return()
  }
  
  
  draw_seed<-round(runif(1,1,1000000))
  #define chunks
  set.seed(draw_seed)

  #define output_name   
  if(str_detect(data, "cases")){
    case_string<-"cases"
  }
  else{
    case_string<-"controls"
  }

  output_prefix<-file.path(dir_out_imputation,paste(run.signifier, "_",
                                                    case_string, "_",
                                                    "imputed",
                                                    sep=""))
  
  output_data<-paste(output_prefix,".vcf.gz", sep="")
  
  #use Beagle5 for imputation
  system2(dir_exec_java,
          c("-jar",
            exec.beagle,
            paste0("ref=", reference.bref3),
            paste0("gt=",instance ),
            paste0("out=", output_prefix),
            paste0("map=", map.plink),
            paste0("chrom=", 19),
            paste0("seed=", draw_seed),
            "nthreads=1 gp=true"))
  
  #return path to output file and simulated file for comparison
  system2(exec.bcftools, c("index -f", paste(output_prefix, ".vcf.gz", sep="")))
  
  
  #index annotated file
  system2(exec.bcftools, c("index -f", output_data))
  return(output_data)
}





##phase with beagle5.1
#INPUT
#input.data: filepath string to dataset in vcf.gz format
#exec.beagle: filepath string to beagle executable
#map.plink: filepath string to genetic map data in plink format
#output.data: string filepath to desired output data in vcf.gz format and location
#seed.no: integer to set the seed to
#OUTPUT
#returns filepath string to output data in vcf.gz format

beagle_phase <- function(input.data,
                         exec.beagle,
                         map.plink,
                         output.data,
                         seed.no) {
  
  out_prefix<-gsub(".vcf.gz","",output.data)
  #execute phasing with Beagle 5.1  
  system2(dir_exec_java,
          c("-jar",
            exec.beagle,
            paste("gt=", input.data, sep = ""),
            paste("out=", out_prefix, sep = ""),
            paste("map=", map.plink, sep = ""),
            paste("seed=", seed.no, sep = ""),
            "nthreads=1 impute=false"))
  

  
  return(output.data)
}
