#new reference functions

reference_wrap_no_load<-function(job.ids.merge, 
                         exec.bcftools, run.signifier){
  
  job_ids_simulation<-unlist(job.ids.merge)
  #initiate the file paths for the reference panels
  dir_reference_cases <- file.path(dir_out_reference,
                                   paste(run.signifier, "_cases_reference.vcf.gz", sep=""))
  dir_reference_controls <- file.path(dir_out_reference,
                                      paste(run.signifier, "_controls_reference.vcf.gz", sep=""))
  
  #divide results into cases and controls
  all_cases <- as.vector(sapply(job_ids_simulation, get_simulation_result_path, run.name = run.signifier, case.flag=1))
  all_controls <- as.vector(sapply(job_ids_simulation, get_simulation_result_path, run.name = run.signifier, case.flag=1))
  
  #check run.signifier
  if(!grepl(run.signifier,all_cases[1])){
    print("Error: run.signifier does not match previous run!")
    return()
  }
  
  dir_reference_cases_bref<-file.path(dir_out_reference,
            paste(run.signifier, "_cases_reference.bref3", sep=""))
  
  dir_reference_controls_bref<-file.path(dir_out_reference,
            paste(run.signifier, "_controls_reference.bref3", sep=""))
  
  
  
  if(!file.exists(dir_reference_cases_bref)){
    #merge files
     system2(exec.bcftools, c("merge", "-O z -o", dir_reference_cases, all_cases))
    #index files
    system2(exec.bcftools, c("index -f", dir_reference_cases))
    #formating into bref 3
    dir_reference_cases_bref<-vcf_to_bref3(dir_reference_cases, dir_exec_java, dir_exec_bref, dir_reference_cases_bref)
    #delete vcf files
    unlink(dir_reference_cases)
    unlink(paste(dir_reference_cases,".csi", sep=""))
  }

  
  if(!file.exists(dir_reference_controls_bref)){
  #merge files
  system2(exec.bcftools, c("merge", "-O z -o", dir_reference_controls, all_controls))
  #index files
  system2(exec.bcftools, c("index -f", dir_reference_controls))
  #formating into bref 3
  dir_reference_controls_bref<-vcf_to_bref3(dir_reference_controls, dir_exec_java, dir_exec_bref,dir_reference_controls_bref)
  #delete vcf files
  unlink(dir_reference_controls)
  unlink(paste(dir_reference_controls, ".csi", sep=""))
  }
  
  is_small_ref<-ifelse(length(job_ids_simulation)>5, 0,1)
  compare_size<-ifelse(is_small_ref, 200000000,300000000)
  
  if(!file.exists(dir_reference_cases_bref)) return(FALSE)
  if(file.info(dir_reference_cases_bref)$size>compare_size){
    list_delete_files<-sapply(1:length(job_ids_simulation)+1,
                              get_simulation_result_path, 
                              run.name=run.signifier, case.flag=1)
    
    list_delete_files<-c(list_delete_files, paste0(list_delete_files,".csi"))
    unlink(list_delete_files)
  }
  
  if(!file.exists(dir_reference_controls_bref)) return(FALSE)
  if(file.info(dir_reference_controls_bref)$size>compare_size){
    
    list_delete_files<-sapply(1:length(job_ids_simulation)+1,
           get_simulation_result_path, 
           run.name=run.signifier, case.flag=0)
    
    list_delete_files<-c(list_delete_files, paste0(list_delete_files,".csi"))
    unlink(list_delete_files)
  }
  return(c(dir_reference_cases_bref,dir_reference_controls_bref))
}
