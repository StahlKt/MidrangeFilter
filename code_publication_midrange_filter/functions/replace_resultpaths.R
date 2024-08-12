#functions to replace loadResult, if the result is a path to a file.


get_results_imputation_path<-function(run.name, case.flag, imputation.flag){
  path_stem<-ifelse(imputation.flag, dir_out_imputation, dir_out_hapgen)
  suffix<-ifelse(imputation.flag, "_imputed.vcf.gz", "_annotated.vcf.gz")
  case_signifier<-ifelse(case.flag, "cases", "controls")
  separator<-ifelse(imputation.flag, "_", "_1.")
  file_path<-file.path(path_stem, paste0(run.name, separator,case_signifier,suffix))
  return(file_path)
}


get_reference_panel_path<- function(run.name, case.flag){
  path_reference<-ifelse(case.flag, file.path(dir_out_reference,
                                      paste0(run.name, "_cases_reference.bref3")),
                                        file.path(dir_out_reference,
                                         paste0(run.name, "_controls_reference.bref3")))
  
  return(path_reference)
}


get_simulation_result_path<-function(run.name, case.flag, job.id){
 
  suffix<-ifelse(case.flag, ".cases_annotated.vcf.gz", ".controls_annotated.vcf.gz")
  file_path<-file.path(dir_out_hapgen, paste0(run.name, "_", job.id, suffix))
  return(file_path)
}

