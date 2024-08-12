library(R.utils)
#run chartable_no_chunks_one_job for char tables

#flag.extended.table kann auch auch c(0,1) gesetzt werden. 
chartable_no_chunks_one_job<-function(run.nr, run.name.stem=run_name_stem){
  
  run_name<-paste(run.name.stem,
                  run.nr,
                  sep = "_")
  
 
  dir_saved_file_char<-file.path(dir_out_char, paste0(run_name, "_acc_char.RDS"))
  
  
  if(all(file.exists(dir_saved_file_char))){
    print(paste(run_name, "characteristics files already exist"))
    return(NA)
  }
  
  
  dir_out_association<-"/home/uni08/stahl11/adjusting_pvalue/data/output/association"
  
  

  pv_table<-readRDS(file.path(dir_out_association, paste0(run_name, "_pv.RDS")))
  
  
  
  #list_pos as all SNPs that are genomewide significant
  list_pos<-pv_table[(FORMAT=="DOSAGE" | FORMAT=="BEST_GUESS") & P_VALUE<5*10^(-8),unique(POS)]
  
  
  rm(pv_table)
  #if there's no signal in the run, abort
  if(length(list_pos)<1){
    
    print(paste("No genomewide significant SNPs in run", run_name))
    return()
  }
  
  #if there are significant snps, then get characteristics
  if(length(list_pos)>0){
    
    region_file_name<-file.path(dir_out_char, paste0(run_name, "_char_regions.text"))
    #set file to read in
    
    set_file<-data.frame(CHROM=19, POS=list_pos)

    write.table(x=set_file,file=region_file_name, quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
    
    
    
    
   get_characteristics_wrap(dir_exec_bcftools,
                            dir_exec_python_bashscript, 
                            region_file_name,
                            run_name)
              

   dir_saved_file_char<-file.path(dir_out_char, paste0(run_name, "_acc_char.RDS"))
   
   return(all(file.exists(dir_saved_file_char)))


  }
  }

    
    #define function to get registries
get_characteristics_wrap<-function( exec.bcftools,
                                    exec.python.script , 
                                    list.pos,
                                    run.signifier,
                                    flag.local=FALSE){
   
        #get file paths of the (whole) imputed files.
        result_case<-get_results_imputation_path(run.signifier,1,1)
        result_control<-get_results_imputation_path(run.signifier,0,1)
        
        
        if(!file.exists(paste0(result_case, ".csi"))){
         c
        }
        
        if(!file.exists(paste0(result_control, ".csi"))){
          system2(exec.bcftools, c("index -f", result_control))
        }
        
        
        #check run.signifier
        if(!grepl(run.signifier, result_case)){
          print("Error: run.signifier does not match previous run!")
          return()
        }
        
        
        
        
        
        
        #read in the imputed snps
        char_table_cases<-read_table_char_snps(result_case,
                                               exec.bcftools,1,
                                               list.pos)
        
        
        
        
        n_sample_cases<-char_table_cases[, length(unique(SAMPLE))]
        
        char_table_cases_MAF<-char_table_cases[, (sum(GP_1)+sum(GP_2)*2)/(2*n_sample_cases), by=.(POS)]
        setnames(char_table_cases_MAF, "V1", "MAF")
        
        #if the alternate allel is not the rarer one, recalculate MAF
        char_table_cases_MAF[MAF>0.5, MAF:=1-MAF]
        
        
        
        char_table_controls<-read_table_char_snps(result_control,
                                                  exec.bcftools,0,
                                                  list.pos)
        
        n_sample_controls<-char_table_controls[, length(unique(SAMPLE))]
        
        char_table_controls_MAF<-char_table_controls[, (sum(GP_1)+sum(GP_2)*2)/(2*n_sample_controls), by=.(POS)]
        setnames(char_table_controls_MAF, "V1", "MAF")
        char_table_controls_MAF[MAF>0.5, MAF:=1-MAF]
        
        
        #combine tables
        char_table<-rbind(char_table_controls_MAF[, CASE:=0],
                          char_table_cases_MAF[, CASE:=1])
        
        #use imputeAccure to get imputation characteristics
        accur_table<-get_table_imputeAccure(result_case, result_control, exec.bcftools,
                                            exec.python.script, list.pos)
        #quality measures with imputeAccure
        char_table<-merge(accur_table,
                          char_table, by=c("POS","CASE"))
        
        dir_saved_file_char<-file.path(dir_out_char, paste0(run.signifier, "_acc_char.RDS"))
        
        saveRDS(char_table, dir_saved_file_char)
        
      
        
        unlink(list.pos)
        
        return(file.exists(dir_saved_file_char))
      }
    
    

                       

