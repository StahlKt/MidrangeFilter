
get_all_mr2<-function(x){
  IA_table<-readRDS(list_runs$PATH[x])[,.(POS, MAF, r2_MACH, CASE)]
  IA_table[, RUN:=list_runs[x,RUN]]
  IA_table[, NAME:= list_runs[x,NAME]]
  
  IA_table<-IA_table[,.(POS, MAF, r2_MACH, RUN,NAME,CASE)]
  IA_table<-IA_table[MAF>0.5, MAF:=1-MAF]
  
  combined_set<-pos_ref_alt[IA_table, on=.(POS)]
  
  
  combined_set<-merge(combined_set, new_pos, by.x="POS", by.y="OLD_POS")
  
  if(length(combined_set$POS)!= length(IA_table$POS)){
    return(list_runs$PATH[x])}
  
  combined_set[,POS:=NULL]
  setnames(combined_set, "NEW_POS", "POS")
  combined_set<-combined_set[,c(1,9,2,3,4,5,6,7,8)]
  setnames(combined_set, "r2_MACH", "Rsq")
  
  
  combined_set_common<-combined_set[MAF>0.05]
  combined_set_lowfreq<-combined_set[MAF>0.005 & MAF<=0.05]
  combined_set_rare<-combined_set[MAF<=0.005]
  
  
  
  if(length(combined_set_common$POS)!=0){
    get_magicalrsq(x, combined_set_common, "common")
    
  }
  
  if(length(combined_set_rare$POS)!=0){
    
    get_magicalrsq(x, combined_set_rare,  "rare")
    
  }
  
  if(length(combined_set_lowfreq$POS)!=0){
    get_magicalrsq(x, combined_set_lowfreq, "lowfreq")
    
    
  }
  
  return(list.files(dir_out_magicr, paste0(list_runs[x,NAME],"_",list_runs[x,RUN])))
  
}





#use magicalrsq functionos as described in the MagicalRsq instructions
get_magicalrsq<-function(i,
                         input.data,
                         how.rare){
  
  dir_start_magicr<-file.path(dir_out_magicr,
                              paste0("starting_point_",
                                     i,".txt.gz"))
  
  write.csv(input.data,
            dir_start_magicr,
            row.names = F) 
  
  #
  dir_mid_magicr<-file.path(dir_out_magicr,
                            paste0("mid_step",
                                   i, ".txt.gz"))
  #
  dir_mid_magicr_no_AF<-file.path(dir_out_magicr,
                                  paste0("mid_step_no_AF",
                                         i, ".txt.gz"))
  
  
  integrated_set_no_AF <- data_integrate(dir_start_magicr, 
                                         chr = 19,
                                         pos_col = 2,
                                         SHIC_dir = file.path(dir_magicalrsq_model, "SHIC/"),
                                         #AF_dir = file.path(dir_magicalrsq_model, "AF/"),
                                         outfile = dir_mid_magicr_no_AF)
  
  
  
  #raum ist hier SHIC ausdequoted aber oben nicht?
  integrated_set <- data_integrate(dir_mid_magicr_no_AF, 
                                   chr = 19,
                                   pos_col = 2,
                                   SHIC_dir = file.path(dir_magicalrsq_model, "SHIC/"),
                                   AF_dir = file.path(dir_magicalrsq_model, "AF/"),
                                   outfile = dir_mid_magicr)
  
  
  
  
  #calculate magicr2
  dir_finished_magicr<-file.path(dir_out_magicr,
                                 paste0("magicr_out",
                                        i, ".txt.gz"))
  dir_finished_magicr_no_AF<-file.path(dir_out_magicr,
                                       paste0("magicr_out_no_AF", 
                                              i, ".txt.gz"))
  
  
  sapply(1:4, function(j){
    dir_model<-file.path(dir_magicalrsq_model,
                         paste0(models[j],"_", how.rare))
    
    
    
    result_table<-as.data.table(calc_MagicalRsq(file = dir_mid_magicr_no_AF, 
                                                model = dir_model,
                                                FeatureCols = c(5,6,10:75), 
                                                keptCols = 1:9,
                                                outfile = dir_finished_magicr_no_AF))
    
    
    result_table[,AF_FLAG:=0]
    
    if(length(integrated_set$POS)!=0){
      out_table<-as.data.table(calc_MagicalRsq(file = dir_mid_magicr, 
                                               model = dir_model,
                                               FeatureCols = c(5,6,10:86), 
                                               keptCols = 1:9,
                                               outfile = dir_finished_magicr))
      
      out_table[,AF_FLAG:=1]
      result_table<-rbind(out_table,result_table)
      
      file.remove(dir_finished_magicr)
    }
    dir_out_table<-paste0(dir_out_magicr, "/",
                          result_table[1,NAME], "_",result_table[1,RUN],
                          "_mr2_", models[j], "_", how.rare, ".RDS")
    saveRDS(result_table, dir_out_table)
    
  })
  
  
  
  file.remove(dir_finished_magicr_no_AF,
              dir_mid_magicr,
              dir_mid_magicr_no_AF,
              dir_start_magicr)
  
  
  return()
}






