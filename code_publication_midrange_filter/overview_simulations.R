#overview runs

library(batchtools)
library(data.table)
#set working directory to midrange filter folder
source(file.path(getwd(),"file_paths.R"))

#sets according to paper/supplement
# b = balanced, so 1000 cases, 1000 controls
# lc = low cases, so 333 cases, 1000 controls
# rd = random deletion
# o5 = illumina omni 5
# o25 = illumina omni 2.5
# oe = illumina express
# sr = simulation with small reference panel, so 5000 instead of 10000
# ipr = improper fit of the reference panel: specific subpopulations in the test data are not present in the refence


run_name_stem<-"set3_b_rd"
disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
run_seed<-10052023



run_name_stem<-"set3_b_o5" 
disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
run_seed<-05072023



run_name_stem<-"set1_b_rd" 

disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))

run_seed<-05072023



run_name_stem<-"set3_lc_rd" 
disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
run_seed<-20072023

run_name_stem<-"set1_lc_rd" 
disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
run_seed<-20072023

run_name_stem<-"set3_sr_b_rd" 
disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
run_seed<-20072023

run_name_stem<-"set1_sr_b_o5" 
disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
run_seed<-20072023



run_name_stem<-"set3_b_o25"
disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
run_seed<-11122023


run_name_stem<-"set3_b_oe" 
  disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
  run_seed<-11122023

run_name_stem<-"set1_b_o25"
  disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
  run_seed<-11122023

run_name_stem<-"set1_b_oe"
disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
  run_seed<-11122023

  
run_name_stem<-"set2_b_o5"
disease_list<-readRDS(file.path(dir_rds, "set2_disease_loci.RDS"))
run_seed<-11122023

run_name_stem<-"set2_b_o25" 
  disease_list<-readRDS(file.path(dir_rds, "set2_disease_loci.RDS"))
  run_seed<-11122023

  
run_name_stem<-"set2_b_oe"
  disease_list<-readRDS(file.path(dir_rds, "set2_disease_loci.RDS"))
  run_seed<-11122023



run_name_stem<-"set4_b_o5" 
  disease_list<-readRDS(file.path(dir_rds, "set4_disease_loci.RDS"))
  run_seed<-11122023

run_name_stem<-"set4_b_o25" 
  disease_list<-readRDS(file.path(dir_rds, "set4_disease_loci.RDS"))
  run_seed<-11122023

run_name_stem<-"set4_b_oe" 
  disease_list<-readRDS(file.path(dir_rds, "set4_disease_loci.RDS"))
  run_seed<-11122023 


  run_name_stem<-"set1_sr_b_o25" 
  disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
    run_seed<-11122023

    run_name_stem<-"set1_sr_b_oe" 
    disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
    run_seed<-11122023
    
  
    run_name_stem<-"set4_sr_b_o5"
    disease_list<-readRDS(file.path(dir_rds, "set4_disease_loci.RDS"))
    run_seed<-11122023 
    

    run_name_stem<- "set3_lc_o5"
    disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
    run_seed<-27022024
    
    run_name_stem<- "set3_lc_o25"
    disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
    run_seed<-27022024
    
    run_name_stem<-"set1_lc_oe"
    disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
    run_seed<-27022024
  
    run_name_stem<-"set2_lc_o5"
    disease_list<-readRDS(file.path(dir_rds, "set2_disease_loci.RDS"))
    run_seed<-27022024
    

    run_name_stem<-"set2_lc_o25"
    disease_list<-readRDS(file.path(dir_rds, "set2_disease_loci.RDS"))
    run_seed<-27022024
    
    
    run_name_stem<-"set2_lc_oe"
    disease_list<-readRDS(file.path(dir_rds, "set2_disease_loci.RDS"))
    run_seed<-27022024
   
    run_name_stem<-"set4_lc_o5"
    disease_list<-readRDS(file.path(dir_rds, "set4_disease_loci.RDS"))
    run_seed<-27022024
    
    run_name_stem<-"set4_lc_o25"
    disease_list<-readRDS(file.path(dir_rds, "set4_disease_loci.RDS"))
    run_seed<-27022024



#added in review
    
run_name_stem<-"set1_ipr_b_o5" 
disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
run_seed<-23102024


run_name_stem<-"set2_ipr_b_o25" 
disease_list<-readRDS(file.path(dir_rds, "set2_disease_loci.RDS"))
run_seed<-23102024


run_name_stem<-"set1_ipr_b_oe" 
disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
run_seed<-23102024


run_name_stem<-"set4_ipr_b_o5" 
disease_list<-readRDS(file.path(dir_rds, "set4_disease_loci.RDS"))
run_seed<-12112024

run_name_stem<-"set3_ipr_b_o5" 
disease_list<-readRDS(file.path(dir_rds,"set3_disease_loci.RDS"))
run_seed<-12112024

run_name_stem<-"set1_ipr_lc_o5" 
disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
run_seed<-12112024


