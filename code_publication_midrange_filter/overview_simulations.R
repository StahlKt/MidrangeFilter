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



run_name_stem<-"set3_b_rd"
disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
run_seed<-10052023



run_name_stem<-"set3_b_o5" 
disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
run_seed<-05072023




# mit  POS grid, sonst normal

run_name_stem<-"set1_b_rd" 

disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))

run_seed<-05072023




#2 mit wenig fällen
run_name_stem<-"set3_lc_rd" 
disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
run_seed<-20072023

run_name_stem<-"set1_lc_rd" 
disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
run_seed<-20072023




#2 mit weniger referenz
run_name_stem<-"set3_sr_b_rd" 
disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
run_seed<-20072023

run_name_stem<-"set1_sr_b_o5" 
disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
run_seed<-20072023






#new snp_arrays mit weniger SNPS

run_name_stem<-"set3_b_o25"
#abgeschickt ref, ANDERER SNP CHIP #done ref
disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
run_seed<-11122023


run_name_stem<-"set3_b_oe" #started simulation(redoing expired)
  disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
  run_seed<-11122023

run_name_stem<-"set1_b_o25"#started simulation
  disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
  run_seed<-11122023

run_name_stem<-"set1_b_oe"
disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
  run_seed<-11122023

  


#validate set 1 
run_name_stem<-"set2_b_o5"
#abgeschickt ref, normaler snp chip #ref done
disease_list<-readRDS(file.path(dir_rds, "set2_disease_loci.RDS"))
run_seed<-11122023

run_name_stem<-"set2_b_o25" 
#started_simulation (redoing expired)
  disease_list<-readRDS(file.path(dir_rds, "set2_disease_loci.RDS"))
  run_seed<-11122023

  
run_name_stem<-"set2_b_oe"
#simulation started
  disease_list<-readRDS(file.path(dir_rds, "set2_disease_loci.RDS"))
  run_seed<-11122023


#set4 disease SNPs
  #kurz Angst, dass alle Char vorher nochmal laufen müssen
run_name_stem<-"set4_b_o5" #done ref
#abgeschickt Ref, normaler snp chip   
  disease_list<-readRDS(file.path(dir_rds, "set4_disease_loci.RDS"))
  run_seed<-11122023

run_name_stem<-"set4_b_o25" 
#started_simulation,, ANDERER SNP CHIP 
  disease_list<-readRDS(file.path(dir_rds, "set4_disease_loci.RDS"))
  run_seed<-11122023

run_name_stem<-"set4_b_oe" #started simulation
# ANDERER SNP CHIP 
  disease_list<-readRDS(file.path(dir_rds, "set4_disease_loci.RDS"))
  run_seed<-11122023 

  
  
  
  
  
    
  #low ref!!! hier muss funktion mehr angepasst werden

#low ref less snp density
  run_name_stem<-"set1_sr_b_o25" #done
#asso started
  disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
    run_seed<-11122023

    run_name_stem<-"set1_sr_b_oe" 
    #imputation started
    disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
    run_seed<-11122023
    
    
    
    
    run_name_stem<-"set4_sr_b_o5"
#reference started
    disease_list<-readRDS(file.path(dir_rds, "set4_disease_loci.RDS"))
    run_seed<-11122023 
    
    

    
    
    
    
    #1.Teil Low Cases  simulation abgeschickt
    run_name_stem<- "set3_lc_o5"
    disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
    run_seed<-27022024
    
    run_name_stem<- "set3_lc_o25"
    disease_list<-readRDS(file.path(dir_rds, "set3_disease_loci.RDS"))
    run_seed<-27022024
    
    run_name_stem<-"set1_lc_oe"
    disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
    run_seed<-27022024
    # 2. Teil low cases
    run_name_stem<-"set2_lc_o5"
    disease_list<-readRDS(file.path(dir_rds, "set2_disease_loci.RDS"))
    run_seed<-27022024
    

    run_name_stem<-"set2_lc_o25"
    disease_list<-readRDS(file.path(dir_rds, "set2_disease_loci.RDS"))
    run_seed<-27022024
    
    
    run_name_stem<-"set2_lc_oe"
    disease_list<-readRDS(file.path(dir_rds, "set2_disease_loci.RDS"))
    run_seed<-27022024
    
   # 3. teil Low cases
    run_name_stem<-"set4_lc_o5"
    disease_list<-readRDS(file.path(dir_rds, "set4_disease_loci.RDS"))
    run_seed<-27022024
    
    run_name_stem<-"set4_lc_o25"
    disease_list<-readRDS(file.path(dir_rds, "set4_disease_loci.RDS"))
    run_seed<-27022024
    
