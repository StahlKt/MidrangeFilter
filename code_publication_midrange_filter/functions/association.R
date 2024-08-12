#functions for running association with linear model 

#runs association test
#INPUT
#which.snp: string of SNP name 
#input.table: data.table containing genetic data in long format
#             colums: SNP, SAMPLE, GT (for simulated genotypes) or DOSAGE or BEST_GUESS
#flag.imputed: Boolean or 0,1, indicating if input.table is imputed(1) or complete (0)
#             for 1 association is run on DOSAGE and BEST_GUESS
#             for 0 association is run on GT
#run.signifier: string that indicates current run, needs to be unique for each run
#OUTPUT
#data.table with association results
#columns SNP, FORMAT, P_VALUE, RUN (run.signifier)
snp_association<-function(which.pos, input.table, flag.imputed, run.signifier){
 library(stats)
  snp_table<-input.table[POS==which.pos]
   if(flag.imputed==0){
    #reduce data 
    
  
    #check if SNP is monomorphic 
    if(snp_table[,length(unique(GT))]==1){
      results<-data.table(
        POS = which.pos,
        FORMAT = "GT",
        P_VALUE = 1,
        RUN = run.signifier
        )
      return(results)
    }
    
    
    #make linear model for SNP and case control Status
    glm_snp<- glm(formula = CASE~GT, family = "binomial", data = snp_table)
  
    #returns p-value
    results<-data.table(
      POS = which.pos,
      FORMAT = "GT",
      P_VALUE = summary(glm_snp)$coefficients["GT","Pr(>|z|)"],
      RUN = run.signifier
    )
    return(results)
  }
  
  else if(flag.imputed==1){
    
    
    #check for monomorphic SNPs first, then lm for association (SNP and case control status)
    if(snp_table[,length(unique(BEST_GUESS))]==1){
      
      results_bg<-data.table(
        POS = which.pos,
        FORMAT = "BEST_GUESS",
        P_VALUE = 1,
        RUN = run.signifier
      )
    }
    else{
      lm_best_guess<- glm(CASE~BEST_GUESS, family = "binomial", data=snp_table)
      results_bg<-data.table(
        POS = which.pos,
        FORMAT = "BEST_GUESS",
        P_VALUE = summary(lm_best_guess)$coefficients["BEST_GUESS", "Pr(>|z|)"],
        RUN = run.signifier
      )
    }

    if(snp_table[,length(unique(DOSAGE))]==1){
      
      results_dosage<-data.table(
        POS = which.pos,
        FORMAT = "DOSAGE",
        P_VALUE = 1,
        RUN = run.signifier
      )
    }
    
    else{
    lm_dosage<- glm(CASE~DOSAGE, family = "binomial", data=snp_table)
    results_dosage<-data.table(
      POS = which.pos,
      FORMAT = "DOSAGE",
      P_VALUE = summary(lm_dosage)$coefficients["DOSAGE","Pr(>|z|)"],
      RUN = run.signifier
      )
    }
    return(rbind(results_dosage, results_bg))
  }
}



# 
# 
# split_vcf<-function(run.nr, which.half, 
#                     exec.bcftools){
#   
#   run_name<-paste(run_name_stem,
#                   run.nr,
#                   sep = "_")
#   
#   
#   reg_imputation<-loadRegistry(get_registries(run_name)[3],
#                                writeable = TRUE, make.default = FALSE)
#   
#   
#   #n_snps<-499659
#   
#   start_2nd<-21835126
#   end_2nd<- 59118783
#   end_1st<-22020888 
#   
#   
#   
#   
#   result_case<-loadResult(1,reg=reg_imputation)[1]
#   
#   result_control<-loadResult(4,reg=reg_imputation)[1]
#   
#   
#   if(which.half==1){
#     out_file_case<-paste0(substr(result_case,
#                                  1,
#                                  nchar(result_case)-7),"_1st.vcf.gz")
#     
#     
#     
#     
#     system2(exec.bcftools, paste0("view --regions 19:1-", 
#                                   end_1st, " -Oz -o ", 
#                                   out_file_case, " ",
#                                   result_case ))
#     
#     #index
#     system2(exec.bcftools, c("index -f", out_file_case))
#     
#     
#     
#     
#     out_file_control<-paste0(substr(result_control,
#                                     1,
#                                     nchar(result_control)-7),"_1st.vcf.gz")
#     
#     
#     
#     
#     system2(exec.bcftools, paste0("view --regions 19:1-", 
#                                   end_1st, " -Oz -o ", 
#                                   out_file_control, " ",
#                                   result_control ))
#     
#     #index
#     system2(exec.bcftools, c("index -f", out_file_control))
#     
#   }
#   if(which.half==2){
#     out_file_case<-paste0(substr(result_case,
#                                  1,
#                                  nchar(result_case)-7),"_2nd.vcf.gz")
#     
#     
#     
#     
#     system2(exec.bcftools, paste0("view --regions 19:", 
#                                   start_2nd, "-", end_2nd, " -Oz -o ", 
#                                   out_file_case, " ",
#                                   result_case ))
#     
#     #index
#     system2(exec.bcftools, c("index -f", out_file_case))
#     
#     
#     
#     
#     out_file_control<-paste0(substr(result_control,
#                                     1,
#                                     nchar(result_control)-7),"_2nd.vcf.gz")
#     
#     
#     
#     
#     system2(exec.bcftools, paste0("view --regions 19:", 
#                                   start_2nd, "-", end_2nd, " -Oz -o ", 
#                                   out_file_control, " ",
#                                   result_control))
#     
#     #index
#     system2(exec.bcftools, c("index -f", out_file_control))
#     
#   }
#   
#   
#   return(c(out_file_case,out_file_control))
# }

