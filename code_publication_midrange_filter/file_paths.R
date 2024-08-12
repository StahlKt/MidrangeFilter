##initialize all variables and functions

library("tools")
library("batchtools")
library("data.table")
library("stringr")

`%notin%` <- Negate(`%in%`)



# home directory/working dir and other folders

dir_home<-"~"

dir_project<-file.path(dir_home, "midrange_filter")
setwd(file.path(dir_project))


dir_data <-file.path(dir_project, "data")
dir_functions <- file.path(dir_project, "functions")

dir_input<- file.path(dir_data,"input")
dir_output<-file.path(dir_data,"output")

dir_rds<-file.path(dir_project,"RDS_files_disease_loci")





#directory for the registries that will be created by batchtools
dir_reg <- file.path(dir_project, "registries")
dir_reg_simulation <- file.path(dir_reg, "simulate_data")
dir_reg_reference <- file.path(dir_reg, "generate_reference")
dir_reg_imputation <- file.path(dir_reg, "phasing_imputation")
dir_reg_association <- file.path(dir_reg, "association")
dir_reg_char <- file.path(dir_reg, "characteristics")




#directory for folder with executables
dir_exec <- file.path(dir_home,"software")



#directory for batchtools configuration file
dir_config<-"./data/input/batchtools.conf.R"

### Directories for required data
#input data for imputation
dir_chr19_raw_vcf <- file.path(dir_input,"chr19.1kg.phase3.v5a.vcf.gz")

dir_chr19_annotated_vcf_original <-file.path(dir_input,"chr19_annotated.vcf.gz")
dir_chr19_prefix_annotated_hap_legend_original<-file.path(dir_input, "chr19_original_annotated")
dir_chr19_annotated_hap_original<-paste(dir_chr19_prefix_annotated_hap_legend_original, ".hap.gz",sep="")
dir_chr19_annotated_legend_original<-paste(dir_chr19_prefix_annotated_hap_legend_original,".legend",sep="")

dir_random_half_txt<- file.path(dir_input, "random_samples.txt")
dir_chr19_annotated_vcf_original_half1<- file.path(dir_input,"chr19_annotated_half1.vcf.gz")
dir_chr19_annotated_vcf_original_half2<- file.path(dir_input,"chr19_annotated_half2.vcf.gz")

dir_chr19_prefix_annotated_hap_legend_original_half1<-paste0(dir_chr19_prefix_annotated_hap_legend_original,"_half1")
dir_chr19_prefix_annotated_hap_legend_original_half2<-paste0(dir_chr19_prefix_annotated_hap_legend_original,"_half2")
dir_chr19_annotated_hap_original_half1<-paste0(dir_chr19_prefix_annotated_hap_legend_original_half1, ".hap.gz")
dir_chr19_annotated_hap_original_half2<-paste0(dir_chr19_prefix_annotated_hap_legend_original_half2, ".hap.gz")
dir_chr19_annotated_legend_original_half1<-paste0(dir_chr19_prefix_annotated_hap_legend_original_half1, ".legend")
dir_chr19_annotated_legend_original_half2<-paste0(dir_chr19_prefix_annotated_hap_legend_original_half2, ".legend")


dir_snp_list_ill_omni5<-file.path(dir_input, "chr19_snp_array_list_ill_omni5.txt")
dir_snp_list_ill_omni25<-file.path(dir_input, "chr19_snp_array_omni25.txt") 
dir_snp_list_ill_omniexpress<-file.path(dir_input, "chr19_snp_array_omniexpress.txt")
  
#list of all SNPs
dir_SNPs_complete<-file.path(dir_input, "snps_complete.txt")


# map files 
dir_map_19_plink <-file.path(dir_input,"plink.chr19.GRCh37.map")
dir_map_19_impute<-file.path(dir_input,"genetic_map_chr19_combined_b37.txt")

#files for magical r2 and liftover
dir_magicalrsq_model<-file.path(dir_input, "model_magicalrsq")
dir_pos_ref_alt<-file.path(dir_input, "reference_pos_alt_ref.txt")
dir_liftover_file<-file.path(dir_input, "liftover_pos_list.RDS")

#executables
dir_exec_beagle <- file.path(dir_exec, "beagle.28Jun21.220.jar")
dir_exec_bcftools <-file.path(dir_exec,"bcftools-1.9","bcftools")
dir_exec_hapgen<- file.path(dir_exec, "hapgen2")
dir_exec_plink<- file.path(dir_exec, "plink")
dir_exec_bref<- file.path(dir_exec, "bref3.19Apr22.7c0.jar")
dir_exec_java <- "java"
dir_exec_imputeAccure<- file.path(dir_exec,"ImputeAccure.py")
dir_exec_python<-"/usr/bin/python"
dir_exec_python_bashscript<- file.path(dir_code, "python_script.sh")


### directiries for output folders
#if these folders don't exist, they will be created

dir_out_imputation <- file.path(dir_output, "imputed")
dir_out_phasing <- file.path(dir_output, "phased")
dir_out_hapgen <- file.path(dir_output, "simulated")
dir_out_reference<-file.path(dir_output, "simulated_reference")
dir_out_association<-file.path(dir_output,"association")
dir_out_char<-file.path(dir_output,"characteristics")
dir_out_magicr<-file.path(dir_output, "magicr")



#source functions folder
function_sources<-list.files(dir_functions, full.names = TRUE)
sapply(function_sources, source)


#global variable
list_chunk_ends_19<-c(1, 17297726, 38851930, 69118783)
