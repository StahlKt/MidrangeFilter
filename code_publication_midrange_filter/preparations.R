#run preparations:

library(batchtools)
library(data.table)



#source all relevant file paths and functions saved in file_paths.R

source(file.path(getwd(), "file_paths.R"))


###### Annotation Data set #####

#annotates test data, which is the basis for simulation

if(!file.exists(dir_chr19_annotated_vcf_original)){

annotate_data(dir_chr19_raw_vcf, dir_exec_bcftools,dir_chr19_annotated_vcf_original)
}

dir_chr19_annotated_vcf_original
###### Convert from vcf to hap/legend ####

#convert original data set into hap/legend format for HapGen2 

if(file.exists(paste(dir_chr19_prefix_annotated_hap_legend_original, ".legend", sep=""))){

vcf_to_hapleg(dir_chr19_annotated_vcf_original, dir_exec_bcftools,dir_chr19_prefix_annotated_hap_legend_original)
}


dir_chr19_prefix_annotated_hap_legend_original



#convert two halves of the dataset
vcf_to_hapleg(dir_chr19_annotated_vcf_original_half1,
              dir_exec_bcftools,
              dir_chr19_prefix_annotated_hap_legend_original_half1)

vcf_to_hapleg(dir_chr19_annotated_vcf_original_half2,
              dir_exec_bcftools,
              dir_chr19_prefix_annotated_hap_legend_original_half2)


#unzip legend files manually, if you have to
file.exists(dir_chr19_annotated_hap_original_half2)
file.exists(dir_chr19_annotated_legend_original_half2)
