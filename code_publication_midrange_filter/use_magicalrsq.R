source(file.path(getwd(),"file_paths.R"))


new_pos<-readRDS(dir_liftover_file)

#list of all runs with significant SNPs
list_runs<-unique(all_sig[, .(NAME, RUN)])
#file.paths to characteristics files
list_runs[,PATH:=paste0(dir_out_char,"/",NAME,"_",RUN,"_acc_char.RDS")]
list_runs<-list_runs[file.exists(PATH)==TRUE]

models<-c("CF_2k_1000G", "CF_2k_TOPMed",
          "UKB_AFR_1k_1000G",
          "UKB_AFR_1k_TOPMed")

sapply(1:401,get_all_mr2)



list_files_magicr<-list.files(dir_out_magicr, full.names = TRUE)
list_CF_1KG<-list_files_magicr[which(str_detect(list_files_magicr, "CF_2k_1000G"))]
list_UKB_1KG<-list_files_magicr[which(str_detect(list_files_magicr, "UKB_AFR_1k_1000G"))]
list_CF_TM<-list_files_magicr[which(str_detect(list_files_magicr, "CF_2k_TOPMed"))]
list_UKB_TM<-list_files_magicr[which(str_detect(list_files_magicr, "UKB_AFR_1k_TOPMed"))]


table_CF_1KG<-rbindlist(lapply(list_CF_1KG, readRDS))[,c(2,6:11)]
table_CF_1KG[, MODEL:="CF_1KG"]
table_UKB_1KG<-rbindlist(lapply(list_UKB_1KG, readRDS))[,c(2,6:11)]
table_UKB_1KG[, MODEL:="UKB_1KG"]
table_CF_TM<-rbindlist(lapply(list_CF_TM, readRDS))[,c(2,6:11)]
table_CF_TM[, MODEL:="CF_TM"]
table_UKB_TM<-rbindlist(lapply(list_UKB_TM, readRDS))[,c(2,6:11)]
table_UKB_TM[, MODEL:="UKB_TM"]

mr2_table<-rbind(table_CF_1KG, table_UKB_1KG, table_CF_TM, table_UKB_TM)

mr2_table<-mr2_table[,.(MIN_MRSQ=min(MagicalRsq), MEAN_MRSQ=mean(MagicalRsq)), by=.(POS, RUN, NAME, AF_FLAG, MODEL)]

setnames(mr2_table, "POS", "NEW_POS")

new_pos<-readRDS(dir_liftover_file)
mr2_table<-new_pos[mr2_table, on="NEW_POS"]

setnames(mr2_table, "OLD_POS", "POS")
mr2_table[,NEW_POS:=NULL]
