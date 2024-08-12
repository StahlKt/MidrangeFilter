#read in the results so we can run the association


#reads in vcf file into a data.table for imputed data for association
#INPUT
#data.imputed: filepath string to imputed data set in vcf.gz fornat
#flag.cases. boolean or 0/1 to indicate if data.imputed is a file of cases (1) or controls(0)
#OUTPUT
#returns data.table with clumns SAMPLE, SNP, POS, GP_0, GP_1, GP_2,
# BEST_GUESS, MAX_GP, DOSAGE, STATUS
read_table_imp_snps<- function(data.imputed, exec.bcftools, flag.cases, 
                               pos.start=0, pos.end=0){

  if(file.exists(data.imputed)){
  
  if((pos.start+pos.end)==0){
  #read in table with samples,chrom pos, ref, alt, GP(genotype probabilities in 3 columns)
  imp_table <- fread(
    cmd = sprintf(
      fmt = "%s %s %s", exec.bcftools, "query -f \'[%SAMPLE,%POS,%GP\\n]\'",
      data.imputed), 
    sep = ","
  )
  }
    else{
      imp_table <- fread(
        cmd = sprintf(
          fmt = "%s %s %s",
          exec.bcftools, paste0("query -f \'[%SAMPLE,%POS,%GP\\n]\' -r 19:", pos.start, 
                                "-", pos.end), data.imputed), 
        sep = ","
      )
    }
  #set names for compability
  setnames(imp_table,c("SAMPLE","POS","0","1", "2"))
  imp_table[,POS:=as.numeric(POS)]
  #add best_guess und max_GP
  imp_table[, BEST_GUESS := as.integer(colnames(.SD)[max.col(.SD, ties.method = "first")]), .SDcols = c("0", "1", "2")]
  #imp_table[, MAX_GP := get(colnames(.SD)[max.col(.SD, ties.method = "first")]), .SDcols = c("0", "1", "2"),  by=.(SAMPLE, SNP,POS)]
  setnames(imp_table, c("0", "1", "2"), c("GP_0", "GP_1", "GP_2"))
  imp_table[, DOSAGE := GP_1 + 2*GP_2, by = .(SAMPLE, POS)]
  imp_table[, CASE:= flag.cases]
  return(imp_table[,.(SAMPLE, POS, BEST_GUESS, DOSAGE, CASE, GP_1, GP_2)])
  
  
  
  } else {
    stop("imputed data not found!")
  }
}




read_table_char_snps<- function(data.imputed, exec.bcftools, flag.cases, 
                               dir.list.char.pos){
  

  
  if(file.exists(data.imputed)){
   
      imp_table <- fread(
        cmd = sprintf(
          fmt = "%s %s %s",
          exec.bcftools, paste0("query -f \'[%SAMPLE,%POS,%GP\\n]\' -R ", dir.list.char.pos), data.imputed), 
        sep = ","
      )
    
    #set names for compability
    setnames(imp_table,c("SAMPLE","POS","0","1", "2"))
    imp_table[,POS:=as.numeric(POS)]
    #add best_guess und max_GP
    imp_table[, BEST_GUESS := as.integer(colnames(.SD)[max.col(.SD, ties.method = "first")]), .SDcols = c("0", "1", "2")]
    #imp_table[, MAX_GP := get(colnames(.SD)[max.col(.SD, ties.method = "first")]), .SDcols = c("0", "1", "2"),  by=.(SAMPLE, SNP,POS)]
    setnames(imp_table, c("0", "1", "2"), c("GP_0", "GP_1", "GP_2"))
    imp_table[, DOSAGE := GP_1 + 2*GP_2, by = .(SAMPLE, POS)]
    imp_table[, CASE:= flag.cases]
    return(imp_table[,.(SAMPLE, POS, BEST_GUESS, DOSAGE, CASE, GP_1, GP_2)])
    
    
    
  } else {
    stop("imputed data not found!")
  }
}





#reads in vcf file into a data.table for simulated data for association
#INPUT
#data.simulated: filepath string to simulated data set in vcf.gz fornat
#flag.cases. boolean or 0/1 to indicate if data.imputed is a file of cases (1) or controls(0)
#OUTPUT
#returns data.table with clumns SAMPLE, SNP, POS, GT, STATUS
read_table_sim_snps<- function(data.simulated, flag.cases){
  #read in simulated data (no imputation) for comparison
  if(file.exists(data.simulated)){
    sim_table<- fread(
          data.simulated)
    
    cols_to_keep<-c(which(names(sim_table) %in% "POS"), grep("id", names(sim_table)))
    sim_table<-sim_table[,..cols_to_keep]
    
    sim_table<- melt(sim_table, id.vars="POS", value.name="GT",
            variable.name = "SAMPLE")
    #set names
    
    
    #swap GT for coded GT
    sim_table[, GT := as.numeric(factor(GT, levels = c("0/0", "0/1", "1/0", "1/1", "0|0", "0|1", "1|0", "1|1"), 
                                        labels = c(0, 1, 1, 2, 0, 1, 1, 2)))-1]
    
    sim_table[, CASE:= flag.cases]
    
    return(sim_table)
    
    
  } else {
    stop("simulated data not found!")
  }
}




##reads in vcf file into a data.table for contructing LD in imputed data
#INPUT
#data.imputed: filepath string to imputed data set in vcf.gz fornat
#OUTPUT
#returns data.table with columns SAMPLE, SNP, POS, GP_0, GP_1, GT_BG
read_table_imp_LD<-function(data.imputed, exec.bcftools){
imp_table<- fread(
  cmd = sprintf(
    fmt = "%s %s %s",
    exec.bcftools, "query -f \'[%SAMPLE,%ID,%POS,%GT,%GP\\n]\'",
    data.imputed), 
  sep = ","
)

setnames(imp_table,c("SAMPLE","SNP","POS","GT_BG","0","1", "2"))
imp_table[, BEST_GUESS := as.integer(colnames(.SD)[max.col(.SD, ties.method = "first")]),
          .SDcols = c("0", "1", "2")]
#set names for compability

setnames(imp_table,c("0","1","2"), c("GP_0","GP_1","GP_2"))

imp_table[, GP_2:=NULL]

return(imp_table)
}