#functions to use ImputeAccure

## get the content from imputed vcf files, melt them and write them into a file fir ImputeAccure
#INPUT
#input.data.cases: filepath string to inputed Case data in vcf.gz
#input.data.controls: filepath string to inputed control data in vcf.gz
#exec.bcftools: filepath string to bcftools executable
#output.file: filepath string to desired location and input file for the ImputeAccure
#OUTPUT 
#returns filepath string to newly written file (output.file)
set_IA_input_data<-function(input.data, 
                            exec.bcftools, output.file, dir.list.char.pos){

 
  

  case_table <- fread(
    cmd = sprintf(
      fmt = "%s %s %s",
      exec.bcftools, paste0( "query -f \'[%SAMPLE,%ID,%POS,%GP\\n]\' -R ", dir.list.char.pos), input.data),
    sep = ","
  )

 

  #set names for compability
  setnames(case_table,c("SAMPLE","SNP","POS","GP_0","GP_1", "GP_2"))
  
  

  case_table[, GP_joint := paste(GP_0, GP_1, GP_2), by=.(SAMPLE, SNP, POS)]
  
  
  case_table[, GP_0:= NULL]
  case_table[, GP_1:= NULL]
  case_table[, GP_2:= NULL]
  
  case_table<-dcast(case_table, SNP + POS ~ SAMPLE, value.var = "GP_joint")
  

  write.table(as.data.frame(case_table), file=output.file, col.names = FALSE, 
              quote = FALSE)
  

  return(output.file)
}




#
#INPUT
#input.data.cases: filepath string to inputed Case data in vcf.gz
#input.data.controls: filepath string to inputed control data in vcf.gz
#exec.imputeaccure: filepath string to imputeAccure executable
#exec.python: filepath string to python executable
#list.pos: vector of integers of positions to inlcude
#OUTPUT 
#returns data.table with SNP, the ImputeAccure Measures and R^2 as accuracy measures for imputation
get_table_imputeAccure<-function(input.data.imputed.cases, input.data.imputed.controls,
                                 exec.bcftools,
                                 exec.python.script,
                                 dir.list.char.pos){
  
  

 
  #input and output for imputeAccure
 
  dir_output_IA_cases<-gsub(".vcf.gz", "_IA.accuracy", input.data.imputed.cases)
  dir_input_IA_cases<-gsub(".vcf.gz", "_IA.imputed", input.data.imputed.cases)
 
  dir_output_IA_controls<-gsub(".vcf.gz","_IA.accuracy", input.data.imputed.controls)
  dir_input_IA_controls<-gsub(".vcf.gz", "_IA.imputed", input.data.imputed.controls)
  
  #write table
 set_IA_input_data(input.data.imputed.cases,
                    exec.bcftools, dir_input_IA_cases, dir.list.char.pos)             
  
  set_IA_input_data(input.data.imputed.controls,
                    exec.bcftools, dir_input_IA_controls, dir.list.char.pos)  
  
  
  system2("bash", c(exec.python.script, dir_input_IA_cases))
  system2("bash", c(exec.python.script, dir_input_IA_controls))
  
  
  # system2("module", "load python")
  # system2(exec.python, c(exec.imputeAccure, "-i", dir_input_IA, "-l 3 -n No,SNP,POS"))
 
 
  # unlink(dir_input_IA)

  acc_table_cases<- fread(dir_output_IA_cases, drop = c(1, 4, 5))
  
  acc_table_controls<- fread(dir_output_IA_controls, drop = c(1, 4, 5))
  
  acc_table<-rbind(acc_table_controls[,CASE:=0], acc_table_cases[, CASE:=1])
  setnames(acc_table, "position", "POS")
  
  #clean_up IA files
  
  unlink(dir_input_IA_cases)
  unlink(dir_input_IA_controls)
  
  if("MAF" %in% names(acc_table)) acc_table[,MAF:=NULL]
  
  return(acc_table)
}


