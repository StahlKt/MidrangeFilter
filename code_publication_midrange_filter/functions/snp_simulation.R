#functions for simulating genetic data with HapGen2


#invokes hapgen2 to simulate genetic data
#INPUT
# input.hap: filepath string to genetic data set in .hap format (needs matching .legend file)
# input.legend: filepath string to genetic data set in .legend format (needs matching .hap file)
# exec.hapgen: filepath string to hapgen2 executable
# map.impute: filepath string to genetic map file in format used by IMPUTE2 and similar
# string.risk.alleles: string containing numbers to indicate risk alleles (see hapgen handbook for -dl option)
# number.controls: integer indicating how many control samples should be generated
# number.cases: integer indicating how many case samples should be generated
# output.prefix: filepath string to location and name prefix for output data in hap/legend format
#OUTPUT
#returns output.prefix, but several data sets are generated, one .hap with .cases attached, one .hap with .controls
#and a joint .legend file (this is accounted for in the wrapper file)

simulate_hapgen2<-function(input.hap,
                           input.legend,
                           exec.hapgen,
                           map.impute, 
                           string.risk.alleles,
                           number.controls, number.cases,
                           output.prefix){
  
  
  system2(exec.hapgen, c("-h", input.hap,
                         "-l", input.legend,
                         "-m", map.impute,
                         "-dl", string.risk.alleles,
                         "-n", number.controls, number.cases,
                         "-no_gens_output",
                         "-o", paste(output.prefix,".gz",sep="")))
  return(output.prefix)
}



#wrap for simulation for batchmap
#returns both paths to both simulated cases and controls in vcf in alphabetical oder
#input: hap/legend file with suffix and path, map data in impute format
# string rist alleles: as a string (using paste eg) position, which allele is the disease risk, 
# the heterozygout and homozygout effekt (see hapgen2's documentation for details)


#INPUT
#rep.no: integer functioning as a counter, since several simulated dataset are later merged to a reference data panel
# input.hap.prefix: filepath string prefix to genetic data set in .hap/legend format
# exec.hapgen: filepath string to hapgen2 executable
# exec.bcftools: filepath string to bcftools executable
# map.impute: filepath string to genetic map file in format used by IMPUTE2 and similar
# string.risk.alleles: string containing numbers to indicate risk alleles (see hapgen handbook for -dl option)
# number.controls: integer indicating how many control samples should be generated
# number.cases: integer indicating how many case samples should be generated
#run.signifier: string that indicates current run, needs to be unique for each run
simulation_wrap<-function(rep.no,
                          input.hap.prefix,
                          exec.hapgen,
                          exec.bcftools,
                          map.impute, string.risk.alleles,
                          number.controls, number.cases,
                          run.signifier){

output_prefix<-paste(dir_out_hapgen,"/",run.signifier,"_",rep.no,sep="")




#unzipp input file if necesssary
if(!file.exists(paste(input.hap.prefix, ".legend", sep=""))){
system2("gunzip", paste(input.hap.prefix, ".legend.gz",sep=""))
system2("gunzip", paste(input.hap.prefix, ".hap.gz",sep=""))
}


#to force hapgen to set different seeds for generating the data set
Sys.sleep(2*rep.no)
#simulate by calling this function
simulate_hapgen2(paste(input.hap.prefix, ".hap",sep=""),
                 paste(input.hap.prefix, ".legend",sep=""),
                 exec.hapgen,
                 map.impute, string.risk.alleles,
                 number.controls, number.cases,
                 output_prefix)

#timestamp for log file
print(paste("simulation complete", Sys.time()))
#save prefix paths to cases and controls
prefix_simulated_cases<-paste(output_prefix,".cases",sep="")  
prefix_simulated_controls<-paste(output_prefix,".controls", sep="")

new_controls<-paste(prefix_simulated_controls,".new.hap", sep="")
new_cases<-paste(prefix_simulated_cases,".new.hap", sep="")

control_samples<-paste(prefix_simulated_controls,".sample", sep="")
case_samples<-paste(prefix_simulated_cases,".sample", sep="")

new_controls_samples<-paste(prefix_simulated_controls,".new.sample", sep="")
new_cases_samples<-paste(prefix_simulated_cases,".new.sample", sep="")


system2("gunzip", paste(prefix_simulated_cases,".haps.gz", sep=""))
system2("gunzip", paste(prefix_simulated_controls,".haps.gz", sep=""))

#delete obsolete blank space in hap file after simulation
system2("sed", c("'s/.$//'", paste(prefix_simulated_controls,".haps", sep=""),">",
                 new_controls))
system2("sed", c("'s/.$//'", paste(prefix_simulated_cases,".haps", sep=""),">",
                 new_cases))

#delete second row of sample (obsolete and hinders convertion) and added run number to ids for merging
system2("sed", c("'2d'", paste(prefix_simulated_cases, ".sample", sep=""), ">",
                 new_cases_samples))
system2("awk", c(paste("'{sub(\"$\", \"",paste("R", rep.no, sep=""), "\", $1)}; 1'",sep = ""), new_cases_samples,
                 ">",case_samples))
system2("awk", c(paste("'{sub(\"$\", \"",paste("R", rep.no, sep=""), "\", $2)}; 2'",sep = ""), case_samples,
                 ">",new_cases_samples))

system2("sed", c("'2d'", paste(prefix_simulated_controls, ".sample", sep=""), ">", new_controls_samples))
system2("awk", c(paste("'{sub(\"$\", \"",paste("R", rep.no, sep=""), "\", $1)}; 1'",sep = ""), new_controls_samples,
                 ">",control_samples))
system2("awk", c(paste("'{sub(\"$\", \"", paste("R", rep.no, sep=""), "\", $2)}; 2'",sep = ""), control_samples,
                 ">",new_controls_samples))

print(paste("conversion prep complete", Sys.time()))
#delete old versions
unlink(paste(prefix_simulated_controls,".haps", sep=""))
unlink(paste(prefix_simulated_cases,".haps", sep=""))
unlink(paste(prefix_simulated_controls,".sample", sep=""))
unlink(paste(prefix_simulated_cases,".sample", sep=""))


######convert data from hap/legend to vcf
#initialize paths
dir_vcf_cases<-paste(prefix_simulated_cases,".vcf.gz",sep="")
dir_vcf_controls<-paste(prefix_simulated_controls,".vcf.gz",sep="")

dir_annotated_vcf_cases<-paste(prefix_simulated_cases,"_annotated.vcf.gz",sep="")
dir_annotated_vcf_controls<-paste(prefix_simulated_controls,"_annotated.vcf.gz",sep="")

#convert
hapleg_to_vcf(prefix_simulated_cases, exec.bcftools, dir_vcf_cases)
hapleg_to_vcf(prefix_simulated_controls, exec.bcftools, dir_vcf_controls)
print(paste("convertion complete", Sys.time()))



#index new files
system2(exec.bcftools, c("index -f", dir_vcf_cases))
system2(exec.bcftools, c("index -f", dir_vcf_controls))

#annotate new data, because Hapgen loses the SNP ID
system2(exec.bcftools, c("annotate --set-id '%CHROM\\:%POS\\_%REF\\_%FIRST_ALT'",
                         "-o", dir_annotated_vcf_cases, "-O z", dir_vcf_cases))

system2(exec.bcftools, c("annotate --set-id '%CHROM\\:%POS\\_%REF\\_%FIRST_ALT'",
                         "-o", dir_annotated_vcf_controls, "-O z", dir_vcf_controls))

print(paste("annotation complete", Sys.time()))
#index new files
system2(exec.bcftools, c("index -f", dir_annotated_vcf_cases))
system2(exec.bcftools, c("index -f", dir_annotated_vcf_controls))


#delete unannotated files
unlink(paste(dir_vcf_cases, sep=""))
unlink(paste(dir_vcf_cases, ".csi", sep=""))
unlink(paste(dir_vcf_controls, sep=""))
unlink(paste(dir_vcf_controls, ".csi", sep=""))



#delete data in old hap/legend format

unlink(new_controls_samples)
unlink(new_cases_samples)
unlink(new_cases)
unlink(new_controls)
unlink(paste(output_prefix,".legend", sep=""))

#returns both cases and controls in vcf
return(c(dir_annotated_vcf_cases, dir_annotated_vcf_controls))

}
