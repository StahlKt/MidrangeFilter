#example script of running the functions for the Midrange Filter
#for specifics on the input parameter for each function, 
#please see the comments in functionsMidrangeFilter.R

#load data.table package
library(data.table)

#Please set the working directory and/or edit the paths to the files and scripts to fit your set-up

# source the script with the functions
source("functionsMidrangeFilter.R")

#and load the data sets in the example file folder:
#this is a data table with imputed SNPs
imputed_snps<-readRDS("example_data/imputed_snps.RDS")
#check properties
names(imputed_snps)
summary(imputed_snps)


#this is a data table listing SNPs which were not imputed but typed by the SNP array
snps_scaffold<-readRDS("example_data/snps_scaffold.RDS")
#check properties
names(snps_scaffold)
summary(snps_scaffold)


#there is an additional file in the example_data folder which reveals the association results,
#if the data set was completely sequenced and not imputed. This is added for comparison. 

#determine the spikes in the data set with the determine_spikes function from the other script.
#the example files have names matching the functions expectation,
#but all column names may be specified as an input parameter
spike_table<-determine_spikes(imputed_snps)
#result table hold one line per detected SNPs
spike_table

#check best guess spikes to determine in which case they fall into:
#are they present in dosage as well or ar they specific to best guess
spike_table_validated<-validate_best_guess_spikes(spike_table)
#this new table has an extra column that specifies the type of spike
spike_table_validated

#remove spikes that contain significant SNPs, that were not imputed
#these spikes are presumed to be true spikes and do not need further filtering
spike_table_validated_checked<-check_scaffold(spike_table_validated, snps_scaffold, imputed_snps)
spike_table_validated_checked

#see spikes, which contain typed SNPs, which are not discarded by the MidrangeFilter as a default.
spikes_typed<-check_scaffold(spike_table_validated, snps_scaffold, imputed_snps,remove.typed.spikes=FALSE)
spikes_typed


#run the Midrange Filter to keep or discard imputed spikes:

#decisions listed for each SNP
result<-midrange_filter(imputed_snps, spike_table_validated_checked, quality.col = "BR2_MIN")
result

#decisions listed for each spike
result_condensed<-midrange_filter(imputed_snps, spike_table_validated_checked, quality.col = "BR2_MIN", condensed.return.table = TRUE)
result_condensed

