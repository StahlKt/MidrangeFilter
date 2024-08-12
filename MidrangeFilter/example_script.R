library(data.table)



# source the docs with the functions
#and load the toy data set
#inluding the list of positions in the scaffold

test_snp_table
snps_scaffold

names(test_snp_table)
summary(test_snp_table)

names(snps_scaffold)
summary(snps_scaffold)


#determine the spikes in the data set
spike_table<-determine_spikes(test_snp_table)
spike_table

#validate best guess spikes in dosage to determine the type of spike:
#best guess specific or not
spike_table_validated<-validate_bg(spike_table)
spike_table_validated

#remove spikes that contain significant SNPs, that were not imputed
spike_table_validated_checked<-check_scaffold(spike_table_validated, snps_scaffold, test_snp_table)
spike_table_validated_checked

#run the Midrange Filter to keep or discard spikes
#decisions listed for each SNP
result<-midrange_filter(test_snp_table, spike_table_validated_checked, quality.col = "IMP_INFO_MIN")
result
#decisions listed for the spikes
result_condensed<-midrange_filter(test_snp_table, spike_table_validated_checked, quality.col = "IMP_INFO_MIN", condensed.return.table = TRUE)
result_condensed
