library(data.table)



# source the docs with the functions
#and load the toy data set
#inluding the list of positions in the scaffold

imputed_snps<-readRDS("example_data/imputed_snps.RDS")
snps_scaffold<-readRDS("example_data/snps_scaffold.RDS")

names(imputed_snps)
summary(imputed_snps)

names(snps_scaffold)
summary(snps_scaffold)


#determine the spikes in the data set
spike_table<-determine_spikes(imputed_snps)
spike_table

#validate best guess spikes in dosage to determine the type of spike:
#best guess specific or not
spike_table_validated<-validate_bg(spike_table)
spike_table_validated

#remove spikes that contain significant SNPs, that were not imputed
spike_table_validated_checked<-check_scaffold(spike_table_validated, snps_scaffold, imputed_snps)
spike_table_validated_checked

#run the Midrange Filter to keep or discard spikes
#decisions listed for each SNP
result<-midrange_filter(test_snp_table, spike_table_validated_checked, quality.col = "BR2_MIN")
result
#decisions listed for each spike
result_condensed<-midrange_filter(test_snp_table, spike_table_validated_checked, quality.col = "BR2_MIN", condensed.return.table = TRUE)
result_condensed
