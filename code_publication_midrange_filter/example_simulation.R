#example on how to run the simulation with the functions specified 
#in full_simulation_applyable.R and the functions folder.
#file_paths.R and configuration files need to be specified to fit the current system
run_name_stem<-"set1_ipr_b_o5" 
disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
run_seed<-23102024

#simulation
batch_simulation(1,
                 dir.ref.basis = dir_chr19_prefix_annotated_hap_legend_seppop_refdata,
                 dir.testdata.basis = dir_chr19_prefix_annotated_hap_legend_seppop_testdata)

sapply(1:length(disease_list),
       batch_simulation, 
       dir.ref.basis = dir_chr19_prefix_annotated_hap_legend_seppop_refdata,
       dir.testdata.basis = dir_chr19_prefix_annotated_hap_legend_seppop_testdata)


#construct reference panel
batch_reference_panel(run_name_stem)


#imputation
batch_phasing_imputation_no_chunks_ver2(1, dir_snp_list_ill_omni5)
sapply(1:length(disease_list), batch_phasing_imputation_no_chunks_ver2, which.snp.array= dir_snp_list_ill_omni5)



#association step 1
batch_association_no_chunks_ver3(1, n.parts=3)
sapply(1:length(disease_list), batch_association_no_chunks_ver3, n.parts=3)
#association step 2
batch_association_follow_ver2(1)
sapply(2:length(disease_list), batch_association_follow_ver2)


#characteristics
batch_char_one_job(run_name_stem)
