#there are hard-coded variables needed, such as shown in overview_simulations.R
#run.signifier is used as the first part of all data sets involved in one run and should be chosen distinctively
#also the run.signifier will be set in the following functions with run_name and should be matching the previous jobs of one iteration


#the functions are meant to be run one after the other with one disease loci set with the same simulation settings. 
#if run.nr is needed, then they can be used with sapply, e.g.

#simulation with set1 and low cases, deletion with illumina omni 5
#run_name_stem<-"set1_lc_o5" 
#disease_list<-readRDS(file.path(dir_rds, "set1_disease_loci.RDS"))
#run_seed<-10101010

#sapply(1:length(disease_list), batch_simulations, number.cases.b=333)

#if run.name.stem is needed, use the function directly, such as

#batch_reference_panel(run_name_stem)


#simulate data sets with hapgen adn batchtools
#run.nr is the index of the specific disease loci and effect sizes specified in disease_list
#seting flag.small.ref to 1 will reduce the individuals in the reference panel to 5*number.reference instead of 10*number.reference
#number.cases.b and number.controls.b are number of individuals in the case and control group respectively

#the simulated data set is needed until the association is done
batch_simulation<-function(run.nr, flag.small.ref=0, 
                           number.cases.b=1000, number.controls.b=1000, number.reference=1000){
  
  ifelse(flag.small.ref,
         number_jobs<-1:6,
         number_jobs<-1:11)
  
  #specify run_name as signifier
  disease_allels<-disease_list[run.nr]
  run_name<-paste(run_name_stem,
                  run.nr,
                  sep = "_")
  
  #check for existing data sets
  bol_ref<-file.exists(sapply(c(0,1),get_reference_panel_path,run.name=run_name))
  bol_test<-file.exists(sapply(c(0,1),get_results_imputation_path,run.name=run_name, imputation.flag=0))
  
  if(all(bol_ref, bol_test)){
    print("test data and reference panel already exists. return")
    return()
  }
  
  
  if(!file.exists(get_registries(run_name)[1])){
    
    
    #create registry for each simulation or clear existing registry
    reg_simulation <- makeRegistry(
      file.dir = get_registries(run_name)[1],
      make.default = FALSE,
      source = file.path(dir_source, "file_paths.R"),
      seed = run_seed
    )
    
    reg_simulation$cluster.functions = makeClusterFunctionsSlurm(template=".batchtools.SCC.tmpl")
  } else{
    
    reg_simulation<-loadRegistry(get_registries(run_name)[1], writeable = TRUE, make.default = FALSE)
    reg_simulation$cluster.functions = makeClusterFunctionsSlurm(template=".batchtools.SCC.tmpl")
    clearRegistry(reg=reg_simulation)
    
  }
  
  #simulation_wrap is the function that simulates data sets
  #first is the test data set, the rest is merged to the reference panel in the following functions
  #reference and test data are based on different parts of the 1000 Genomes Project data set

  batchMap(
    simulation_wrap,
    rep.no=number_jobs,
    number.controls = c(number.controls.b, rep(number.reference,(length(number_jobs)-1))), 
    number.cases = c(number.cases.b, rep(number.reference,(length(number_jobs)-1))),
    input.hap.prefix = c(dir_chr19_prefix_annotated_hap_legend_original_half1,
                         rep(dir_chr19_prefix_annotated_hap_legend_original_half2, (length(number_jobs)-1))),
    reg=reg_simulation,
    more.args = list(
      exec.hapgen = dir_exec_hapgen,
      exec.bcftools = dir_exec_bcftools,
      map.impute = dir_map_19_impute,
      string.risk.alleles = disease_allels,
      run.signifier = run_name
    )
  )
  
  
  
  
  if(any(bol_test, bol_ref)){
    if(bol_ref){
      number_jobs<-1
    } else{
      find_cases<-sapply(number_jobs, get_simulation_result_path,
                         case.flag=1, run.name=run_name)
      find_controls<-sapply(number_jobs, get_simulation_result_path,
                            case.flag=0, run.name=run_name)
      
      cases_simulated<-which(!file.exists(find_cases))
      controls_simulated<-which(!file.exists(find_controls))
      
      number_jobs<-unique(c(cases_simulated, controls_simulated))
      
      
    }
    
  }
  
  
  #submit jobs
  invisible(sapply(number_jobs, submit_jobs_wrap, registry=reg_simulation, mem=80000,
                   wall.time=60, partition.name="medium"))
  
  return()
}





#generate refence panel
#this creates one registry for all simulations of the same setting (run_name)
#jobs.not.ready is set, if some simulations of the setting need to be held back by their run number
#flag.small.ref as above
#reference panels are needed until the imputation step is done
batch_reference_panel<-function(run.name.stem,
                                number.runs=length(disease_list),
                                flag.small.ref=0,
                                jobs.not.ready=0){
  
  ifelse(flag.small.ref,number_jobs<-list(2:6),number_jobs<-list(2:11))
  
  list_runs<-paste0(run.name.stem,"_", 1:number.runs)
  
  
  reg_reference <- makeRegistry(
    file.dir = get_registries(run.name.stem)[2], 
    make.default = FALSE,
    source = file.path(dir_source, "file_paths.R"),
  )
  
  
  reg_reference$cluster.functions = makeClusterFunctionsSlurm(template=".batchtools.SCC.tmpl")
  
  ids_reference<- batchMap(
    reference_wrap_no_load, 
    run.signifier=list_runs,
    reg=reg_reference,
    more.args = list(
      exec.bcftools = dir_exec_bcftools,
      job.ids.merge = number_jobs
    )
  )
  ids_to_run<-1:number.runs
  
  ids_to_run<-ids_to_run[!file.exists(paste0(dir_out_reference, "/", list_runs, "_cases_reference.bref3"))]
  
  if(!jobs.not.ready==0) ids_to_run<-ids_to_run[ids_to_run%notin% jobs.not.ready]
  
  
  invisible(sapply(ids_to_run, submit_jobs_wrap, registry=reg_reference, mem=80000,
                   wall.time=150, partition.name="medium"))
  
}




#phasing and imputation, 
#which.snp.array expects path to a list with SNPs that remain in the data.set
#if flag.random.deletion.b is set to 1, SNPs are not deleted according to the list, but accordning to a fraction. 
#the fraction is hardcoded to 20% remaining in this function, but can be adjusted easily in two lines
#these two lines are indicated by a comment or search "fraction.to.keep"
#the imputed data set is needed until the gathering of characteristics is done
batch_phasing_imputation_no_chunks_ver2<-function(run.nr, which.snp.array=NULL, flag.random.deletion.b=0){
  
  
  
  run_name<-paste(run_name_stem,
                  run.nr,
                  sep = "_")
  
  
  
  bol_case<-file.exists(get_results_imputation_path(run_name, 1,1))
  bol_control<-file.exists(get_results_imputation_path(run_name, 0,1))
  
  if(all(bol_case, bol_control)){
    print("imputed files exist already")
    return()
  }
   
  #make experiment registry for phasing and Imputation
  reg_imputation <- makeExperimentRegistry(
    file.dir = get_registries(run_name)[3], 
    make.default = FALSE, 
    seed = run_seed,
    source = file.path(dir_source,"file_paths.R")
  )
  reg_imputation$cluster.functions = makeClusterFunctionsSlurm(template=".batchtools.SCC.tmpl")
  
  
  
  ##### define problem as deleting data and phasing #####
  
  addProblem("phasing_cases",
             data = get_simulation_result_path(run_name,1,1),
             fun = problem_wrap,
             seed = run_seed,
             cache = TRUE, reg = reg_imputation)
  
  addProblem("phasing_controls",
             data = get_simulation_result_path(run_name,0,1),
             fun = problem_wrap,
             seed = run_seed,
             cache = TRUE, reg = reg_imputation)
  
  
  arguments_problem <- list(
    phasing_cases = data.table(
      exec.bcftools = dir_exec_bcftools,
      exec.plink = dir_exec_plink,
      exec.beagle = dir_exec_beagle,
      map.plink = dir_map_19_plink,
      flag.random.deletion = flag.random.deletion.b,
      #change 0.2 to the fraction to keep, if 20% is not desired
      fraction.to.keep = 0.2,
      #snp arrray
      list.snps.to.keep = which.snp.array,
      run.signifier = run_name
      
    ),
    phasing_controls = data.table(
      exec.bcftools = dir_exec_bcftools,
      exec.plink = dir_exec_plink,
      exec.beagle = dir_exec_beagle,
      map.plink = dir_map_19_plink,
      flag.random.deletion = flag.random.deletion.b,
      #change 0.2 to the fraction to keep, if 20% is not desired
      fraction.to.keep = 0.2,
      #snp array
      list.snps.to.keep = which.snp.array,
      run.signifier = run_name
    )
  )
  
  
  
  ######## define algorithm as imputation #######
  
  
  
  addAlgorithm("imputation_cases",
               fun = beagle_impute_no_chunks,
               reg = reg_imputation)
  addAlgorithm("imputation_controls",
               fun= beagle_impute_no_chunks,
               reg = reg_imputation)
  
  
  arguments_algorithm <-list(
    imputation_cases = data.table(
      flag.cases.algorithm = 1,
      exec.beagle = dir_exec_beagle,
      map.plink = dir_map_19_plink,
      reference.bref3 = get_reference_panel_path(run_name, 1),
      exec.bcftools = dir_exec_bcftools,
      run.signifier = run_name
    ),
    imputation_controls = data.table(
      flag.cases.algorithm = 0,
      exec.beagle = dir_exec_beagle,
      map.plink = dir_map_19_plink,
      reference.bref3 = get_reference_panel_path(run_name, 0),
      exec.bcftools = dir_exec_bcftools,
      run.signifier = run_name  
    )
  )
  
  
  
  ###### combine problem and algorithm to Experiments #####
  
  #Add Experiments
  addExperiments(prob.designs = arguments_problem,
                 algo.designs = arguments_algorithm,
                 reg = reg_imputation)
  
  
  
  #remove crossover experiments:
  job_table_imp<-getJobTable(reg=reg_imputation)
  
  ids_to_remove<-job_table_imp[(problem=="phasing_cases" & algorithm=="imputation_controls") |
                                 (problem=="phasing_controls" & algorithm=="imputation_cases"), job.id]
  
  removeExperiments(ids_to_remove, reg=reg_imputation)
  

  ###### submit jobs ######
  invisible(sapply(findJobs(reg=reg_imputation)$job.id, submit_jobs_wrap, registry=reg_imputation, mem=80000,
                   wall.time=120, partition.name="medium"))
}


#run association on imputed and simulated data. 
#this uses another global variable: list_chunk_ends_19<-c(1, 17297726, 38851930, 69118783)
#defined in file_paths, because reading in too many imputed SNPs at once was not time efficient
#as is, the complete simulated data set does not get chunked and the imputed set is divided into 3 parts

#the association registries are needed until the association_follow function is done.
batch_association_no_chunks_ver3<-function(run.nr, n.parts=3){
  
  run_name<-paste(run_name_stem,
                  run.nr,
                  sep = "_")
  
  
  reg_asso<- makeRegistry(
    file.dir = get_registries(run_name)[4],
    make.default = FALSE,
    source = file.path(dir_source,"file_paths.R")
  )
  
  reg_asso$cluster.functions = makeClusterFunctionsSlurm(template=".batchtools.SCC.tmpl")
  
  
  ids_asso<-batchMap(
    fun=function(sim.imp.flag, exec.bcftools, run.signifier){

      #for imputed
      if(sim.imp.flag>0){
        
        #read in matching cases and controls
        result_case<-get_results_imputation_path(run.signifier,1,1)
        
        
        result_control<-get_results_imputation_path(run.signifier,0,1)
        
        
        
        if(!file.exists(paste0(result_case, ".csi"))){
          system2(exec.bcftools, c("index -f", result_case))
        }
        
        if(!file.exists(paste0(result_control, ".csi"))){
          system2(exec.bcftools, c("index -f", result_control))
        }
        
        
        
        pos_start<-list_chunk_ends_19[sim.imp.flag]
        pos_end<-list_chunk_ends_19[(sim.imp.flag+1)]
        
        
        
        combined_table<-rbind(read_table_imp_snps(result_case,exec.bcftools,1, 
                                                  pos_start, pos_end),
                              read_table_imp_snps(result_control,exec.bcftools,0, 
                                                  pos_start, pos_end))
        
        
        
        
        list_pos<-combined_table[,unique(POS)]
        
        #run association
        return(rbindlist(lapply(list_pos, 
                                snp_association,
                                input.table=combined_table,
                                flag.imputed=1,
                                run.signifier=run.signifier)
        ))
      }
      
      #for simulated, not imputed
      if(sim.imp.flag==0){
        simulated_cases<-get_results_imputation_path(run.signifier,1,0)
        simulated_controls<-get_results_imputation_path(run.signifier,0,0)
        
        
        
        if(!file.exists(paste0(simulated_cases, ".csi"))){
          system2(exec.bcftools, c("index -f", simulated_cases))
        }
        
        if(!file.exists(paste0(simulated_controls, ".csi"))){
          system2(exec.bcftools, c("index -f", simulated_controls))
        }
        
        combined_table<-rbind(read_table_sim_snps(simulated_cases,1),
                              read_table_sim_snps(simulated_controls,0))
        
        
        #run association
        return(rbindlist(lapply(combined_table[,unique(POS)], 
                                snp_association,
                                input.table=combined_table,
                                flag.imputed=0,
                                run.signifier=run.signifier)
        )
        ) 
      }
    },
    sim.imp.flag=((1:(n.parts+1))-1), 
    
    more.args = list(
      exec.bcftools = dir_exec_bcftools,
      run.signifier = run_name
    ),
    reg=reg_asso
  )
  
  #send off jobs
  invisible(sapply(1:(n.parts+1), submit_jobs_wrap, registry=reg_asso, mem=80000,
                   wall.time=180, partition.name="medium"))
}

#after the association with batchtools is complete, this is run localy to gather the association results and saves them
#in the association output folder. 
batch_association_follow_ver2<-function(run.nr){
  
  run_name<-paste(run_name_stem,
                  run.nr,
                  sep = "_")
  
  
  reg_asso<-loadRegistry(get_registries(run_name)[4],
                         writeable = TRUE, make.default = FALSE)
  
  reg_asso$cluster.functions = makeClusterFunctionsSlurm(template=".batchtools.SCC.tmpl")
  
  
  
  table_p_values<-rbindlist(reduceResultsList(ids = findDone(reg=reg_asso)$job.id,
                                              reg = reg_asso),
                            idcol = "job.id")
  table_p_values[,job.id:=NULL]
  table_p_values<-unique(table_p_values)
  
  saveRDS(table_p_values, file.path(dir_out_association, paste(run_name,"_pv.RDS", sep="")))
  
  return()
  
}


#gathers characteristics such as MAF and imputation quality measures
#tables with characteristics of significant SNPs are saved as RDS files
#as the reference panel, one registry for the whole setting.
batch_char_no_chunks_ver4_one_job<-function(run.name.stem, number.runs=length(disease_list)){
  
  
  
  reg_char <- makeRegistry(
    file.dir = get_registries(run.name.stem)[5], 
    make.default = FALSE,
    source = file.path(dir_source, "file_paths.R"),
  )
  
  
  reg_char$cluster.functions = makeClusterFunctionsSlurm(template=".batchtools.SCC.tmpl")
  
  #run function to get characteristics table
  ids_char<-batchMap(function(run.name.stem, number.runs){
    sapply(1:number.runs,
           chartable_no_chunks_one_job, 
           run.name.stem=run.name.stem)
    },
    run.name.stem=run.name.stem,
    number.runs=number.runs,
    reg=reg_char
  )
  
  
  #send off jobs
  invisible(sapply(ids_char$job.id, submit_jobs_wrap, registry=reg_char, mem=80000,
                   wall.time=120, partition.name="medium"))
  
  return()
  
}






