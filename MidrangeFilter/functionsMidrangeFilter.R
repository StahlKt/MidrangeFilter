#functions to run the MidrangeFilter as a post imputation quality method

############ function to sort SNPs into spikes######################
#the input parameter snp.table needs to have at least one columns containing positions
#three additional (optional) columns:
# chromosome  
#the genotype format which the association is based on (dosage or best guess)
#and p-value
#to use the MidrangeFilter in a later function, another column containing imputation quality measures 
#such as IMPUTE info or Beagle R^2 is needed, but could be added after this function.
#if p-value is missing, all SNPs are considered equally significant
#if the chromosome column is missing, the function assumes all SNPs are on the same chromosome
#if the format column is missing, the function assumes the dosage format
#all column names can be specified by the respective inputparameter, if they differ from the set default names
#theshold.p.significance sets the significant threshold for the association analysis
#threshold.p.grouping sets a threshold on the p-value for the aggregation algorithm
#to help determining the location of spikes.
#We recommend to set the grouping threshold less stringend than the signifance threshold
#the max.distance parameter determines the maximum distance between SNPs assigned to the same spike
#the function returns an output table with one line per found spike
determine_spikes<-function(snp.table, 
                           pos.col="POS", chrom.col="CHROM", pvalue.col="P_VALUE", format.col="FORMAT", 
                           threshold.p.grouping=5e-7, threshold.p.significance=5e-8, max.distance=2e6){
  
  #check which columns are supplied in the data table
  col_names<-names(snp.table)
  which_cols_are_set<-c(pos.col, chrom.col, pvalue.col, format.col) %in% col_names 
  #stop function, if no positions are supplied
  if(!which_cols_are_set[1]) stop("Positions are not found in the snp table. Check the variables 'snp.table' and 'pos.col'.")
#list given columns
  print(paste("Columns given:", c(pos.col, chrom.col, pvalue.col,format.col)[which_cols_are_set]))
  snp_table<-copy(snp.table)
  #if format or chromosome columns are not supplied, add them as a default
  if(!which_cols_are_set[2])  snp_table[ , eval(chrom.col) := 0]
  if(!which_cols_are_set[4])  snp_table[ , eval(format.col) := 0]

  
  #determine how many separate spike calculation are needed
  #Each chromosome is calculated separately 
  #Each genotype format is calculated separately
  if(which_cols_are_set[3]){
    #if p-value column is supplied, significant SNPs in the table are counted for each format
    n_lists<-length(unique(snp_table[get(pvalue.col)<threshold.p.significance, paste(get(chrom.col), get(format.col))]))
    which_combinations<-unique(snp_table[get(pvalue.col)<threshold.p.significance, .(get(chrom.col), get(format.col))])
    setnames(which_combinations, c("V1", "V2"), c(chrom.col, format.col))
    
  } else{ 
    #if p-value is not supplied, all SNPs in the table are counted for each format
    n_lists<-length(unique(snp_table[, paste(get(chrom.col), get(format.col))]))
    which_combinations<-unique(snp_table[, .(get(chrom.col), get(format.col))])
    setnames(which_combinations, c("V1", "V2"), c(chrom.col, format.col))
  }
  
  #set spike table in several steps, using the number of separate spike calculations from above
  spike_table<-rbindlist(lapply(1:n_lists, function(x){

    #reduce table for specific chromosome and format combination
    small_table<-snp_table[get(chrom.col)== which_combinations[x,get(chrom.col)]]
    small_table<-small_table[get(format.col)== which_combinations[x,get(format.col)]]
    
    if(which_cols_are_set[3]){ 
      #if p-value column is supplied, all SNPs below grouping threshold are used to determine spikes
      pos_list<-sort(as.vector(unlist(small_table[get(pvalue.col)<=threshold.p.grouping, ..pos.col])))
    }else {
      #if p-value column is not supplied, all SNPs in the table are used to determine spikes
      pos_list<-sort(as.vector(unlist(small_table[, ..pos.col])))
      
    }

    #determine spikes with limits based on the positions and the specified distance
    spikes<-as.data.table(matrix(get_spike_limits(pos_list, max.distance),ncol=3 ,byrow = TRUE))
    setnames(spikes, names(spikes), c("SPIKE_START", "SPIKE_END", "SPIKE_CENTER"))

    #add back chromosome column and format colum
    spikes[, eval(chrom.col) :=which_combinations[x, get(chrom.col)]]
    spikes[, eval(format.col) :=which_combinations[x, get(format.col)]]
    
    
    #check, if spikes are present, which only contain grouping SNPs and no significant SNPs
    #and remove them
    if(which_cols_are_set[3]){
      for(i in 1:spikes[,.N]){
        low_pos<-spikes[i, SPIKE_START]
        up_pos<-spikes[i, SPIKE_END]
        ifelse(small_table[between(get(pos.col), low_pos, up_pos), min(get(pvalue.col))]<threshold.p.significance,
               spikes[i, REMOVE:=0], spikes[i, REMOVE:=1])
      }
      spikes<-spikes[REMOVE==0]
      spikes[, REMOVE:=NULL]
    }

    
  }))
  #delete unneccesary colums
  if(!all(which_cols_are_set)){
    added_colls<-c(pos.col, chrom.col,pvalue.col,format.col)[!which_cols_are_set]
    for(i in 1:length(added_colls)) spike_table[, added_colls[i]:=NULL] 
  }
  #return table with spikes
  return(spike_table)
}



#recursive function called by determine_spikes to assing positions to spikes
#pos.list is a numeric vector containing positions on the same chromosome and in the same format
#max.distance determines the maximum distance between SNPs assigned to the same spike as in determine-spikes
get_spike_limits <- function(pos.list, max.distance=2e6){
  #initial range of spike: 
  #start with lower end (smallest position number) and add the distance to consider SNPs as part
  #of the same spike. All SNPs within this distance are added to the current spike.
  pos_in_spike_range <- pos.list[pos.list <(min(pos.list)+max.distance)]
  
  
  
  #set new limits of spike around included positions
  #given the way this function is set up, there are no SNPs with a smaller position within range
  #so only the upper range is extended by the speficied distance
  spike_start <- min(pos_in_spike_range)
  spike_end <- max(pos_in_spike_range)+max.distance
  
  #check, if new positions are added with the new upper limit
  pos_in_spike_range_new<-pos.list[pos.list>=spike_start &pos.list<=spike_end]
  
  #repeat with new positions added to the spike, until no further positions are added
  while(length(pos_in_spike_range_new)>length(pos_in_spike_range)){
    pos_in_spike_range<-pos_in_spike_range_new
    spike_end<-max(pos_in_spike_range)+max.distance
    pos_in_spike_range_new<-pos.list[pos.list>=spike_start &pos.list<=spike_end]
  }
  
  #update upper limit to final size of the spike
  spike_end<-max(pos_in_spike_range)
  
  #if all positions from the pos.list input parameter are assigned to this spike,
  #return the center and limits of the spike
  if(all(pos.list %in% pos_in_spike_range)){
    return(c(spike_start, spike_end, round(mean(pos_in_spike_range))))
  } else{
    #if positions from the pos.list input parameter remain unassinged,
    #call the function again with the unassuinged positions as input.
    #this initializes a new spike separate from the SNPs grouped in the current iteration.
    return(c(spike_start, spike_end, round(mean(pos_in_spike_range)),
             get_spike_limits(pos.list[pos.list %notin% pos_in_spike_range])))
  }
}


############# optional function for spikes in the best guess columns ######################
#the resulting table determines, if best guess spikes are specific to the format or 
#if they are present in dosage format as well. 
#this function expects a spike table such as the output of determine_spikes
#else a supply data.table with at least following cols:
#lower limit of spike, upper limit of spike, and format in which the spike is present
#optional column for chromosome
#if the chromosome column is missing, the function assumes all SNPs are on the same chromosome
#all column names can be specified by the respective inputparameter, if they differ from the set default names
#the best.guess input parameter specifies, which entry in the format column signifies a spike in best guess format
#outer.margin determines if spikes are physically close enough to be considered the same spike in either formats,
#if significant positions differ between formats. Analogue to max.distance in determine_spikes
#the function returns a spike table with an added "TYPE" column, 
#which is needed for the midrange filter for best guess spikes 
validate_best_guess_spikes<-function(spike.table,
                      format.col="FORMAT",chrom.col="CHROM",
                      spike.start.col="SPIKE_START", spike.end.col="SPIKE_END",
                      best.guess = "BEST_GUESS", outer.margin=2e6){
  
  #check if chromosome column is present
  remove_chrom<-!chrom.col %in% names(spike.table)
  spike_table<-copy(spike.table)
  if(remove_chrom){
    print(paste("The column" ,chrom.col, "was not found in both tables. All spikes are treated as they are located on the same chromosome."))
    spike_table[,eval(chrom.col):=0]
  }
  #initialize TYPE column with dosage
  spike_table[,TYPE:="DOSAGE"]
  #separate spike table into spikes in dosage or best guess
  bg_table<-spike_table[get(format.col)==best.guess]
  dos_table<-spike_table[get(format.col)!=best.guess]

  #check, if best guess spikes and dosage spikes have an overlap
  for (i in 1:bg_table[,.N]){
    chr<-bg_table[i, get(chrom.col)]
    low <- bg_table[i, get(spike.start.col)]-outer.margin
    up <- bg_table[i, get(spike.end.col)]+outer.margin
    
    #which line in dosage contains an overlapping spike
    check_in_dosage<-dos_table[get(chrom.col)==chr & (up>=get(spike.start.col) &up<=get(spike.end.col) |
                                                        low>=get(spike.start.col) &low<=get(spike.end.col) |
                                                        low<=get(spike.start.col) &up>=get(spike.end.col)),
                               which=TRUE]
    
    if(length(check_in_dosage)>0){
      #if at least one overlap is found between a best guess spike and a dosage spike,
      #adjust the limits of that overlapping spike to encompasse both spikes in the dos table
      dos_table[check_in_dosage, eval(spike.end.col):= max(get(spike.end.col), bg_table[i, get(spike.end.col)])]
      dos_table[check_in_dosage, eval(spike.start.col):= min(bg_table[i, get(spike.start.col)], get(spike.start.col))]
      #set spike type as mixed
      dos_table[check_in_dosage, TYPE:="MIXED"]
    } else{ #if not, mark spike as unique to best guess in the best guess spike table
      bg_table[i, TYPE:= "BEST GUESS ONLY"]}
  }
  #add back the best guess specific spikes to the dosage table
  combine_back<-rbind(dos_table, bg_table[TYPE=="BEST GUESS ONLY"])
  if(remove_chrom) combine_back[, eval(chrom.col):=NULL]
  return(combine_back)
}





###################### function to check if significant SNPs within a spike contain typed SNPs #############
#these spikes are not filtered further by the midrange filter
#and are generally of high quality. 
#input parameter expect a spike.table such as the output of determine_spikes or validate_best_guess_spikes
#a positions.scaffold file which contains the location of typed SNPs. 
#the paramenter snp.table expects the same table as determine_spikes.
#all column names can be specified by the respective inputparameter, if they differ from the set default names
#the input parameter remove.spikes determines, if the output table contains completely imputed spikes (TRUE) or
#if the output table contains only spikes including typed significant SNPs (FALSE). 
check_scaffold<-function(spike.table, positions.scaffold, snp.table,
                         threshold.p.significance=5e-8,
                         pvalue.col="P_VALUE", pos.col="POS",
                         spike.start.col="SPIKE_START", spike.end.col="SPIKE_END", 
                         chrom.col="CHROM", remove.typed.spikes=TRUE){

  #check if chromosome column is present for all input paramenters
  remove_chrom<-!all(chrom.col %in% names(spike.table),chrom.col %in% names(snp.table), chrom.col %in% names(positions.scaffold))

  #if chromosome column is not present, all SNPs are treated as on the same chromosome
  if(remove_chrom){
    print(paste("The column" ,chrom.col,
                " was not found in all tables. All spikes are treated as they are located on the same chromosome."))
    
    #find singificant typed SNPs
    list_sig_snps_scaffold<-snp.table[threshold.p.significance>get(pvalue.col), unique(get(pos.col))]
    list_sig_snps_scaffold<-list_sig_snps_scaffold[list_sig_snps_scaffold %in% positions.scaffold[,get(pos.col)]]
    
    #determine, which spikes the typed SNPs are located in
    to_remove<-spike.table[,between(list_sig_snps_scaffold,
                                    get(spike.start.col),
                                    get(spike.end.col)),
                           by=get(spike.start.col)] [V1==TRUE, unique(get)]

    #return either spike table with only imputed spikes or spike table only with spikes including typed SNPs
    if(remove.typed.spikes){
    return(spike.table[(get(spike.start.col)%in% to_remove)==FALSE])
      }
     return(spike.table[(get(spike.start.col)%in% to_remove)==TRUE])
  } else{
    #if chromosome column is supplied

    #list significant SNPs of scaffold for each chromosome
    list_sig_snps_scaffold<-snp.table[threshold.p.significance>get(pvalue.col), .(get(pos.col), get(chrom.col))]
    setnames(list_sig_snps_scaffold, c("V1", "V2"), c(pos.col, chrom.col))

    #chromosomes need to be treated separately, so different entries in the column are counted
    list_chrom<-spike.table[, unique(get(chrom.col))]
    to_remove<-data.table()
    #for each chromosome determine, if spikes contain significant typed SNPs
    for(i in length(list_chrom)){
      
      list_sig_snps_scaffold_chrom<-list_sig_snps_scaffold[get(chrom.col)==list_chrom[i], get(pos.col)]
      
      start_cols<-spike.table[get(chrom.col)==list_chrom[i],
                  between(list_sig_snps_scaffold_chrom,
                          get(spike.start.col),
                          get(spike.end.col)),
                  by=get(spike.start.col)] [V1==TRUE, unique(get)]
      
      add_to_remove<-data.table(spike.start.col=start_cols,
                                chrom.col=list_chrom[i])
      
    
      setnames(add_to_remove, c("spike.start.col", "chrom.col"), c(spike.start.col, chrom.col))
     
      to_remove<-rbind(to_remove, add_to_remove)
   
    }
    rows_to_remove<-spike.table[to_remove, which=TRUE, on=c(spike.start.col,chrom.col)]

   #return either spike table with only imputed spikes or spike table only with spikes including typed SNPs
    if(remove.typed.spikes){
      return(spike.table[-rows_to_remove])
    } 
    return(spike.table[-rows_to_remove])
  }
}




################# function to run the MidrangeFilter as a post-imputation quality control method ###################
# input parameters snp.table and spike.table are described for the functions above
#and are meant to be the same as the input parameters as the previous functions
#quality.col expects the name of the column in snp.table containing imputation quality measures 
#such as IMPUTE info or Beagle R^2 for each SNP 
#best.guess.type expects the entry in the type column (type.col), which indicates a best guess specific spike, 
#the default matches the output of validate_best_guess_spike
#threshold.quality.low and threshold.quality.high sets thresholds as lowest acceptable imputation 
#quality and high confidence imputation quality respectively as described in the manuscript
#all column names can be specified by the respective inputparameter, if they differ from the set default names
#theshold.p.significance sets the significant threshold for the association analysis
#condensed.return.table indicates, if the return table displayes the decision to keep or discard 
#for each SNP, i.e. one row per SNP (FALSE) or if the decision is displayed for each spike, 
#i.e. one row per spike (TRUE)
#the decision is given as a column in the return table. DISCARD==TRUE means the MidrangeFilter
#considers this spike as likely false and would discard the signal.
midrange_filter<-function(snp.table, spike.table,
                          quality.col, pvalue.col="P_VALUE",
                          pos.col="POS", chrom.col="CHROM", best.guess.type="BEST GUESS ONLY",
                          spike.start.col="SPIKE_START", spike.end.col= "SPIKE_END",
                          threshold.quality.low= 0.3, threshold.quality.high = 0.8, type.col = "TYPE",
                          threshold.p.significance=5e-8, condensed.return.table=FALSE){
  #check which columns are provided
  remove_chrom<-!all(chrom.col %in% names(snp.table), chrom.col %in% names(spike.table))
  no_type_col<-!type.col %in% names(spike.table)
  spike_table<-copy(spike.table)
  snp_table<-copy(snp.table)
  #number spikes
  spike_table[, SPIKE_NUMBER:= 1:.N]
  
  
  if(remove_chrom){
    #warn if chromosome column is not found
    print(paste("The column" ,chrom.col, " was not found in both tables. All spikes are treated as they are located on the same chromosome."))
    spike_table[, eval(chrom.col):=0]
    snp_table[,eval(chrom.col):=0]
  } 
  
  if(no_type_col){
    #warn if no type column is found
    print(paste("The column" ,type.col, " was not found in spike.table. All spikes are treated as Mixed/Dosage types."))
    spike_table[,eval(type.col):="DOSAGE/MIXED"]
  }

#update number spikes, if different chromosomes are present
  for (i in 1:spike_table[,.N]){
    
    chr<-spike_table[i,get(chrom.col)]
    
    low<-as.numeric(spike_table[i, get(spike.start.col)])
    up<-as.numeric(spike_table[i, get( spike.end.col)])
    snp_table[get(chrom.col)==chr & between(get(pos.col), low, up), SPIKE_NUMBER:=i]
  }
  
#if p-value col is not provided, all SNPs in snp.table are considered significant
  if(pvalue.col %in% names(snp_table)){
    full_table<-snp_table[get(pvalue.col)<threshold.p.significance]
  } else{
    full_table<-snp_table[!is.na(SPIKE_NUMBER)]
  }

#add spike.table information to snp.table to make decision
  full_table<-spike_table[, .(SPIKE_NUMBER, get(type.col), get(spike.start.col), get(spike.end.col))][full_table, on="SPIKE_NUMBER"]
  setnames(full_table, c("V2", "V3", "V4"), c(type.col,spike.start.col, spike.end.col))
  
  #determine best and worst imputation quality within each spike
  full_table[, MIN_Q:=min(get(quality.col), na.rm = TRUE), SPIKE_NUMBER]
  full_table[, MAX_Q:=max(get(quality.col), na.rm = TRUE), SPIKE_NUMBER]

  #initialize discard column. FALSE means not to discard the spike.
  full_table[,DISCARD:=FALSE]
  #decision MidrangeFilter:
  #discard all spikes, where the worst SNP cannot meet the lower imputation quality threshold
  full_table[MIN_Q<threshold.quality.low, DISCARD:=TRUE]
  #discard all best guess specific spikes, where the best SNP cannot meet the higher imputation quality threshold
  full_table[get(type.col)==best.guess.type & MAX_Q<threshold.quality.high, DISCARD:=TRUE]

  #optional remove default chromosome column
  if(remove_chrom)full_table[,eval(chrom.col):=NULL]


  #condense output table to one row for each spike, if indicated by input parameter
  if(condensed.return.table){
    full_table[,eval(pos.col):=NULL]
    full_table[,eval(pvalue.col):=NULL]
    full_table[,FORMAT:=NULL]
    full_table[,eval(quality.col):=NULL]
    
    return(unique(full_table))
  } 
  return(full_table)
} 



