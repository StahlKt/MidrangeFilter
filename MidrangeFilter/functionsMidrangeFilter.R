determine_spikes<-function(snp.table,pos.col="POS", chrom.col="CHROM", pvalue.col="P_VALUE", format.col="FORMAT", 
                           threshold.p.grouping=5e-7, threshold.p.significance=5e-8, max.distance=2e6){
  
  
  col_names<-names(snp.table)
  which_cols_are_set<-c(pos.col, chrom.col, pvalue.col, format.col) %in% col_names 
  
  if(!which_cols_are_set[1]) stop("Positions are not found in the snp table. Check the variables 'snp.table' and 'pos.col'.")
  
  
  print(paste("Columns given:", c(pos.col, chrom.col, pvalue.col,format.col)[which_cols_are_set]))
  snp_table<-copy(snp.table)
  if(!which_cols_are_set[2])  snp_table[ , eval(chrom.col) := 0]
  if(!which_cols_are_set[4])  snp_table[ , eval(format.col) := 0]
  #how many separate lists do we need
  
  #if p-value column is supplied
  if(which_cols_are_set[3]){
    
    
    
    n_lists<-length(unique(snp_table[get(pvalue.col)<threshold.p.significance, paste(get(chrom.col), get(format.col))]))
    which_combinations<-unique(snp_table[get(pvalue.col)<threshold.p.significance, .(get(chrom.col), get(format.col))])
    setnames(which_combinations, c("V1", "V2"), c(chrom.col, format.col))
    
  } else{
    n_lists<-length(unique(snp_table[, paste(get(chrom.col), get(format.col))]))
    which_combinations<-unique(snp_table[, .(get(chrom.col), get(format.col))])
    setnames(which_combinations, c("V1", "V2"), c(chrom.col, format.col))
  }
  
  #for each combination, send positions to the position only function
  spike_table<-rbindlist(lapply(1:n_lists, function(x){
    
    small_table<-snp_table[get(chrom.col)== which_combinations[x,get(chrom.col)]]
    
    small_table<-small_table[get(format.col)== which_combinations[x,get(format.col)]]
    
    if(which_cols_are_set[3]){  # eval(pos.col)
      pos_list<-sort(as.vector(unlist(small_table[get(pvalue.col)<=threshold.p.grouping, ..pos.col])))
    }else {
      pos_list<-sort(as.vector(unlist(small_table[, ..pos.col])))
      
    }
    
    spikes<-as.data.table(matrix(get_spike_limits(pos_list, max.distance),ncol=3 ,byrow = TRUE))
    setnames(spikes, names(spikes), c("SPIKE_START", "SPIKE_END", "SPIKE_CENTER"))
    
    spikes[, eval(chrom.col) :=which_combinations[x, get(chrom.col)]]
    spikes[, eval(format.col) :=which_combinations[x, get(format.col)]]
    
    
    #check, if we have spikes that stay under significance and remove them 
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
  
  return(spike_table)
}



#recursive function to assing positions to spikes

get_spike_limits <- function(pos.list, max.distance=2e6){
  #initial range of spike
  
  pos_in_spike_range <- pos.list[pos.list <(min(pos.list)+max.distance)]
  
  
  
  #set limits of spike around included positions, lower limit can't be below 1
  spike_start <- min(pos_in_spike_range)
  spike_end <- max(pos_in_spike_range)+max.distance
  
  #check, if new positions are added with the new limits
  pos_in_spike_range_new<-pos.list[pos.list>=spike_start &pos.list<=spike_end]
  
  #repeat with new positions added to the spike, until no further positions are added
  while(length(pos_in_spike_range_new)>length(pos_in_spike_range)){
    pos_in_spike_range<-pos_in_spike_range_new
    spike_end<-max(pos_in_spike_range)+max.distance
    pos_in_spike_range_new<-pos.list[pos.list>=spike_start &pos.list<=spike_end]
  }
  
  #update limits to final size of the spike
  spike_end<-max(pos_in_spike_range)
  
  #if all positions are assigned, return the center and limits of the spike
  if(all(pos.list %in% pos_in_spike_range)){
    return(c(spike_start, spike_end, round(mean(pos_in_spike_range))))
  } else{
    #if positions remain unassinged, call the function again with the unassuinged positions as input. 
    return(c(spike_start, spike_end, round(mean(pos_in_spike_range)),
             get_spike_limits(pos.list[pos.list %notin% pos_in_spike_range])))
  }
}

#optional function.
#supply data table with at least following cols:
#lower limit of spike, 
#upper limit of spike,
#format

validate_bg<-function(spike.table, format.col="FORMAT", best.guess = "BEST_GUESS", 
                      spike.start.col="SPIKE_START", spike.end.col="SPIKE_END", chrom.col="CHROM",
                      outer.margin=2e6){
  
  
  remove_chrom<-!chrom.col %in% names(spike.table)
  spike_table<-copy(spike.table)
  if(remove_chrom){
    print(paste("The column" ,chrom.col, "was not found in both tables. All spikes are treated as they are located on the same chromosome."))
    spike_table[,eval(chrom.col):=0]
  }
  spike_table[,TYPE:="DOSAGE"]
  bg_table<-spike_table[get(format.col)==best.guess]
  
  dos_table<-spike_table[get(format.col)!=best.guess]
  
  for (i in 1:bg_table[,.N]){
    chr<-bg_table[i, get(chrom.col)]
    low <- bg_table[i, get(spike.start.col)]-outer.margin
    up <- bg_table[i, get(spike.end.col)]+outer.margin
    
    #which line in dosagoe contains an overlapping spike
    check_in_dosage<-dos_table[get(chrom.col)==chr & (up>=get(spike.start.col) &up<=get(spike.end.col) |
                                                        low>=get(spike.start.col) &low<=get(spike.end.col) |
                                                        low<=get(spike.start.col) &up>=get(spike.end.col)),
                               which=TRUE]
    #if one overlap is found, adjust the limits of that overlapping spike to encompasse both spikes
    if(length(check_in_dosage)>0){
      dos_table[check_in_dosage, eval(spike.end.col):= max(get(spike.end.col), bg_table[i, get(spike.end.col)])]
      dos_table[check_in_dosage, eval(spike.start.col):= min(bg_table[i, get(spike.start.col)], get(spike.start.col))]
      dos_table[check_in_dosage, TYPE:="MIXED"]
    } else{bg_table[i, TYPE:= "BEST GUESS ONLY"]} #if not, mark spike as unique to best guess
  }
  combine_back<-rbind(dos_table, bg_table[TYPE=="BEST GUESS ONLY"])
  if(remove_chrom) combine_back[, eval(chrom.col):=NULL]
  return(combine_back)
}



midrange_filter<-function(snp.table, spike.table, quality.col, pvalue.col="P_VALUE",
                          pos.col="POS", chrom.col="CHROM", best.guess.type="BEST GUESS ONLY",
                          spike.start.col="SPIKE_START", spike.end.col= "SPIKE_END",
                          threshold.quality.low= 0.3, threshold.quality.high = 0.8, type.col = "TYPE",
                          threshold.p.significance=5e-8, condensed.return.table=FALSE){
  
  remove_chrom<-!all(chrom.col %in% names(snp.table), chrom.col %in% names(spike.table))
  no_type_col<-!type.col %in% names(spike.table)
  spike_table<-copy(spike.table)
  snp_table<-copy(snp.table)
  spike_table[, SPIKE_NUMBER:= 1:.N]
  
  
  if(remove_chrom){
    print(paste("The column" ,chrom.col, " was not found in both tables. All spikes are treated as they are located on the same chromosome."))
    spike_table[, eval(chrom.col):=0]
    snp_table[,eval(chrom.col):=0]
  } 
  
  if(no_type_col){
    print(paste("The column" ,type.col, " was not found in spike.table. All spikes are treated as Mixed/Dosage types."))
    spike_table[,eval(type.col):="DOSAGE/MIXED"]
  }
  
  for (i in 1:spike_table[,.N]){
    
    chr<-spike_table[i,get(chrom.col)]
    
    low<-as.numeric(spike_table[i, get(spike.start.col)])
    up<-as.numeric(spike_table[i, get( spike.end.col)])
    snp_table[get(chrom.col)==chr & between(get(pos.col), low, up), SPIKE_NUMBER:=i]
  }
  
  
  if(pvalue.col %in% names(snp_table)){
    full_table<-snp_table[get(pvalue.col)<threshold.p.significance]
  } else{
    full_table<-snp_table[!is.na(SPIKE_NUMBER)& get(pvalue.col)<threshold.p.significance]
  }
  
  full_table<-spike_table[, .(SPIKE_NUMBER, get(type.col), get(spike.start.col), get(spike.end.col))][full_table, on="SPIKE_NUMBER"]
  setnames(full_table, c("V2", "V3", "V4"), c(type.col,spike.start.col, spike.end.col))
  full_table[, MIN_Q:=min(get(quality.col), na.rm = TRUE), SPIKE_NUMBER]
  full_table[, MAX_Q:=max(get(quality.col), na.rm = TRUE), SPIKE_NUMBER]
  
  full_table[,DISCARD:=FALSE]
  
  full_table[get(type.col)==best.guess.type & MAX_Q<threshold.quality.high, DISCARD:=TRUE]
  full_table[get(type.col)!=best.guess.type & MIN_Q<threshold.quality.low, DISCARD:=TRUE]
  
  if(remove_chrom)full_table[,eval(chrom.col):=NULL]
  if(condensed.return.table){
    full_table[,eval(pos.col):=NULL]
    full_table[,eval(pvalue.col):=NULL]
    full_table[,FORMAT:=NULL]
    full_table[,eval(quality.col):=NULL]
    
    return(unique(full_table))
  } 
  return(full_table)
} 



