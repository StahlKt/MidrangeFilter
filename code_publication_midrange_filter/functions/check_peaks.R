#this is the unclean basis of the spike aggregation function
#it works the same, only some parameters are hardcoded
#please refer to the functions in the midrange filter folder for a cleaner and commented version.


check_peak<-function(run.name, run.number){
  dir_pv<-file.path(dir_out_association, paste0(run.name, "_", run.number, "_pv.RDS"))
  
  if(!file.exists(dir_pv)){
    dir_pv<-file.path(dir_out_association, paste0(run.name, "_", run.number, "_pv_NEW.RDS"))
  }
  
  
  asso_table<-readRDS(dir_pv)

  
  pos_list_bg<-asso_table[FORMAT=="BEST_GUESS" &P_VALUE<5e-8, unique(POS)]
  pos_list_dos<-asso_table[FORMAT=="DOSAGE" &P_VALUE<5e-8, unique(POS)]
  
  if(length(pos_list_bg)==0 &length(pos_list_dos)==0) return(NA)
  pos_list_bg_wide<-asso_table[FORMAT=="BEST_GUESS" &P_VALUE<5e-7, unique(POS)]
  pos_list_dos_wide<-asso_table[FORMAT=="DOSAGE" &P_VALUE<5e-7, unique(POS)]
  
  
  pos_length<-c(length(pos_list_bg), length(pos_list_dos))
  
  
  
  
  
  
  if(pos_length[1]!=0){
    peaks_bg<-as.data.table(matrix(get_number_peaks(pos_list_bg_wide),ncol=3 ,byrow = TRUE))
    setnames(peaks_bg, names(peaks_bg), c("PEAK_MEAN", "UP_LIM", "LOW_LIM"))
    peaks_bg[,FORMAT:="BEST_GUESS"]
    peaks_bg[,PEAK_IN_FORMAT:=1]
    
     peaks_bg<-as.data.table(peaks_bg %>% rowwise() %>% filter(any(between(pos_list_bg, LOW_LIM, UP_LIM))))
  }
  
  if(pos_length[2]!=0){
    peaks_dos<-as.data.table(matrix(get_number_peaks(pos_list_dos_wide),ncol=3 ,byrow = TRUE))
    setnames(peaks_dos, names(peaks_dos), c("PEAK_MEAN", "UP_LIM", "LOW_LIM"))
    peaks_dos[,FORMAT:="DOSAGE"]
    peaks_dos[,PEAK_IN_FORMAT:=1]
     peaks_dos<-as.data.table(peaks_dos %>% rowwise() %>% filter(any(between(pos_list_dos, LOW_LIM, UP_LIM))))
  }
  
  
  # peaks_bg
  # peaks_dos
 
 
  # 
  
 
  if(pos_length[1]==0){
    peaks_bg<-as.data.table(matrix(get_number_peaks(pos_list_dos),ncol=3 ,byrow = TRUE))
    setnames(peaks_bg, names(peaks_bg), c("PEAK_MEAN", "UP_LIM", "LOW_LIM"))
    peaks_bg[,FORMAT:="BEST_GUESS"]
    peaks_bg[,PEAK_IN_FORMAT:=0]
  }
  
  
  if(pos_length[2]==0){
    peaks_dos<-as.data.table(matrix(get_number_peaks(pos_list_bg),ncol=3 ,byrow = TRUE))
    setnames(peaks_dos, names(peaks_dos), c("PEAK_MEAN", "UP_LIM", "LOW_LIM"))
    peaks_dos[,FORMAT:="DOSAGE"]
    peaks_dos[,PEAK_IN_FORMAT:=0]
  }
  
  peaks_bg
  peaks_dos
  if(pos_length[2]!=0 &pos_length[1]!=0){
    
    #chekc for overlap
    dos_limits_with_bg_peaks<-rbind(as.data.table(peaks_dos %>% rowwise() %>% filter(any(between(peaks_bg$LOW_LIM, LOW_LIM, UP_LIM)))),
                                    as.data.table(peaks_dos %>% rowwise() %>% filter(any(between(peaks_bg$UP_LIM, LOW_LIM, UP_LIM)))))
    bg_limits_with_dos_peaks<-rbind(as.data.table(peaks_bg %>% rowwise() %>%  filter(any(between(peaks_dos$LOW_LIM, LOW_LIM, UP_LIM)))),
                                    as.data.table(peaks_bg %>% rowwise() %>%  filter(any(between(peaks_dos$UP_LIM, LOW_LIM, UP_LIM)))))
    
    dos_limits_with_bg_peaks<-unique(dos_limits_with_bg_peaks)
    bg_limits_with_dos_peaks<-unique(bg_limits_with_dos_peaks)
    
    
    dos_to_add_to_bg<-peaks_dos[LOW_LIM %notin% dos_limits_with_bg_peaks$LOW_LIM]
    
    
    
    
    dos_to_add_to_bg[,PEAK_IN_FORMAT:=0]
    dos_to_add_to_bg[,FORMAT:="BEST_GUESS"]
    
    bg_to_add_to_dos<-peaks_bg[LOW_LIM %notin% bg_limits_with_dos_peaks$LOW_LIM]
    
    bg_to_add_to_dos[,PEAK_IN_FORMAT:=0]
    bg_to_add_to_dos[,FORMAT:="DOSAGE"]
    peaks_bg<-rbind(peaks_bg, dos_to_add_to_bg)
    peaks_dos<-rbind(peaks_dos, bg_to_add_to_dos)
    
    
    
    if(length(peaks_bg$PEAK_MEAN)!=length(peaks_dos$PEAK_MEAN)){
      ifelse(length(peaks_bg$PEAK_MEAN)>length(peaks_dos$PEAK_MEAN), 
             peaks_to_reduce<-peaks_bg,
             peaks_to_reduce<-peaks_dos)
      
      list_to_test<-peaks_to_reduce[PEAK_IN_FORMAT==0]$PEAK_MEAN
      
      for(i in 1:length(list_to_test)){
        res<-peaks_to_reduce[list_to_test[i]>LOW_LIM& list_to_test[i]<UP_LIM &PEAK_IN_FORMAT==1, PEAK_MEAN]
        
        if(length(res)!=0) peaks_to_reduce<-peaks_to_reduce[PEAK_MEAN%notin%list_to_test[i]]
      }
      
      ifelse(length(peaks_bg$PEAK_MEAN)>length(peaks_dos$PEAK_MEAN), 
             peaks_bg<-peaks_to_reduce,
             peaks_dos<-peaks_to_reduce)
      
    }
    
    
    rm(bg_to_add_to_dos)
    rm(dos_to_add_to_bg)
    rm(bg_limits_with_dos_peaks)
    rm(dos_limits_with_bg_peaks)
    
    
  }
  
 
  setkey(peaks_dos, PEAK_MEAN)
  setkey(peaks_bg, PEAK_MEAN)
  peaks_bg[,PEAK_NUMBER:=1:length(PEAK_MEAN)]
  peaks_dos[,PEAK_NUMBER:=1:length(PEAK_MEAN)]
  
  
  
  peak_table<-rbind(peaks_dos,peaks_bg)
  peak_table
  number_peaks<-peak_table[,max(PEAK_NUMBER)]
  rm(peaks_dos)
  rm(peaks_bg)
  
  peak_table[LOW_LIM<0, LOW_LIM:=0]
  
  
  
  peak_chars<-rbindlist(lapply(1:number_peaks, function(x){
    low_lim_bg<- peak_table[PEAK_NUMBER==x & FORMAT=="BEST_GUESS",LOW_LIM]
    up_lim_bg<-peak_table[PEAK_NUMBER==x & FORMAT=="BEST_GUESS", UP_LIM]
    low_lim_dos<-peak_table[PEAK_NUMBER==x & FORMAT=="DOSAGE",LOW_LIM]
    up_lim_dos<-peak_table[PEAK_NUMBER==x & FORMAT=="DOSAGE",UP_LIM]
    
    
    
    
    data.table(PEAK_NUMBER=x,
               FORMAT =c("BEST_GUESS", "DOSAGE"),
               PEAK_N.SIG = c(asso_table[FORMAT=="BEST_GUESS" & POS>low_lim_bg &POS<up_lim_bg & P_VALUE<5e-8, 
                                         .N],
                              asso_table[FORMAT=="DOSAGE" & POS>low_lim_dos &POS<up_lim_dos & P_VALUE<5e-8, 
                                         .N]),
               PEAK_N.TOP=c(asso_table[FORMAT=="BEST_GUESS" & POS>low_lim_bg &POS<up_lim_bg & P_VALUE<5e-6, 
                                       .N],asso_table[FORMAT=="DOSAGE" & POS>low_lim_dos &POS<up_lim_dos & P_VALUE<5e-6, 
                                                      .N]),
               PEAK_STRENGTH_TOP=c(asso_table[FORMAT=="BEST_GUESS" & POS>low_lim_bg &POS<up_lim_bg & P_VALUE<5e-6, 
                                              sum(-log10(P_VALUE))],
                                   asso_table[FORMAT=="DOSAGE" & POS>low_lim_dos &POS<up_lim_dos & P_VALUE<5e-6, 
                                              sum(-log10(P_VALUE))]),
               PEAK_STRENGTH_ALL=c(asso_table[FORMAT=="BEST_GUESS" & POS>low_lim_bg &POS<up_lim_bg, 
                                              sum(-log10(P_VALUE))],
                                   asso_table[FORMAT=="DOSAGE" & POS>low_lim_dos &POS<up_lim_dos, 
                                              sum(-log10(P_VALUE))]))
    
  }))
  
  
  
  
  peak_table<-merge(peak_table,peak_chars)
  

 peak_table[,PEAK_BREADTH:= UP_LIM-LOW_LIM]              
 peak_table[,RUN:=run.number]
 peak_table[,NAME:=run.name]
 
 return(peak_table)
}



#recursive function to check for number of peaks

get_number_peaks<-function(pos.list){

      

     pos_in_peak_range<-pos.list[pos.list <(min(pos.list)+2e6)]
     
     peak_lower_limit<- min(pos_in_peak_range)-2e6
     peak_upper_limit<- max(pos_in_peak_range)+2e6
     
     pos_in_peak_range_new<-pos.list[pos.list>peak_lower_limit &pos.list<peak_upper_limit]
     
     while(length(pos_in_peak_range_new)>length(pos_in_peak_range)){
       pos_in_peak_range<-pos_in_peak_range_new
       peak_lower_limit<-min(pos_in_peak_range)-2e6
       peak_upper_limit<-max(pos_in_peak_range)+2e6
       pos_in_peak_range_new<-pos.list[pos.list>peak_lower_limit &pos.list<peak_upper_limit]
     }
     
     peak_lower_limit<-min(pos_in_peak_range)-5e5
     peak_upper_limit<-max(pos_in_peak_range)+5e5
    if(all(pos.list %in% pos_in_peak_range)){
      return(c(round(mean(pos_in_peak_range)), peak_upper_limit, peak_lower_limit))
    } else{
      
    
   return(c(round(mean(pos_in_peak_range)), peak_upper_limit, peak_lower_limit,
            get_number_peaks(pos.list[pos.list %notin% pos_in_peak_range])))
    }
}

