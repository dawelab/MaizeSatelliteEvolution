#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)
rep<-args[1]
dir=args[2]

library('dplyr')
library('reshape2')
library('tidyr')
library('stringr')
library('rpart')
library('data.table')

setwd(dir)
fit2_test<- readRDS("/scratch/rdp22327/Dawe/consensus/lda_model.rds")
`%!in%` = Negate(`%in%`)

path =getwd()

##smoothing classifications for an array

spur_smooth<- function(k) {
  bins<- k
  bins$bin_start<- as.numeric(str_split(bins$read, pattern="_", simplify=T)[,1])
  bins$bin_end<- bins$bin_start+10000
  bins$array<- str_split(bins$read, pattern="_", simplify=T)[,3]
  #bins$tree_pred<- as.character(bins$tree_pred)
  #bins$tree_pred<- as.character( bins$tree_pred)
  #bins$comp<- bins$tree_pred == bins$lda_pred
  bins$lda_pred_smooth<-NA
  bins_smooth<- as.data.frame(matrix( nrow=0, ncol=ncol(bins)))
  colnames(bins_smooth)<- colnames(bins)
  
  #if a single bin is off, compared to the 2 preceeding and following bins, change the bin to match
  for( i in 1:length(unique( bins$array))){
    sub<- bins[bins$array == unique(bins$array)[i],] %>% arrange(bin_start)
    sub$num<- 1:nrow(sub)
    for(j in 1:nrow( sub)){
      if( j > 2 & j < (nrow(sub)-2)){
        if( sub$lda_pred[j] != sub$lda_pred[j-1] & sub$lda_pred[j] != sub$lda_pred[j+1] & sub$lda_pred[j+1] == sub$lda_pred[j-1]  & sub$clust_val[j+1] == sub$clust_val[j-1] ){
          sub$lda_pred_smooth[j] <- as.character(sub$lda_pred[j+1])
          sub$clust_val[j] <- sub$clust_val[j+1]
          sub$lda_pred_pos[j] <-sub$lda_pred_pos[j+1]
        }
      }
      if( is.na(sub$lda_pred_smooth[j])){ sub$lda_pred_smooth[j]<-as.character(sub$lda_pred[j])}
    }
    bins_smooth<- rbind(bins_smooth, sub)
  }
  return(bins_smooth)
}

# Combine consecutive bins with same class/threshold level
comb_over<- function(fil){
  bins_smooth<- fil
  bins_smooth$tar<-NA
  bins_smooth_comb<- as.data.frame(matrix(nrow=0, ncol=ncol(bins_smooth)))
  colnames( bins_smooth_comb)<- colnames( bins_smooth)
  for( i in 1:length(unique( bins_smooth$array))){
    
    #sub one array at a time
    sub<- bins_smooth[bins_smooth$array == unique(bins_smooth$array)[i],] %>% arrange(bin_start)
    
    #merge bins with less than 5 monomers with the preceding bin
    change<- 1
    while( change > 0){
      change<- 0
      for(r in 1:nrow(sub)){
      #unless it's the last row, combine with bin below
        if(sub$num_monomer[r] < 5 & r < nrow(sub) & sub$bin_start[r+1] != sub$bin_start[r]){ 
          sub$bin_start[r+1] <- sub$bin_start[r]
          sub$tar[r]<- "remove"
          sub$num_monomer[r+1]<- sub$num_monomer[r]+ sub$num_monomer[r+1]
          change<- change+1
        }else {
        #IF it's the last row, combine with bin above
          if(sub$num_monomer[r] < 5 & i == nrow(sub) & sub$tar[r] %!in% "remove" ){ 
            sub$end[r-1] <- sub$end[r]
            sub$tar[r]<- "remove"
            sub$num_monomer[r-1]<- sub$num_monomer[r]+ sub$num_monomer[r-1]
            change<- change+1
          }
        }
      } }
    sub<- sub[is.na(sub$tar),]
    sub$num<- 1:nrow(sub)
   
    j<- 2
    tar<- 1
    sub$tar[tar]<-"YES"
    #group similar regions together-- consecutive, same class, & same threshold
    #here, I'll add all the data to the first (tar) entry of the time
    while(j <= nrow(sub)){
      if( sub$lda_pred_smooth[j] == sub$lda_pred_smooth[tar] & sub$bin_start[j] <= sub$bin_end[tar] & sub$bin_end[j] >= sub$bin_start[tar] & sub$clust_val[j] == sub$clust_val[tar]){
        sub$bin_end[tar] <- sub$bin_end[j]
      }else{
       tar<- sub$num[j]
        sub$tar[tar]<-"YES"
      }
      j<-j+1
    }
    #filter sub to the combined regions
    sub2<- sub[sub$tar %in% "YES",]
    
    #now loops remove overlaps.
    #multi loops used so we can give preference to certain classes
    if( nrow(sub2)>1){
    
      #second, HOR
      for( k in 1:nrow(sub2)){
        if( sub2$lda_pred_smooth[k] == "HOR"){
          if(k == 1){sub2$start[k+1] <-  sub2$end[k] }else{
            if( k == nrow(sub2)){ sub2$end[k-1] <-  sub2$start[k]}else{
              sub2$start[k+1] <-  sub2$end[k] 
              sub2$end[k-1] <-  sub2$start[k]
            }
          }
          
        }
      }
      

      #third, Order
      for( k in 1:nrow(sub2)){
        if( sub2$lda_pred_smooth[k] == "Order"){
          if(k == 1){sub2$start[k+1] <-  sub2$end[k] }else{
            if( k == nrow(sub2)){ sub2$end[k-1] <-  sub2$start[k]}else{
              sub2$start[k+1] <-  sub2$end[k] 
              sub2$end[k-1] <-  sub2$start[k]
            }
          }
          
        }
      }
      
      #fourth, Fix disorder bins
      if( sub2$lda_pred_smooth[1] == "Disorder"){sub2$end[1]<- sub2$start[2] }
      for( k in 1:nrow(sub2)){
        if( sub2$lda_pred_smooth[k] == "Disorder"){
          if(k == 1){sub2$start[k+1] <-  sub2$end[k] }else{
            if( k == nrow(sub2)){ sub2$end[k-1] <-  sub2$start[k]}else{
              sub2$start[k+1] <-  sub2$end[k] 
              sub2$end[k-1] <-  sub2$start[k]
            }
          }
          
        }
      }

    
      #again, group similar bins, if shifts occurred above
      j<- 2
      tar<- 1
      sub2$tar<-NA
      sub2<- sub2[sub2$bin_end > sub2$bin_start,]
      sub2$tar[tar]<-"YES"
      sub2$num<- 1:nrow(sub2)
      # regroup similar regions together-- consecutive and same class --> the other shifts so far may have removed separating bins 
      # here, I'll add all the data to the first (tar) entry of the time
      while(j <= nrow(sub2)){
        if( sub2$lda_pred_smooth[j] == sub2$lda_pred_smooth[tar] & sub2$bin_start[j] <= sub2$bin_end[tar] & sub2$bin_end[j] >= sub2$bin_start[tar] & sub$clust_val[j] == sub$clust_val[tar] ){
          sub2$end[tar] <- sub2$end[j]
        }else{
          tar<- sub2$num[j]
          sub2$tar[tar]<-"YES"
        }
        j<-j+1
      }
      #filter sub to the combined regions
      sub3<- sub2[sub2$tar %in% "YES",]
      
      #remove any remaining overlaps simply 
      
      if( nrow(sub3)>1){
        for(f in 1:(nrow(sub3)-1)){
          if( sub3$bin_end[f] != sub3$bin_start[f+1] ){ sub3$bin_start[f+1]<- sub3$bin_end[f] }
        }
      }
    }else{
      sub3<- sub2
    }
    
    bins_smooth_comb<- rbind(bins_smooth_comb, sub3 )
  }
  
  return(bins_smooth_comb)
}


fil=paste(rep,"_reads_class.csv", sep="")
combined_Multi_bins<- read.csv( fil)
# Now back to main data and applying those functions
combined_Multi_bins2<- spur_smooth(combined_Multi_bins) 
combined_Multi_bins3<- comb_over(combined_Multi_bins2) 
##outputting data
write.table(combined_Multi_bins3,  paste(rep, "_fin_bins_combined.txt", sep=""), quote=F, row.names=F)
