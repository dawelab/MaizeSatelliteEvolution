#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)
dir= args[1]

#setting up other data needed
setwd(dir)

data_files<- dir(getwd(), pattern=c( ".blat_bins_nam_groups.csv"), full.names=FALSE)
`%!in%` = Negate(`%in%`)
clust_nams<- c(LETTERS[1:25], 1:9, letters[1:25] ) ##set characters here for the pattern strings. monomers from private clusters are all called
 "Z"
priv_lab<- "Z"

library("stringr")
library("dplyr")


pattern_string_v2<- function(dat){

  dat$array<- paste(str_split(dat$nam, pattern = ":", simplify=T)[,1], str_split(dat$nam, pattern = ":", simplify=T)[,2], sep=":")
  clust_uni<-dat %>% group_by(clust) %>% dplyr::summarize(count=n_distinct(array))
  singles<- clust_uni[clust_uni$count==1,]$clust
  
  if( length(singles) > 0){
    dat[dat$clust %in% singles ,]$clust<- max(singles)
    dat_clust_rename<- as.data.frame(sort(unique(dat$clust)))
    dat_clust_rename$letters<- clust_nams[1:nrow(dat_clust_rename)]
    dat_clust_rename[dat_clust_rename[,1] == max(singles), ]$letters<- priv_lab
  } else{ 
    dat_clust_rename<- as.data.frame(sort(unique(dat$clust)))
    dat_clust_rename$letters<- clust_nams[1:nrow(dat_clust_rename)]
    }
  dat_sub<- dat[,c("nam", "clust")]
   dat_sub$clust<- as.factor(dat_sub$clust)
  
   dat_sub$letter<- NA
  for(i in 1:nrow(dat_sub)){
     dat_sub$letter[i]<-dat_clust_rename[dat_clust_rename$`sort(unique(dat$clust))` == dat_sub$clust[i], ]$letters
  }

 dat_sub$dir<- str_split(dat_sub$nam, pattern = ":", simplify=T)[,1]
  dat_sub$bin<- str_split(dat_sub$nam, pattern = ":", simplify=T)[,2]
  dat_sub$orig<- str_split(dat_sub$nam, pattern = ":", simplify=T)[,3]
 return( dat_sub)
}



kmer_counts <- function(seq, kval){
  if( nchar(seq) < kval){
    return(NA)
  }else{
    start<- 1
    end<- kval
    mers <- vector()
    while(end<=nchar(seq)){
      mers<-append(mers, substr(seq, start, end))
      start<-start+1
      end<-end+1
    }
    return(mers)
  }
}

filtered_kmer_counts<- function(sequence){
  ##first, calling kmer counting for every value of k from k=2 until k when all kmers only have a count of 1
  check<-1
  k<-2
  while(check>=1){
    if(k==2){
      out<- kmer_counts(sequence,k)
       k_out<- as.data.frame(table(out))
      check<- nrow(k_out[k_out$Freq >1,])
      dat_out<- k_out
    }else{
      out<- kmer_counts(sequence,k)
      k_out<- as.data.frame(table(out))
      check<- nrow(k_out[k_out$Freq >1,])
      if(check>0){
        dat_out<- rbind(dat_out, k_out)
      }
    }
    k<-k+1
  }
 #second type of check-- when looking for counts of a kmer, finding substrings does not allow overlaps (mean, if the seq is GGGG, this will make sure n=2, not n=3)
  for(i in 1:nrow(dat_out)){
    dat_out$Freq[i]<- length(unlist(gregexpr(dat_out$out[i], sequence))) 
  }
  
  ##then, remove kmers that only occur once and kmers that contain a private monomer
  n<-  as.character(dat_out$out)
  dat_out_Z<- if(length(grep("Z", dat_out$out))>0){dat_out[-grep("Z", dat_out$out) ,]}else{ dat_out}
  dat_out_Z_sub <- dat_out_Z[dat_out_Z$Freq >1,]
  
  ##then, remove kmers that are all one letter (i.e. CCCCC or CC) or are more than half one letter (CCAC)
  sum_uniq_lets<- function(x){
    co<- str_count(x,LETTERS)
    m<- co
    length(co[co>0])
    return(sum(length(co[co>0])))
  }

max_char_p<- function(x){
    max(str_count(x,LETTERS)/nchar(x))
  }
  
dat_out_Z_sub$charuni<- as.vector(sapply(as.character(dat_out_Z_sub$out),sum_uniq_lets ))
  dat_out_Z_sub$char_prop_m<- as.vector(sapply(as.character(dat_out_Z_sub$out), max_char_p ))
  dat_out_Z_sub$char_num<- nchar(as.character(dat_out_Z_sub$out))
  dat_out_Z_sub_char<- dat_out_Z_sub[dat_out_Z_sub$charuni>1 & dat_out_Z_sub$char_prop_m <=.50 & dat_out_Z_sub$char_num > 2,]
  return(dat_out_Z_sub_char)
}

filtered_kmer_counts_2<- function(dat_out_Z_sub_char){
 if(nrow(dat_out_Z_sub_char)>1){
    dat_out_Z_sub_char$filt<-NA
    dat_out_Z_sub_char$count<- NA
    dat_out_Z_sub_char$count_self<- NA
    #then, remove larger kmers that contain smaller kmers that occur MORE often
    
    for(i in 1:nrow(dat_out_Z_sub_char)){
      if( dat_out_Z_sub_char$filt[i] %!in% "filt"){
      tar_freq<- dat_out_Z_sub_char$Freq[i]
      tar_num<- dat_out_Z_sub_char$char_num[i]
      val<- dat_out_Z_sub_char$out[i]
      dat_out_Z_sub_char<- dat_out_Z_sub_char %>% mutate(filt = ifelse( char_num>tar_num & grepl( val, out) & Freq < tar_freq, "filt", filt))
      }
    }

fin<- dat_out_Z_sub_char[dat_out_Z_sub_char$filt %!in% "filt",]
    
    #then, remove smaller kmers that are contained in larger kmers and occur AS often
    dat_out_Z_sub2<- as.data.frame(matrix(ncol=ncol(fin), nrow=0))
    
    for(i in 1:(nrow(fin))){
      tar_freq<- fin$Freq[i]
      tar_num<- fin$char_num[i]
      val<- fin$out[i]
      match_tar<-fin %>% filter(char_num >= tar_num) %>% filter(grepl( val, out)) %>% filter(Freq == tar_freq) %>% arrange(-char_num) %>% slice(1)
      dat_out_Z_sub2 <-unique(rbind(dat_out_Z_sub2, match_tar))
    }
    
     return(dat_out_Z_sub2)
  }else{
    return(dat_out_Z_sub_char)}
  
} 


###
for(j in 1:length(data_files)){
  dat<- read.csv(data_files[j])
  pattern_out<-pattern_string_v2(dat)
  
  HOR_patterns_all<- as.data.frame(matrix(nrow=0, ncol=7))
  
  for( p in unique(pattern_out$dir)){
        sub<- pattern_out[pattern_out$dir %in% p,]
        for(b in unique(sub$bin)){
                sub2<- sub[sub$bin %in% b,]
                dir_bin=paste("/scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/",p,sep="")
                data_fil<- dir(dir_bin, pattern=c(b), full.names=TRUE)
                for( d in data_fil){ #match with total monomer list
                        orig_dat<- read.table(d)
                        colnames(orig_dat)<- c("V1", "V2")
                        orig_dat_M<- merge(orig_dat,sub2,by.x="V2", by.y="orig", all=T, fill=NA )
                        orig_dat_M_sub<- select(orig_dat_M, c("V1", "letter"))
                        orig_dat_M_sub$letter[is.na(orig_dat_M_sub$letter)] <- "Z"
                        spl_start<- str_split(orig_dat_M_sub$V1, pattern=":", simplify=T)[,2]
                        orig_dat_M_sub$start<-  as.numeric(str_split(spl_start, pattern="-", simplify=T)[,1])
                        pattern<- paste(orig_dat_M_sub[order(orig_dat_M_sub$start),]$letter, collapse="")
                        HOR_patterns<- filtered_kmer_counts(pattern)
                
                        if(nrow(HOR_patterns)>0){
                                HOR_patterns$dir<- p
                                HOR_patterns$bin<- b
                                HOR_patterns_all<- rbind(HOR_patterns_all,HOR_patterns )
                                }
                }
          }
    }
    if(nrow(HOR_patterns_all) > 0 ){
    check<-HOR_patterns_all %>% group_by(out) %>% dplyr::summarize(count=n_distinct(bin)) %>% filter(count >1)
    
    if(nrow(check) > 0 ){
    HOR_patterns_all_sub<- HOR_patterns_all[HOR_patterns_all$out %in% check$out, ]
    HOR_patterns_all_sub_M <- HOR_patterns_all_sub %>% group_by(out, charuni, char_prop_m, char_num) %>% summarise(Freq=sum(Freq)) %>% as.data.frame()
    
    HOR_patterns_all_sub_M_filt<- filtered_kmer_counts_2(HOR_patterns_all_sub_M)
    
    nam_split<- str_split(data_files[j], pattern="[.]|_", simplify=T)
        write.table(HOR_patterns_all_sub, paste(nam_split[1,4], nam_split[1,6], "SHARED_ALL_HOR.out", sep="_") , col.names=F, row.names=F, quote=F)
        
    HOR_patterns_all_sub2<-  HOR_patterns_all_sub[ HOR_patterns_all_sub$out %in% HOR_patterns_all_sub_M_filt$out, ]
    write.table(HOR_patterns_all_sub2, paste(nam_split[1,4], nam_split[1,6], "SHARED_FILT_HOR.out", sep="_") , col.names=F, row.names=F, quote=F)
    }
 }
}


