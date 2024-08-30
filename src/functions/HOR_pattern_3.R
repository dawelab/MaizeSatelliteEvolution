#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)
bins_fil= args[1]
dir= args[2]

setwd(dir)
reval<- read.csv("HOR_re-eval_clusters.csv")
bins<- read.table(bins_fil)

#setting up other data needed
data_files<- dir(getwd(), pattern=c( ".blat_bins_nam_groups.csv"), full.names=FALSE)
`%!in%` = Negate(`%in%`)
clust_nams<- c(LETTERS[1:25], 1:9, letters[1:25] ) ##set characters here for the pattern strings. monomers from private clusters are all called "Z"
priv_lab<- "Z"

#library("kmer")
library("stringr")
library("dplyr")

###FUNCTIONS
##converts monomer grouping and cluster information into a pattern string
pattern_string<- function(dat){
  clust_uni<- dat %>% group_by(clust) %>% dplyr::summarize(count=n())
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
  dat_sub$start<- as.numeric(str_split(dat_sub$nam, pattern = ":|-", simplify=T)[,2])
  dat_sub_or<- dat_sub[order(dat_sub$start),]
  dat_sub_or$num<-1:nrow(dat_sub_or)
  dat_sub_or$clust<- as.factor(dat_sub_or$clust)
  nam<-str_split(data_files[j], pattern=".csv", simplify=T)[,1]
  
  dat_sub_or$letter<- NA
  for(i in 1:nrow(dat_sub_or)){
    dat_sub_or$letter[i]<-dat_clust_rename[dat_clust_rename$`sort(unique(dat$clust))` == dat_sub_or$clust[i], ]$letters
  }
 return( select(dat_sub_or, c("nam", "letter")))   
 # return( paste(dat_sub_or$letter, collapse=""))
}

##k-mer counts in a pattern
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

##filter the k-mer set down
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
    
      
    return(dat_out_Z_sub2[,1:2])
  }else{
    return(dat_out_Z_sub_char[,1:2])}
  
}

bed_convert<- function(bed_fil, patt){
  as_bed<- as.data.frame(matrix(nrow=0, ncol=3))
  for(i in 1:nrow(bed_fil)){ 
    coords<- unlist(gregexpr(bed_fil[i,1], patt))
    co_sub<- as.data.frame((matrix(nrow=length(coords), ncol=3)))
    co_sub$V1<- as.character(bed_fil[i,1])
    co_sub$V2<- coords
    co_sub$V3<- co_sub$V2+nchar(co_sub$V1)-1
    as_bed<- rbind(as_bed, co_sub)
  }
  return(as_bed)
}


HOR_patt_check<-function(bed_fil, patt){
  as_bed<- bed_convert(bed_fil, patt)
  as_bed_or_original<-  as_bed[order(as_bed$V2),]
  as_bed_or<- as_bed[order(as_bed$V2),]
  
  changes<-1
  while(changes > 0){
    changes<- 0
    if(nrow(as_bed_or)<=2){
      end<-2
    }else{end<- nrow(as_bed_or)-1}
    for(i in 2:end){
      if(as_bed_or[i,2] <= as_bed_or[i-1,3]){ as_bed_or[i,2]<- as_bed_or[i-1,3]+1; changes<- changes+1}
    }
   as_bed_or<- as_bed_or[as_bed_or$V3-as_bed_or$V2 >0,]
  }
  
  as_bed_or$chr_num<- nchar(as_bed_or$V1)
  as_bed_or$patt_coord<-as_bed_or$V3-as_bed_or$V2+1
  as_bed_or2<- as_bed_or[as_bed_or$patt_coord >= as_bed_or$chr_num*.5,] %>% group_by(V1) %>% filter(n() >= 1) %>% dplyr::count(V1)
  
  as_bed_2<- bed_convert(as_bed_or2, patt) 
  as_bed_2_or<- as_bed_2[order(as_bed_2$V2),]
  
  return(as_bed_2_or)
}

purity_calc<- function(as_bed_2_or, patt){
  as_bed_int<- as.data.frame(matrix(nrow=0, ncol=2))
  as_bed_2_or[,2]<- as.numeric( as_bed_2_or[,2])
  as_bed_2_or[,3]<- as.numeric( as_bed_2_or[,3])
  as_bed_int[1,1]<- as_bed_2_or[1,2]
  as_bed_int[1,2]<- as_bed_2_or[1,3]
  start<-1
  for(i in 2:nrow(as_bed_2_or)){
    if(as_bed_2_or[i,2] <= as_bed_int[start,2] ){
      as_bed_int[start,2] <- as_bed_2_or[i,3]
    }else{
      start<- start+1
      as_bed_int[start,1]<-as_bed_2_or[i,2]
      as_bed_int[start,2]<-as_bed_2_or[i,3]
      }
  }
  as_bed_int$V3<- as.numeric( as_bed_int$V2)- as.numeric(as_bed_int$V1)+1
 len<-  sum(as_bed_int$V3)
 return(len/nchar(patt))
}
#####


string_out<- as.data.frame(matrix(nrow=length(data_files), ncol=7))
colnames(string_out)<-c("chr","bin","clust_val","start", "rep","pattern", "purity")
HOR_bed_filt_all<- as.data.frame(matrix(nrow=0, ncol=4))
for(j in 1:length(data_files)){
  dat<- read.csv(data_files[j]) #brining in that clustering info
  nam_split<- str_split(data_files[j], pattern = "_|[.]", simplify=T)
  string_out$bin[j]<-paste(nam_split[1,3], nam_split[1,7], sep="_") #monomer name contains cluster and coordinate info that we need
  string_out$chr[j]<-nam_split[1,3]    
#if(string_out$chr[j] == "chr10") { #for AB10_HALF_AB10
  string_out$clust_val[j]<-reval[reval$bin %in% string_out$bin[j], ]$clust_val #labeleding monomer by cluster
  string_out$start[j]<-reval[reval$bin %in% string_out$bin[j], ]$start #labeling monomer by coordinate
  string_out$rep[j]<-as.character(reval[reval$bin %in% string_out$bin[j], ]$rep) #make sure repeat type is the sam
  pattern_out<-pattern_string(dat) #this converts those monomer cluster labels to letters and puts the monomers in order based on coordinate. Now we have a character string!
  write.table(pattern_out, paste(nam_split[1,3], nam_split[1,7], string_out$rep[j], ".out", sep="_") , col.names=F, row.names=F, quote=F) #output those strings
  string_out$pattern[j]<-paste(pattern_out$letter, collapse="")
    HOR_patterns<- filtered_kmer_counts(string_out$pattern[j]) #reports all the kmers
    if(nrow(HOR_patterns) > 0 ){
      HOR_bed_filt<- HOR_patt_check(HOR_patterns, string_out$pattern[j]) #filtering kmers and checking
    }else{
      HOR_bed_filt<- HOR_patterns
    }  
    if(nrow( HOR_bed_filt) >0 ){
      string_out$purity[j]<-purity_calc( HOR_bed_filt,  string_out$pattern[j]) #purity is # mono in any patter/total
      HOR_bed_filt$bin<- paste(string_out$chr[j],"_", string_out$start[j], sep="")
      HOR_bed_filt<-cbind(HOR_bed_filt$bin, HOR_bed_filt$V1, HOR_bed_filt$V2, HOR_bed_filt$V3)
      colnames(HOR_bed_filt)<-c("bin", "Pattern", "Start", "End")
      HOR_bed_filt_all<- rbind(HOR_bed_filt_all, HOR_bed_filt)
    }else{ 
      string_out$purity[j]<-0}
}
#}
write.csv( HOR_bed_filt_all, "HOR_bed.csv", quote=F)
write.csv( string_out,"string_out.csv", quote=F)

