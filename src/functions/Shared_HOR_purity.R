#!/usr/bin/Rscript

library("dplyr")

path="/scratch/rdp22327/Dawe/scaffolding/Mo17_scaff"
setwd(path)
strings<- read.table("ALL_SHARED_Pattern_String.out", fill=TRUE)
patts<- read.table("/scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/HOR_newpatt/ALL_SHARED_ALL_HOR.bed", fill=TRUE)
colnames(patts)<- c("pattern", "start", "end","line", "bin", "group", "level")

#fix two issue lines from Tx779 that are missing column 7 (repeat type)
hold<- strings[strings$V10 %in% "",7:9]
strings[strings$V10 %in% "",8:10]<- hold
strings[strings$V10 %in% "",7]<- "NA"

teo_lines<- c("TIL01.3cell.HiFi_half_Mo17", "TIL11.2cell.HiFi_half_Mo17", "Zea-mays-ssp-mexicana-TIL25_4cell-hifi_half_Mo17")
maize_lines<- c("AB10_HALF_AB10", "AB10_half_Mo17", "B73_half_Mo17", "CG108_half_Mo17", "CG119_fix",
		"CG44_half_Mo17", "CML442_half_Mo17", "K64_half_Mo17", "Mo17_half_Mo17", "Tx777_half_Mo17", "Tx779_half_Mo17")
		
bin_count<- patts %>% group_by(pattern, group, level ) %>% summarise(n_bins=n_distinct(bin))
teo_count<- patts  %>% group_by(pattern, group, level) %>% filter(line %in% teo_lines) %>% summarise(n_teo_lines=n_distinct(line))
maize_count<- patts  %>% group_by(pattern, group, level) %>% filter(line %in% maize_lines) %>% summarise(n_maize_lines=n_distinct(line))

m1<- merge(bin_count[bin_count$n_bins>1,], teo_count, by=c("pattern", "group", "level"))
m2<- merge(m1, maize_count, by=c("pattern", "group","level"))
patts_counts<- merge(m2, patts, by=c("pattern", "group","level"))

#only Mo17 shared within arrays
bin_multi_Mo17<- patts %>% group_by(pattern, group, level ) %>% filter(line %in% "Mo17_half_Mo17")  %>% summarise(n_bins_Mo17=n_distinct(bin)) %>% filter(n_bins_Mo17>1)
patts_counts_Mo17_multibins<- merge(patts_counts[patts_counts$line %in% "Mo17_half_Mo17",], bin_multi_Mo17, by=c("pattern", "group", "level"))
write.table(patts_counts_Mo17_multibins,"patts_counts_Mo17_multibins", quote=F, col.names=F, row.names=F )

#only CG108 shared within arrays
bin_multi_CG108<- patts %>% group_by(pattern, group, level ) %>% filter(line %in% "CG108_half_Mo17")  %>% summarise(n_bins_CG108=n_distinct(bin)) %>% filter(n_bins_CG108>1)
patts_counts_CG108_multibins<- merge(patts_counts[patts_counts$line %in% "CG108_half_Mo17",], bin_multi_CG108, by=c("pattern", "group", "level"))
write.table(patts_counts_CG108_multibins,"patts_counts_CG108_multibins", quote=F, col.names=F, row.names=F )

Mo17_string<- strings[strings$V2 %in% "Mo17_half_Mo17",] #Mo17 strings
Mo17_string$shared_purity<-0
Mo17_string$teo_purity<-0
Mo17_string$maiz_purity<-0

CG108_string<- strings[strings$V2 %in% "CG108_half_Mo17",] #CG108 strings
CG108_string$shared_purity<-0
CG108_string$teo_purity<-0
CG108_string$maiz_purity<-0

patts_counts_Mo17<- patts_counts[patts_counts$line %in% "Mo17_half_Mo17", ]
patts_counts_CG108<- patts_counts[patts_counts$line %in% "CG108_half_Mo17", ]
for(i in 1:nrow(Mo17_string)){
	new_patt<- Mo17_string[i,10]
	all_shared_HOR<- patts_counts_Mo17[patts_counts_Mo17$bin %in% Mo17_string$V4[i],] %>% select("line", "start", "end")
	all_shared_teo<- patts_counts_Mo17[patts_counts_Mo17$bin %in% Mo17_string$V4[i] & patts_counts_Mo17$n_teo_lines > 0,] %>% select("line", "start", "end")
	all_shared_maize<- patts_counts_Mo17[patts_counts_Mo17$bin %in% Mo17_string$V4[i] & patts_counts_Mo17$n_teo_lines == 0 & patts_counts_Mo17$n_maize_lines >1 ,] %>% select("line", "start", "end")


	if(nrow(all_shared_HOR)>0){
		all_shared_HOR2<- HOR_patt_check(all_shared_HOR, new_patt)
		Mo17_string$shared_purity[i]<- purity_calc(all_shared_HOR2, new_patt)
	}
	if(nrow(all_shared_teo)>0){
		all_shared_teo2<- HOR_patt_check(all_shared_teo, new_patt)
		Mo17_string$teo_purity[i]<- purity_calc(all_shared_teo2, new_patt)
	}
	if(nrow(all_shared_maize)>0){
		all_shared_maize2<- HOR_patt_check(all_shared_maize, new_patt)
		Mo17_string$maiz_purity[i]<- purity_calc(all_shared_maize2, new_patt)
	}
}

write.table(Mo17_string,"Mo17_string", quote=F, col.names=F, row.names=F )



for(i in 1:nrow(CG108_string)){
	new_patt<- CG108_string[i,10]
	all_shared_HOR<- patts_counts_CG108[patts_counts_CG108$bin %in% CG108_string$V4[i],] %>% select("line", "start", "end")
	all_shared_teo<- patts_counts_CG108[patts_counts_CG108$bin %in% CG108_string$V4[i] & patts_counts_CG108$n_teo_lines > 0,] %>% select("line", "start", "end")
	all_shared_maize<- patts_counts_CG108[patts_counts_CG108$bin %in% CG108_string$V4[i] & patts_counts_CG108$n_teo_lines == 0 & patts_counts_CG108$n_maize_lines >1 ,] %>% select("line", "start", "end")


	if(nrow(all_shared_HOR)>0){
		all_shared_HOR2<- HOR_patt_check(all_shared_HOR, new_patt)
		CG108_string$shared_purity[i]<- purity_calc(all_shared_HOR2, new_patt)
	}
	if(nrow(all_shared_teo)>0){
		all_shared_teo2<- HOR_patt_check(all_shared_teo, new_patt)
		CG108_string$teo_purity[i]<- purity_calc(all_shared_teo2, new_patt)
	}
	if(nrow(all_shared_maize)>0){
		all_shared_maize2<- HOR_patt_check(all_shared_maize, new_patt)
		CG108_string$maiz_purity[i]<- purity_calc(all_shared_maize2, new_patt)
	}
}

write.table(CG108_string,"CG108_string", quote=F, col.names=F, row.names=F )


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


HOR_patt_check<-function(as_bed, patt){
  as_bed_or_original<-  as_bed[order(as_bed$start),]
  as_bed_or<- as_bed[order(as_bed$start),]
  
  changes<-1
  while(changes > 0){
    changes<- 0
    if(nrow(as_bed_or)<=2){
      end<-2
    }else{end<- nrow(as_bed_or)-1}
    for(i in 2:end){
      if(as_bed_or[i,2] <= as_bed_or[i-1,3]){ as_bed_or[i,2]<- as_bed_or[i-1,3]+1; changes<- changes+1}
    }
   as_bed_or<- as_bed_or[as_bed_or$end-as_bed_or$start >0,]
  }
  
  return(as_bed_or)
}

