
#shared HOR 
struct_all<- read.table("~/Desktop/all_structure_conserved_wMo17_v2_mod")
colnames(struct_all)<- c("line", "element", "chr", "ar_start","bin_start", "bin_end", "group", "threshold","mon_num", "structure")
struct_all_sub<- struct_all[struct_all$line %in% c("Mo17_half_Mo17","CG108_half_Mo17","AB10_HALF_AB10") , ]
struct_all_sub_HOR<- struct_all_sub[struct_all_sub$structure %in% "HOR",]
struct_all_sub_HOR$bin<- paste(struct_all_sub_HOR$chr, struct_all_sub_HOR$bin_start, sep="_")

#String self & shared
HOR_allStrings<- read.table("~/Downloads/ALL_SHARED_Pattern_String.out", fill=NA)
HOR_allStrings_sub<- HOR_allStrings[HOR_allStrings$V2 %in% c("Mo17_half_Mo17", "CG108_half_Mo17", "AB10_HALF_AB10"),]
colnames(HOR_allStrings_sub)<- c("line_bin","line", "chr", "bin", "thresh", "bin_start","rep", "String", "String_purity","String_shared")
HOR_allStrings_sub$threshold<- HOR_allStrings_sub$thresh*100

teo_ALL<- read.table("~/Downloads/teosinte_ALL_SHARED_ALL_HOR.bed")
mai_ALL<- read.table("~/Downloads/maize_ALL_SHARED_ALL_HOR.bed")
ALL_HOR<- rbind(teo_ALL, mai_ALL)
colnames(ALL_HOR)<-c("patt", "start", "end","line","bin","threshold", "group")

ALL_HOR_grouped<- ALL_HOR %>% group_by(threshold, group, patt, line) %>% summarise()
ALL_HOR_grouped_all<- ALL_HOR_grouped %>% group_by(threshold, group, patt) %>% summarise(n=n()) %>% 
  filter(n>1)

ALL_HOR_grouped_teo<- ALL_HOR_grouped %>% group_by(threshold, group, patt) %>% 
  filter(line %in% c("TIL11.2cell.HiFi_half_Mo17", "TIL01.3cell.HiFi_half_Mo17", "Zea-mays-ssp-mexicana-TIL25_4cell-hifi_half_Mo17")) %>%
  summarise(n_teo=n())

ALL_HOR_grouped_mai<- ALL_HOR_grouped %>% group_by(threshold, group, patt) %>% 
  filter(line %!in% c("TIL11.2cell.HiFi_half_Mo17", "TIL01.3cell.HiFi_half_Mo17", "Zea-mays-ssp-mexicana-TIL25_4cell-hifi_half_Mo17")) %>%
  summarise(n_mai=n())

struct_all_CentC<- struct_all %>% filter(line %in% "Mo17_half_Mo17" & element %in% "CentC")
share_counts <- mai_ALL %>% filter(V4 %in% "Mo17_half_Mo17" & V7 %in% struct_all_CentC$group) %>% select(V1, V5, V6, V7) %>% 
  unique() %>% group_by(V1, V6, V7) %>% summarise(n_shared=n()) %>% filter(n_shared>1)
M_mai_ALL<- merge(mai_ALL,share_counts, by=c("V1", "V6", "V7" ) ) %>% filter(V4 %in% "Mo17_half_Mo17") %>%
  select(V5) %>% unique()

nrow(M_mai_ALL)/nrow(struct_all_CentC)

ALL_HOR_grouped_M1<- merge(ALL_HOR_grouped_all, ALL_HOR_grouped_teo, by=c("threshold", "group", "patt"), all=T)
ALL_HOR_grouped_M2<- merge(ALL_HOR_grouped_M1,ALL_HOR_grouped_mai,by=c("threshold", "group", "patt"), all=T  )
ALL_HOR_counts<- merge(ALL_HOR, ALL_HOR_grouped_M2,by=c("threshold", "group", "patt") , all=T )
ALL_HOR_counts[is.na(ALL_HOR_counts)] <- 0

adjusted_total<- function(bed, nchar){
  bed_sort<- bed[order(bed$start),]
  for(j in 2:nrow(bed_sort)){
    if(bed_sort$start[j]<=bed_sort$end[j-1]){
      bed_sort$start[j]<-bed_sort$end[j-1]+1
    }
  }
  bed_sort$tot<- bed_sort$end- bed_sort$start+1
  return(sum(bed_sort[bed_sort$tot>0,]$tot)/nchar)
}

HOR_allStrings_sub$teo_shared_purity<- NA
HOR_allStrings_sub$mai_shared_purity<- NA
HOR_allStrings_sub$all_shared_purity<- NA
HOR_allStrings_sub$group<- NA
for(i in 1:nrow(HOR_allStrings_sub)){
  patts<- ALL_HOR_counts[ALL_HOR_counts$line %in% HOR_allStrings_sub$line[i] & ALL_HOR_counts$bin %in% HOR_allStrings_sub$bin[i],]
  if(nrow(patts)>0){
    HOR_allStrings_sub$group[i]<- patts$group[1]
    
    HOR_allStrings_sub$all_shared_purity[i]<- adjusted_total(patts[order(patts$start),], nchar(HOR_allStrings_sub$String_shared[i]))
    if(nrow(patts[order(patts$start),] %>% filter(n_mai>0)) >0 ){
      HOR_allStrings_sub$mai_shared_purity[i]<- adjusted_total(patts[order(patts$start),] %>% filter(n_mai>0), nchar(HOR_allStrings_sub$String_shared[i]))
    }
    if(nrow(patts[order(patts$start),] %>% filter(n_teo>0)) >0 ){
      HOR_allStrings_sub$teo_shared_purity[i]<- adjusted_total(patts[order(patts$start),] %>% filter(n_teo>0), nchar(HOR_allStrings_sub$String_shared[i]))
    }
  }
}

HOR_allStrings_sub$String_purity<- as.numeric(HOR_allStrings_sub$String_purity)
for(i in 1:nrow(HOR_allStrings_sub)){
  if(HOR_allStrings_sub$String_purity[i]>1 ){
    HOR_allStrings_sub$String_purity[i]<-1
  }
}

for(i in 1:nrow(HOR_allStrings_sub)){
  if(HOR_allStrings_sub$all_shared_purity[i] > 0  & !is.na(HOR_allStrings_sub$all_shared_purity[i] )){
    
  if(HOR_allStrings_sub$all_shared_purity[i] > HOR_allStrings_sub$String_purity[i] ){
    HOR_allStrings_sub$all_shared_purity[i] <-  HOR_allStrings_sub$String_purity[i]
  }
  }
  
}


###Mega Knob
val<-"123"
dat<- read.table("~/Desktop/array_comps_Mo17/BigChr7_TE.bed")
CG108_HOR<- read.csv("~/Downloads/CG108_string_out.csv")
arrays<- read.table(paste("~/Desktop/new_groups/all_nam_consensus_mono_123.fasta.fa.blat.sub", sep="") )
arrays<- arrays %>% group_by( V1, V3) %>% slice(which.max(V2)) %>% filter(V3 %!in% "Cent4") %>% as.data.frame()
arrays$chr<- str_split(arrays$V1, pattern="_", simplify=T)[,1]
arrays$line<- str_split(arrays$V1, pattern="_|:", simplify=T)[,2]

array_start <- str_split(arrays$V1, pattern=":", simplify=T)[,2]
mono_start <- str_split(arrays$V1, pattern=":", simplify=T)[,3]
arrays$array_start<- str_split(array_start, pattern = "-", simplify=T)[,1] %>% as.numeric()
arrays$mono_start<- as.numeric(str_split(mono_start, pattern="-", simplify=T)[,1])
arrays$start<- arrays$array_start + arrays$mono_start
arrays$structure<- NA
arrays$threshold<- NA
arrays[arrays$line %in% "Zea-mays-ssp-mexicana-TIL25",]$line<-"TIL25"
arrays$line<- as.factor(arrays$line)
arrays$line<- factor(arrays$line, levels=rev(c( "B73","AB10" , "Mo17", "CG108", "CG119","CG44","Tx777", "Tx779",
                                                "CML442", "K64","TIL01","TIL11","TIL25")))
##########NEED N gaps######
N_lin_sub<- N_lin[N_lin$line %in% unique(arrays$line) & N_lin$chr %in% "chr7",]
arrays_coords<- arrays %>% group_by(line) %>% summarise(min=min(start), max=max(start), n=n())
arrays_coords$len<- arrays_coords$max-arrays_coords$min
arrays_coords


N_lin_sub$incl<- "NO"
for(i in 1:nrow(N_lin_sub)){
  arrays_coords_2<- arrays_coords[arrays_coords$line %in% N_lin_sub$line[i] & arrays_coords$min <  N_lin_sub$V3[i] & arrays_coords$max>  N_lin_sub$V3[i], ]
  if(nrow(arrays_coords_2)>0){
    N_lin_sub$incl[i]<- "YES"
  }
}
N_lin_sub2<-  N_lin_sub[ N_lin_sub$incl %in% "YES",]
N_lin_sub2$line<- as.factor(N_lin_sub2$line)
N_lin_sub2$line<- factor(N_lin_sub2$line, levels=rev(c( "B73","AB10" , "Mo17", "CG108", "CG119","CG44","Tx777", "Tx779",
                                                        "CML442", "K64","TIL01","TIL11","TIL25")))
#6b
ggplot()+geom_point(data=arrays, size=.1, aes(x=start, y=V2, color=V3))+
  geom_point(data=N_lin_sub2, size=1, aes(x=V3, y=1), color="black")+
  facet_wrap(~line, nrow=13, strip.position = "right") +# ,scales="free_x"
  theme_classic() +ggtitle(val)+ylab("Jacc to consensus")+xlab("Array Position") +scale_colour_manual(values = c("#882255", "#CC6677"))


arrays<- read.table(paste("~/Desktop/new_groups/all_nam_consensus_mono_123.fasta.fa.blat.sub", sep="") )
arrays2<- read.table(paste("~/Desktop/new_groups/all_nam_consensus_mono_132.fasta.fa.blat.sub", sep="") )
arrays3<- read.table(paste("~/Desktop/new_groups/all_nam_consensus_mono_117.fasta.fa.blat.sub", sep="") )
arrays<- rbind(arrays, arrays2, arrays3)

arrays<- arrays %>% group_by( V1, V3) %>% slice(which.max(V2)) %>% filter(V3 %!in% "Cent4") %>% as.data.frame()
arrays$chr<- str_split(arrays$V1, pattern="_", simplify=T)[,1]
arrays$line<- str_split(arrays$V1, pattern="_|:", simplify=T)[,2]

arrays <- arrays %>% filter(line %in% "CG108")

array_start <- str_split(arrays$V1, pattern=":", simplify=T)[,2]
mono_start <- str_split(arrays$V1, pattern=":", simplify=T)[,3]
arrays$array_start<- str_split(array_start, pattern = "-", simplify=T)[,1] %>% as.numeric()
arrays$mono_start<- as.numeric(str_split(mono_start, pattern="-", simplify=T)[,1])
arrays$start<- arrays$array_start + arrays$mono_start

#7b
ggplot()+geom_point(data=arrays, size=.1, aes(x=start, y=V2, color=V3))+
  facet_wrap(~chr, nrow=3, strip.position = "right") +
  theme_classic() +ggtitle(val)+ylab("Jacc to consensus")+xlab("Array Position") +scale_colour_manual(values = c("#882255", "#CC6677"))


####
mai_ALL_knob<- mai_ALL %>% filter(V4  %in% "CG108_half_Mo17") %>% filter(V7 %in% c(123))

HOR_allStrings_sub_knob<- HOR_allStrings_sub %>% filter(line  %in% "CG108_half_Mo17") %>% filter(group %in% c(123))

mai_ALL_shared_patts<- mai_ALL_knob%>% group_by(V1, V6) %>% select(V1, V5, V6, V7) %>%  unique() %>% summarise(n_bin=n()) %>% filter(n_bin >1)
mai_ALL_shared_bins<-mai_ALL_knob %>% filter(V1 %in% mai_ALL_shared_patts$V1) %>% select(V5) %>% unique()
#nrow(mai_ALL_shared_bins) = 234
simi_bins<- data.frame(matrix(nrow=0,ncol=3))

for (i in 1:234){
  for (j in (i+1):234){
    sub<- mai_ALL_knob[mai_ALL_knob$V5 %in% c(mai_ALL_shared_bins[i,1], mai_ALL_shared_bins[j,1]),]
    sub_counts<- sub %>% select(V1, V6, V5) %>% unique() %>% group_by(V1, V6) %>% summarise(n_bins=n())
    if(nrow(sub_counts)>0){
      b<- data.frame(matrix(nrow=1,ncol=3))
      b[1,1]<- mai_ALL_shared_bins[i,1]
      b[1,2]<- mai_ALL_shared_bins[j,1]
      b[1,3]<- nrow( sub_counts[ sub_counts$n_bins >1,])/nrow( sub_counts)
      simi_bins<- rbind(simi_bins, b)
    }
  }
}

simi_bins$b1<- as.numeric(str_split(simi_bins$X1, pattern="_", simplify = T)[,2])
simi_bins$b2<- as.numeric(str_split(simi_bins$X2, pattern="_", simplify = T)[,2])

#6ci
ggplot()+geom_point(data=simi_bins, size=.1, aes(x=b1, y=b2, alpha=X3))+
  geom_point(data=simi_bins, size=.1, aes(x=b2, y=b1, alpha=X3))+
  theme_classic() 

##
mai_ALL_shared_patts<- mai_ALL_knob%>% group_by(V1, V6) %>% select(V1, V5, V6, V7) %>%  unique() %>% summarise(n_bin=n()) %>% filter(n_bin >1)
mai_ALL_shared_patts_allcount<- mai_ALL_knob%>% group_by(V1, V6) %>% summarise(n_bin=n()) %>% filter(n_bin >1)

dat<- read.table("~/Desktop/CG108_TE.bed")
dat$chr<- str_split(dat$V1, pattern="_", simplify = T)[,1]

arrays_chr7<- arrays[arrays$chr %in% "chr7",]
HOR_allStrings_sub_knob$String_purity <- as.numeric(HOR_allStrings_sub_knob$String_purity)
HOR_allStrings_sub_knob$all_shared_purity<- as.numeric(HOR_allStrings_sub_knob$all_shared_purity)

#6civ
ggplot()+
  geom_point(data=arrays_chr7, size=.1, aes(x=start, y=V2, color=V3))+
  geom_rect(data= dat[dat$chr %in% "chr7" & dat$V3-dat$V2>=1000,], aes(xmin=V2 , xmax=V3, ymax=0, ymin=-.25), fill="darkblue" )+
  geom_rect(data=CG108_HOR[CG108_HOR$chr %in% "chr7" & CG108_HOR$purity>0,], aes(xmin=start, xmax=start+10000, ymax=-.25, ymin=-.5,alpha=purity), fill="darkred")+
  geom_rect(data= HOR_allStrings_sub_knob[HOR_allStrings_sub_knob$all_shared_purity > 0 ,], aes(xmin=bin_start , xmax=bin_start+10000, ymax=-.5, ymin=-.75, alpha=all_shared_purity), fill="darkred" )+
  xlim(min(arrays_chr7$start),180000000)+
  scale_colour_manual(values = c("#882255", "#CC6677"))+
  theme_classic() +ggtitle(val)+ylab("Jacc to consensus")+xlab("Array Position") 


HOR_knob_coord<- HOR_knob %>% filter(V1 %in% c("ADBA", "AADB", "BAAD", "ADB", "ADJ", "DBA", "DBR","DBEL", "BEL","DBE")) %>%
  group_by(V1, V6) %>% select(V1, chr, V5, V6, V7, start) %>%  unique()

HOR_knob_coord$V1<- as.factor(HOR_knob_coord$V1)
HOR_knob_coord$V1<- factor(HOR_knob_coord$V1, 
                           levels=c("DBE", "BEL", "DBEL","DBR", "DBA","ADJ","ADB","BAAD","AADB" ,"ADBA"),
                           ordered=T)
#6cii / 7b
ggplot()+
  geom_point(data=HOR_knob_coord, size=1, aes(x=start, y=V1))+
  facet_wrap(~chr, strip.position = "right", nrow = 3)+ #, scales="free_x"
  theme_classic() 


