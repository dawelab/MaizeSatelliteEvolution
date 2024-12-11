##cent

HOR_allStrings_sub_MoCent<- HOR_allStrings_sub[HOR_allStrings_sub$line %in% "Mo17_half_Mo17" & HOR_allStrings_sub$rep %in% "CentC",]
MoCent_patt<- merge(HOR_allStrings_sub_MoCent, mai_ALL , by.x=c("line","bin"), by.y=c("V4", "V5"))
dat<- read.table("~/Desktop/Mo17_TE.bed")
dat$chr<- str_split(dat$V1, pattern="_", simplify = T)[,1]

Cent_arrays<- position_info_plot %>% filter(rep %in% "CentC" & line %in% "Mo17") %>% select(chr,distinct) %>% as.data.frame()
#c(25,  54, 66, 78, 79,  104, 113, 121, 129, 143, 144, 44) %>% as.data.frame()
Chip<- read.table("SRR21509778_SRR21509776_1000_rpkmNORM.bw")
Chip$chr<- str_split(Chip$V1, pattern="_", simplify = T)[,1]
Chip2<- read.table("SRR21509779_SRR21509777_1000_rpkmNORM.bw")
Chip2$chr<- str_split(Chip2$V1, pattern="_", simplify = T)[,1]

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
 
mai_ALL<- read.table("~/Downloads/maize_ALL_SHARED_ALL_HOR.bed")
mai_ALL_coord<- merge(mai_ALL, HOR_allStrings_sub, by.x=c("V5"), by.y=c("bin"))

C_100kb<- as.data.frame(matrix(nrow=0, ncol=5))
for(c in unique(Chip2$chr)){
  C_sub<- Chip[Chip$chr %in% c,]
  C_sub$V2<- as.numeric(C_sub$V2)
  C_sub$V3<- as.numeric(C_sub$V3)
  C2_sub<- Chip2[Chip2$chr %in% c,]
  C2_sub$V2<- as.numeric(C2_sub$V2)
  C2_sub$V3<- as.numeric(C2_sub$V3)
  
  min<-1
  max<- min+100000
  while(min<max(C_sub$V3)){
    C_sub_bin<- C_sub[C_sub$V2>min & C_sub$V3<max,]
    C2_sub_bin<- C2_sub[C2_sub$V2>min & C2_sub$V3<max,]
    
    C_100kb_ROW<- as.data.frame(matrix(nrow=1, ncol=5))
    C_100kb_ROW[1,1]<-c
    C_100kb_ROW[1,2]<-min
    C_100kb_ROW[1,3]<-max
    C_100kb_ROW[1,4]<-mean(C_sub_bin$V4)
    C_100kb_ROW[1,5]<-mean(C2_sub_bin$V4)
    C_100kb<- rbind(C_100kb_ROW, C_100kb)
    
    min<- min+100000
    max<- max+100000
  }

}

for(C in unique(Cent_arrays$chr)){
  Cent_arrays_vals<- Cent_arrays[Cent_arrays$chr %in% C, ]$distinct
  arrays<- as.data.frame(matrix(nrow=0, ncol=8))
  for(i in Cent_arrays_vals){
    arrays_sub<- read.table(paste("~/Desktop/new_groups/all_nam_consensus_mono_", i,".fasta.fa.blat.sub", sep="") )
    arrays_sub<- arrays_sub %>% group_by( V1, V3) %>% slice(which.max(V2)) %>% as.data.frame()
    arrays_sub$chr<- str_split(arrays_sub$V1, pattern="_", simplify=T)[,1]
    arrays_sub$line<- str_split(arrays_sub$V1, pattern="_|:", simplify=T)[,2]
    arrays_sub<- arrays_sub[arrays_sub$line %in% "Mo17",]
    array_start <- str_split(arrays_sub$V1, pattern=":", simplify=T)[,2]
    mono_start <- str_split(arrays_sub$V1, pattern=":", simplify=T)[,3]
    arrays_sub$array_start<- str_split(array_start, pattern = "-", simplify=T)[,1] %>% as.numeric()
    arrays_sub$mono_start<- as.numeric(str_split(mono_start, pattern="-", simplify=T)[,1])
    arrays_sub$start<- arrays_sub$array_start + arrays_sub$mono_start
    arrays<- rbind(arrays, arrays_sub)
  }
  if(C %in% "chr1"){
    start<- 137500000
    end<-139500000
  }else{
    if(C %in% "chr2"){
      start<- 98300000
      end<-98935632
    }else{
      if(C %in% "chr3"){
        start<-87500000
        end<-90000000
      }else{
        if(C %in% "chr4"){
          start<-109000000 
          end<- 112000000
        }else{
          if(C %in% "chr5"){
            start<-102500000
            end<-107500000
          }else{
            if(C %in% "chr6"){
              start<-53500000
              end<-55500000
            }else{
              if(C %in% "chr7"){
                start<- 54000000
                end<- 60000000
              }else{
                if(C %in% "chr8"){
                  start<- 50000000
                  end<- 55000000
                }else{
                  if(C %in% "chr9"){
                    start<- 59500000
                    end<- 62500000
                  }else{
                    if(C %in% "chr10"){
                      start<- 48000000
                      end<-53000000
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  #fig 5, S4, S5
  ggplot()+
    geom_point(data=arrays[arrays$V3 %in% "CentC",], size=.1, aes(x=start, y=V2, color=V3))+
    geom_rect(data= C_100kb[C_100kb$V1 %in% C & C_100kb$V3>=start & C_100kb$V3<=end,], aes(xmin=V2 , xmax=V3, ymax=0, ymin=-.25, alpha=V4/max(C_100kb[C_100kb$V1 %in% C & C_100kb$V3>=start & C_100kb$V4<=end,]$V4, na.rm=T)), fill="darkgreen" )+
    geom_rect(data=  C_100kb[C_100kb$V1 %in% C & C_100kb$V3>=start & C_100kb$V3<=end,], aes(xmin=V2 , xmax=V3, ymax=-.25, ymin=-.5, alpha=V5/max(C_100kb[C_100kb$V1 %in% C & C_100kb$V3>=start & C_100kb$V4<=end,]$V5,  na.rm=T)), fill="darkgreen" )+
    geom_rect(data= dat[dat$chr %in% C,], aes(xmin=V2 , xmax=V3, ymax=-.5, ymin=-.75), fill="darkblue" )+
    geom_rect(data= HOR_allStrings_sub_MoCent[HOR_allStrings_sub_MoCent$chr %in% C & HOR_allStrings_sub_MoCent$String_purity > 0 & !is.na(HOR_allStrings_sub_MoCent$String_purity),], aes(xmin=bin_start , xmax=bin_start+10000, ymax=-.75, ymin=-1, alpha=String_purity), fill="darkred" )+
    geom_rect(data= HOR_allStrings_sub_MoCent[HOR_allStrings_sub_MoCent$chr %in% C & HOR_allStrings_sub_MoCent$all_shared_purity > 0 & !is.na(HOR_allStrings_sub_MoCent$all_shared_purity),], aes(xmin=bin_start , xmax=bin_start+10000, ymax=-1, ymin=-1.25, alpha=all_shared_purity), fill="darkred" )+
    xlim(start,end)+
    scale_colour_manual(values = c("#6699CC"))+
    theme_classic() +ggtitle(C)+ylab("Jacc to consensus")+xlab("Array Position") + theme(legend.position="none")
  
    
  patts<- mai_ALL_coord[mai_ALL_coord$group %in% Cent_arrays_vals & mai_ALL_coord$V4 %in% "Mo17_half_Mo17", ] %>% select(V1,  bin_start)  %>% 
    unique() %>% group_by(V1) %>% summarise(n=n()) %>% filter(n>=2)
  patt_plot<- mai_ALL_coord[mai_ALL_coord$group %in% Cent_arrays_vals & mai_ALL_coord$V4 %in% "Mo17_half_Mo17" & mai_ALL_coord$V1 %in% patts$V1, ]
 
   ggplot() + 
    geom_point(data=patt_plot, size=.1, aes(x=bin_start, y=V1))+ xlim(start,end)+
     theme_classic()

}




