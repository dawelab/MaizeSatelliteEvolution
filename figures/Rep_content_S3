#R plotting

library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)

setwd("~/Desktop")
main_chr<- paste("chr",1:10,sep = "")


#
setwd("~/Desktop/old_rep_content")
files<- as.array(dir(path = ".",pattern="output") ) 

dat<- as.data.frame(matrix(nrow=0, ncol=5))
for( i in files){
  nam<- str_split(i, pattern = "_", simplify = T)[,2]
  d<- read.table(i)
  d$nam<- nam
  dat<- rbind(dat, d)
}

dat$prop<- dat$V3/dat$V4

dat$lab<- dat$V1
dat[dat$V1 %in% c("HiFi_Assem_Chr", "Final_Assem"),]$lab<-"HiFi_Assembly_Chr"
dat[dat$V1 %in% c("HiFi_Assem", "HiFi_Assem_Chr", "HiFi_Contigs"),]$lab<-"HiFi_Contigs"
dat[dat$V1 %in% c("Assem", "assem"),]$lab<-"Published_Assembly"
Mo17_chr<- dat[dat$V1 %in% "assem" & dat$nam %in% "Mo17", ]
Mo17_chr$lab<- "Assem_Chr"
dat<- rbind(dat, Mo17_chr)
dat_sub<- dat[dat$lab %in% c("IL", "ONT", "PB", "HiFi", "Published_Assembly", "Assem_Chr"),]
dat_sub$lab<- as.factor(dat_sub$lab)
dat_sub$lab<- factor(dat_sub$lab, levels = c("IL", "ONT", "PB", "HiFi", "Published_Assembly", "Assem_Chr"))

#S3a
ggplot()+geom_bar(data=dat_sub,aes(x=nam, y=prop,fill=lab), stat = "identity",position = "dodge")+
  facet_wrap(~V2, nrow=4, scales="free_y", strip.position = "right")+theme_classic()+
  scale_fill_manual(values = c("#009E73","#F0E442","#56B4E9", "#0072B2", "grey", "gray40"))

fil<-"B73_assem.bp.p_ctg.noseq.gfa_Filt_repvals_unanchoredEnds"
fil<-"MO17_assem.bp.p_ctg.noseq.gfa_Filt_repvals_unanchoredEnds"
lin<- fil 
print(lin)
rep_cont<- read.table("~/Downloads/B73_assem.bp.p_ctg.noseq.gfa_Filt_repvals_unanchoredEnds")
rep_cont<- read.table("~/Downloads/MO17_assem.bp.p_ctg.noseq.gfa_Filt_repvals_unanchoredEnds")

rep_cont %>% group_by( V5) %>% count()
rep_cont$V3<- str_split(rep_cont$V3, pattern=":", simplify=T)[,3] %>% as.numeric()
rep_cont$V2<- str_split(rep_cont$V2, pattern=":", simplify=T)[,3] %>% as.numeric()
rep_cont_tot_rep<- rep_cont %>% group_by(V1) %>% dplyr::summarise(rep_val=sum(V4)) %>% as.data.frame()
rep_cont_max_rep<- rep_cont %>% group_by(V1) %>% dplyr::top_n(1, V4 )%>% select( c(V1,V5 )) %>%  as.data.frame() %>% dplyr::rename("max_rep" = "V5") 
rep_cont_2<- merge(rep_cont, rep_cont_tot_rep)
rep_cont_info<- merge(rep_cont_2, rep_cont_max_rep)
rep_cont_info_lim<- rep_cont_info %>% dplyr::select(c(V1, V2, V3, rep_val, max_rep)) %>% unique()
rep_cont_info_lim$color<- NA
rep_cont_info_lim[rep_cont_info_lim$max_rep %in% "0",]$color<- "gray"
rep_cont_info_lim[rep_cont_info_lim$max_rep %in% "CentC",]$color<- "#6699CC" 
rep_cont_info_lim[rep_cont_info_lim$max_rep %in% "TR1",]$color<- "#CC6677" 
rep_cont_info_lim[rep_cont_info_lim$max_rep %in% "knob180",]$color<-"#882255"
if(nrow(rep_cont_info_lim[is.na(rep_cont_info_lim$color),])>0){
  rep_cont_info_lim[is.na(rep_cont_info_lim$color),]$color<- "#88CCEE"}
rep_cont_info_lim$len_adj<- NA
for(i in 1:nrow(rep_cont_info_lim)){
  if(rep_cont_info_lim$V2[i]<100000){
    rep_cont_info_lim$len_adj[i]<- 1
  }else{
    if( rep_cont_info_lim$V2[i]>100000 & rep_cont_info_lim$V2[i]<250000 ){
      rep_cont_info_lim$len_adj[i]<- 2
    }else{
      if( rep_cont_info_lim$V2[i]>250000 & rep_cont_info_lim$V2[i]<500000 ){
        rep_cont_info_lim$len_adj[i]<- 3
      }else{
        if( rep_cont_info_lim$V2[i]>500000 ){
          rep_cont_info_lim$len_adj[i]<- 4
        }
      }
    }
  }
}

#S3b&c
#B73 contig plots
ggplot() + geom_jitter(data=rep_cont_info_lim, aes(x= rep_val/V2, y=V3, size=len_adj), colour = rep_cont_info_lim$color) +
  ylab("reads")+xlab("length")+theme_classic() +geom_hline(yintercept=20, linetype="dashed", color="grey")+
  theme(legend.position = "none")+
  ggtitle(paste(lin,"contigs"))+ylim(0,30)+xlim(0,1) +xlab("Proportion ETR") +ylab("Read Depth")
##
#Mo17 contig plots
ggplot() + geom_jitter(data=rep_cont_info_lim, aes(x= rep_val/V2, y=V3, size=len_adj), colour = rep_cont_info_lim$color) +
  ylab("reads")+xlab("length")+theme_classic() +geom_hline(yintercept=66, linetype="dashed", color="grey")+
  theme(legend.position = "none")+
  ggtitle(paste(lin,"contigs"))+ylim(0,80)+xlim(0,1) +xlab("Proportion ETR") +ylab("Read Depth")
##

#S3d
data_files<-grep(dir(path = "~/Desktop",pattern="Filt_repvals_unanchoredEnds$"), pattern='chr_', invert=TRUE, value=TRUE)
Mo17_scaff<- read.table("~/Desktop/scaff_info_Mo172")
Mo17_scaff$cont<- str_split( Mo17_scaff$V2, pattern="_", simplify = T)[,1]
Mo17_scaff_sub<- Mo17_scaff[Mo17_scaff$cont %in% chrs, ]
B73_scaff<- read.table("~/Desktop/scaff_info_B732")
B73_scaff$cont<- str_split( B73_scaff$V2, pattern="_", simplify = T)[,1]
B73_scaff_sub<- B73_scaff[B73_scaff$cont %in% chrs, ]

count_conts<- as.data.frame(matrix(nrow=0, ncol=5))
rep_conts<- as.data.frame(matrix(nrow=0, ncol=4))
for(i in data_files){
  dat<- read.table(i)
  lin<- dat[1,1]
  dat$V3<- str_split( dat$V3, pattern=":", simplify=T)[,3] %>% as.numeric()
  dat$V2<- str_split( dat$V2, pattern=":", simplify=T)[,3] %>% as.numeric()
  lin<- str_split(i, pattern="_", simplify = T)[,1]
  Mo17_s1<- Mo17_scaff_sub[Mo17_scaff_sub$V1 %in% lin,]
  M_Mo17<- merge(dat, Mo17_s1, by.x=c("V1"), by.y=c("V3"))
  M_Mo17_sum<- M_Mo17 %>% group_by(V5) %>% summarise(sum=sum(V4)) %>% as.data.frame()
  M_Mo17_sum$lin<- lin
  M_Mo17_sum$ref<- "Mo17"
  M_Mo17_count<- M_Mo17 %>% select(c("V1", "V2.x")) %>% unique() %>% summarise(n=n(), sum=sum( V2.x))
  M_Mo17_count_unanchor<- M_Mo17 %>% select(c("V1","V6")) %>% unique() %>% filter(V6 ==2) %>% summarise(n=n())
  M_Mo17_count$lin<- lin
  M_Mo17_count$unanchor<- M_Mo17_count_unanchor[1,1]
  M_Mo17_count$ref<- "Mo17"
  
  B73_s1<- B73_scaff_sub[B73_scaff_sub$V1 %in% lin,]
  M_B73<- merge(dat, B73_s1, by.x=c("V1"), by.y=c("V3"))
  M_B73_sum<- M_B73 %>% group_by(V5) %>% summarise(sum=sum(V4)) %>% as.data.frame()
  M_B73_sum$lin<- lin
  M_B73_sum$ref<- "B73"
  M_B73_count<- M_B73 %>% select(c("V1", "V2.x")) %>% unique() %>% summarise(n=n(), sum=sum( V2.x))
  M_B73_count_unanchor<- M_B73 %>% select(c("V1","V6")) %>% unique() %>% filter(V6 ==2) %>% summarise(n=n())
  M_B73_count$lin<- lin
  M_B73_count$unanchor<- M_B73_count_unanchor[1,1]
  M_B73_count$ref<- "B73"
  
  count_conts<- rbind(count_conts,   M_Mo17_count,   M_B73_count)
  rep_conts<- rbind(rep_conts,  M_B73_sum,  M_Mo17_sum)
}

rep_len_conts<- merge(count_conts, rep_conts, by=c("lin", "ref"))
rep_conts_sub<-  rep_conts[  rep_conts$V5 %!in% "0", ]


setwd("~/Downloads/ETR_in_HiFi-main/Additional_assemblies/repeat_content")
path=getwd()
data<- rownames(file.info(dir(path, pattern="output.vals", full.names=TRUE)))
data_all<- as.data.frame(matrix(nrow=0, ncol=5))
for(dat in 1:length(data)){
  i<- data[dat]
  nam<- str_split(i, pattern="_", simplify=T)[,6]
  d<- read.table(i)
  d$nam<- nam
  data_all<- rbind(data_all,d)
}
data_allsub<-data_all[data_all$V1 %in% c("HiFi", "Contigs", "HiFi_Contigs", "HiFi_Assem", "Assem", "assem") & data_all$V2 %in% c("knob180", "CentC", "Cent4","TR1"),]
data_allsub[data_allsub$nam %in% "Ab10",]$nam<- "AB10"
cov<-read.csv("~/Desktop/HiFi_ReadCoverage.csv", header=F)
cov$V2<- as.numeric(str_split(cov$V2, pattern="x", simplify=T)[,1])
cov[cov$V1 %in% "TIL11.2cell.HiFi",]$V1<- "TIL11"
cov[cov$V1 %in% "TIL01.3cell.HiFi",]$V1<- "TIL01"
cov[cov$V1 %in% "Zea-mays-ssp-mexicana-TIL25_4cell-hifi",]$V1<- "TIL25"


m_data_allsub<- merge(data_allsub, cov, by.y="V1",  by.x="nam")
m_data_allsub[m_data_allsub$V1 %in% c("HiFi_Assem","Contigs", "HiFi_Contigs",  "Assem", "assem"),]$V2.y<-1

rep_len_conts_sel<- rep_len_conts %>% filter(V5 %in% c("CentC", "Cent4", "TR1", "knob180")) %>%
  select(c("lin", "ref", "V5", "sum.y", "sum.x","n"))
rep_len_conts_sel$n<- 1
rep_len_conts_sel$ref<- paste(rep_len_conts_sel$ref, "Assem", sep="_")
colnames(rep_len_conts_sel)<- c("nam", "V1", "V2.x", "V3", "V4", "V2.y")
rep_len_conts_sel[rep_len_conts_sel$nam %in% "TIL01.3cell.HiFi" , ]$nam<- "TIL01"
rep_len_conts_sel[rep_len_conts_sel$nam %in% "TIL11.2cell.HiFi" , ]$nam<- "TIL11"
rep_len_conts_sel[rep_len_conts_sel$nam %in% "Zea-mays-ssp-mexicana-TIL25" , ]$nam<- "TIL25"

info<- rbind(m_data_allsub, rep_len_conts_sel)
info[info$V1 %in% "HiFi_Assem", ]$V1<- "Contigs"
info[info$V1 %in% "HiFi_Contigs", ]$V1<- "Contigs"

info$nam<- as.factor(info$nam)
info$nam<- factor(info$nam, levels=c( "B73","AB10" , "Mo17", "CG108", "CG119","CG44","Tx777", "Tx779",
                                      "CML442", "K64","TIL01","TIL11","TIL25"))

info[info$V1 %in% "Assem",]$V1<- "Published Assembly"
info[info$V1 %in% "assem",]$V1<- "Published Assembly"

info[info$V1 %in% "Published Assembly" & info$nam %in% "AB10",]$V3<- info[info$V1 %in% "Published Assembly" & info$nam %in% "AB10",]$V3*1000000

info$V1<- as.factor(info$V1)
info$V1<- factor(info$V1, levels=c( "HiFi" , "Contigs", "B73_Assem","Mo17_Assem", "Published Assembly"))

ggplot()+geom_bar(data=info, stat="identity",position="dodge", aes(x=nam, y=V3/V2.y/1000000, fill=V1))+
  facet_wrap(~V2.x, nrow=4, scales="free_y", strip.position = "right") +
  theme_classic()+scale_fill_manual(values=c("#0072B2", "lightblue", "lightgrey", "darkgrey", "black"))


