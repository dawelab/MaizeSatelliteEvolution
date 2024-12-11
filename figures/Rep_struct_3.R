library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)

'%!in%'<- Negate('%in%')

###HOR class summaries

all_norep<- read.table("~/Desktop/All_Raw_norep2")
cov<- read.csv("~/Desktop/HiFi_ReadCoverage.csv", header=F)
cov$V2<- as.numeric(str_split(cov$V2, pattern="x", simplify=T)[,1])
cov[cov$V1 %in% "Zea-mays-ssp-mexicana-TIL25_4cell-hifi",]$V1<- "TIL25"
cov[cov$V1 %in% "TIL11.2cell.HiFi",]$V1<- "TIL11"
cov[cov$V1 %in% "TIL01.3cell.HiFi",]$V1<- "TIL01"
cov[cov$V1 %in% "AB10",]$V1<- "Ab10B73"
all_norep_cov<- merge(all_norep, cov, by="V1")
all_norep_cov[all_norep_cov$V1 %in% "Ab10B73",]$V1<-"AB10"
all_norep_cov$V3_mod<- as.numeric(all_norep_cov$V4)/as.numeric(all_norep_cov$V2.y)
all_norep_cov_sub<- all_norep_cov %>% select("V1", "V2.x", "V3", "V3_mod")
all_norep_cov_sub$data<- "Reads"
all_array<- read.table("~/Downloads/All_Arrays_Total.txt")
all_array$data<- "Chromosomes"
colnames(all_norep_cov_sub)<- colnames(all_array)
all_counts_class<- rbind(all_norep_cov_sub, all_array)

all_counts_class$V1- as.factor(all_counts_class$V1)
all_counts_class$V1<- factor(all_counts_class$V1, levels=c( "B73","AB10" , "Mo17", "CG108", "CG119","CG44","Tx777", "Tx779",
                                                          "CML442", "K64","TIL01","TIL11","TIL25"))
all_counts_class$V3_2<- all_counts_class$V3
all_counts_class[all_counts_class$V3_2 %!in% "HOR", ]$V3_2<- "not HOR"
all_counts_class$V3_2<- as.factor(all_counts_class$V3_2)
all_counts_class$V3_2<- factor(all_counts_class$V3_2, levels=c( "not HOR", "HOR"))

all_counts_class$V3<- as.factor(all_counts_class$V3)
all_counts_class$V3<- factor(all_counts_class$V3, levels=c( "Disorder","Order","HOR"))

#3a
ggplot() +geom_bar(data=all_counts_class[all_counts_class$data %in% "Chromosomes",], aes(x=V1, y=V4, fill=V3), color="black",stat="identity",position = "stack")+
  facet_wrap(~V2, scales="free_y", nrow=4, strip.position = "right")+theme_classic()+ 
  #theme(axis.text.x = element_text(angle=45, hjust=1))+
  scale_fill_manual(values = c( "white","darkgray","black"))

###HOR Pattern Length and Freq
setwd("~/Desktop")
data_files<-dir(path = "~/Desktop", pattern="_HOR_bed.csv")
HOR_all_local<- as.data.frame(matrix(nrow=0, ncol=3))
for(i in data_files){
  d<- read.csv(i)
  d_nam<- str_split(i, pattern="_", simplify = T)[,1]
  d$X<- d_nam
  HOR_all_local<- rbind(HOR_all_local, d[,1:3])
}
HOR_all_local$n_char<- nchar( HOR_all_local$Pattern)
HOR_all_local_patt_freq <- HOR_all_local %>% group_by(bin, Pattern, n_char) %>% summarise(n=n())

#combine with satellite type
data_files<-dir(path = "~/Desktop", pattern="_string_out.csv")
string_all_local<- as.data.frame(matrix(nrow=0, ncol=8))
for(i in data_files){
  d<- read.csv(i)
  d_nam<- str_split(i, pattern="_", simplify = T)[,1]
  d$X<- d_nam
  string_all_local<- rbind(string_all_local, d)
}
string_all_local_sub<- select(string_all_local, c("bin", "rep", "X")) 

HOR_all_local_patt_freq_info<- merge(HOR_all_local_patt_freq , string_all_local_sub, by=c("bin"))
HOR_all_local_patt_freq_info<-HOR_all_local_patt_freq_info[HOR_all_local_patt_freq_info$rep %in% c("CentC", "Cent4", "knob180", "TR1"),]
HOR_all_local_patt_freq_info_mean<- HOR_all_local_patt_freq_info %>% 
  group_by(rep) %>%summarise(mean_count=mean(n_char))
HOR_all_local_patt_freq_info_mean2<- HOR_all_local_patt_freq_info %>% 
  group_by(rep) %>%summarise(mean_n=mean(n))

#3b
ggplot()+geom_bar( data=HOR_all_local_patt_freq_info, aes(x=n_char, fill=rep))+
  geom_vline(data= HOR_all_local_patt_freq_info_mean, aes(xintercept = mean_count), linetype="dashed") + 
  facet_wrap(~rep, nrow=4, scales = "free_y", strip.position = "right") +theme_classic()+
  scale_fill_manual(values = c( "#88CCEE", "#6699CC" ,"#882255", "#CC6677"))+xlim(2,15)

#3c
ggplot()+geom_bar( data=HOR_all_local_patt_freq_info[HOR_all_local_patt_freq_info$n>2,], aes(x=n, fill=rep))+
  geom_vline(data= HOR_all_local_patt_freq_info_mean2, aes(xintercept = mean_n), linetype="dashed") + 
  facet_wrap(~rep, nrow=4, scales = "free_y", strip.position = "right") +theme_classic()+
 c+xlim(2,20)

#group comparisons

group_comp<- read.table("~/Downloads/comp_all_combined.txt")

position_info_plot_n <- position_info_plot %>% 
  select(distinct, line) %>% unique()  %>% group_by(distinct) %>% summarise(n=n())
position_info_plot_ngpaless <- position_info_plot %>% 
  select(distinct, line,  N) %>% unique()  %>% filter(N %in% "NO") %>%
  group_by(distinct) %>% summarise(ngapless=n())

position_info_counts_gapless<- merge(position_info_plot_n , position_info_plot_ngpaless, by="distinct")
group_comp_n<- merge(position_info_counts_gapless, group_comp, by.x="distinct", by.y="gene_up")
cor.test(group_comp_n$max_monolen, group_comp_n$max_HOR)
#group_comp
reg<-lm(formula =  max_HOR ~ max_monolen, 
        data=group_comp_n)                       
#get intercept and slope value 
coeff<-coefficients(reg)           
intercept<-coeff[1] 
slope<- coeff[2] 
#3d
ggplot()+  geom_abline(intercept = intercept, slope = slope, color="gray",  
                       linetype="dashed", size=1)+
geom_point( data=group_comp_n, aes(x=max_monolen,y=  max_HOR, size=max_len, colour=element, alpha=n))+
  theme_classic()+
  scale_colour_manual(values = c( "#88CCEE", "#6699CC" ,"#882255", "#CC6677"))


group_comp_n2<-group_comp_n[!is.na(group_comp_n$matNoGap_mean) & group_comp_n$ngapless>1,]
group_comp_n2[group_comp_n2$matNoGap_mean>1,]$matNoGap_mean<- 1
cor.test(group_comp_n2$max_HOR, group_comp_n2$matNoGap_mean)
#group_comp
reg2<-lm(formula =  jacc_mono_av~max_HOR , 
          data=group_comp_n2)                       
#get intercept and slope value 
coeff2<-coefficients(reg2)           
intercept2<-coeff2[1] 
slope2<- coeff2[2] 
#3e
ggplot()+geom_point( data=group_comp_n2, aes(y=matNoGap_mean,x=  max_HOR, size=max_len, colour=element, alpha=n))+
  theme_classic()+
  geom_abline(intercept = intercept2, slope = slope2, color="gray",  
              linetype="dashed", size=1)+
scale_colour_manual(values = c( "#88CCEE", "#6699CC" ,"#882255", "#CC6677"))
