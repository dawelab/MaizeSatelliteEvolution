

'%!in%'<- Negate('%in%')

#Fig 1d
#PANgenome repeat arrays
N<- read.table("~/Desktop/all_100_gaps")
lines<- c("AB10_half_Mo17","B73_half_Mo17","CG108_half_Mo17","CG119_fix","CG44_half_Mo17", "AB10_half_AB10",
          "CML442_half_Mo17","K64_half_Mo17", "Mo17_half_Mo17", "Tx777_half_Mo17","Tx779_half_Mo17",
          "TIL01.3cell.HiFi_half_Mo17", "TIL11.2cell.HiFi_half_Mo17", "Zea-mays-ssp-mexicana-TIL25_4cell-hifi_half_Mo17")
N$line<- str_split(N$V1, pattern="_", simplify = T)[,1]
N$chr<- str_split(N$V2, pattern="_", simplify = T)[,1]
N_lin<- N[N$V1 %in% lines, ]

comp_all_combined<- read.table("~/Desktop/new_groups/comp_all_combined_V9")
position_info<- read.table("~/Desktop/grouped_arrays_pos_member", header = T)
position_infoAb10<- read.table("~/Desktop/Ab10_arrays_tail_grouped", header = T)
position_info<-rbind(position_info, position_infoAb10)
position_info$len<- position_info$array_end- position_info$array_start
position_info[position_info$chr %in% "chr10_RagTag",]$chr<-"chr10"
position_info$line<- str_split(position_info$dat, pattern="_", simplify = T)[,1]

position_info[position_info$line %in% "Zea-mays-ssp-mexicana-TIL25",]$line<- "TIL25"
position_info[position_info$line %in% "TIL11.2cell.HiFi",]$line<- "TIL11"
position_info[position_info$line %in% "TIL01.3cell.HiFi",]$line<- "TIL01"

position_info$line<- as.factor(position_info$line)
position_info$line<- factor(position_info$line, levels=c( "B73","AB10" , "Mo17", "CG108", "CG119","CG44","Tx777", "Tx779",
                                                          "CML442", "K64","TIL01","TIL11","TIL25"))

position_info$filt<- NA
position_info[position_info$dat %in% "AB10_half_Mo17" & position_info$chr %in% "chr10" & position_info$array_start > 150000000,]$filt<- "YES"
position_info_plot<- position_info[is.na(position_info$filt),]

position_info_plot$N<- "NO"
for(i in 1:nrow(position_info_plot)){
  N_sub<- N_lin[N_lin$V1 %in% position_info_plot$dat[i] & N_lin$chr %in% position_info_plot$chr[i] & N_lin$V3 > position_info_plot$array_start[i] & N_lin$V3< position_info_plot$array_end[i],]
  if(nrow(N_sub>0)){
    position_info_plot$N[i]<- "YES"
  }
}

position_info_plot %>% select(c("gene_up", "dat", "N")) %>% unique() %>% filter(N =="YES") %>% group_by(dat) %>% summarise(n=n())

#chr lens
chr_len<- read.table("~/Desktop/Mo17_t2t.fna.fai") 
chr_len$chr<- str_split(chr_len$V1, pattern="_", simplify = T)[,1]

m_coord<- max(chr_len$V2)
for (i in main_chr){
  position_info_plot_chr<- position_info_plot[position_info_plot$chr %in% i, ]
  len<- as.numeric(chr_len[chr_len$chr %in% i, ]$V2)
  if(nrow(position_info_plot_chr[position_info_plot_chr$rep %in% "Cent4" & position_info_plot_chr$chr %!in% "chr4" , ])>0 ) {
    position_info_plot_chr[position_info_plot_chr$rep %in% "Cent4" & position_info_plot_chr$chr %!in% "chr4" , ]$rep<- "knob180"}
  if(nrow(position_info_plot_chr[is.na(position_info_plot_chr$ref_start), ])>0 ) {
    position_info_plot_chr[is.na(position_info_plot_chr$ref_start), ]$ref_start<- position_info_plot_chr[is.na(position_info_plot_chr$ref_start), ]$ref_end}
  if(length(unique(position_info_plot_chr$rep)) < 4){
  g<- ggplot() +
  geom_line(data=position_info_plot_chr, aes(x = ref_start, y =line, colour = rep, group= gene_up))+
  geom_segment(aes(x =0, xend = len, y =as.factor("B73"), yend = as.factor("B73")), size=1)+
  geom_point(data=position_info_plot_chr, aes(x = ref_start, y =line,  colour = rep, size=len, shape=N))+
  scale_shape_manual(values=c(19, 1))+ xlim(0,m_coord )+
    geom_point(data=position_info_plot_chr[position_info_plot_chr$rep %in% "knob180" & position_info_plot_chr$distinct %!in% "123",], aes(x = ref_start, y =line, size=len, shape=N),color="#882255" )+
    geom_line(data=position_info_plot_chr[position_info_plot_chr$rep %in% "knob180" & position_info_plot_chr$distinct %!in% "123",], aes(x = ref_start, y =line,  group= gene_up), color="#882255")+
  geom_point(data=position_info_plot_chr[position_info_plot_chr$rep %in% "knob180" & position_info_plot_chr$distinct %!in% "121",], aes(x = ref_start, y =line, size=len, shape=N),color="#882255" )+
  geom_line(data=position_info_plot_chr[position_info_plot_chr$rep %in% "knob180" & position_info_plot_chr$distinct %!in% "121",], aes(x = ref_start, y =line,  group= gene_up), color="#882255")+
  theme_classic()+  theme(legend.position = "none")+
  scale_color_manual(values = c(  "#6699CC" ,"#882255", "#CC6677" ))
  }else{
    g<- ggplot() +
      geom_line(data=position_info_plot_chr, aes(x = ref_start, y =line, colour = rep, group= gene_up))+
      geom_segment(aes(x =0, xend = len, y =as.factor("B73"), yend = as.factor("B73")), size=1)+
      geom_point(data=position_info_plot_chr, aes(x = ref_start, y =line,  colour = rep, size=len, shape=N))+
      scale_shape_manual(values=c(19, 1))+ xlim(0,m_coord )+
      geom_point(data=position_info_plot_chr[position_info_plot_chr$rep %in% "knob180" & position_info_plot_chr$distinct %!in% "123",], aes(x = ref_start, y =line, size=len, shape=N),color="#882255" )+
      geom_line(data=position_info_plot_chr[position_info_plot_chr$rep %in% "knob180" & position_info_plot_chr$distinct %!in% "123",], aes(x = ref_start, y =line,  group= gene_up), color="#882255")+
      geom_point(data=position_info_plot_chr[position_info_plot_chr$rep %in% "knob180" & position_info_plot_chr$distinct %!in% "121",], aes(x = ref_start, y =line, size=len, shape=N),color="#882255" )+
      geom_line(data=position_info_plot_chr[position_info_plot_chr$rep %in% "knob180" & position_info_plot_chr$distinct %!in% "121",], aes(x = ref_start, y =line,  group= gene_up), color="#882255")+
      theme_classic()+  theme(legend.position = "none")+
      scale_color_manual(values = c(  "#88CCEE", "#6699CC" ,"#882255", "#CC6677" ))
  }
  ggsave(g, file=paste(i, "positions.png", sep="_"), width=10, height=2.5)
}

#
AFS_count<- position_info_plot %>% group_by(distinct, rep) %>% summarise(n_lines=n_distinct(dat))
AFS_count_2<- AFS_count %>% group_by(n_lines,rep) %>% summarise(N=n())
AFS_count_all<- position_info_plot %>% group_by(distinct) %>% summarise(n_lines=n_distinct(dat))
AFS_count_2_all<- AFS_count %>% group_by(n_lines) %>% summarise(N=n())
AFS_count_2_all$rep<- "all"
AFS_count_2_plusall<-rbind(AFS_count_2, AFS_count_2_all)

line_count<- position_info_plot %>% group_by(line, rep) %>% summarise(n_arrays=n())
line_count_all<- position_info_plot %>% group_by(line) %>% summarise(n_arrays=n())
line_count_all$rep<- "all"
line_count_plusall<-rbind(line_count, line_count_all)

#1b
ggplot()+geom_bar(data=line_count_plusall[line_count_plusall$rep %in% c("Cent4", "CentC", "knob180", "TR1"),], aes(x=line,y=n_arrays, fill=rep), position="stack", stat = "identity") +
  #facet_wrap(~rep, nrow=5, strip.position = "right")+
  theme_classic()+
  scale_fill_manual(values = c( "#88CCEE", "#6699CC" ,"#882255", "#CC6677" ))

#1c
ggplot()+geom_bar(data=AFS_count_2_plusall[AFS_count_2_plusall$rep %!in% "all",], aes(x=n_lines, y=N, fill=rep), stat = "identity") +
  theme_classic()+
  scale_fill_manual(values = c(  "#88CCEE", "#6699CC" ,"#882255", "#CC6677"))

chrs<- paste("chr", 1:10,sep="")

setwd("~/Desktop")

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

#1a
ggplot() + geom_bar( data=rep_conts[rep_conts$ref %in% "Mo17" & rep_conts$V5 %in% c("CentC", "Cent4", "TR1", "knob180"),], aes(x=lin, y=sum/1000000,fill=V5), stat = "identity")+
  theme_classic()+  scale_fill_manual(values = c( "#88CCEE", "#6699CC" ,"#882255", "#CC6677"))
