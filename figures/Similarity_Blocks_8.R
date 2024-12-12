#block comps

arrays<- read.table(paste("~/Desktop/new_groups/all_nam_consensus_mono_123.fasta.fa.blat.sub", sep="") )
arrays2<- read.table(paste("~/Desktop/new_groups/all_nam_consensus_mono_132.fasta.fa.blat.sub", sep="") )
#arrays3<- read.table(paste("~/Desktop/new_groups/all_nam_consensus_mono_117.fasta.fa.blat.sub", sep="") )
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

ggplot()+geom_point(data=arrays, size=.1, aes(x=start, y=V2, color=V3))+
  facet_wrap(~chr, nrow=3, strip.position = "right") +# ,scales="free_x"
  theme_classic() +ggtitle(val)+ylab("Jacc to consensus")+xlab("Array Position") +scale_colour_manual(values = c("#882255", "#CC6677"))
#####
#bin comps
ggplot()+geom_point(data=arrays[arrays$chr %in% "chr7" & arrays$V3 %in% c("knob180", "TR1"),], size=.1, aes(x=start, y=V2, color=V3))+
  geom_rect(mapping=aes(xmin=158500000, xmax=159500000, ymin=0, ymax=1), colour="black", alpha=0)+
  geom_rect(mapping=aes(xmin=162000000, xmax=163000000, ymin=0, ymax=1), colour="black", alpha=0)+  
  geom_rect(mapping=aes(xmin=173200000, xmax=174200000, ymin=0, ymax=1), colour="black", alpha=0)+          
  theme_classic() +ylab("Jacc to consensus")+xlab("Array Position") +scale_colour_manual(values = c("#882255", "#CC6677"))


ggplot()+geom_point(data=arrays[arrays$chr %in% "chr8" & arrays$V3 %in% c("knob180", "TR1"),], size=.1, aes(x=start, y=V2, color=V3))+
  geom_rect(mapping=aes(xmin=187000000, xmax=188000000, ymin=0, ymax=1), colour="black", alpha=0)+
  geom_rect(mapping=aes(xmin=180750000, xmax=181750000, ymin=0, ymax=1), colour="black", alpha=0)+  
  geom_rect(mapping=aes(xmin=168000000, xmax=169000000, ymin=0, ymax=1), colour="black", alpha=0)+          
  theme_classic() +ylab("Jacc to consensus")+xlab("Array Position") +
  scale_colour_manual(values = c("#882255", "#CC6677"))+xlim(160000000, 195000000)


##Bin Comps
library("viridis")
Chr7_12<- read.table("chr7_bin_12.blat.sub")
Chr7_12$V1<- as.numeric(str_split(Chr7_12$V1, pattern=":|-", simplify=T)[,2])
Chr7_12$V2<- as.numeric(str_split(Chr7_12$V2, pattern=":|-", simplify=T)[,2])
Chr7_12$V3<- as.numeric(Chr7_12$V3)
ggplot()+geom_point(data=Chr7_12, size=.01, aes(x=V1, y=V2, alpha=V3, colour=V3))+
  geom_point(data=Chr7_12[Chr7_12$V3 ==1,], size=.01, aes(x=V1, y=V2, alpha=V3), colour="#7A0403FF")+
  theme_classic() + scale_colour_viridis(option="turbo")+coord_fixed()

Chr7_13<- read.table("chr7_bin_13.blat.sub")
Chr7_13$V1<- as.numeric(str_split(Chr7_13$V1, pattern=":|-", simplify=T)[,2])
Chr7_13$V2<- as.numeric(str_split(Chr7_13$V2, pattern=":|-", simplify=T)[,2])
Chr7_13$V3<- as.numeric(Chr7_13$V3)
ggplot()+geom_point(data=Chr7_13, size=.01, aes(x=V1, y=V2, alpha=V3, colour=V3))+
  geom_point(data=Chr7_13[Chr7_13$V3 ==1,],  size=.01, aes(x=V1, y=V2, alpha=V3), colour="#7A0403FF")+
  theme_classic() + scale_colour_viridis(option="turbo")+coord_fixed()

Chr78_1<- read.table("chr78_bin_1.blat.sub")
Chr78_1$V1<- as.numeric(str_split(Chr78_1$V1, pattern=":|-", simplify=T)[,2])
Chr78_1$V2<- as.numeric(str_split(Chr78_1$V2, pattern=":|-", simplify=T)[,2])
Chr78_1$V3<- as.numeric(Chr78_1$V3)
ggplot()+geom_point(data=Chr78_1, size=.001, aes(x=V1, y=V2, alpha=V3, colour=V3))+
  geom_point(data=Chr78_1[Chr78_1$V3 ==1,],  size=.01, aes(x=V1, y=V2, alpha=V3), colour="#7A0403FF")+
  theme_classic() + scale_colour_viridis(option="turbo")+coord_fixed()

Chr8_12<- read.table("chr8_bin_12.blat.sub")
Chr8_12$V1<- as.numeric(str_split(Chr8_12$V1, pattern=":|-", simplify=T)[,2])
Chr8_12$V2<- as.numeric(str_split(Chr8_12$V2, pattern=":|-", simplify=T)[,2])
Chr8_12$V3<- as.numeric(Chr8_12$V3)
ggplot()+geom_point(data=Chr8_12, size=.01, aes(x=V1, y=V2, alpha=V3, colour=V3))+
  geom_point(data=Chr8_12[Chr8_12$V3 ==1,],  size=.01, aes(x=V1, y=V2, alpha=V3), colour="#7A0403FF")+
  theme_classic() + scale_colour_viridis(option="turbo")+coord_fixed()

Chr8_13<- read.table("chr8_bin_13.blat.sub")
Chr8_13$V1<- as.numeric(str_split(Chr8_13$V1, pattern=":|-", simplify=T)[,2])
Chr8_13$V2<- as.numeric(str_split(Chr8_13$V2, pattern=":|-", simplify=T)[,2])
Chr8_13$V3<- as.numeric(Chr8_13$V3)
ggplot()+geom_point(data=Chr8_13, size=.01, aes(x=V1, y=V2, alpha=V3, colour=V3))+
  geom_point(data=Chr8_13[Chr8_13$V3 ==1,], size=.01, aes(x=V1, y=V2, alpha=V3), colour="#7A0403FF")+
  theme_classic() + scale_colour_viridis(option="turbo")+coord_fixed()
