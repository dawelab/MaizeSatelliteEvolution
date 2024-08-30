##Conserved_Rep_Arrays  w/teosinte

setwd("~/Desktop")

library("dplyr")
library("ggplot2")
library("stringr")
library(tidyverse)
library(igraph)

#reading in and prepping df
"%!in%" <- Negate("%in%")
chr_lens<- read.table("Mo17_t2t.fna.fai")
chr_lens$chr<- chr_lens$V1
chr_order<- paste("chr", 1:10, sep="")
gaps<- read.table("all_100_gaps")

arrays<- read.table("Mo17_scaff_array_gen_cords2")
arrays$V2<- str_split(arrays$V2, pattern="_", simplify=T)[,1]
colnames(arrays)<- c("dat", "chr","array_start", "array_end", "array_len", "element",
                  "gene_up_end", "gene_up","gene_down_start", "gene_down")

arrays$dat<- as.factor(arrays$dat)
arrays$dat<-factor(arrays$dat, levels=c("Mo17_half_Mo17","AB10_half_Mo17", "B73_half_Mo17", "CG108_half_Mo17", 
                    "CG119_half_Mo17", "CG119_fix","CG44_half_Mo17", "CML442_half_Mo17", "K64_half_Mo17",
                    "Tx777_half_Mo17", "Tx779_half_Mo17", "TIL01.3cell.HiFi_half_Mo17", 
                    "TIL11.2cell.HiFi_half_Mo17", "Zea-mays-ssp-mexicana-TIL25_4cell-hifi_half_Mo17"))
arrays$gene_up<- str_split(arrays$gene_up, pattern = ",", simplify=T)[,1] 
arrays$gene_down<- str_split(arrays$gene_down, pattern = ",", simplify=T)[,1]

# merging list of conserved genes with coordinates with arrays to get plooting coordaintes projected onto the reference
ref_pos<- read.table("Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.gff3.txt.bed")
colnames(ref_pos)<- c("chr", "ref_start", "ref_end", "gene")
ref_pos_up<- ref_pos[,c("ref_end", "gene")]
ref_pos_down<- ref_pos[,c("ref_start", "gene")]
arrays_pos_up<- merge(arrays, ref_pos_up, by.x="gene_down", by.y="gene", all.x=T)
arrays_pos<- merge(arrays_pos_up, ref_pos_down, by.x="gene_up", by.y="gene", all.x=T)
arrays_pos<- arrays_pos[arrays_pos$dat %!in% "CG119_half_Mo17", ]

arrays_pos$group<- paste(arrays_pos$gene_up, arrays_pos$gene_down, arrays_pos$chr, sep="_")
arrays_pos %>% group_by(gene_up) %>% summarise(n=n()) #152 #1794 empty
arrays_pos %>% group_by(gene_down) %>% summarise(n=n()) #155 #529 empty
groups<- arrays_pos %>% group_by(gene_up,gene_down, chr) %>% summarise(n=n()) #173
arrays_pos %>% group_by(group) %>% summarise(n=n()) #173

arrays_pos[arrays_pos$gene_up %in% ".", ]$gene_up<- paste(arrays_pos[arrays_pos$gene_up %in% ".", ]$gene_down, "edit", sep="_")
arrays_pos[arrays_pos$gene_down %in% ".", ]$gene_down<- paste(arrays_pos[arrays_pos$gene_down %in% ".", ]$gene_up, "edit", sep="_")

#edge list for arrays. if arrays share up or down gene, give edge. This allows us the group by shared up OR shared down genes
arrays_pos$unique_name<- paste(arrays_pos$dat, arrays_pos$array_start, arrays_pos$chr, sep="_" )
arrays_pos_chr_edges<- as.data.frame(matrix(nrow=0, ncol=2))
for(i in 1:nrow(arrays_pos)){
  sub<- arrays_pos[arrays_pos$gene_up %in% arrays_pos$gene_up[i] | arrays_pos$gene_down %in% arrays_pos$gene_down[i],]
  new_edges<- as.data.frame(cbind(arrays_pos$unique_name[i], sub$unique_name))
  arrays_pos_chr_edges<- rbind(arrays_pos_chr_edges, new_edges)
  if(i%%1000 ==0){
    print(i)
}}

#cluster, group based on cluster membership. 
g1<- graph_from_edgelist(as.matrix(arrays_pos_chr_edges))
#plot(g1)
V(g1)$membership<- components(g1)$membership
membership<- as.data.frame(components(g1)$membership)
membership$unique_name<- row.names(membership)
arrays_pos_member<- merge(arrays_pos, membership, by="unique_name")
arrays_pos_member$label<- arrays_pos_member$"components(g1)$membership" ##150 groups
chr_arrays_pos_member<- arrays_pos_member %>% group_by(label) %>% summarise(n=length(unique(chr)))

##now group arrays within the same region (i.e. if two arrays in Mo17 sit between the same two conserved genes, merge into one)
arrays_pos_member$rep<- arrays_pos_member$element
grouped_arrays_pos_member<-arrays_pos_member %>% group_by(gene_up, gene_down, dat, chr) %>% 
  summarise(array_start=min(array_start), array_end=max(array_end),
  rep=unique(rep), gene_up_end=min(gene_up_end), 
  gene_down_start=min(gene_down_start), ref_end=min(ref_end),
  ref_start=min(ref_start), distinct=unique(label )) %>% as.data.frame()
write.table(arrays_pos_member, "arrays_pos_member", col.names=T, row.names=F, quote=F)
write.table(grouped_arrays_pos_member, "grouped_arrays_pos_member", col.names=T, row.names=F, quote=F)

#####

#repeat with B73
B73_arrays<- read.table("~/Desktop/B73_scaff_array_gen_cords")
B73_arrays$V2<- str_split(B73_arrays$V2, pattern="_", simplify=T)[,1]
colnames(B73_arrays)<- c("dat", "chr","array_start", "array_end", "array_len", "element",
                     "gene_up_end", "gene_up","gene_down_start", "gene_down")

B73_arrays$dat<- as.factor(B73_arrays$dat)
B73_arrays$dat<-factor(B73_arrays$dat, levels=c("Mo17_out_HALF","AB10_out_HALF","AB10_HALF_AB10", "B73_out_HALF", "CG108_out_HALF", 
                                        "CG119_out_HALF","CG44_out_HALF", "CML442_out_HALF", "K64_out_HALF",
                                        "Tx777_out_HALF", "Tx779_out_HALF", "TIL01.3cell.HiFi_out_HALF", 
                                        "TIL11.2cell.HiFi_out_HALF", "Zea-mays-ssp-mexicana-TIL25_4cell-hifi_out_HALF"))
B73_arrays$gene_up<- str_split(B73_arrays$gene_up, pattern = ",", simplify=T)[,1] 
B73_arrays$gene_down<- str_split(B73_arrays$gene_down, pattern = ",", simplify=T)[,1]


B73_ref_pos<- read.table("~/Desktop/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.txt.bed")
colnames(B73_ref_pos)<- c("chr", "ref_start", "ref_end", "gene")
B73_ref_pos_up<- B73_ref_pos[,c("ref_end", "gene")]
B73_ref_pos_down<- B73_ref_pos[,c("ref_start", "gene")]
B73_arrays_pos_up<- merge(B73_arrays, B73_ref_pos_up, by.x="gene_down", by.y="gene", all.x=T)
B73_arrays_pos<- merge(B73_arrays_pos_up, B73_ref_pos_down, by.x="gene_up", by.y="gene", all.x=T)
B73_arrays_pos<- B73_arrays_pos[B73_arrays_pos$dat %!in% "AB10_out_HALF", ]

B73_arrays_pos$group<- paste(B73_arrays_pos$gene_up, B73_arrays_pos$gene_down, B73_arrays_pos$chr, sep="_")
B73_arrays_pos %>% group_by(gene_up) %>% summarise(n=n()) #157
B73_arrays_pos %>% group_by(gene_down) %>% summarise(n=n()) #156
groups<- B73_arrays_pos %>% group_by(gene_up,gene_down, chr) %>% summarise(n=n()) #181
B73_arrays_pos %>% group_by(group) %>% summarise(n=n()) #181

B73_arrays_pos[B73_arrays_pos$gene_up %in% ".", ]$gene_up<- paste(B73_arrays_pos[B73_arrays_pos$gene_up %in% ".", ]$gene_down, "edit", sep="_")
B73_arrays_pos[B73_arrays_pos$gene_down %in% ".", ]$gene_down<- paste(B73_arrays_pos[B73_arrays_pos$gene_down %in% ".", ]$gene_up, "edit", sep="_")

B73_arrays_pos$unique_name<- paste(B73_arrays_pos$dat, B73_arrays_pos$array_start, B73_arrays_pos$chr, sep="_" )
B73_arrays_pos_chr_edges<- as.data.frame(matrix(nrow=0, ncol=2))
for(i in 1:nrow(B73_arrays_pos)){
  sub<- B73_arrays_pos[B73_arrays_pos$gene_up %in% B73_arrays_pos$gene_up[i] | B73_arrays_pos$gene_down %in% B73_arrays_pos$gene_down[i],]
   new_edges<- as.data.frame(cbind(B73_arrays_pos$unique_name[i], sub$unique_name))
   B73_arrays_pos_chr_edges<- rbind(B73_arrays_pos_chr_edges, new_edges)
  if(i%%1000 ==0){
    print(i)
}}


g2<- graph_from_edgelist(as.matrix(B73_arrays_pos_chr_edges))
V(g2)$membership<- components(g2)$membership
membership2<- as.data.frame(components(g2)$membership)
membership2$unique_name<- row.names(membership2)
B73_arrays_pos_member<- merge(B73_arrays_pos, membership2, by="unique_name")
B73_arrays_pos_member$label<- B73_arrays_pos_member$"components(g2)$membership" ##150 groups
chr_B73_arrays_pos_member<- B73_arrays_pos_member %>% group_by(label) %>% summarise(n=length(unique(chr))) #151

B73_arrays_pos_member$rep<- B73_arrays_pos_member$element

grouped_B73_arrays_pos_member<-B73_arrays_pos_member %>% group_by(gene_up, gene_down, dat, chr) %>% 
  summarise(array_start=min(array_start), array_end=max(array_end),
            rep=unique(rep), gene_up_end=min(gene_up_end), 
            gene_down_start=min(gene_down_start), ref_end=min(ref_end),
            ref_start=min(ref_start), distinct=unique(label )) %>% as.data.frame()
write.table(B73_arrays_pos_member, "B73_arrays_pos_member", col.names=T, row.names=F, quote=F)
write.table(grouped_B73_arrays_pos_member, "grouped_B73_arrays_pos_member", col.names=T, row.names=F, quote=F)


#########################################################################################################################
#plotting
####plotting
#fill in missing information
grouped_arrays_pos_member$ref_start[is.na(grouped_arrays_pos_member$ref_start)]<- grouped_arrays_pos_member[is.na(grouped_arrays_pos_member$ref_start),]$ref_end 
grouped_arrays_pos_member$ref_end[is.na(grouped_arrays_pos_member$ref_end)]<- grouped_arrays_pos_member[is.na(grouped_arrays_pos_member$ref_end),]$ref_start
grouped_arrays_pos_member$dat2<- str_split(grouped_arrays_pos_member$dat, pattern="_", simplify=T)[,1]
grouped_arrays_pos_member$dat2<- str_split(grouped_arrays_pos_member$dat2, pattern="[.]", simplify=T)[,1]
grouped_arrays_pos_member[grouped_arrays_pos_member$dat2 %in% "Zea-mays-ssp-mexicana-TIL25", ]$dat2 <- "TIL25"
grouped_arrays_pos_member$dat2<- as.factor(grouped_arrays_pos_member$dat2)
grouped_arrays_pos_member$dat2<-factor(grouped_arrays_pos_member$dat2, levels=c("Mo17","AB10", "B73", "CG108", 
                                                                                "CG119", "CG44", "CML442", "K64",
                                                                                "Tx777", "Tx779", "TIL01", 
                                                                                "TIL11", "TIL25"))
grouped_arrays_pos_member$chr<-factor(grouped_arrays_pos_member$chr, levels=c("chr1","chr2","chr3","chr4",
                                                                              "chr5","chr6","chr7","chr8",
                                                                              "chr9","chr10"))
p<- ggplot(data=grouped_arrays_pos_member, aes(x = ref_start, y =dat2, label=distinct)) +
  geom_point( aes(colour = rep, size=array_end-array_start))+
  geom_line(aes(colour = rep, group=distinct))+
  #geom_text()+
  theme_classic()+facet_wrap(~chr, strip.position = "right", nrow=5)+
  scale_color_manual(values = c( "darkcyan", "darkturquoise" ,"deeppink", "purple"))
ggsave("Mo17scaff_ArrayPositionPlot.png", plot = p, width = 12, height = 6, units = "in", dpi = 300)

###B73
grouped_B73_arrays_pos_member$ref_start[is.na(grouped_B73_arrays_pos_member$ref_start)]<- grouped_B73_arrays_pos_member[is.na(grouped_B73_arrays_pos_member$ref_start),]$ref_end 
grouped_B73_arrays_pos_member$ref_end[is.na(grouped_B73_arrays_pos_member$ref_end)]<- grouped_B73_arrays_pos_member[is.na(grouped_B73_arrays_pos_member$ref_end),]$ref_start
grouped_B73_arrays_pos_member$dat2<- str_split(grouped_B73_arrays_pos_member$dat, pattern="_", simplify=T)[,1]
grouped_B73_arrays_pos_member$dat2<- str_split(grouped_B73_arrays_pos_member$dat2, pattern="[.]", simplify=T)[,1]
grouped_B73_arrays_pos_member[grouped_B73_arrays_pos_member$dat2 %in% "Zea-mays-ssp-mexicana-TIL25", ]$dat2 <- "TIL25"
grouped_B73_arrays_pos_member$dat2<- as.factor(grouped_B73_arrays_pos_member$dat2)
grouped_B73_arrays_pos_member$dat2<-factor(grouped_B73_arrays_pos_member$dat2, levels=c("Mo17","AB10", "B73", "CG108", 
                                                                                "CG119", "CG44", "CML442", "K64",
                                                                                "Tx777", "Tx779", "TIL01", 
                                                                                "TIL11", "TIL25"))
grouped_B73_arrays_pos_member$chr<-factor(grouped_B73_arrays_pos_member$chr, levels=c("chr1","chr2","chr3","chr4",
                                                                              "chr5","chr6","chr7","chr8",
                                                                              "chr9","chr10"))
p2<- ggplot(data=grouped_B73_arrays_pos_member, aes(x = ref_start, y =dat2, label=distinct)) +
  geom_point( aes(colour = rep, size=array_end-array_start))+
  geom_line(aes(colour = rep, group=distinct))+
  #geom_text()+
  theme_classic()+facet_wrap(~chr, strip.position = "right", nrow=5)+
  scale_color_manual(values = c( "darkcyan", "darkturquoise" ,"deeppink", "purple"))
ggsave("B73scaff_ArrayPositionPlot.png", plot = p2, width = 12, height = 6, units = "in", dpi = 300)

####
#n positions
group_freq<- grouped_arrays_pos_member %>% group_by(distinct, rep) %>% 
  summarise(n=n_distinct(dat)) %>% group_by(rep, n) %>% summarise(n_pos=n())
afs<-ggplot() +geom_bar( data=group_freq, aes(x=n, y=n_pos,fill=rep),stat = "identity" ) +theme_classic()+
  #geom_point( data=group_freq, aes(x=n, y=n_pos*n)) +theme_classic()+
  facet_wrap(~rep, strip.position = "right", nrow=4)+
  scale_fill_manual(values = c( "darkcyan", "darkturquoise" ,"deeppink", "purple"))+
  scale_y_continuous( name = "Number of Positions")
ggsave("Mo17scaff_AFS.png", plot = afs, width = 12, height = 3, units = "in", dpi = 300)

ar_counts<- grouped_arrays_pos_member %>% group_by(dat2, rep) %>% 
  summarise(n=n()) 
ArCount<-ggplot() +geom_bar( data=ar_counts, aes(x=dat2, y=n,fill=rep),stat = "identity" ) +theme_classic()+
  facet_wrap(~rep, strip.position = "right", nrow=5)+
  scale_fill_manual(values = c( "darkcyan", "darkturquoise" ,"deeppink", "purple"))+
  scale_y_continuous( name = "Number of Positions")
ggsave("Mo17scaff_ArrayCount.png", plot = ArCount, width = 12, height = 3, units = "in", dpi = 300)



group_freq_B73<- grouped_B73_arrays_pos_member %>% group_by(distinct, rep) %>% 
  summarise(n=n_distinct(dat)) %>% group_by(rep, n) %>% summarise(n_pos=n())
afs2<-ggplot() +geom_bar( data=group_freq_B73, aes(x=n, y=n_pos,fill=rep),stat = "identity" ) +theme_classic()+
  #geom_point( data=group_freq_B73, aes(x=n, y=n_pos*n)) +theme_classic()+
  facet_wrap(~rep, strip.position = "right", nrow=4)+
  scale_fill_manual(values = c( "darkcyan", "darkturquoise" ,"deeppink", "purple"))+
  scale_y_continuous( name = "Number of Positions")
ggsave("B73scaff_AFS.png", plot = afs2, width = 12, height = 3, units = "in", dpi = 300)
ar_counts_B73<-grouped_B73_arrays_pos_member %>% group_by(dat2, rep) %>% 
  summarise(n=n()) 
ArCount2<-ggplot() +geom_bar( data=ar_counts_B73, aes(x=dat2, y=n,fill=rep),stat = "identity" ) +theme_classic()+
  facet_wrap(~rep, strip.position = "right", nrow=5)+
  scale_fill_manual(values = c( "darkcyan", "darkturquoise" ,"deeppink", "purple"))+
  scale_y_continuous( name = "Number of Positions")
ggsave("B73scaff_ArrayCount.png", plot = ArCount2, width = 12, height = 3, units = "in", dpi = 300)






