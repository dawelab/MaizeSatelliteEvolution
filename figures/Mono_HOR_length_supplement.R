#mono_coords
mono_coords<-read.table("~/Desktop/assem/mono_HOR_coords")
mono_coords1<-read.table("~/Desktop/assem/AB10_HALF_AB10_sub_mono_info_lens.bed")
mono_coords<- rbind(mono_coords, mono_coords1)
bins_interest<- unique(HOR_all_local$bin)
HOR_all_local_coords<- as.data.frame(matrix(nrow=0, ncol=8))
for( b in 41685:length(bins_interest)){
  HOR_all_local_sub<- HOR_all_local[HOR_all_local$bin %in% bins_interest[b], ]
  mono_coords_sub<- mono_coords[mono_coords$V1 %in% bins_interest[b] , ] %>% unique() %>% arrange(V3)
  if(nrow(mono_coords_sub)>0){
  mono_coords_sub$num<- 1:nrow(mono_coords_sub)
  HOR_all_local_sub$Start_coord<-NA
  HOR_all_local_sub$End_coord<- NA
  
  for(i in 1:nrow(HOR_all_local_sub)){
    HOR_all_local_sub$Start_coord[i]<- mono_coords_sub[mono_coords_sub$num %in% HOR_all_local_sub$Start[i],]$V3
    HOR_all_local_sub$End_coord[i]<- mono_coords_sub[mono_coords_sub$num %in% HOR_all_local_sub$End[i],]$V4
  }
  
HOR_all_local_sub$len<- HOR_all_local_sub$End_coord -  HOR_all_local_sub$Start_coord
HOR_all_local_coords<- rbind(HOR_all_local_coords, HOR_all_local_sub)
}
}

write.table(HOR_all_local_coords, file="HOR_all_local_coords", quote=F, row.names = F)
HOR_all_local_coord_rep_len_mean<- HOR_all_local_coords %>% group_by(X, bin, Pattern) %>% summarise(len_mean=mean(len), n=n()) %>% filter(n>=3)
mono_coords_uni<- mono_coords %>% group_by(V1, V5) %>% summarise(n=n())
HOR_all_local_coord_rep<- merge(HOR_all_local_coord_rep_len_mean, mono_coords_uni, by.x="bin", by.y="V1")
write.table(HOR_all_local_coord_rep, file="HOR_all_local_coord_rep", quote=F, row.names = F)

HOR_all_local_coord_rep_mean<- HOR_all_local_coord_rep %>% 
  group_by(V5) %>%summarise(mean_count=mean(len_mean))

HOR_all_local_coord_rep$n_char<- nchar( HOR_all_local_coord_rep$Pattern)

ggplot()+geom_histogram( data=HOR_all_local_coord_rep, aes(x=len_mean, fill=V5), bins=50)+
  geom_vline(data= HOR_all_local_coord_rep_mean, aes(xintercept = mean_count), linetype="dashed") + 
  facet_wrap(~V5, nrow=4, scales = "free_y", strip.position = "right") +theme_classic()+
  scale_fill_manual(values = c( "#88CCEE", "#6699CC" ,"#882255", "#CC6677"))+xlim(0,3000)

ggplot()+geom_density( data=HOR_all_local_coord_rep, aes(x=len_mean, fill=V5))+
  geom_vline(data= HOR_all_local_coord_rep_mean, aes(xintercept = mean_count), linetype="dashed") + 
  facet_wrap(~V5, nrow=4, scales = "free_y", strip.position = "right") +theme_classic()+
  scale_fill_manual(values = c( "#88CCEE", "#6699CC" ,"#882255", "#CC6677"))+xlim(0,3000)

ggplot()+geom_point( data=HOR_all_local_coord_rep, aes(x=len_mean, y=n_char,  colour=V5, size=n.x))+
  #geom_vline(data= HOR_all_local_coord_rep_mean, aes(xintercept = mean_count), linetype="dashed") + 
  #facet_wrap(~V5, nrow=4, scales = "free_y", strip.position = "right") +theme_classic()+
  scale_colour_manual(values = c( "#88CCEE", "#6699CC" ,"#882255", "#CC6677"))+xlim(0,3000)


Cent4_len<- read.table("mono_Cent4_ragtag.scaffold.fasta_arrays.fasta_lens")
Cent4_len$rep<- "Cent4"
CentC_len<- read.table("mono_CentC_ragtag.scaffold.fasta_arrays.fasta_lens")
CentC_len$rep<- "CentC"
TR1_len<- read.table("mono_TR1_ragtag.scaffold.fasta_arrays.fasta_lens")
TR1_len$rep<- "TR1"
knob180_len<- read.table("mono_knob180_ragtag.scaffold.fasta_arrays.fasta_lens")
knob180_len$rep<- "knob180"

mono_lens<- rbind( Cent4_len, CentC_len, TR1_len, knob180_len)

mono_lens_mean<- mono_lens %>% 
  group_by(rep) %>%summarise(mean_count=mean(V1))

ggplot()+geom_histogram( data=mono_lens, aes(x=V1, fill=rep))+
  geom_vline(data= mono_lens_mean, aes(xintercept = mean_count), linetype="dashed") + 
  facet_wrap(~rep, nrow=4, scales = "free", strip.position = "right") +theme_classic()+
  scale_fill_manual(values = c( "#88CCEE", "#6699CC" ,"#882255", "#CC6677"))

