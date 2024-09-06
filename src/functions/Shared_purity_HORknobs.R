##purity

library('stringr')

strings<-string_out
patts<- read.table("/scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/HOR_newpatt/HORknobs_ALL_SHARED_ALL_HOR.bed", fill=TRUE)
colnames(patts)<- c("pattern", "start", "end","line", "bin", "group", "level")
patts$start_coord<- as.numeric(str_split(patts$bin, pattern="_", simplify=T)[,2])
patts$chr<- str_split(patts$bin, pattern="_", simplify=T)[,1]

struct<- read.table("/scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/all_structure_conserved_wMo17_v2_mod")
struct_shared<- struct[struct$V1 %in% "CG108_half_Mo17",]
struct_shared$bin<- paste(struct_shared$V3, struct_shared$V6, sep="_")

patts_struct<- merge(patts, struct_shared, by.x=c("bin"), by.y=c("bin"))

bin_count<- patts_struct %>% group_by(pattern, group, level ) %>% summarise(n_bins=n_distinct(V7)) %>% filter(n_bins>1)
patts_multi<- merge(patts, bin_count, by=c("pattern","group"))

strings$purity<- 0
for(i in 1:nrow(strings)){
	new_patt<- strings[i,3]
	all_shared_HOR<- patts_multi[patts_multi$bin %in% strings$bin[i],] %>% select("line", "start", "end")
		if(nrow(all_shared_HOR)>0){
			all_shared_HOR2<- HOR_patt_check(all_shared_HOR, new_patt)
			strings$purity[i]<- purity_calc(all_shared_HOR2, new_patt)
		}
}

write.table(strings,"HORknobs_strings", quote=F, col.names=F, row.names=F )

