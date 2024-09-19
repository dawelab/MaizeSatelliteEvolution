#####graph plot and cluster for data set with Mo17

library("igraph")
library("ggplot2")
library("dplyr")
library("stringr")
library(reshape2)

setwd("~/Desktop/new_groups")
path="~/Desktop/new_groups"

struct<- read.table("~/Desktop/all_structure_conserved_wMo17_v2_mod")
colnames(struct)<- c("line", "element", "chr", "ar_start","bin_start", "bin_end", "group", "threshold","mon_num", "structure")
struct$len<- struct$bin_end-struct$bin_start

struct_sum<- struct %>% group_by(line, element, chr, ar_start, group, structure) %>% 
  summarise(n=sum(len), mono_num=sum(mon_num))

 #Grouping structure to get the total # of monomers, and number of monomers in each class
struct_sum_group<- struct_sum[,1:5] %>% unique()
struct_sum_group_HOR<- struct_sum %>% group_by(line, element, chr, ar_start, group) %>% filter(structure == "HOR") %>% as.data.frame()
struct_sum_group_Order<- struct_sum %>% group_by(line, element, chr, ar_start, group) %>% filter(structure == "Order")%>% as.data.frame()
struct_sum_group_Disorder<- struct_sum %>% group_by(line, element, chr, ar_start, group) %>% filter(structure == "Disorder")%>% as.data.frame()
struct_sum_group1<- merge(struct_sum_group, struct_sum_group_HOR, c("line", "element", "chr", "ar_start", "group"), all.x=T)
struct_sum_group2<- merge(struct_sum_group1, struct_sum_group_Order, c("line", "element", "chr", "ar_start", "group"), all.x=T)
struct_sum_group3<- merge(struct_sum_group2, struct_sum_group_Disorder, c("line", "element", "chr", "ar_start", "group"), all.x=T)
struct_sum_group3_sel<- dplyr::select(struct_sum_group3, c("line", "element", "chr", "ar_start", 
                                                    "group", "n.x", "mono_num.x", "n.y", "mono_num.y","n", "mono_num" ))
colnames(struct_sum_group3_sel)<- c("line", "element", "chr", "ar_start","group", "HOR_n", "HOR_mono_num",
                                    "Order_n", "Order_mono_num", "Disorder_n", "Disorder_mono_num")
struct_sum_group3_sel[is.na(struct_sum_group3_sel)]<-0
struct_sum_group3_sel$line<- str_split(struct_sum_group3_sel$line, pattern="_", simplify = T)[,1]

#N gaps
N<- read.table("~/Desktop/all_100_gaps")
lines<- c("AB10_half_Mo17","B73_half_Mo17","CG108_half_Mo17","CG119_fix","CG44_half_Mo17", 
          "CML442_half_Mo17","K64_half_Mo17", "Mo17_half_Mo17", "Tx777_half_Mo17","Tx779_half_Mo17",
          "TIL01.3cell.HiFi_half_Mo17", "TIL11.2cell.HiFi_half_Mo17", "Zea-mays-ssp-mexicana-TIL25_4cell-hifi_half_Mo17")
N$line<- str_split(N$V1, pattern="_", simplify = T)[,1]
N$chr<- str_split(N$V2, pattern="_", simplify = T)[,1]
N_lin<- N[N$V1 %in% lines, ]

#Maize and teo lines
maize_lines<- c("AB10","B73","CG108","CG119","CG44", "CML442","K64", "Mo17", "Tx777","Tx779")
teosinte_lines<- c("TIL01.3cell.HiFi", "TIL11.2cell.HiFi", "Zea-mays-ssp-mexicana-TIL25_4cell-hifi")

data_files<-rownames(file.info(dir(path, pattern="_arrays_otherMASK.blast.sub_test2", full.names=TRUE)))
vals<- str_split(data_files, pattern="/|_", simplify=T)[,7] %>% as.data.frame()
comp_mono_vals_struct_all<- as.data.frame(matrix(nrow=0, ncol=16))
comp_all_combined<- as.data.frame(matrix(nrow=0, ncol=20))

for(v in 1:length(vals[,1])){
  i<-vals[v,1]

  #get the data for array comp-- make sure all the lines are good, no malformed lines
  comp_mono<- read.table(paste(i,"_arrays_otherMASK.blast.sub_test2", sep=""))
  
  comp_mono$line1<-  str_split( comp_mono$V1, pattern="_", simplify = T)[,2]
  comp_mono$line2<-  str_split( comp_mono$V1, pattern="_", simplify = T)[,6]
  comp_mono$chr1<-  str_split( comp_mono$V1, pattern="_", simplify = T)[,1]
  comp_mono$chr2<-  str_split( comp_mono$V1, pattern="_", simplify = T)[,5]
  comp_mono<- comp_mono[ comp_mono$line1 %in% lines2 &  comp_mono$line2 %in% lines2 & comp_mono$chr1 %in% comp_mono$chr2,]
  if(nrow(  comp_mono)>1){
  comp_mono$data1<- str_split(comp_mono$V1, pattern="_chr", simplify=T)[,1]
  comp_mono$data2<- paste("chr", str_split(comp_mono$V1, pattern="_chr", simplify=T)[,2], sep="")
  comp_mono$arr_coords1<- str_split(comp_mono$data1, pattern=":", simplify=T)[,2]
  comp_mono$arr_coords2<- str_split(comp_mono$data2, pattern=":", simplify=T)[,2]
  comp_mono$len1<- as.numeric(str_split(comp_mono$arr_coords1, pattern="-", simplify=T)[,2])-as.numeric(str_split(comp_mono$arr_coords1, pattern="-", simplify=T)[,1])
  comp_mono$len2<- as.numeric(str_split(comp_mono$arr_coords2, pattern="-", simplify=T)[,2])-as.numeric(str_split(comp_mono$arr_coords2, pattern="-", simplify=T)[,1])
  
  comp_mono<- comp_mono[comp_mono$len1 > 0 & comp_mono$len2 > 0 , ]
  
  comp_mono$V4<- as.numeric(comp_mono$V4)
  comp_mono<- comp_mono[!is.na(comp_mono$V4),]
  comp_mono$V5<- as.numeric(comp_mono$V5)
  comp_mono<- comp_mono[!is.na(comp_mono$V5),]

  #top hits
  comp_mono_top<- comp_mono %>%filter(startsWith(as.character(V1), 'chr')) %>% 
    dplyr::filter(grepl('_chr', V1)) %>% group_by(V1, data1, len1, data2, len2) %>% filter(!is.na(as.numeric(V2)))%>% 
    filter(!is.na(as.numeric(V3))) %>% 
    summarise(V4=sum(V4), V5=mean(V5)) %>% as.data.frame()
  comp_mono_vals1<- dplyr::select(comp_mono_top, c("data1", "len1") ) 
  colnames(comp_mono_vals1)<- c("data", "len")
  comp_mono_vals2<- dplyr::select(comp_mono_top, c("data2", "len2") ) 
  colnames(comp_mono_vals2)<- c("data", "len")
  comp_mono_vals<- rbind(comp_mono_vals1, comp_mono_vals2) %>% unique()
  comp_mono_vals$gene_up<- i
  comp_mono_vals$line<- str_split(comp_mono_vals$data, pattern="_", simplify=T)[,2]
  if(as.numeric(i) >=151 ){
    comp_mono_vals$line<-"AB10" 
    }
  comp_mono_vals$coords<- str_split(comp_mono_vals$data, pattern=":", simplify=T)[,2]
  comp_mono_vals$ar_start<- as.numeric(str_split(comp_mono_vals$coords, pattern="-", simplify=T)[,1])
  comp_mono_vals$ar_end<- as.numeric(str_split(comp_mono_vals$coords, pattern="-", simplify=T)[,2])
  comp_mono_vals$chr<- str_split(comp_mono_vals$data, pattern="_", simplify = T)[,1]
  comp_mono_vals$group<- comp_mono_vals$gene_up
  comp_mono_vals$N<- 0
  for(k in 1:nrow(comp_mono_vals)){
    comp_mono_vals$N[k]<- nrow(N_lin[N_lin$line %in% comp_mono_vals$line[k] & N_lin$V3 > comp_mono_vals$ar_start[k] &  N_lin$V3 < comp_mono_vals$ar_end[k]  & N_lin$chr %in%  comp_mono_vals$chr[k], ])
  }
  no_gap_lines<- comp_mono_vals[comp_mono_vals$N==0,]$line

    #grab monomer to consensus information
  if(file.exists(paste("all_nam_consensus_mono_",i,".fasta.fa.blat.sub", sep="")) ){
    arrays<- read.table(paste("all_nam_consensus_mono_",i,".fasta.fa.blat.sub", sep="") )
    arrays<- arrays %>% group_by( V1, V3) %>% slice(which.max(V2)) %>% as.data.frame()
    if(nrow(arrays)>5){ #if at least 5 monomers in a single array, keep going 
    
  struct_sum_sub<-struct_sum_group3_sel[struct_sum_group3_sel$group %in% i, ] #struct_sum[struct_sum$group %in% i,]  #struct_sumHOR[struct_sumHOR$group %in% i,]
  comp_mono_vals_struct<- merge(comp_mono_vals,  struct_sum_sub, 
                                by=c("line", "ar_start", "group","chr" ), all.x=T)
  comp_mono_vals_struct[is.na(comp_mono_vals_struct)] <- 0

  #adding majority repeat to array info
  for( j in 1:nrow(comp_mono_vals_struct)){
    if(comp_mono_vals_struct$element[j]==0){
    comp_mono_vals_struct$element[j]<- struct_sum[struct_sum$group %in% i & struct_sum$line %in% comp_mono_vals_struct$line[j] & struct_sum$ar_start %in% comp_mono_vals_struct$ar_start[j],]$element[1]
    }
  }

  #monomer structure info
  comp_mono_vals_struct$HOR_prop<- (comp_mono_vals_struct$HOR_mono_num)/(comp_mono_vals_struct$HOR_mono_num+comp_mono_vals_struct$Order_mono_num + comp_mono_vals_struct$Disorder_mono_num)
  comp_mono_vals_struct_s<- comp_mono_vals_struct
  comp_mono_vals_struct_all<- rbind( comp_mono_vals_struct_all, comp_mono_vals_struct)

  #bed file of non-monomer space, so I can calculate how much  of the array is monomer
  if(file.exists(paste("non_mono_",i,"_arrays.fasta.bed", sep="")) ){
    mono_len<- read.table(paste("non_mono_",i,"_arrays.fasta.bed", sep=""))
    mono_len$len<-mono_len$V3-mono_len$V2
    mono_len_sum<- mono_len %>% group_by(V1) %>% summarise(non_mono_len=sum(len))
    comp_mono_vals_struct_s_len<- merge(comp_mono_vals_struct_s, mono_len_sum, by.x="data", by.y="V1")
    comp_mono_vals_struct_s_len$mono_len<- comp_mono_vals_struct_s_len$len-comp_mono_vals_struct_s_len$non_mono_len
  }else{ 
    comp_mono_vals_struct_s_len$mono_len<- comp_mono_vals_struct_s_len$len
    comp_mono_vals_struct_s_len$non_mono_len<-0
  }

  #all grouped
  comp_mono_vals_struct_s_len_grouped<- comp_mono_vals_struct_s_len %>% group_by(line, ar_start, group) %>% 
    summarise(sum_mono_len=max(mono_len), sum_len=max(len), HOR_mono_num=sum(HOR_mono_num) ,  
              Disorder_mono_num=sum(Disorder_mono_num), Order_mono_num=sum(Order_mono_num)) %>%
    group_by(line) %>% summarise(sum_mono_len=sum(sum_mono_len), sum_len=sum(sum_len),  prop_HOR=sum(HOR_mono_num)/sum(HOR_mono_num+ Order_mono_num + Disorder_mono_num))

  
  comp_mono_top$line1<- str_split(comp_mono_top$data1, pattern="_|:", simplify=T)[,2]
  comp_mono_top$line2<- str_split(comp_mono_top$data2, pattern="_|:", simplify=T)[,2]
  if(as.numeric(i) >=151  ){
    comp_mono_top$line1<-"AB10"
    comp_mono_top$line2<-"AB10" 
  }
  comp_mono_top_sub_sum<- comp_mono_top %>% group_by(line1, line2)%>% summarise(V4=sum(V4))

  comp_mono_top2<- merge( comp_mono_top_sub_sum,comp_mono_vals_struct_s_len_grouped, by.x="line1", by.y="line" )
   
  comp_mono_top2$jacc<- comp_mono_top2$V4/comp_mono_top2$sum_mono_len 

  comp_mono_top_sub<- dplyr::select(comp_mono_top2, c("line1", "line2","jacc"))
  comp_mono_top_sub[comp_mono_top_sub$jacc>1,]$jacc<- 1

  mat<- acast(comp_mono_top_sub, line1 ~ line2, value.var = 'jacc')
  
  # Get the lower triangular indices
  lower_tri_indices <- lower.tri(mat, diag = FALSE)
  # Create an empty matrix of the same dimensions
  lower_tri_max <- matrix(NA, nrow(mat), ncol(mat))
  lower_tri_min <- matrix(NA, nrow(mat), ncol(mat))
  # Populate the lower triangular part using vectorized pmax
  lower_tri_max[lower_tri_indices] <- pmax(mat[lower_tri_indices], t(mat)[lower_tri_indices])
  lower_tri_min[lower_tri_indices] <- pmin(mat[lower_tri_indices], t(mat)[lower_tri_indices])
  # Print the result
  mean(lower_tri_max, na.rm=T)
  mean(lower_tri_min, na.rm=T)
  diag(mat)=NA
  
  if(length(comp_mono_top_sub$line1[comp_mono_top_sub$line1 %in% maize_lines])>1){
    matMaize<- acast(comp_mono_top_sub[comp_mono_top_sub$line1 %in% maize_lines & comp_mono_top_sub$line2 %in% maize_lines, ], line1 ~ line2, value.var = 'jacc')
    diag( matMaize)=NA
    matMaize_mean<-  mean( matMaize, na.rm=T)
    }else{
      matMaize_mean<- NA
    }
  
  if(length(comp_mono_top_sub$line1[comp_mono_top_sub$line1 %in% teosinte_lines])>1){
    matTeo<- acast(comp_mono_top_sub[comp_mono_top_sub$line1 %in% teosinte_lines & comp_mono_top_sub$line2 %in% teosinte_lines, ], line1 ~ line2, value.var = 'jacc')
    diag( matTeo)=NA
    matTeo_mean<-  mean( matTeo, na.rm=T)
  }else{
    matTeo_mean<- NA
  }
  
  if(length(comp_mono_top_sub$line1[comp_mono_top_sub$line1 %in% no_gap_lines])>1){
    matNoGap<- acast(comp_mono_top_sub[comp_mono_top_sub$line1 %in% no_gap_lines & comp_mono_top_sub$line2 %in% no_gap_lines, ], line1 ~ line2, value.var = 'jacc')
    diag( matNoGap)=NA
    matNoGap_mean<-  mean( matNoGap, na.rm=T)
  }else{
    matNoGap_mean<- NA
  }
  
  if( dim(mat)[1] == 1){
    mat<- acast(comp_mono_top_sub, line1 ~ line2, value.var = 'jacc')
    lower_tri_min<- mat
    lower_tri_max<- mat
    
  }

  arrays$chr<- str_split(arrays$V1, pattern="_", simplify=T)[,1]
  arrays$line<- str_split(arrays$V1, pattern="_|:", simplify=T)[,2]
  
  #aary information
  comp_mono_top_sub_summary<-as.data.frame(matrix(nrow=1, ncol=0))
  comp_mono_top_sub_summary$gene_up<- i
  comp_mono_top_sub_summary$element<- comp_mono_vals_struct_s_len$element[1]
  comp_mono_top_sub_summary$n_arrays<- length(unique( arrays$line) )#length(V(comp_top_sub75_m_g))
  comp_mono_top_sub_summary$n_arraysNoGap<-length(no_gap_lines)
  comp_mono_top_sub_summary$jacc_mono_av<- mean(mat, na.rm=T)
  comp_mono_top_sub_summary$jacc_mono_med<- median(mat, na.rm=T)
  comp_mono_top_sub_summary$jacc_mono_max<- mean(lower_tri_max, na.rm=T)
  comp_mono_top_sub_summary$jacc_mono_maxMean<- max(lower_tri_max, na.rm=T)
  comp_mono_top_sub_summary$jacc_mono_min<- mean(lower_tri_min, na.rm=T)
  comp_mono_top_sub_summary$jacc_mono_minMean<- min(lower_tri_min, na.rm=T)
  comp_mono_top_sub_summary$matMaize_mean<- matMaize_mean
  comp_mono_top_sub_summary$matTeo_mean<- matTeo_mean
  comp_mono_top_sub_summary$matNoGap_mean<- matNoGap_mean
  comp_mono_top_sub_summary$min_len<- min(comp_mono_top2$sum_len)
  comp_mono_top_sub_summary$max_len<- max(comp_mono_top2$sum_len)
  comp_mono_top_sub_summary$max_monolen<- max(comp_mono_top2$sum_mono_len)
  comp_mono_top_sub_summary$max_HOR<- max(comp_mono_vals_struct$HOR_prop, na.rm=T)
  comp_mono_top_sub_summary$av_HOR<- mean(comp_mono_vals_struct$HOR_prop, na.rm=T)
  comp_mono_top_sub_summary$av_HOR_NoGap<- mean(comp_mono_vals_struct[comp_mono_vals_struct$line %in% no_gap_lines , ]$HOR_prop, na.rm=T)
  comp_all_combined<- rbind(comp_all_combined, comp_mono_top_sub_summary)

    }
    array_start <- str_split(arrays$V1, pattern=":", simplify=T)[,2]
    mono_start <- str_split(arrays$V1, pattern=":", simplify=T)[,3]
    arrays$array_start<- str_split(array_start, pattern = "-", simplify=T)[,1] %>% as.numeric()
    arrays$mono_start<- as.numeric(str_split(mono_start, pattern="-", simplify=T)[,1])
    arrays$start<- arrays$array_start + arrays$mono_start
    arrays$structure<- NA
    arrays$threshold<- NA
    
    struct$line<- str_split(struct$line, pattern="_", simplify = T)[,1]
    sub<-struct[struct$group %in% i,]


    #now, plotting monomer information
    for(j in 1:nrow(sub)){
      if(nrow(arrays[arrays$line %in% sub$line[j] & arrays$chr %in% sub$chr[j] & arrays$start>=sub$bin_start[j] & arrays$start<=sub$bin_end[j], ])){
        arrays[arrays$line %in% sub$line[j] & arrays$chr %in% sub$chr[j] & arrays$start>=sub$bin_start[j] & arrays$start<=sub$bin_end[j], ]$structure<- sub$structure[j]
        arrays[arrays$line %in% sub$line[j] & arrays$chr %in% sub$chr[j] & arrays$start>=sub$bin_start[j] & arrays$start<=sub$bin_end[j], ]$threshold<- sub$threshold[j]
      }
    }
    if( nrow( arrays[is.na(arrays$structure),]) >0) {arrays[is.na(arrays$structure),]$structure<- "Disorder"}
    ggplot()+geom_point(data=arrays, aes(x=mono_start, y=V2, shape=V3, color=structure), size=1)+
      facet_wrap(~line, nrow=13, scales="free_x", strip.position = "right") +
      theme_classic() +ggtitle(i)+ylab("Jacc to consensus")+xlab("Array Position")
    ggsave(paste(i,"mono","png",sep="."), width = 10, height = 15)
    
    ar<- str_split(arrays$V1, pattern=":", simplify=T)
    arrays$array<- paste(ar[,1], ar[,2], sep=":")
    arrays_or<- as.data.frame(matrix(nrow=0, ncol=ncol(arrays)+1))
    for( lin in unique(arrays$array)){
      sub<- arrays[arrays$array %in% lin, ]
      sub_or<- sub[order(sub$mono_start),]
      sub_or$num<-1:nrow(sub_or)
      arrays_or<- rbind(arrays_or, sub_or)
    }
    
    ggplot()+geom_point(data=arrays_or, aes(x=num, y=V2, shape=V3, color=structure), size=1)+
      facet_wrap(~line, nrow=13, scales="free_x", strip.position = "right") +
      theme_classic() +ggtitle(i)+ylab("Jacc to consensus")+xlab("Monomer Order")
    ggsave(paste(i,"mono_order","png",sep="."), width = 10, height = 15)
    
  }
  }
  print(i)
  
}

#write.table( comp_all_combined, "comp_all_combined", quote=F)
