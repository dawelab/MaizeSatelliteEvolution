#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)
bins_fil= args[1]
dir= args[2]

library('dplyr')
library('reshape2')
library('stringr')
library('igraph')
library('caret')

setwd(dir)
bins<- read.table(bins_fil)
`%!in%` = Negate(`%in%`)

data_files<- file.info(dir(dir, pattern=c( "sub_HOR"), full.names=FALSE)) ##all the comparison files
data_files2<- row.names(data_files[data_files$size !=0, ])
data_files3<-grep(glob2rx("sub_HOR*.blat"), data_files2, value = TRUE)

## First section of code reads in the pairwise monomer comparison scores for th_HOR_bins.bede bin, turns them into a similarity matrix, and creates networks for similarity thresholds, as extracted from file name
# Summarizing information from those networks are then saved in data frames

All <- lapply(data_files3,function(i){
	x<- read.table(i)
	y<- as.data.frame(cbind( x$V1, as.character(x$V10), x$V11, as.character(x$V14), x$V15))
})

All2<- list()
All_networks<- list()
network_summary<- as.data.frame(matrix(nrow=length(data_files3), ncol=16))

for( i in 1:length(All) ){
    colnames(All[[i]])<- c("match", "nam1", "len1", "nam2", "len2")
    All[[i]]$match<- as.numeric(as.character(All[[i]]$match))
    All[[i]]$len1<- as.numeric(as.character(All[[i]]$len1))
    All[[i]]$len2<- as.numeric(as.character(All[[i]]$len2))
    All[[i]]$jacc<- round(All[[i]]$match/( All[[i]]$len1 + All[[i]]$len2 - All[[i]]$match ), 6)
    
    nam_val<- str_split(data_files3[i], pattern="_|[.]", simplify=T)
    nam_chr<-  nam_val[,3]
    ar_start<- as.numeric(nam_val[,7])
    rep<- nam_val[,8]
    directory<-paste(nam_val[,4],nam_val[,5], nam_val[,6], sep="_" )
    info<- bins[bins$V8 %in% nam_chr & bins$V9 %in% directory & bins$V2 %in% ar_start ,]
    c<- info$V5
    

    dat2<- aggregate(jacc ~ nam1 + nam2, All[[i]], max)
    rev_dat2<-dat2
    colnames(rev_dat2)<- c("nam2",  "nam1", "jacc")
    dat2_d<- rbind( dat2, rev_dat2)
    dat2_mat<- acast(dat2_d, nam1~nam2, value.var="jacc", fun.aggregate = max, fill=0.0000)
    
    for(j in 1:nrow( dat2_mat)){
      dat2_mat[j,j]<-1
    }
    All2[[i]]<- dat2_mat
    
    if( nrow(dat2_mat) > 1){
      dup<- findCorrelation(dat2_mat,cutoff = .9999999, names=TRUE) ##concat identical 
      keep<- colnames(dat2_mat)[colnames(dat2_mat) %!in% dup]
      nams<- as.data.frame(keep)

    }else{
      nams<- as.data.frame(dat2_mat)
      keep<- colnames(nams)
    }
    
    nams$count<- NA
    
    for( j in 1:nrow(nams)){
      row<-as.data.frame(dat2_mat[row.names(dat2_mat) %in% nams[j,1],])
      colnames(row)<- c("val")
      row$nam<- row.names(row)
      sub<-row[row$val==1,]
      nams$count[j]<-nrow(sub)
      max_ident<- max( nams$count)
    }
    
    
    count_nam<- nams$count
    colnames( nams)<- c("V1", "V3")
    

	if( length(keep) > 1) {
    	dat2_mat_sub<- dat2_mat[keep, keep]
    	dat2_mat_sub[ dat2_mat_sub < c ] <- 0
    	dat2_mat_sub[dat2_mat_sub >= c ] <- 1
        
    	network2 <- graph_from_adjacency_matrix(dat2_mat_sub,  mode="undirected", diag=F, weighted = T)
    	group_counts<-as.data.frame(components(network2)$membership) 
   	 	colnames( group_counts)<- c("V2")
    	group_counts$V1<- row.names(group_counts)
    	group_counts_or_names <- merge( group_counts, nams, by="V1") #V1 is the monomer, V2 is group, V3 is number of monomers
    	count2<- group_counts_or_names %>% group_by(V2) %>% summarise(V3=sum(V3))
    	count2_or<- count2[order(-count2$V3),]
    	V(network2)$label <- count_nam
    
    	network_summary[i,1]<- paste(nam_chr, ar_start, sep="_")
    	network_summary[i,2]<- nam_chr
    	network_summary[i,3]<- ar_start
    	network_summary[i,4]<- info$V3
    	network_summary[i,5]<- c
    	network_summary[i,6]<- length(unique(All[[i]]$nam1)) #number of monomers
    	network_summary[i,7]<- (length(dup)+ nrow( nams[nams$V3 > 1,]))/network_summary[i,6] #number of duplicate monomers, prop of total
    	network_summary[i,8]<- length(unique(components(network2)$membership))/network_summary[i,6] #number of separated clusters, prop of total
    	network_summary[i,9]<- count2_or$V3[1]/network_summary[i,6] #bigggest cluster, prop of total
    	network_summary[i,10]<- count2_or$V3[2]/network_summary[i,6]
    	network_summary[i,11]<- count2_or$V3[3]/network_summary[i,6]
    	network_summary[i,12]<- modularity(network2, components(network2)$membership) #modularity 
    	network_summary[i,13]<- mean_distance(network2, directed=F) #mean distance between nodes
    	network_summary[i,14]<- mean(dat2_d$jacc)
    	network_summary[i,15]<- max_ident/network_summary[i,4]
    	network_summary[i,16]<- rep
     }else{
    	network_summary[i,1]<- paste(nam_chr, ar_start, sep="_")
    	network_summary[i,2]<- nam_chr
    	network_summary[i,3]<- ar_start
    	network_summary[i,4]<- info$V3
    	network_summary[i,5]<- c
        network_summary[i,6]<- length(unique(All[[i]]$nam1)) #number of monomers
        network_summary[i,7]<-(length(dup)+ nrow( nams[nams$V3 > 1,]))/network_summary[i,6] #number of duplicate monomers
        network_summary[i,8]<- 1 #number of separated clusters
        network_summary[i,9]<- 1
        network_summary[i,10]<- 0
        network_summary[i,11]<- 0
        network_summary[i,12]<- 0 #modularity 
        network_summary[i,13]<- 0 #mean distance between nodes
        network_summary[i,14]<- 1
        network_summary[i,15]<- 1
        network_summary[i,16]<- rep
        }
        
    nams_group<- as.data.frame(matrix(nrow=0, ncol=2))
      if( length(keep) > 1) {
      for( k in 1:nrow(nams)){
  		row<-as.data.frame(dat2_mat[as.character(nams[k,1]),])
  		colnames(row)<- c("val")
  		row$nam<- row.names(row)
  		sub<-row[row$val==1,]
  		sub$group<- as.character(nams[k,1])
  		sub$clust<- group_counts[nams$V1[k],]$V2
  		nams$count[k]<-nrow(sub)
  		nams_group<- rbind(nams_group, sub)
	}

	nams_group$bin<-data_files3[i]
	
    write.csv(nams_group, paste( data_files3[i], "_bins_nam_groups.csv", sep=""), quote=F)
    }
}
colnames(network_summary)<- c("bin","chr",  "start", "end", "clust_val", "num_monomer", "dup_monomers", "uncon_clust", "prop_clust1", "prop_clust2","prop_clust3", "modularity", "mean_dist", "avg_jacc", "max_contract", "rep")
write.csv( network_summary,"HOR_re-eval_clusters.csv", quote=F)

