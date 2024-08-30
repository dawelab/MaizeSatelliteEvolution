#!/usr/bin/Rscript

##Rscript --vanilla network_sliding_smooth_reads.R <repeat> <directory>
## Example: Rscript --vanilla network_sliding_smooth_reads.R CentC /scratch/rdp22327/Dawe/mono/read_out
#module load R/4.3.0-foss-2020b

args = commandArgs(trailingOnly = TRUE)
rep<-args[1]
dir=args[2]

library('dplyr')
library('reshape2')
library('stringr')
library('igraph')
library('caret')
library('MASS')
library('rpart')
library('data.table')

setwd(dir)
fit2_test<- readRDS("lda_model.rds")
`%!in%` = Negate(`%in%`)

path =getwd()
data_files<- file.info(dir(path, pattern=".sub", full.names=TRUE))
data_files<- row.names(data_files[data_files$size !=0, ])
clust_val <- list(.90, .91, .92, .93, .94, .95, .96, .97, .98, .99)

All <- lapply(data_files,function(i){
    x<- read.table(i)
    y<- as.data.frame(cbind( x$V1, as.character(x$V10), x$V11, as.character(x$V14), x$V15))
  })

  All2<- list()
  All_networks<- list()
  network_summary_list<- list()

  for( i in 1:length(clust_val)){
    network_summary_list[[i]]<- as.data.frame(matrix(nrow=length(All), ncol = 11))
  }

for( i in 1:length(All) ){
        nam<- str_split( data_files[i], pattern="/", simplify=T)[,7]
    colnames(All[[i]])<- c("match", "nam1", "len1", "nam2", "len2")
    All[[i]]$match<- as.numeric(as.character(All[[i]]$match))
    All[[i]]$len1<- as.numeric(as.character(All[[i]]$len1))
    All[[i]]$len2<- as.numeric(as.character(All[[i]]$len2))
    All[[i]]$jacc<- round(All[[i]]$match/( All[[i]]$len1 + All[[i]]$len2 - All[[i]]$match ), 6)
    
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

    for( c in 1:length(clust_val)){
                dat2_mat_sub<- dat2_mat[keep, keep]
        dat2_mat_sub[ dat2_mat_sub < clust_val[c] ] <- 0
        dat2_mat_sub[dat2_mat_sub >= clust_val[c] ] <- 1

                network2 <- graph_from_adjacency_matrix(dat2_mat_sub,  mode="undirected", diag=F, weighted = T)
        group_counts<-as.data.frame(components(network2)$membership)
        colnames( group_counts)<- c("V2")
        group_counts$V1<- row.names(group_counts)
        group_counts_or_names <- merge( group_counts, nams, by="V1") #V1 is the monomer, V2 is group, V3 is number of monomers
        count2<- group_counts_or_names %>% group_by(V2) %>% summarise(V3=sum(V3))
        count2_or<- count2[order(-count2$V3),]
        V(network2)$label <- count_nam
        network_summary_list[[c]][i,1]<- nam
        network_summary_list[[c]][i,2]<- clust_val[c]
                network_summary_list[[c]][i,3]<- length(unique(All[[i]]$nam1)) #number of monomers
        network_summary_list[[c]][i,4]<- (length(dup)+ nrow( nams[nams$V3 > 1,]))/network_summary_list[[c]][i,3] #prop of duplicate monomers
        network_summary_list[[c]][i,5]<- length(unique(components(network2)$membership))/network_summary_list[[c]][i,3]
        network_summary_list[[c]][i,6]<- count2_or$V3[1]/network_summary_list[[c]][i,3] 
        network_summary_list[[c]][i,7]<- count2_or$V3[2]/network_summary_list[[c]][i,3] #count in 2nd largest cluster
        network_summary_list[[c]][i,8]<- modularity(network2, components(network2)$membership) #modularity
        network_summary_list[[c]][i,9]<- mean_distance(network2, directed=F) #mean distance between nodes
        network_summary_list[[c]][i,10]<- mean(dat2_d$jacc)
        network_summary_list[[c]][i,11]<- max_ident/network_summary_list[[c]][i,3]
    }
  }
      

network_summary_or_list<- list()

  # Network summaries are then used to predict structure using LDA model
  for( i in 1:length(clust_val)){
    network_summary<- network_summary_list[[i]]
    colnames(network_summary)<- c("read", "clust_val","num_monomer", "dup_monomers", "uncon_clust", "prop_clust1", "prop_clust2", "modularity", "mean_dist", "avg_jacc", "max_contract")
        network_summary<-network_summary[network_summary$num_monomer >=5,]
        
    if(nrow(network_summary[is.na(network_summary$mean_dist),])>0){ network_summary[is.na(network_summary$mean_dist),]$mean_dist<- 10}
    if(nrow(network_summary[is.na(network_summary$prop_clust2),])>0){network_summary[is.na(network_summary$prop_clust2),]$prop_clust2<- 0}
    if(nrow(network_summary[is.na(network_summary$modularity),])>0){network_summary[is.na(network_summary$modularity),]$modularity<- 0}
    
    network_summary$lda_pred<-predict(fit2_test, newdata=network_summary, type = 'class')$class
    network_summary$lda_pred_pos<-apply(predict(fit2_test, newdata=network_summary, type = 'class')$posterior,1,max)
    
    network_summary_or_list[[i]]<- network_summary
  }
  
 combined_Multi_bins<- as.data.frame(matrix(nrow=0, ncol=ncol(network_summary_or_list[[1]] )))
 for( i in 1:nrow(network_summary_or_list[[1]])){
    comp<-rbindlist(lapply(network_summary_or_list, function(x)(x[i,]))) ##pull all the same bins from the different clustering files to compare
    comp<- comp[comp$lda_pred_pos >= .80,]
    if( "Exp" %in% comp$lda_pred ){ comp2<- comp[comp$lda_pred  %in% "Exp",] ##sub down to lines with the highest-order classification
    }else if ( "HOR" %in% comp$lda_pred){ comp2<- comp[comp$lda_pred  %in% "HOR",] 
    }else if ( "Order" %in% comp$lda_pred){ comp2<- comp[comp$lda_pred  %in% "Order",] 
    } else { comp2<- comp[comp$lda_pred  %in% "Disorder",]}
    
    best<- as.data.frame(comp2[which.max(comp2$lda_pred_pos),] )##pull line with higest posterior probability
    combined_Multi_bins<- rbind(combined_Multi_bins,best)
    if(nrow(best) > 1 ){ print(paste("Line is issue: ", i, places$V1[p]))}
  }  

write.csv(combined_Multi_bins, paste(rep,"_reads_class.csv"), quote=F, row.names=F)
