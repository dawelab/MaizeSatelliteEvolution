#!/usr/bin/Rscript

##Rscript --vanilla network_sliding_smooth_reads.R <repeat> <directory>
## Example: Rscript --vanilla network_sliding_smooth_reads.R CentC /scratch/rdp22327/Dawe/mono/read_out
#module load R/4.3.0-foss-2020b

args = commandArgs(trailingOnly = TRUE)
rep<-args[1]
dir=args[2]


library('dplyr')
library('reshape2')
library('tidyr')
library('stringr')
library('igraph')
library('caret')
library('MASS')
library('rpart')
library('data.table')

setwd(dir)
fit2_test<- readRDS("/scratch/rdp22327/Dawe/consensus/lda_model.rds")
`%!in%` = Negate(`%in%`)


path =getwd()
data_files<- file.info(dir(path, pattern=".sub", full.names=TRUE))
data_files<- row.names(data_files[data_files$size !=0, ])
clust_val <- list(.90, .91, .92, .93, .94, .95, .96, .97, .98, .99)
network_summary_list<- list()

#create a list of DF that will contain summary information for all similarity threshold networks
  for( i in 1:length(clust_val)){
    network_summary_list[[i]]<- as.data.frame(matrix(nrow=length(data_files), ncol = 11))
  }

  
  for(i in 1:length(data_files)){
  print(data_files[i])
    y<- as.data.frame(read.table(data_files[i]))
    nam<- str_split( data_files[i], pattern="/", simplify=T)[,9] #column for name
    colnames(y)<- c("match", "nam1", "len1", "nam2", "len2") #DF from should have match score, mono name 1, mono 1 len, mono name 2, mono 2 len
    y$match<- as.numeric(as.character(y$match))
    y$len1<- as.numeric(as.character(y$len1))
    y$len2<- as.numeric(as.character(y$len2))
    y$jacc<- round(y$match/( y$len1 + y$len2 - y$match ), 6)
    y_sub<- y[y$jacc>=.90, ] #filterpairwise similarities to at least the minimum threshold to cut down on data size
    dat2<- aggregate(jacc ~ nam1 + nam2, y_sub, max) # find max value for each pair
    rev_dat2<-dat2 # get opposite direction version, just to make sure the max val is symmetric in DF
    colnames(rev_dat2)<- c("nam2",  "nam1", "jacc")
    dat2_d<- rbind( dat2, rev_dat2)
    dat2_mat<- acast(dat2_d, nam1~nam2, value.var="jacc", fun.aggregate = max, fill=0.0000) # turn paurewuse list to a matrix!
  for(j in 1:nrow( dat2_mat)){
      dat2_mat[j,j]<-1
    } #for similarities of self-self, set as 1
   if( nrow(dat2_mat) > 1){
      dup<- findCorrelation(dat2_mat,cutoff = .9999999, names=TRUE) ##concat identical
      keep<- colnames(dat2_mat)[colnames(dat2_mat) %!in% dup] ##dedup identical
      nams<- as.data.frame(keep)
    }else{
      nams<- as.data.frame(dat2_mat)
      keep<- colnames(nams)
    }
     nams$count<- NA 

 for( j in 1:nrow(nams)){ ##count number of identical monomers for each deduped subtype
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
 	dat2_mat_sub<- dat2_mat  # now, finally, set edged => cut off as 1, remove is similarity value is less. Here, I make a copy of the data, just in case
        dat2_mat_sub[ dat2_mat_sub < clust_val[c] ] <- 0
        dat2_mat_sub[dat2_mat_sub >= clust_val[c] ] <- 1

       network2 <- graph_from_adjacency_matrix(dat2_mat_sub,  mode="undirected", diag=F, weighted = T) #make into a network
        group_counts<-as.data.frame(components(network2)$membership)
        colnames( group_counts)<- c("V2")
        group_counts$V1<- row.names(group_counts)
        group_counts_or_names <- merge( group_counts, nams, by="V1") #V1 is the monomer, V2 is group, V3 is number of monomers (from max_ident loop above)
        count2<- group_counts_or_names %>% group_by(V2) %>% summarise(V3=sum(V3))
        count2_or<- count2[order(-count2$V3),] 

        #finally, putting network summary information in the list of DFs
        network_summary_list[[c]][i,1]<- nam
        network_summary_list[[c]][i,2]<- clust_val[c]
	network_summary_list[[c]][i,3]<- length(unique(y$nam1)) #number of monomers
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
 	network_summary<-network_summary[network_summary$num_monomer >=5,] #minimum monomers in a region to get a class is 5

   #setting some default values to empty rows
    if(nrow(network_summary[is.na(network_summary$mean_dist),])>0){ network_summary[is.na(network_summary$mean_dist),]$mean_dist<- 10}  #setting high distance rather than NA if average distance empty
    if(nrow(network_summary[is.na(network_summary$prop_clust2),])>0){network_summary[is.na(network_summary$prop_clust2),]$prop_clust2<- 0} #setting zero rather than NA if second cluster empty
    if(nrow(network_summary[is.na(network_summary$modularity),])>0){network_summary[is.na(network_summary$modularity),]$modularity<- 0} #setting zero rather than NA if modularity empty
    
    network_summary$lda_pred<-predict(fit2_test, newdata=network_summary, type = 'class')$class #predict!
    network_summary$lda_pred_pos<-apply(predict(fit2_test, newdata=network_summary, type = 'class')$posterior,1,max) #get score on that prediction
    
    network_summary_or_list[[i]]<- network_summary #add value and class to original DF
  }

##now, combine all those values, find most likely class and optimal threshold
  combined_Multi_bins<- as.data.frame(matrix(nrow=0, ncol=ncol(network_summary_or_list[[1]] )))
 for( i in 1:nrow(network_summary_or_list[[1]])){
    comp<-rbindlist(lapply(network_summary_or_list, function(x)(x[i,]))) ##pull all the same bins from the different clustering files to compare
    comp<- comp[comp$lda_pred_pos >= .80,] #set minimal probability to consider the classification.
    if( "Exp" %in% comp$lda_pred ){ comp2<- comp[comp$lda_pred  %in% "Exp",] ##sub down to lines with the highest-order classification. Defunct class, rereverted to "Order" later on
    }else if ( "HOR" %in% comp$lda_pred){ comp2<- comp[comp$lda_pred  %in% "HOR",] 
    }else if ( "Order" %in% comp$lda_pred){ comp2<- comp[comp$lda_pred  %in% "Order",] 
    } else { comp2<- comp[comp$lda_pred  %in% "Disorder",]} ##if the region is not given a pred
    
    best<- as.data.frame(comp2[which.max(comp2$lda_pred_pos),] )##pull line with higest posterior probability
    combined_Multi_bins<- rbind(combined_Multi_bins,best)
    if(nrow(best) > 1 ){ print(paste("Line is issue: ", i, places$V1[p]))}
  }  

write.csv(combined_Multi_bins, paste(rep,"_reads_class.csv", sep=""), quote=F, row.names=F)
  

