#1 Remove RNA from survival results
#library(clusterProfiler)
cancer_name <- read.csv('/path/cancer_name.csv',header = F)
cancer_name <- as.character(cancer_name[,1])
gene_protein <- read.csv('/path/gene_protein_ID.txt',header = T)
gene_protein <- gene_protein[gene_protein[,3]!='',]
surv_data_clean <- list()
surv_data <- list()
for (i in 1:15) {
  surv_data[[i]] <- list()
  names(surv_data)[i] <- cancer_name[i]
  surv_data[[i]][[1]] <- read.csv(paste('/path/',cancer_name[i],'/surv_t.csv',sep = ''),row.names = 1)
  surv_data[[i]][[2]] <- read.csv(paste('/path/',cancer_name[i],'/surv_up.csv',sep = ''),row.names = 1)
  surv_data[[i]][[3]] <- read.csv(paste('/path/',cancer_name[i],'/surv_down.csv',sep = ''),row.names = 1)
  names(surv_data[[i]]) <- c('surv_t','surv_up','surv_down')
  surv_data_clean[[i]] <- list()
  for (j in 1:3) {
    surv_data_clean[[i]][[j]] <- surv_data[[i]][[j]]
    surv_data_clean[[i]][[j]] <- surv_data_clean[[i]][[j]][rownames(surv_data_clean[[i]][[j]])%in%unique(gene_protein$Gene.name),]
  }
  names(surv_data_clean[[i]]) <- c('surv_t','surv_up','surv_down')
}
#names(surv_data_clean) <- cancer_name

#Calculate significant genes

surv_data_V2 <- list()
surv_data_sig_V2 <- list()
for (j in 1:15) {
  a <- apply(surv_data_clean[[j]][[1]], 2, max)
  b <- surv_data_clean[[j]][[1]][,which.max(a[a!=0])]
  surv_data_sig_V2[[j]] <- surv_data_clean[[j]]
  surv_data_sig_V2[[j]][[1]] <- surv_data_sig_V2[[j]][[1]][which(surv_data_sig_V2[[j]][[1]][,which.max(a[a!=0])] >= 0.8),]
  surv_data_sig_V2[[j]][[2]] <- surv_data_sig_V2[[j]][[2]][which(surv_data_sig_V2[[j]][[2]][,which.max(a[a!=0])] >= min(surv_data_sig_V2[[j]][[1]][,which.max(a[a!=0])])),]
  surv_data_sig_V2[[j]][[3]] <- surv_data_sig_V2[[j]][[3]][which(surv_data_sig_V2[[j]][[3]][,which.max(a[a!=0])] >= min(surv_data_sig_V2[[j]][[1]][,which.max(a[a!=0])])),]
}


#Construct network 
surv_raw_data_V2 <- list()
result <- list()
surv_sig_inter_V2 <- list()
surv_sig_dgree_V2 <- list()

similar_value <- function(value){length(which(which(value==1)%in%which(a1==1)))/((length(which(a1==1))+length(which(value==1)))-2*length(which(which(a1==1)%in%which(value==1)))+1)}

length(cancer_name)
for (i in 1:15) {
  file_name <- dir(paste('/sibcb1/jihongbinlab1/caixinlei/data/database/project/survival_analysis/analysis/network_build/',cancer_name[i],'/raw_data',sep = ''))
  a <- which.min(which(apply(surv_data_sig_V2[[i]][[1]], 2, min)>0.1))
  surv_raw_data_V2[[i]] <- read.csv(paste('/path/',cancer_name[i],'/raw_data/',file_name[[a]],sep = ''),row.names = 1,header = T)
  surv_raw_data_V2[[i]] <- surv_raw_data_V2[[i]][rownames(surv_raw_data_V2[[i]])%in%unique(gene_protein$Gene.name),]
  surv_raw_data_V2[[i]] <- as.matrix(surv_raw_data_V2[[i]])
  surv_raw_data_V2[[i]][surv_raw_data_V2[[i]]>=0.05] <- 0
  surv_raw_data_V2[[i]][surv_raw_data_V2[[i]]!=0] <- 1
  a <- surv_raw_data_V2[[i]][rownames(surv_raw_data_V2[[i]])%in%rownames(surv_data_sig_V2[[i]][[1]]),]
  result[[i]] <- data.frame(row.names = rownames(surv_raw_data_V2[[i]]))
  for (x in 1:nrow(a)) {
    a1 <- a[x,]
    a2 <- surv_raw_data_V2[[i]]
    result[[i]][,x] <- apply(a2, 1, similar_value)
  }

  colnames(result[[i]]) <- rownames(a)
  result_matrix <- as.matrix(result[[i]])
  a <- seq(0.01,1,0.01)
  b <- data.frame()
  for (x in 1:100) {
    b[x,1] <- a[x]
    b[x,2] <- length(which(result_matrix>a[x]))-ncol(result_matrix)
  }
  surv_sig_inter_V2[[i]] <- data.frame()
  for (x in 1:ncol(result[[i]])) {
    a <- 1:length(result[[i]][,x])
    n <- nrow(surv_sig_inter_V2[[i]])
    surv_sig_inter_V2[[i]][(n+1):(n+length(a)),1] <- rep(colnames(result[[i]])[x],length(a))
    surv_sig_inter_V2[[i]][(n+1):(n+length(a)),2] <- result[[i]][a,x]
    surv_sig_inter_V2[[i]][(n+1):(n+length(a)),3] <- rownames(result[[i]])[a]
  }
  a <- character()
  b <- character()
  paste_all <- function(ob){paste0(ob,collapse = '')}
  a <- apply(surv_sig_inter_V2[[i]][,1:3], 1, paste_all)
  b <- apply(surv_sig_inter_V2[[i]][,3:1], 1, paste_all)
  c <- surv_sig_inter_V2[[i]][!b%in%a,]
  d <- surv_sig_inter_V2[[i]][b%in%a,]
  e <- unique(d$V1)
  f <- as.data.frame(t(data.frame(row.names = colnames(c))))
  colnames(f) <- colnames(c)
  for (x in 1:length(e)) {
    f <- rbind(f,d[d[,1]==e[x],])
    d <- d[d[,3]!=e[x],]
  }
  c <- rbind(c,f)
  surv_sig_inter_V2[[i]] <- c
  surv_sig_inter_V2[[i]][surv_sig_inter_V2[[i]][,1]==surv_sig_inter_V2[[i]][,3],2] <- 0
  surv_sig_inter_V2[[i]] <- surv_sig_inter_V2[[i]][surv_sig_inter_V2[[i]][,2]!=0,]
  surv_sig_inter_V2[[i]] <- surv_sig_inter_V2[[i]][order(surv_sig_inter_V2[[i]][,2],decreasing = T),]
  write.csv(surv_sig_inter_V2[[i]],paste('/path/result/',cancer_name[i],'_surv_sig_inter_all.csv',sep = ''),quote = F,row.names = F)
  
  surv_sig_inter_V2[[i]] <- surv_sig_inter_V2[[i]][surv_sig_inter_V2[[i]][,2]>=surv_sig_inter_V2[[i]][15000,2],]
  b <- surv_sig_inter_V2[[i]][,3:1]
  colnames(b) <- colnames(surv_sig_inter_V2[[i]])
  b <- rbind(surv_sig_inter_V2[[i]][,1:3],b)
  a <- unique(c(surv_sig_inter_V2[[i]][,1],surv_sig_inter_V2[[i]][,3]))
  surv_sig_dgree_V2[[i]] <- data.frame(a)
  sum_dgree <- function(ob){sum(b[,1]==ob)}
  surv_sig_dgree_V2[[i]][,2] <- apply(surv_sig_dgree_V2[[i]], 1, sum_dgree)
  surv_sig_dgree_V2[[i]] <- surv_sig_dgree_V2[[i]][order(surv_sig_dgree_V2[[i]][,2],decreasing = T),]
  
  write.csv(surv_sig_inter_V2[[i]],paste('/path/result/',cancer_name[i],'_surv_sig_inter.csv',sep = ''),quote = F,row.names = F)
  write.csv(surv_sig_dgree_V2[[i]],paste('/path/result/',cancer_name[i],'_surv_sig_dgree.csv',sep = ''),quote = F,row.names = F)
  write.csv(result[[i]],paste('/path/result/',cancer_name[i],'_surv_sig_result.csv',sep = ''),quote = F)
}
