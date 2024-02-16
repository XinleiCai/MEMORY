#Calculating significant genes
library(grid)
cancer_name <- read.csv('cancer_name.csv',header = F)
cancer_name <- as.character(cancer_name[,1])
surv_data_V2 <- list()
surv_data_sig_V2 <- list()
for (j in 1:length(cancer_name)) {
  # surv_data[[j]] <- list()
  # names(surv_data)[j] <- cancer_name[j]
  # surv_data[[j]][[1]] <- read.csv(paste('result/',cancer_name[j],'/data/surv_t.csv',sep = ''),row.names = 1)
  # surv_data[[j]][[2]] <- read.csv(paste('result/',cancer_name[j],'/data/surv_up.csv',sep = ''),row.names = 1)
  # surv_data[[j]][[3]] <- read.csv(paste('result/',cancer_name[j],'/data/surv_down.csv',sep = ''),row.names = 1)
  # names(surv_data[[j]]) <- c('surv_t','surv_up','surv_down')
  a <- apply(surv_data_clean[[j]][[1]], 2, max)
  b <- surv_data_clean[[j]][[1]][,which.max(a[a!=0])]
  surv_data_sig_V2[[j]] <- surv_data_clean[[j]]
  surv_data_sig_V2[[j]][[1]] <- surv_data_sig_V2[[j]][[1]][which(surv_data_sig_V2[[j]][[1]][,which.max(a[a!=0])] >= 0.8),]
  surv_data_sig_V2[[j]][[2]] <- surv_data_sig_V2[[j]][[2]][which(surv_data_sig_V2[[j]][[2]][,which.max(a[a!=0])] >= min(surv_data_sig_V2[[j]][[1]][,which.max(a[a!=0])])),]
  surv_data_sig_V2[[j]][[3]] <- surv_data_sig_V2[[j]][[3]][which(surv_data_sig_V2[[j]][[3]][,which.max(a[a!=0])] >= min(surv_data_sig_V2[[j]][[1]][,which.max(a[a!=0])])),]
}##读取数据以及计算sig_gene
names(surv_data_sig_V2) <- cancer_name



#Survival gene network was constructed
gene_protein <- read.csv('/database/gene_name_id_data/gene_protein_ID.txt',header = T)
gene_protein <- gene_protein[gene_protein[,3]!='',]

surv_raw_data_V2 <- list()
result <- list()
surv_sig_inter_V2 <- list()
surv_sig_dgree_V2 <- list()
# for (i in 1:length(cancer_name)) {
#   file_name <- dir(paste('result/',cancer_name[i],'/raw_data',sep = ''))
#   a <- which.min(which(apply(surv_data_sig_V2[[i]][[1]], 2, min)>0.1))
#   print(file_name[[a]])
# }
similar_value <- function(value){length(which(which(value==1)%in%which(a1==1)))/((length(which(a1==1))+length(which(value==1)))-2*length(which(which(a1==1)%in%which(value==1)))+1)}
# similar_value_all <- function(ob){
#     a1 <- ob
#     a2 <- surv_raw_data_V2[[i]]
#     apply(a2, 1, similar_value)
# }


length(cancer_name)
for (i in c(1:15)) {
  file_name <- dir(paste('result/',cancer_name[i],'/raw_data',sep = ''))
  a <- which.min(which(apply(surv_data_sig_V2[[i]][[1]], 2, min)>0.1))
  surv_raw_data_V2[[i]] <- read.csv(paste('result/',cancer_name[i],'/raw_data/',file_name[[a]],sep = ''),row.names = 1,header = T)
  surv_raw_data_V2[[i]] <- surv_raw_data_V2[[i]][rownames(surv_raw_data_V2[[i]])%in%unique(gene_protein$Gene.name),]
  surv_raw_data_V2[[i]] <- as.matrix(surv_raw_data_V2[[i]])
  surv_raw_data_V2[[i]][surv_raw_data_V2[[i]]>=0.05] <- 0
  surv_raw_data_V2[[i]][surv_raw_data_V2[[i]]!=0] <- 1
  a <- surv_raw_data_V2[[i]][rownames(surv_raw_data_V2[[i]])%in%rownames(surv_data_sig_V2[[i]][[1]]),]
  result[[i]] <- data.frame(row.names = rownames(surv_raw_data_V2[[i]]))
  # c <- apply(a[1:2,], 1, similar_value_all)
  for (x in 1:nrow(a)) {
    a1 <- a[x,]
    a2 <- surv_raw_data_V2[[i]]
    result[[i]][,x] <- apply(a2, 1, similar_value)
  }
  #   for (j in 1:nrow(surv_raw_data_V2[[i]])) {
  #     result[[i]][j,x] <- length(which(which(surv_raw_data_V2[[i]][j,]==1)%in%which(a[x,]==1)))/((length(which(a[x,]==1))+length(which(surv_raw_data_V2[[i]][j,]==1)))-2*length(which(which(a[x,]==1)%in%which(surv_raw_data_V2[[i]][j,]==1)))+1)
  #   }
  # }
  # for (x in 1:nrow(a)) {
  #   a1 <- a[x,]
  #   a2 <- surv_raw_data_V2[[i]]
  #   apply(array, margin, ...)
  #   for (j in 1:nrow(surv_raw_data_V2[[i]])) {
  #     result[[i]][j,x] <- length(which(which(surv_raw_data_V2[[i]][j,]==1)%in%which(a[x,]==1)))/((length(which(a[x,]==1))+length(which(surv_raw_data_V2[[i]][j,]==1)))-2*length(which(which(a[x,]==1)%in%which(surv_raw_data_V2[[i]][j,]==1)))+1)
  #   }
  # }
  # rownames(result[[i]]) <- rownames(surv_raw_data_V2[[i]])
  colnames(result[[i]]) <- rownames(a)
  result_matrix <- as.matrix(result[[i]])
  a <- seq(0.01,1,0.01)
  b <- data.frame()
  for (x in 1:100) {
    b[x,1] <- a[x]
    b[x,2] <- length(which(result_matrix>a[x]))-ncol(result_matrix)
  }
  plot(b)
  surv_sig_inter_V2[[i]] <- data.frame()
  for (x in 1:ncol(result[[i]])) {
    # a <- which(result[[i]][,x]>=max(b[which(b[,2]>=2000),1]))##此处还需更改，应当需要一个合适的数字
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
  # for (x in 1:nrow(surv_sig_inter_V2[[i]])) {
  #   a[x] <- paste0(surv_sig_inter_V2[[i]][x,1:3],collapse = '')
  #   b[x] <- paste0(surv_sig_inter_V2[[i]][x,3:1],collapse = '')
  # }
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
  # for (x in 1:nrow(surv_sig_inter_V2[[i]])) {
  #   ifelse(surv_sig_inter_V2[[i]][x,1]==surv_sig_inter_V2[[i]][x,3],surv_sig_inter_V2[[i]][x,2] <- 0,0)
  # }
  surv_sig_inter_V2[[i]] <- surv_sig_inter_V2[[i]][surv_sig_inter_V2[[i]][,2]!=0,]
  surv_sig_inter_V2[[i]] <- surv_sig_inter_V2[[i]][order(surv_sig_inter_V2[[i]][,2],decreasing = T),]
  surv_sig_inter_V2[[i]] <- surv_sig_inter_V2[[i]][surv_sig_inter_V2[[i]][,2]>=surv_sig_inter_V2[[i]][15000,2],]
  b <- surv_sig_inter_V2[[i]][,3:1]
  colnames(b) <- colnames(surv_sig_inter_V2[[i]])
  b <- rbind(surv_sig_inter_V2[[i]][,1:3],b)
  
  a <- unique(c(surv_sig_inter_V2[[i]][,1],surv_sig_inter_V2[[i]][,3]))
  surv_sig_dgree_V2[[i]] <- data.frame(a)
  sum_dgree <- function(ob){sum(b[,1]==ob)}
  surv_sig_dgree_V2[[i]][,2] <- apply(surv_sig_dgree_V2[[i]], 1, sum_dgree)
  # for (x in 1:length(a)) {
  #   surv_sig_dgree_V2[[i]][x,1] <- a[x]
  #   surv_sig_dgree_V2[[i]][x,2] <- sum(b[,1]==a[x])
  # }
  surv_sig_dgree_V2[[i]] <- surv_sig_dgree_V2[[i]][order(surv_sig_dgree_V2[[i]][,2],decreasing = T),]
  
  write.csv(surv_sig_inter_V2[[i]],paste('result/network/raw_data/',cancer_name[i],'_surv_sig_inter.csv',sep = ''),quote = F,row.names = F)
  write.csv(surv_sig_dgree_V2[[i]],paste('result/network/raw_data/',cancer_name[i],'_surv_sig_dgree.csv',sep = ''),quote = F,row.names = F)
  write.csv(result[[i]],paste('result/network/raw_data/',cancer_name[i],'_surv_sig_result.csv',sep = ''),quote = F)
  surv_raw_data_V2[[i]] <- data.frame()
  result[[i]] <- data.frame()
  a <- data.frame()
  b <- data.frame()
  c <- data.frame()
}

names(surv_sig_inter_V2) <- cancer_name
names(surv_sig_dgree_V2) <- cancer_name
names(result) <- cancer_name
