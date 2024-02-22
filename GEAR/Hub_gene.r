surv_sig_inter_V2 <- list()
surv_sig_dgree_V2 <- list()
cancer_name <- read.csv('/path/cancer_name.csv')
for (i in 1:15) {
  surv_sig_inter_V2[[i]] <- read.csv(paste('/path/',cancer_name[i],'_surv_sig_inter.csv',sep = ''))
  surv_sig_dgree_V2[[i]] <- read.csv(paste('/path/',cancer_name[i],'_surv_sig_dgree.csv',sep = ''))
}
names(surv_sig_inter_V2) <- cancer_name
names(surv_sig_dgree_V2) <- cancer_name

#Calculate hub genes
hub_gene <- list()
# hub_type <- read.csv('',)
for (i in 1:15) {
  a <- surv_sig_inter_V2[[i]]
  colnames(a) <- c('gene1','value','gene2')
  b <- a[,3:1]
  colnames(b) <- c('gene1','value','gene2')
  a <- rbind(a,b)
  gene <- unique(a$gene1)
  b <- data.frame(1:length(gene))
  for (j in 1:length(gene)) {
    k1 <- a$gene2[a$gene1%in%gene[j]]
    k2 <- a$gene2[a$gene1%in%k1]
    k3 <- a$gene2[a$gene1%in%k2]
    b[j,2] <- length(unique(c(k1,k2,k3)))
    b[j,3] <- sum(a$gene1%in%gene[j])
  }
  b[,1] <- gene
  colnames(b) <- c('gene','k3','degree')
  b <- b[b$gene%in%rownames(surv_data_sig_V2[[i]][[1]]),]
  b <- b[order(b$k3,decreasing = T),]
  a <- which.min(which(apply(surv_data_sig_V2[[i]][[1]], 2, min)>0.1))
  a <- as.data.frame(cbind(rownames(surv_data_sig_V2[[i]][[1]]),surv_data_sig_V2[[i]][[1]][,a]))
  colnames(a) <- c('gene','surv_sig_p')
  b <- merge(b,a,by='gene',all=F)
  b <- b[order(b$k3,decreasing = T)[1:10],]
  # b <- b[b$k3>=median(b$k3)&b$degree>=median(b$degree),]
  hub_gene[[i]] <- b
}
names(hub_gene) <- cancer_name
