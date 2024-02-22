expression <- read.csv('/path/cancer_tpm.csv',row.names = 1) 
expression <- expression[apply(expression, 1, mean)>1,]
lifetime <- read.csv('/path/cancer_tpm.csv',row.names = 1)
lifetime[,6] <- gsub('-','.',lifetime[,6])
expression <- expression[,which(colnames(expression) %in% lifetime[,6])]
lifetime <- lifetime[which(lifetime[,6] %in% colnames(expression)),]
lifetime <- lifetime[which(!duplicated(lifetime[,6])),]
lifetime$Sample.ID <- factor(lifetime$Sample.ID, levels = colnames(expression))
lifetime <- lifetime[order(lifetime$Sample.ID),]
lifetime[which(lifetime[,2]=='Alive'),2] <- 0
lifetime[which(lifetime[,2]=='Dead'),2] <- 1

#sampling 
grade <- seq(from=0.1,by=0.1,length=1)
load(file="/sibcb1/jihongbinlab1/caixinlei/data/database/project/survival_analysis/3luad/samp_id_all.rdata")

#survical analysis
library(survival)
library(survminer)
for (k in 1) {
  HR <- data.frame()
  surv_p <- data.frame()
  for (j in 1:1000) {
    d <- expression[,as.numeric(samp_id[[k]][j,])]
    for (i in 1:nrow(expression)){tryCatch({
      a <- d[i,]
      b <- lifetime[lifetime[,6] %in% colnames(a),]
      b[b[,6]%in%colnames(a)[which(a <= median(as.numeric(a)))],7] <- 'low'
      b[b[,6]%in%colnames(a)[which(a > median(as.numeric(a)))],7] <- 'high'
      c <- survdiff(Surv(b[,5],as.numeric(b[,2]))~b[,7], data=b)
      surv_p[i,j] <- 1-pchisq(c$chisq,1)
      HR[i,j] <- (c$obs[2]/c$exp[2])/(c$obs[1]/c$exp[1])},warning = function(w){
      print('warning')}, error = function(e){surv_p[i,j] <- 1
                                             HR[i,j] <- 1
                                             print('error')})
    }
    print(paste(grade[k],j))
  }
  rownames(surv_p) <- rownames(expression)
  rownames(HR) <- rownames(expression)
  write.csv(surv_p,paste('/path/',grade[k],'_surv_p.csv',sep = ''))
  write.csv(HR,paste('/path/',grade[k],'_HR.csv',sep = ''))
  print(paste(grade[k],'_surv_ok'))
}

##output
surv_p_t <- data.frame()
for (i in 1:nrow(surv_p)){
     surv_p_t[i,1] <- length(which(surv_p[i,] <= 0.05))/1000
}
colnames(surv_p_t) <- c('sample_size')
rownames(surv_p_t) <- rownames(expression)
surv_p_up <- data.frame()
surv_p_down <- data.frame()
surv_p <- as.matrix(surv_p)
HR <- as.matrix(HR)
for (i in 1:nrow(surv_p)){
     surv_p_up[i,1] <- length(which(surv_p[i,] <= 0.05 & HR[i,] > 1))/1000
     surv_p_down[i,1] <- length(which(surv_p[i,] <= 0.05 & HR[i,] < 1))/1000
}
colnames(surv_p_up) <- c('sample_size')
colnames(surv_p_down) <- c('sample_size')
rownames(surv_p_up) <- rownames(expression)
rownames(surv_p_down) <- rownames(expression)
write.csv(surv_p_t,paste0('/path/',[sample_size],'_surv_p_t.csv'))
write.csv(surv_p_up,paste0('/path/',[sample_size],'_surv_p_up.csv'))
write.csv(surv_p_down,paste0('/path/',[sample_size],'_surv_p_down.csv'))
