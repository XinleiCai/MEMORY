grade <- seq(from=0.1,by=0.1,length=1)
samp_id <-list()
for (k in 1:10) {
  a <- data.frame()
  for(i in 1:1000){
    a[i,1:c(grade[k]*length(total_samples))] <- sample(1:length(total_samples),size=grade[k],replace = FALSE)
  }
  samp_id[[k]] <- a
  print(paste(grade[k],'sample_ok',sep=''))
}
names(samp_id) <- grade
save(samp_id, file="/path/samp_id_all.rdata")
