
source("JTK_CYCLEv3.1.R")

project <- "Per3"

options(stringsAsFactors=FALSE)
annot <- read.delim("Per3 30nM.txt")
data <- read.delim("Per3 30nM.txt")

rownames(data) <- data[,1]
data <- data[,-1]
jtkdist(7,3)

periods <- 0:24
jtk.init(24)

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("JTK",project,"txt",sep="."),row.names=F,col.names=T,quote=F,sep="\t")

