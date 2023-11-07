
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


library(ConsensusClusterPlus)       #引用包
expFile="exp.txt"            #表达输入文件
workDir="E:\\data\\H.jieluoxuan\\阿尔兹海默病\\136Diagnostic\\10.subCluster"     #工作目录
setwd(workDir)       #设置工作目录

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

#聚类
maxK=9
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="pam",
              distance="euclidean",
              seed=123456,
              plot="png")


#输出分型结果
clusterNum=2        #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("m6Acluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$m6Acluster))
cluster$m6Acluster=letter[match(cluster$m6Acluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="subCluster.txt", sep="\t", quote=F, col.names=F)

