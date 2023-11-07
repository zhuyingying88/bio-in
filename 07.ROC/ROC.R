#install.packages("pROC")


library(pROC)                  #引用包
expFile="diffGeneExp.txt"      #表达数据文件
geneFile="gene.txt"      #交集基因列表文件
setwd("E:\\data\\H.jieluoxuan\\阿尔兹海默病\\136Diagnostic\\06.ROC")    #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="con", 0, 1)

#读取基因列表文件
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#对交集基因进行循环，绘制ROC曲线
for(x in as.vector(geneRT[,1])){
	#绘制ROC曲线
	roc1=roc(y, as.numeric(rt[x,]))
	ci1=ci.auc(roc1, method="bootstrap")
	ciVec=as.numeric(ci1)
	pdf(file=paste0("ROC.",x,".pdf"), width=5, height=5)
	plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=x)
	text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
	dev.off()
}



