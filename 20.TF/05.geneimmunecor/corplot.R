#install.packages("corrplot")


library(corrplot)           
inputFile="input.txt"       
setwd("E:\\data\\H.jieluoxuan\\�����Ⱥ�Ĭ��\\136Diagnostic\\16.geneimmunecor") 

rt=read.table(inputFile,sep="\t",header=T,row.names=1)      #读取文件
rt=t(rt)      #数据转置
M=cor(rt)     #相关型矩�?

#绘制相关性图�?
pdf(file="corpot1.pdf",width=25,height=25)
corrplot(M,
         method = "circle",
         order = "hclust", #聚类
         type = "upper",
         col=colorRampPalette(c("green", "white", "red"))(50)
         )
dev.off()

#第二个图
pdf(file="corpot2.pdf",width=20,height=20)
corrplot(M,
         order="original",
         method = "color",
         number.cex = 0.7, #相关系数
         addCoef.col = "black",
         diag = TRUE,
         tl.col="black",
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

