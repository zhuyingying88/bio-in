
#install.packages("vioplot")


library(vioplot)                      #引用包 
inputFile="CIBERSORT-Results.txt"     #输入文件
setwd("C:\\Users\\Administrator\\Desktop\\当前研究任务\\酸谈计划\\解螺旋兼职\\生信-阿尔兹海默症 3+相关项目\\数据\\R代码及其说明\\17.vioplot")     #设置工作目录

#读取免疫细胞浸润文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#对样品分组
con=grepl("_con", rownames(rt), ignore.case=T)
treat=grepl("_treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
rt=rbind(conData,treatData)

#输出小提琴图
outTab=data.frame()
pdf(file="vioplot.pdf", height=8, width=13)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.05),
     main="",xlab="", ylab="Immune Infiltration",
     pch=21,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制vioplot，对照组用蓝色表示，实验组用红色表示
for(i in 1:ncol(rt)){
	  if(sd(rt[1:conNum,i])==0){
	    rt[1,i]=0.00001
	  }
	  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
	    rt[(conNum+1),i]=0.00001
	  }
	  conData=rt[1:conNum,i]
	  treatData=rt[(conNum+1):(conNum+treatNum),i]
	  vioplot(conData,at=3*(i-1),lty=1,add = T,col = 'blue')
	  vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
	  wilcoxTest=wilcox.test(conData,treatData)
	  p=wilcoxTest$p.value
	  if(p<0.05){
	      cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
		  outTab=rbind(outTab,cellPvalue)
	  }
	  mx=max(c(conData,treatData))
	  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
	  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("Control", "Alzheimer's Disease"),
       lwd=3,bty="n",cex=1,
       col=c("blue","red"))
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()

#输出免疫细胞和p值表格文件
write.table(outTab,file="immuneDiff.xls",sep="\t",row.names=F,quote=F)



