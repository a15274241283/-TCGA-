
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")

library(limma)
library(ggpubr)

setwd("C:\\Users\\Administrator\\Desktop\\96glycolysis\\19.scatter")        #设置工作目录
expFile="glyGeneExp.txt"                                          #表达输入文件
riskFile="risk.txt"                                               #风险输入文件
conNum=41                                                         #normal组样品数目
treatNum=473                                                     #tumor组样品数目

#读取表达输入文件
outTab=data.frame()
Group=c(rep("Normal",conNum),rep("Tumor",treatNum))
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#读取风险输入文件
risk=read.table(riskFile,sep="\t",header=T,row.names=1,check.names=F)
modelGene=colnames(risk)[3:(ncol(risk)-2)]

#提取模型基因表达量
data=t(data[modelGene,])
data=cbind(data,Group)
data=rbind(ID=colnames(data),data)
write.table(data,file="data.txt",sep="\t",col.names=F,quote=F)

#绘制图形
data=read.table("data.txt",sep="\t",header=T,row.names=1,check.names=F)
for(gene in colnames(data)[1:(ncol(data)-1)]){
	subData=data[,c(gene,"Group")]
	colnames(subData)=c("gene","Group")
	#设置比较组
	group=levels(factor(subData$Group))
	comp=combn(group,2)
	my_comparisons=list()
    for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	#绘制boxplot
	boxplot=ggboxplot(subData, x="Group", y="gene", color="Group",
	          xlab="",
	          ylab=paste(gene,"expression"),
	          legend.title="",
	          add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	pdf(file=paste0(gene,".pdf"),width=5.5,height=5)
	print(boxplot)
	dev.off()
}



