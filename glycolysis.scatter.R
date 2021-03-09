
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")

library(limma)
library(ggpubr)

setwd("C:\\Users\\Administrator\\Desktop\\96glycolysis\\19.scatter")        #���ù���Ŀ¼
expFile="glyGeneExp.txt"                                          #���������ļ�
riskFile="risk.txt"                                               #���������ļ�
conNum=41                                                         #normal����Ʒ��Ŀ
treatNum=473                                                     #tumor����Ʒ��Ŀ

#��ȡ���������ļ�
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

#��ȡ���������ļ�
risk=read.table(riskFile,sep="\t",header=T,row.names=1,check.names=F)
modelGene=colnames(risk)[3:(ncol(risk)-2)]

#��ȡģ�ͻ��������
data=t(data[modelGene,])
data=cbind(data,Group)
data=rbind(ID=colnames(data),data)
write.table(data,file="data.txt",sep="\t",col.names=F,quote=F)

#����ͼ��
data=read.table("data.txt",sep="\t",header=T,row.names=1,check.names=F)
for(gene in colnames(data)[1:(ncol(data)-1)]){
	subData=data[,c(gene,"Group")]
	colnames(subData)=c("gene","Group")
	#���ñȽ���
	group=levels(factor(subData$Group))
	comp=combn(group,2)
	my_comparisons=list()
    for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	#����boxplot
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


