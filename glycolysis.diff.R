#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")

library("limma")

setwd("C:\\Users\\Administrator\\Desktop\\96glycolysis\\10.diff")           #���ù���Ŀ¼
inputFile="glyGeneExp.txt"                                        #�����ļ�
pvalFilter=0.05                                                   #pֵ�ٽ�ֵ
logFCfilter=0                                                     #logFC�ٽ�ֵ
conNum=41                                                         #normal����Ʒ��Ŀ
treatNum=473                                                      #tumor����Ʒ��Ŀ

#��ȡ�����ļ�
outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#�������
for(i in row.names(data)){
	  geneName=unlist(strsplit(i,"\\|",))[1]
	  geneName=gsub("\\/", "_", geneName)
	  rt=rbind(expression=data[i,],grade=grade)
	  rt=as.matrix(t(rt))
	  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
	  conGeneMeans=mean(data[i,1:conNum])
	  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
	  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
	  pvalue=wilcoxTest$p.value
	  conMed=median(data[i,1:conNum])
	  treatMed=median(data[i,(conNum+1):ncol(data)])
	  diffMed=treatMed-conMed
	  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		  outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	  }
}

#������л���Ĳ������
write.table(outTab,file="all.txt",sep="\t",row.names=F,quote=F)

#����������
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$pValue))<pvalFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)
write.table(outDiff,file="diff.txt",sep="\t",row.names=F,quote=F)

#����������ı����ļ�
heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="diffGeneExp.txt",sep="\t",col.names=F,quote=F)
