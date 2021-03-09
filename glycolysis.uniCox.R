

#install.packages("survival")

library(survival)

setwd("C:\\Users\\Administrator\\Desktop\\96glycolysis\\12.uniCox")                   #工作目录（需修改）
pFilter=0.05                                                                #显著性过滤条件
rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)    #读取输入文件
rt$futime=rt$futime/365
rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)

#读取差异分析结果文件
diffRT=read.table("diff.txt",header=T,sep="\t",check.names=F,row.names=1)       #读取输入文件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
	   if(sd(rt[,gene])<0.01){next}
	   cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
	   coxSummary = summary(cox)
	   coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	   if( ((coxSummary$conf.int[,"exp(coef)"]>1) & (diffRT[gene,"logFC"]>0)) | ((coxSummary$conf.int[,"exp(coef)"]<1) & (diffRT[gene,"logFC"]<0)) ){
		   if(coxP<pFilter){
		   	   group=ifelse(rt[,gene]>median(rt[,gene]),"high","low")
		       diff=survdiff(Surv(futime, fustat) ~group,data = rt)
		       pValue=1-pchisq(diff$chisq,df=1)
		       if(pValue<pFilter){
			       sigGenes=c(sigGenes,gene)
			       outTab=rbind(outTab,
			                    cbind(gene=gene,
			                         #KM=pValue,
			                         HR=coxSummary$conf.int[,"exp(coef)"],
			                         HR.95L=coxSummary$conf.int[,"lower .95"],
			                         HR.95H=coxSummary$conf.int[,"upper .95"],
					                 coxPvalue=coxP) )
				}
		   }
	   }
}
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)    #输出基因和p值表格文件
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)



