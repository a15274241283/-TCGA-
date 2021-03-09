
#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)

setwd("C:\\Users\\Administrator\\Desktop\\96glycolysis\\14.survival")       #���ù���Ŀ¼

rt=read.table("risk.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
tab=table(rt$risk)
labels=paste0(names(tab),"(n=",tab,")")

fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)

#������������
pdf(file="survival.pdf",onefile = FALSE,
       width = 5.5,             #ͼƬ�Ŀ���
       height =5)               #ͼƬ�ĸ߶�
ggsurvplot(fit, 
           data=rt,
           pval=pValue,
           pval.size=6,
           legend.labs=labels,
           legend.title="Risk",
           xlab="Time(years)",
           break.time.by = 1,
           palette=c("red", "blue") )
dev.off()

summary(fit)    #�鿴����������

