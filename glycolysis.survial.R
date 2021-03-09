
#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)

setwd("C:\\Users\\Administrator\\Desktop\\96glycolysis\\14.survival")       #设置工作目录

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

#绘制生存曲线
pdf(file="survival.pdf",onefile = FALSE,
       width = 5.5,             #图片的宽度
       height =5)               #图片的高度
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

summary(fit)    #查看五年生存率


