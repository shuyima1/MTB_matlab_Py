library(ggplot2)
setwd('/Users/Shuyi/Dropbox/InfectiousDiseaseProject/MTB/Metabolism/gemini_files_code_data_0/')

fbootJaa2 = read.csv('fbootstrapgrowthratiosJaa2.txt',header=TRUE,sep='\t')
fbootG = read.csv('fbootstrapgrowthratiosGriffin.txt',header=TRUE,sep='\t')
fTFstats = read.csv('PROMBIG50TFstats.txt',header=TRUE,sep='\t')

#ggplot()+geom_boxplot(data=fbootJaa2,aes(fill=factor(TF),x=factor(TF,levels=fTFstats[with(fTFstats, order(Griffin)),1]),y=BootstrapgrRatio)) + geom_point(data=fTFstats,aes(x=factor(TF,levels=fTFstats[with(fTFstats, order(Griffin)),1]),y=Griffin),colour='red',size = 3) + labs(list(x = "TF", y = "Predicted Growth Ratio")) + theme(text= element_text (size=18),legend.position="none",axis.text.x = element_text(angle=90,vjust=0.5,size=16))
#ggplot()+geom_boxplot(data=fbootJaa2,aes(fill=factor(TF),x=factor(TF,levels=fTFstats[with(fTFstats, order(NumRxns)),1]),y=BootstrapgrRatio)) + geom_point(data=fTFstats,aes(x=factor(TF,levels=fTFstats[with(fTFstats, order(NumRxns)),1]),y=Jaa2grRatio),colour='red',size = 3) + labs(list(x = "TF", y = "Predicted Growth Ratio")) + theme(text= element_text (size=18),legend.position="none",axis.text.x = element_text(angle=90,vjust=0.5,size=16))

ggplot()+geom_boxplot(data=fbootJaa2,aes(fill=factor(Sassetti),x=factor(TF,levels=fTFstats[with(fTFstats, order(factor(Sassetti,levels=c('non-essential','essential')),NumRxns)),1]),y=BootstrapgrRatio)) +  scale_fill_manual(values=c("non-essential" = "light green","growth-defect" = "red","essential" = "red","no-data"="grey")) + labs(list(x = "TF", y = "Predicted Growth Ratio")) + theme(text= element_text (size=18),legend.position="none",axis.text.x = element_text(angle=90,vjust=0.5,size=16))


#ggplot()+geom_boxplot(data=fbootG,aes(fill=factor(TF),x=factor(TF,levels=fTFstats[with(fTFstats, order(NumRxns)),1]),y=BootstrapgrRatioGriffin)) + geom_point(data=fTFstats,aes(x=factor(TF,levels=fTFstats[with(fTFstats, order(NumRxns)),1]),y=GriffingrRatio),colour='red',size = 3) + labs(list(x = "TF", y = "Predicted Growth Ratio")) + theme(text= element_text (size=18),legend.position="none",axis.text.x = element_text(angle=90,vjust=0.5,size=16))
#ggplot()+geom_boxplot(data=fbootG,aes(fill=factor(TF),x=factor(TF,levels=fTFstats[with(fTFstats, order(Griffin)),1]),y=BootstrapgrRatioGriffin)) + geom_point(data=fTFstats,aes(x=factor(TF,levels=fTFstats[with(fTFstats, order(Griffin)),1]),y=Griffin),colour='red',size = 3) + labs(list(x = "TF", y = "Predicted Growth Ratio")) + theme(text= element_text (size=18),legend.position="none",axis.text.x = element_text(angle=90,vjust=0.5,size=16))

ggplot()+geom_boxplot(data=fbootG,aes(fill=factor(fbootJaa2$Griffin >= 0.15),x=factor(TF,levels=fTFstats[with(fTFstats, order(-Griffin,NumRxns)),1]),y=BootstrapgrRatioGriffin)) + scale_fill_manual(values=c("TRUE" = "light green", "FALSE" = "red")) + geom_point(data=fTFstats,aes(x=factor(TF,levels=fTFstats[with(fTFstats, order(-Griffin,NumRxns)),1]),y=Griffin),colour='blue',size = 3) + labs(list(x = "TF", y = "Predicted Growth Ratio")) + theme(text= element_text (size=18),legend.position="none",axis.text.x = element_text(angle=90,vjust=0.5,size=16))