## This code takes sensitivity and specificity data from single gene deletion simulations for iNJ661, GSMN-TB Vitro, and the updated models,
## and plots ROC curves (1-Specificity vs. Sensitivity) using ggplot2

library(ggplot2)
setwd('/Users/Shuyi/Dropbox/InfectiousDiseaseProject/MTB/Metabolism/')

## Reading in the 7H9 performance data
perfJam7H9 = read.csv('perfJam7H9.txt',header=FALSE)
perfBJ7H9 = read.csv('perfBJ7H9.txt',header=FALSE)
perfBV7H9 = read.csv('perfBV7H9.txt',header=FALSE)

## Reading in the Griffin performance data
sensBVGriffin = read.csv('sensBVGriffin.txt',header=FALSE)
specBVGriffin = read.csv('specBVGriffin.txt',header=FALSE)
specBJGriffin = read.csv('specBJGriffin.txt',header=FALSE)
sensBJGriffin = read.csv('sensBJGriffin.txt',header=FALSE)
sensJamGriffin = read.csv('sensJamGriffin.txt',header=FALSE)
specJamGriffin = read.csv('specJamGriffin.txt',header=FALSE)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette3 <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette4 <- c("#CC6666", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#perf7H9 = data.frame(sensitivity=c(perfJam7H9$V1,perfBV7H9$V1,perfBJ7H9$V1,0,1),specificity=c(perfJam7H9$V2,perfBV7H9$V2,perfBJ7H9$V2,1,0),model=c(rep(c('Jamshidi'),times=21),rep(c('BesteVitro'),times=21),rep(c('UpdatedModel'),times=21),rep(c('Random'),times=2)))
#perf7H9$model = factor(perf7H9$model,levels=c("Jamshidi","BesteVitro","UpdatedModel","Random"))
perf7H9 = data.frame(sensitivity=c(perfJam7H9$V1,perfBV7H9$V1,perfBJ7H9$V1),specificity=c(perfJam7H9$V2,perfBV7H9$V2,perfBJ7H9$V2),model=c(rep(c('iNJ661'),times=21),rep(c('GSMN-TB-Vitro'),times=21),rep(c('UpdatedModel'),times=21)))
perf7H9$model = factor(perf7H9$model,levels=c("iNJ661","GSMN-TB-Vitro","UpdatedModel"))
qplot(1-specificity,sensitivity,color = model,data = perf7H9, geom='line',xlab='1-Specificity',ylab='Sensitivity',main='EssentialityPrediction Performance for 7H9 Media')+scale_color_manual(values=cbbPalette4) + geom_line(size=2) + theme(plot.title = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=14), axis.title.x = element_text(size=16),axis.text.x = element_text(size=14),axis.title.y = element_text(size=16),axis.text.y = element_text(size=14))+coord_cartesian(ylim=c(0, 0.8),xlim=c(0,0.8))+geom_abline(intercept=0,slope=1,linetype='dashed')

#perfGriffin2 = data.frame(sensitivity=c(sensJamGriffin$V2,sensBVGriffin$V2,sensBJGriffin$V2,0,1),specificity=c(specJamGriffin$V2,specBVGriffin$V2,specBJGriffin$V2,1,0),model=c(rep(c('Jamshidi'),times=21),rep(c('BesteVitro'),times=21),rep(c('UpdatedModel'),times=21),rep(c('Random'),times=2)))
#perfGriffin2$model = factor(perfGriffin2$model,levels=c("Jamshidi","BesteVitro","UpdatedModel","Random"))
perfGriffin2 = data.frame(sensitivity=c(sensJamGriffin$V2,sensBVGriffin$V2,sensBJGriffin$V2),specificity=c(specJamGriffin$V2,specBVGriffin$V2,specBJGriffin$V2),model=c(rep(c('iNJ661'),times=21),rep(c('GSMN-TB-Vitro'),times=21),rep(c('UpdatedModel'),times=21)))
perfGriffin2$model = factor(perfGriffin2$model,levels=c("iNJ661","GSMN-TB-Vitro","UpdatedModel"))
qplot(1-specificity,sensitivity,color = model,data = perfGriffin2, geom='line',xlab='1-Specificity',ylab='Sensitivity',main='EssentialityPrediction Performance for Griffin Media')+scale_color_manual(values=cbbPalette4) + geom_line(size=2) + theme(plot.title = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=14), axis.title.x = element_text(size=16),axis.text.x = element_text(size=14),axis.title.y = element_text(size=16),axis.text.y = element_text(size=14))+coord_cartesian(ylim=c(0,0.8),xlim=c(0,0.8))+geom_abline(intercept=0,slope=1,linetype='dashed')


#qplot(1-specificity,sensitivity,color = model,data = perf7H9, geom='line',xlab='1-Specificity',ylab='Sensitivity',main='EssentialityPrediction Performance for 7H9 Media')+scale_color_manual(values=cbbPalette2) + geom_line(size=2) + theme(plot.title = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=14), axis.title.x = element_text(size=16),axis.text.x = element_text(size=14),axis.title.y = element_text(size=16),axis.text.y = element_text(size=14))+coord_cartesian(ylim=c(0, 1),xlim=c(0,1))

##Plot only the existing models
#perf7H9 = data.frame(sensitivity=c(perfJam7H9$V1,perfBV7H9$V1,perfBJ7H9$V1,0,1),specificity=c(perfJam7H9$V2,perfBV7H9$V2,perfBJ7H9$V2,1,0),model=c(rep(c('Jamshidi'),times=21),rep(c('BesteVitro'),times=21),rep(c('UpdatedModel'),times=21),rep(c('Random'),times=2)))
#perf7H9$model = factor(perf7H9$model,levels=c("Jamshidi","BesteVitro","UpdatedModel","Random"))
perf7H9b = data.frame(sensitivity=c(perfJam7H9$V1,perfBV7H9$V1),specificity=c(perfJam7H9$V2,perfBV7H9$V2),model=c(rep(c('iNJ661'),times=21),rep(c('GSMN-TB-Vitro'),times=21)))
perf7H9b$model = factor(perf7H9$model,levels=c("iNJ661","GSMN-TB-Vitro","UpdatedModel"))
qplot(1-specificity,sensitivity,color = model,data = perf7H9b, geom='line',xlab='1-Specificity',ylab='Sensitivity',main='EssentialityPrediction Performance for 7H9 Media')+scale_color_manual(values=cbbPalette4) + geom_line(size=2) + theme(plot.title = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=14), axis.title.x = element_text(size=16),axis.text.x = element_text(size=14),axis.title.y = element_text(size=16),axis.text.y = element_text(size=14))+coord_cartesian(ylim=c(0, 0.8),xlim=c(0,0.8))+geom_abline(intercept=0,slope=1,linetype='dashed')

#perfGriffin2 = data.frame(sensitivity=c(sensJamGriffin$V2,sensBVGriffin$V2,sensBJGriffin$V2,0,1),specificity=c(specJamGriffin$V2,specBVGriffin$V2,specBJGriffin$V2,1,0),model=c(rep(c('Jamshidi'),times=21),rep(c('BesteVitro'),times=21),rep(c('UpdatedModel'),times=21),rep(c('Random'),times=2)))
#perfGriffin2$model = factor(perfGriffin2$model,levels=c("Jamshidi","BesteVitro","UpdatedModel","Random"))
perfGriffinb = data.frame(sensitivity=c(sensJamGriffin$V2,sensBVGriffin$V2),specificity=c(specJamGriffin$V2,specBVGriffin$V2),model=c(rep(c('iNJ661'),times=21),rep(c('GSMN-TB-Vitro'),times=21)))
perfGriffinb$model = factor(perfGriffinb$model,levels=c("iNJ661","GSMN-TB-Vitro","UpdatedModel"))
qplot(1-specificity,sensitivity,color = model,data = perfGriffinb, geom='line',xlab='1-Specificity',ylab='Sensitivity',main='EssentialityPrediction Performance for Griffin Media')+scale_color_manual(values=cbbPalette4) + geom_line(size=2) + theme(plot.title = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=14), axis.title.x = element_text(size=16),axis.text.x = element_text(size=14),axis.title.y = element_text(size=16),axis.text.y = element_text(size=14))+coord_cartesian(ylim=c(0,0.8),xlim=c(0,0.8))+geom_abline(intercept=0,slope=1,linetype='dashed')


