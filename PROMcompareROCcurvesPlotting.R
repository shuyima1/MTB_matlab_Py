library(ggplot2)
setwd('/Users/Shuyi/Dropbox/InfectiousDiseaseProject/MTB/Metabolism/')

perfPROMoldSassetti = read.csv('perfPROMoldsassetti.txt',header=FALSE)
sensPROMoldGriffin = read.csv('sensPROMoldgriffin.txt',header=FALSE)
specPROMoldGriffin = read.csv('specPROMoldgriffin.txt',header=FALSE)

perfPROMBJ7aa2Sassetti = read.csv('perfPROMBJ7aa2sassetti.txt',header=FALSE)
sensPROMBJGriffin = read.csv('sensPROMBJgriffin.txt',header=FALSE)
specPROMBJGriffin = read.csv('specPROMBJgriffin.txt',header=FALSE)

perfPROMBJOPSassetti = read.csv('perfPROMBJOPsassetti.txt',header=FALSE)
sensPROMBJOPGriffin = read.csv('sensPROMBJOPgriffin.txt',header=FALSE)
specPROMBJOPGriffin = read.csv('specPROMBJOPgriffin.txt',header=FALSE)

cbbPalette3 <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette4 <- c("#CC6666", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette5 <- c("#E69F00", "#56B4E9", "#009E73",  "#D55E00", "#F0E442","#CC6666", "#000000", "#0072B2", "#CC79A7")

PROMoldPerfMedia = data.frame(sensitivity=c(perfPROMoldSassetti$V1,sensPROMoldGriffin$V2), specificity=c(perfPROMoldSassetti$V2,specPROMoldGriffin$V2), model=c(rep('7H9',times=51),rep('Griffin',times=51)))

PROMoldPerfSassettivsMet = data.frame(sensitivity=c(perfPROMoldSassetti$V1,perfJam7H9$V1), specificity=c(perfPROMoldSassetti$V2,perfJam7H9$V2), model=c(rep('MTBPROM1.0',times=51),rep('iNJ661',times=21)))

qplot(1-specificity,sensitivity,color = model,data = PROMoldPerfSassettivsMet, geom='line',xlab='1-Specificity',ylab='Sensitivity',main='Essentiality Performance in 7H9')+scale_color_manual(values=cbbPalette5) + geom_line(size=2) + theme(legend.position = "bottom", plot.title = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=14), axis.title.x = element_text(size=16),axis.text.x = element_text(size=14),axis.title.y = element_text(size=16),axis.text.y = element_text(size=14))+coord_cartesian(ylim=c(0, 0.8),xlim=c(0,0.8))+geom_abline(intercept=0,slope=1,linetype='dashed')

PROMsassetticompare = data.frame(sensitivity=c(perfPROMoldSassetti$V1,perfPROMBJ7aa2Sassetti$V1,perfPROMBJOPSassetti$V1), specificity=c(perfPROMoldSassetti$V2,perfPROMBJ7aa2Sassetti$V2,perfPROMBJOPSassetti$V2), model=c(rep('MTBPROM1.0',times=51),rep('MTBPROM2.0',times=51),rep('MTBPROM2.0+Operon',times=51)))
PROMsassetticompare2 = data.frame(sensitivity=c(perfPROMoldSassetti$V1,perfPROMBJ7aa2Sassetti$V1,perfJam7H9$V1), specificity=c(perfPROMoldSassetti$V2,perfPROMBJ7aa2Sassetti$V2,perfJam7H9$V2), model=c(rep('MTBPROM1.0',times=51),rep('MTBPROM2.0',times=51),rep('iNJ661',times=21)))

qplot(1-specificity,sensitivity,color = model,data = PROMsassetticompare, geom='line',xlab='1-Specificity',ylab='Sensitivity',main='Essentiality Performance in 7H9')+scale_color_manual(values=cbbPalette5) + geom_line(size=2) + theme(legend.position = "bottom", plot.title = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=14), axis.title.x = element_text(size=16),axis.text.x = element_text(size=14),axis.title.y = element_text(size=16),axis.text.y = element_text(size=14))+coord_cartesian(ylim=c(0, 0.8),xlim=c(0,0.8))+geom_abline(intercept=0,slope=1,linetype='dashed')

PROMgriffincompare = data.frame(sensitivity=c(sensPROMoldGriffin$V2,sensPROMBJGriffin$V2,sensPROMBJOPGriffin$V2), specificity=c(specPROMoldGriffin$V2,specPROMBJGriffin$V2,specPROMBJOPGriffin$V2), model=c(rep('MTBPROM1.0',times=51),rep('MTBPROM2.0',times=51),rep('MTBPROM2.0+Operon',times=51)))
PROMgriffincompare2 = data.frame(sensitivity=c(sensPROMoldGriffin$V2,sensPROMBJGriffin$V2,sensJamGriffin$V2), specificity=c(specPROMoldGriffin$V2,specPROMBJGriffin$V2,specJamGriffin$V2), model=c(rep('MTBPROM1.0',times=51),rep('MTBPROM2.0',times=51),rep('iNJ661',times=21)))
qplot(1-specificity,sensitivity,color = model,data = PROMgriffincompare, geom='line',xlab='1-Specificity',ylab='Sensitivity',main='Essentiality Performance in Griffin')+scale_color_manual(values=cbbPalette5) + geom_line(size=2) + theme(legend.position = "bottom", plot.title = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=14), axis.title.x = element_text(size=16),axis.text.x = element_text(size=14),axis.title.y = element_text(size=16),axis.text.y = element_text(size=14))+coord_cartesian(ylim=c(0, 1.1),xlim=c(0,1.1))+geom_abline(intercept=0,slope=1,linetype='dashed')
qplot(1-specificity,sensitivity,color = model,data = PROMgriffincompare2, geom='line',xlab='1-Specificity',ylab='Sensitivity',main='Essentiality Performance in Griffin')+scale_color_manual(values=cbbPalette5) + geom_line(size=2) + theme(legend.position = "bottom", plot.title = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=14), axis.title.x = element_text(size=16),axis.text.x = element_text(size=14),axis.title.y = element_text(size=16),axis.text.y = element_text(size=14))+coord_cartesian(ylim=c(0, 1.1),xlim=c(0,1.1))+geom_abline(intercept=0,slope=1,linetype='dashed')

