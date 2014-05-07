setwd('/Users/Shuyi/Dropbox/InfectiousDiseaseProject/MTB/Metabolism/gemini_files_code_data_0/')

library(pROC)

predPROM = read.csv('PROMmodJpredictions040714.txt',header=TRUE,sep='\t')
predPROMold = read.csv('PROMoldpredictions.txt',header=TRUE,sep='\t')

perfGriffin = data.frame(Griffin=predPROM$Griffin,PROM=predPROM$PROMGriffin)

perfSassetti=data.frame(Sassetti=predPROM$Sassetti,PROM=predPROM$PROMJaa2)
perfSassetti$Sassetti[perfSassetti$Sassetti == 'growth-defect'] = 'essential'

rocOLD = roc(response=predPROMold$Sassetti,predictor=predPROMold$PROMold,levels=c('non-essential','Essential'),ci=TRUE,plot=TRUE)

rocGriffin = roc(response=c(perfGriffin$Griffin < 0.15),predictor=perfGriffin$PROM,ci=TRUE,plot=TRUE)
rocSassetti = roc(response=perfSassetti$Sassetti,predictor=perfSassetti$PROM,levels=c('non-essential','essential'),ci=TRUE,plot=TRUE)

roc.test(rocSassetti,rocGriffin)
roc.test(rocOLD,rocSassetti)'
