

##Required packages

library(pROC)

####Network expansion from neurodegenerative diseases
astro.combin=readRDS("input_tables\\\neuro_traits\\all_neuro.rds")
####Benchamrk associations from diseases database (Jensen lab) and ChEMBL
EFO_CHM=readRDS("input_tables\\\neuro_traits\\benchmark_file.rds.rds")

temp=astro.combin[,"EFO"]
temp=temp[!duplicated(temp)]

sum(temp%in%EFO_CHM[,"EFO"])

EFO_CHM[EFO_CHM[,"EFO"]%in%temp,c("count.chmbl.combin","AUC.chmbl.combin")]
EFO_CHM[EFO_CHM[,"EFO"]%in%temp,c("count.dis.COMBIN.2","AUC.COMBIN.2")]

efo$name[gsub("_",":",temp[temp%in%EFO_CHM[,"EFO"]])]

######Vamos a ir construyendo ROCs para cada una de las diseases

###ALS score Q3 & ChEMBL

temp.combin=cbind(astro.combin[	astro.combin[,"EFO"]=="EFO_0000253" & astro.combin[,"padj"]=="0",c("ENSG","page.rank")],0)
colnames(temp.combin)[ncol(temp.combin)]="ROC"
temp=unlist(strsplit(EFO_CHM[EFO_CHM[,"EFO"]=="EFO_0000253","unique.dis.cut.4"],";"))
temp.combin[temp.combin[,"ENSG"]%in%temp,"ROC"]=1

als.dis=roc(	as.numeric(temp.combin[,"ROC"]),
             as.numeric(temp.combin[,"page.rank"]))

temp.combin=cbind(astro.combin[	astro.combin[,"EFO"]=="EFO_0000253" & astro.combin[,"padj"]=="0",c("ENSG","page.rank")],0)
colnames(temp.combin)[ncol(temp.combin)]="ROC"
temp=unlist(strsplit(EFO_CHM[EFO_CHM[,"EFO"]=="EFO_0000253","unique.chmbl.combin"],";"))
temp.combin[temp.combin[,"ENSG"]%in%temp,"ROC"]=1

als.chm=roc(	as.numeric(temp.combin[,"ROC"]),
             as.numeric(temp.combin[,"page.rank"]))


###Parkinson score Q3 & ChEMBL

temp.combin=cbind(astro.combin[	astro.combin[,"EFO"]=="EFO_0002508" & astro.combin[,"padj"]=="0",c("ENSG","page.rank")],0)
colnames(temp.combin)[ncol(temp.combin)]="ROC"
temp=unlist(strsplit(EFO_CHM[EFO_CHM[,"EFO"]=="EFO_0002508","unique.dis.cut.4"],";"))
temp.combin[temp.combin[,"ENSG"]%in%temp,"ROC"]=1

par.dis=roc(	as.numeric(temp.combin[,"ROC"]),
             as.numeric(temp.combin[,"page.rank"]))

temp.combin=cbind(astro.combin[	astro.combin[,"EFO"]=="EFO_0002508" & astro.combin[,"padj"]=="0",c("ENSG","page.rank")],0)
colnames(temp.combin)[ncol(temp.combin)]="ROC"
temp=unlist(strsplit(EFO_CHM[EFO_CHM[,"EFO"]=="EFO_0002508","unique.chmbl.combin"],";"))
temp.combin[temp.combin[,"ENSG"]%in%temp,"ROC"]=1

par.chm=roc(	as.numeric(temp.combin[,"ROC"]),
             as.numeric(temp.combin[,"page.rank"]))

###Alzheimers score Q3 & ChEMBL

temp.combin=cbind(astro.combin[	astro.combin[,"EFO"]=="EFO_0000249" & astro.combin[,"padj"]=="0",c("ENSG","page.rank")],0)
colnames(temp.combin)[ncol(temp.combin)]="ROC"
temp=unlist(strsplit(EFO_CHM[EFO_CHM[,"EFO"]=="EFO_0000249","unique.dis.cut.4"],";"))
temp.combin[temp.combin[,"ENSG"]%in%temp,"ROC"]=1

alz.dis=roc(	as.numeric(temp.combin[,"ROC"]),
             as.numeric(temp.combin[,"page.rank"]))

temp.combin=cbind(astro.combin[	astro.combin[,"EFO"]=="EFO_0000249" & astro.combin[,"padj"]=="0",c("ENSG","page.rank")],0)
colnames(temp.combin)[ncol(temp.combin)]="ROC"
temp=unlist(strsplit(EFO_CHM[EFO_CHM[,"EFO"]=="EFO_0000249","unique.chmbl.combin"],";"))
temp.combin[temp.combin[,"ENSG"]%in%temp,"ROC"]=1

alz.chm=roc(	as.numeric(temp.combin[,"ROC"]),
             as.numeric(temp.combin[,"page.rank"]))


####Con otros valores para Disease

###ALS score Q3

temp.combin=cbind(astro.combin[	astro.combin[,"EFO"]=="EFO_0000253" & astro.combin[,"padj"]=="0",c("ENSG","page.rank")],0)
colnames(temp.combin)[ncol(temp.combin)]="ROC"
temp=unlist(strsplit(EFO_CHM[EFO_CHM[,"EFO"]=="EFO_0000253","unique.dis.cut.3"],";"))
temp.combin[temp.combin[,"ENSG"]%in%temp,"ROC"]=1

als.dis=roc(	as.numeric(temp.combin[,"ROC"]),
             as.numeric(temp.combin[,"page.rank"]))

temp.combin=cbind(astro.combin[	astro.combin[,"EFO"]=="EFO_0000253" & astro.combin[,"padj"]=="0",c("ENSG","page.rank")],0)
colnames(temp.combin)[ncol(temp.combin)]="ROC"
temp=unlist(strsplit(EFO_CHM[EFO_CHM[,"EFO"]=="EFO_0000253","unique.chmbl.combin"],";"))
temp.combin[temp.combin[,"ENSG"]%in%temp,"ROC"]=1

als.chm=roc(	as.numeric(temp.combin[,"ROC"]),
             as.numeric(temp.combin[,"page.rank"]))


###Parkinson score Q3

temp.combin=cbind(astro.combin[	astro.combin[,"EFO"]=="EFO_0002508" & astro.combin[,"padj"]=="0",c("ENSG","page.rank")],0)
colnames(temp.combin)[ncol(temp.combin)]="ROC"
temp=unlist(strsplit(EFO_CHM[EFO_CHM[,"EFO"]=="EFO_0002508","unique.dis.cut.3"],";"))
temp.combin[temp.combin[,"ENSG"]%in%temp,"ROC"]=1

par.dis=roc(	as.numeric(temp.combin[,"ROC"]),
             as.numeric(temp.combin[,"page.rank"]))

temp.combin=cbind(astro.combin[	astro.combin[,"EFO"]=="EFO_0002508" & astro.combin[,"padj"]=="0",c("ENSG","page.rank")],0)
colnames(temp.combin)[ncol(temp.combin)]="ROC"
temp=unlist(strsplit(EFO_CHM[EFO_CHM[,"EFO"]=="EFO_0002508","unique.chmbl.combin"],";"))
temp.combin[temp.combin[,"ENSG"]%in%temp,"ROC"]=1

par.chm=roc(	as.numeric(temp.combin[,"ROC"]),
             as.numeric(temp.combin[,"page.rank"]))

###Alzheimers score Q3

temp.combin=cbind(astro.combin[	astro.combin[,"EFO"]=="EFO_0000249" & astro.combin[,"padj"]=="0",c("ENSG","page.rank")],0)
colnames(temp.combin)[ncol(temp.combin)]="ROC"
temp=unlist(strsplit(EFO_CHM[EFO_CHM[,"EFO"]=="EFO_0000249","unique.dis.cut.3"],";"))
temp.combin[temp.combin[,"ENSG"]%in%temp,"ROC"]=1

alz.dis=roc(	as.numeric(temp.combin[,"ROC"]),
             as.numeric(temp.combin[,"page.rank"]))

temp.combin=cbind(astro.combin[	astro.combin[,"EFO"]=="EFO_0000249" & astro.combin[,"padj"]=="0",c("ENSG","page.rank")],0)
colnames(temp.combin)[ncol(temp.combin)]="ROC"
temp=unlist(strsplit(EFO_CHM[EFO_CHM[,"EFO"]=="EFO_0000249","unique.chmbl.combin"],";"))
temp.combin[temp.combin[,"ENSG"]%in%temp,"ROC"]=1

alz.chm=roc(	as.numeric(temp.combin[,"ROC"]),
             as.numeric(temp.combin[,"page.rank"]))


