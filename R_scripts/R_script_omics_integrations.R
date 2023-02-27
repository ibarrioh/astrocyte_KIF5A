


####OMICS INTEGRATIONS
#
##Required packages
library(devtools)
library(igraph)
library(reshape2)

jaccard.IB<-function(set1,set2){
  jac=sum(set1%in%set2)/sum(!duplicated(c(set1,set2)))
  return(jac)
}

###########################
###########################
###Functions to run the network based expansion, "anotations" (to format the tables) and "astro" (expansion per se)
###There are 5 different datasets: proteomics, transcriptomics, GWAS (all ALS), GWAS (only familial), GWAS (only sporadic)
###Finally, we decided to use for the analysis GWAS "all ALS".

anotation<-function(node.list,genes,disease){
  
  genes=genes[genes[,"disease"]%in%disease,]
  
  for (i in 1:nrow(node.list)){
    if(node.list[i,"ENSG"]%in%genes[,"gene"]){
      node.list[i,"padj"]=max(as.numeric(genes[genes[,"gene"]%in%node.list[i,"ENSG"],"padj"]))
    }else{
      node.list[i,"padj"]=0
    }}
  
  return(node.list)
  
}

astro<-function(node.gwas,edge.string,all.nodes){
  
  ##Diffusion
  net=graph_from_data_frame(d=edge.string,vertices=node.gwas,directed=F)
  
  E(net)$weight=as.numeric(as.character(edge.string[,"combined_score"]))
  
  net.clean=igraph::simplify(net,
                             remove.loops = T,
                             remove.multiple = T ,
                             edge.attr.comb = c(weight="max","ignore"))
  
  page.rank=page_rank(net.clean, personalized=as.numeric(node.gwas[,"padj"]), weights=E(net.clean)$weight)
  
  node.gwas=cbind(node.gwas,page.rank$vector)
  colnames(node.gwas)[ncol(node.gwas)]="page.rank"
  
  deg=igraph::degree(net.clean)
  
  ##Network filter
  
  node.filter=node.gwas[as.numeric(node.gwas[,"page.rank"])>quantile(as.numeric(node.gwas[,"page.rank"]))[4],]
  
  colnames(node.filter)[1]="ENSP"
  
  edge.filter=edge.string[	as.character(edge.string[,1])%in%node.filter[,"ENSP"] & 
                             as.character(edge.string[,2])%in%node.filter[,"ENSP"] ,]
  
  node.filter=node.filter[node.filter[,"ENSP"]%in%c(as.character(edge.filter[,1]),as.character(edge.filter[,2])),]
  
  
  edge.string=edge.filter[,1:2]
  
  net=graph_from_data_frame(d=edge.string,vertices=node.filter,directed=F)
  
  E(net)$weight=as.numeric(as.character(edge.filter[,"combined_score"]))
  
  net.clean=igraph::simplify(net,
                             remove.loops = T,
                             remove.multiple = T ,
                             edge.attr.comb = c(weight="max","ignore"))
  
  cwt=cluster_walktrap(	net.clean, 
                        weights = E(net.clean)$weight, 
                        steps = 6,
                        merges = TRUE, 
                        modularity = TRUE, 
                        membership = TRUE)
  
  degree=igraph::degree(net.clean)
  
  node.filter=cbind(node.filter,degree,cwt$membership,cwt$modularity)
  
  colnames(node.filter)[(ncol(node.filter)-2):ncol(node.filter)]=c("node.degree","cluster.walktrap","modularity.walktrap")
  
  clust=as.matrix(as.data.frame(table(cwt$membership)))
  
  ####Recluster
  
  if(sum(as.numeric(clust[,2])>=300)>0){
    
    temp=rbind(c("0","0"),clust[as.numeric(clust[,2])>=300,])
    
    for (i in 1:nrow(temp)){
      
      node.re=node.filter[node.filter[,"cluster.walktrap"]==temp[i,1],]	
      edge.re=edge.filter[	as.character(edge.filter[,1])%in%node.re[,"ENSP"] &
                             as.character(edge.filter[,2])%in%node.re[,"ENSP"],]
      node.re=node.re[	node.re[,"ENSP"]%in%c(as.character(edge.re[,1]),as.character(edge.re[,2])),]
      
      net=graph_from_data_frame(d=as.data.frame(edge.re[,1:2]),vertices=node.re,directed=F)
      
      E(net)$weight=as.numeric(as.character(edge.re[,"combined_score"]))
      
      net.re=igraph::simplify(	net,
                               remove.loops = T,
                               remove.multiple = T ,
                               edge.attr.comb = c(weight="max","ignore"))
      
      cwt.re=cluster_walktrap(net.re, 
                              weights = E(net.re)$weight, 
                              steps = 6,
                              merges = TRUE, 
                              modularity = TRUE, 
                              membership = TRUE)
      
      node.re[,"cluster.walktrap"]=paste(node.re[,"cluster.walktrap"],cwt.re$membership,sep=";")
      node.re[,"modularity.walktrap"]=cwt.re$modularity
      
      node.filter[node.filter[,"cluster.walktrap"]%in%temp[i,1],"cluster.walktrap"]=node.re[,"cluster.walktrap"]
      node.filter[node.filter[,"cluster.walktrap"]%in%temp[i,1],"modularity.walktrap"]=node.re[,"modularity.walktrap"]
      
    }
    
  }
  
  ####Re-Recluster
  
  clust=as.matrix(as.data.frame(table(node.filter[,"cluster.walktrap"])))
  
  if(sum(as.numeric(clust[,2])>=300)>0){
    
    temp=rbind(c("0","0"),clust[as.numeric(clust[,2])>=300,])
    
    for (i in 1:nrow(temp)){
      
      node.re=node.filter[node.filter[,"cluster.walktrap"]==temp[i,1],]	
      edge.re=edge.filter[	as.character(edge.filter[,1])%in%node.re[,"ENSP"] &
                             as.character(edge.filter[,2])%in%node.re[,"ENSP"],]
      node.re=node.re[	node.re[,"ENSP"]%in%c(as.character(edge.re[,1]),as.character(edge.re[,2])),]
      
      net=graph_from_data_frame(d=as.data.frame(edge.re[,1:2]),vertices=node.re,directed=F)
      
      E(net)$weight=as.numeric(as.character(edge.re[,"combined_score"]))
      
      net.re=igraph::simplify(	net,
                               remove.loops = T,
                               remove.multiple = T ,
                               edge.attr.comb = c(weight="max","ignore"))
      
      cwt.re=cluster_walktrap(net.re, 
                              weights = E(net.re)$weight, 
                              steps = 6,
                              merges = TRUE, 
                              modularity = TRUE, 
                              membership = TRUE)
      
      node.re[,"cluster.walktrap"]=paste(node.re[,"cluster.walktrap"],cwt.re$membership,sep=";")
      node.re[,"modularity.walktrap"]=cwt.re$modularity
      
      node.filter[node.filter[,"cluster.walktrap"]%in%temp[i,1],"cluster.walktrap"]=node.re[,"cluster.walktrap"]
      node.filter[node.filter[,"cluster.walktrap"]%in%temp[i,1],"modularity.walktrap"]=node.re[,"modularity.walktrap"]
      
    }
    
  }
  
  ####Re-Re-Recluster
  
  clust=as.matrix(as.data.frame(table(node.filter[,"cluster.walktrap"])))
  
  if(sum(as.numeric(clust[,2])>=300)>0){
    
    temp=rbind(c("0","0"),clust[as.numeric(clust[,2])>=300,])
    
    for (i in 1:nrow(temp)){
      
      node.re=node.filter[node.filter[,"cluster.walktrap"]==temp[i,1],]	
      edge.re=edge.filter[	as.character(edge.filter[,1])%in%node.re[,"ENSP"] &
                             as.character(edge.filter[,2])%in%node.re[,"ENSP"],]
      node.re=node.re[	node.re[,"ENSP"]%in%c(as.character(edge.re[,1]),as.character(edge.re[,2])),]
      
      net=graph_from_data_frame(d=as.data.frame(edge.re[,1:2]),vertices=node.re,directed=F)
      
      E(net)$weight=as.numeric(as.character(edge.re[,"combined_score"]))
      
      net.re=igraph::simplify(	net,
                               remove.loops = T,
                               remove.multiple = T ,
                               edge.attr.comb = c(weight="max","ignore"))
      
      cwt.re=cluster_walktrap(net.re, 
                              weights = E(net.re)$weight, 
                              steps = 6,
                              merges = TRUE, 
                              modularity = TRUE, 
                              membership = TRUE)
      
      node.re[,"cluster.walktrap"]=paste(node.re[,"cluster.walktrap"],cwt.re$membership,sep=";")
      node.re[,"modularity.walktrap"]=cwt.re$modularity
      
      node.filter[node.filter[,"cluster.walktrap"]%in%temp[i,1],"cluster.walktrap"]=node.re[,"cluster.walktrap"]
      node.filter[node.filter[,"cluster.walktrap"]%in%temp[i,1],"modularity.walktrap"]=node.re[,"modularity.walktrap"]
      
    }
    
  }
  
  
  
  if(all.nodes==F){
    
    return(node.filter)	
    
  }else{
    
    node.gwas=cbind(node.gwas,deg)
    colnames(node.gwas)[ncol(node.gwas)]="degree"
    
    temp=node.filter[,c("ENSP","cluster.walktrap")]
    
    node.gwas=as.matrix(merge(node.gwas,temp,by.x="ENSG",by.y="ENSP",all.x=T))
    
    return(node.gwas)
    
  }
}


###########################
###########################
###Loading input tables with the starting signal (all.gene.gwas) and the 
###node information (all.node.gwas) of the interactome (string) 

all.gene.gwas=readRDS("input_tables\\omic_integration\\all_gene_gwas.rds")
all.node.gwas=readRDS("input_tables\\common\\all_node_gwas.rds")
string=readRDS("input_tables\\common\\Combined_STRING_Intact_Biogrid_FILTER.rds")

###Loop for  the network expansion

dis=as.matrix(as.data.frame(table(all.gene.gwas[,"disease"])))
dis[,2]=paste("trait",c(1:nrow(dis)),sep="_")
temp=table(all.gene.gwas[all.gene.gwas[,"gene"]%in%all.node.gwas[,"ENSG"],"disease"])
dis=dis[dis[,1]%in%names(temp[temp>2]),]

for (i in 1:nrow(dis)){
  
  print(i*100/nrow(dis))
  
  node.gwas=anotation(all.node.gwas,all.gene.gwas,disease=dis[i,1])
  
  temp.nodes=astro(node.gwas,as.data.frame(string),all.nodes=T)
  
  path=paste(	"output_tables\\omic_integration\\RDS_astro\\",
              "nodes.",
              dis[i,1],
              ".rds",	sep="")
  
  saveRDS(temp.nodes,path)
  rm(temp.nodes)
  rm(node.gwas)
  gc()
  
}

#############################
####combine the tables in one

setwd("output_tables\\omic_integration\\RDS_astro")

path=paste(	"output_tables\\omic_integration\\RDS_astro\\",
            list.files(pattern= ".rds"),sep="")

########
primer=cbind(readRDS(path[1]),list.files(pattern= ".rds")[1])

for (i in 2:length(path)){
  
  primer=rbind(primer,cbind(readRDS(path[i]),list.files(pattern= ".rds")[i]))
  
}

saveRDS(primer,"output_tables\\omic_integration\\all_together.rds")

astro.all=primer

#######################################################################################################
####Calculating significant modules 

colnames(astro.all)[7]="Trait"
astro.all[,"Trait"]=unlist(strsplit(astro.all[,"Trait"],"\\."))[c(F,T,F)]
astro.all=cbind(	astro.all,0,0,0,1,1)
colnames(astro.all)[(ncol(astro.all)-4):ncol(astro.all)]=c("Selected.cluster","Selected.fisher","Selected.KS","padj.fisher","padj.KS")
dis=astro.all[!duplicated(astro.all[,"Trait"]),"Trait"]

for (i in 1:length(dis)){
  
  print(i)
  
  node=astro.all[astro.all[,"Trait"]%in%dis[i],]
  
  clusters=as.matrix(as.data.frame(table(node[,"cluster.walktrap"])))
  clusters=cbind(clusters,0,0,0)
  colnames(clusters)=c("clust","all.nodes","gwa.nodes","fisher","KS")
  clusters=clusters[as.numeric(clusters[,"all.nodes"])>=10,]
  
  for (j in 1:nrow(clusters)){
    
    clusters[j,"gwa.nodes"]=sum(node[node[,"cluster.walktrap"]%in%clusters[j,"clust"],"padj"]!="0")
    
  }
  
  clusters=clusters[as.numeric(clusters[,"gwa.nodes"])>0,]
  
  if(length(clusters)>5){
    ##KS	
    
    for (j in 1:nrow(clusters)){
      
      x=log10(as.numeric(node[!is.na(node[,"cluster.walktrap"]),"page.rank"]))
      y=log10(as.numeric(node[node[,"cluster.walktrap"]==clusters[j,"clust"],"page.rank"]))
      
      clusters[j,"KS"]=ks.test(x,y,alternative="greater")$p.value
      
    }
    
    ##Fisher
    
    for (j in 1:nrow(clusters)){
      
      matrix=cbind(	c(as.numeric(clusters[j,"gwa.nodes"]),sum(node[,"padj"]!="0")-as.numeric(clusters[j,"gwa.nodes"])),
                    c(as.numeric(clusters[j,"all.nodes"]),sum(!duplicated(node[,1]))-as.numeric(clusters[j,"all.nodes"])))
      
      
      clusters[j,"fisher"]=fisher.test(matrix, alternative="greater")$p.value
    }
    
    clusters[,"fisher"]=p.adjust(as.numeric(clusters[,"fisher"]),method="BH")
    clusters[,"KS"]=p.adjust(as.numeric(clusters[,"KS"]),method="BH")
    
    
    node=cbind(	node,0,0,0,1,1)
    
    colnames(node)[(ncol(node)-4):ncol(node)]=c("Selected.cluster","Selected.fisher","Selected.KS","padj.fisher","padj.KS")
    
    node[node[,"cluster.walktrap"]%in%clusters[,"clust"],"Selected.cluster"]=1
    node[node[,"cluster.walktrap"]%in%clusters[as.numeric(clusters[,"fisher"])<=0.05,"clust"],"Selected.fisher"]=1	
    node[node[,"cluster.walktrap"]%in%clusters[as.numeric(clusters[,"KS"])<=0.05,"clust"],"Selected.KS"]=1
    
    for (j in 1:nrow(clusters)){
      
      node[node[,"cluster.walktrap"]%in%clusters[j,"clust"],"padj.KS"]=clusters[j,"KS"]
      node[node[,"cluster.walktrap"]%in%clusters[j,"clust"],"padj.fisher"]=clusters[j,"fisher"]
      
    }
    
    astro.all[astro.all[,"Trait"]%in%dis[i],c(	"Selected.cluster",
                                               "Selected.fisher",
                                               "Selected.KS",
                                               "padj.fisher",
                                               "padj.KS")]=node[,c(	"Selected.cluster",
                                                                    "Selected.fisher",	
                                                                    "Selected.KS",
                                                                    "padj.fisher",
                                                                    "padj.KS")]
    
    
  }else{	
    
    astro.all[astro.all[,"Trait"]%in%dis[i],c("Selected.cluster",
                                              "Selected.fisher",
                                              "Selected.KS",
                                              "padj.fisher",
                                              "padj.KS")]=NA
    
    
  }}

saveRDS(astro.all,"output_tables\\omic_integration\\all_together_sig.rds")


#######################################################################################################
####Calculating Jaccar index for modudule overlaps for fig2B

####We select GWAS all ALS (called "combo")

sig.KS=table(astro.all[astro.all[,"Selected.KS"]=="1" & 
                         astro.all[,"Trait"]%in%c("Prot","Tran","combo"),"cluster.walktrap"])

jaccard=matrix(0,length(sig.KS),length(sig.KS))

colnames(jaccard)=names(sig.KS)
rownames(jaccard)=names(sig.KS)

for (i in 1:nrow(jaccard)){
  
  for (j in 1:ncol(jaccard)){
    
    A=astro.all[astro.all[,"cluster.walktrap"]==rownames(jaccard)[i],"ENSG"]
    B=astro.all[astro.all[,"cluster.walktrap"]==colnames(jaccard)[j],"ENSG"]
    
    jaccard[i,j]=jaccard.IB(A,B)
    
  }}


###2 column conversion for network figure

ut <- upper.tri(jaccard)
jaccard.table=data.frame(	i = rownames(jaccard)[row(jaccard)[ut]],
                          j = rownames(jaccard)[col(jaccard)[ut]],
                          jaccard=t(jaccard)[ut])

jaccard.table=as.matrix(jaccard.table)

jaccard.table=jaccard.table[order(as.numeric(jaccard.table[,3]),decreasing=T),]

sum(as.numeric(jaccard.table[,3])>0)


node.jac=c(jaccard.table[as.numeric(jaccard.table[,3])>0,1],jaccard.table[as.numeric(jaccard.table[,3])>0,2])

node.jac=cbind(node.jac[!duplicated(node.jac)],0)


edge.jac=jaccard.table[as.numeric(jaccard.table[,3])>0,]

colnames(edge.jac)=c("A","B","jaccard.index")

colnames(node.jac)=c("ID.name","Set")

node.jac[,"Set"]=unlist(strsplit(node.jac[,"ID.name"],"-"))[c(T,F)]

#write.csv(node.jac,"output_tables\\omic_integration\\ThreeSets_Jaccard_node.csv")
#write.csv(edge.jac,"output_tables\\omic_integration\\ThreeSets_Jaccard_edge.csv")

##############################################################################################################################################
##############################################################################################################################################















