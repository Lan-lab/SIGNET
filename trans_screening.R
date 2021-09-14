Copaierd<-function(Outlier){
  copaired <- data.frame()
  num = 1
  for (i in 1:dim(Outlier)[1]) {
    print(i)
    for(j in 1:dim(Outlier)[2]){
      if(Outlier[i,j]!=0){
        copaired[num,1] = rownames(Outlier)[i]
        copaired[num,2] = colnames(Outlier)[j]
        num = num +1
      }
    }
  }
  print(num)
  return(copaired)
}

Genesets<-function(copaired){
  genesets_list <- levels(as.factor(copaired$V1))
  genesets= as.data.frame(matrix(nrow=length(genesets_list),ncol=1))
  for (i in 1:length(genesets_list)){
    genesets[i,1] = paste(copaired[copaired$V1==genesets_list[i],2],collapse = ",")
  }
  rownames(genesets)=genesets_list
  colnames(genesets)="gene_list"
  return (genesets)
}
# add tf gene into genesets
Genesets2<-function(copaired){
  genesets_list <- levels(as.factor(copaired$V1))
  genesets= as.data.frame(matrix(nrow=length(genesets_list),ncol=1))
  for (i in 1:length(genesets_list)){
    genesets[i,1] = paste(copaired[copaired$V1==genesets_list[i],2],collapse = ",")
  }
  rownames(genesets)=genesets_list
  colnames(genesets)="gene_list"
  return (genesets)
}


# gESD.test version 3.0
gESD.test <- function(x,alpha=0.05,iteration){
  DNAME <- deparse1(substitute(x))
  stopifnot(is.numeric(x))
  n <- length(x)
  if (iteration > n) 
    stop("iteration times is bigger than the length of 'x'")
  if (sd(x)==0)
    return(NULL)
  i <- 1
  X <- numeric()
  I <- numeric()
  R <- numeric()
  L <- numeric()
  ares <- abs(x - mean(x))/sd(x)
  num <- c(1:n)
  tdf <- data.frame(num, x, ares)
  for (i in 1:iteration) {
    if (i == 1){
      R[i] <- max(ares)
      I[i] <- which(ares==max(ares))[1]
      X[i] <- x[I[i]]
      p <- 1 - alpha/(2*(n-i+1))
      t <- qt(p,n-i-1)
      L[i] <- t*(n-i) / sqrt((n-i-1+t**2)*(n-i+1))
    }
    else if (i != 1){
      tdf <- tdf[order(tdf$ares, decreasing = T),]
      tdf <- tdf[-1,]
      if(sd(tdf$x)==0)
        return(NULL)
      ares <- abs(tdf$x - mean(tdf$x))/sd(tdf$x)
      tdf$ares <- ares
      R[i] <- max(ares)
      I[i] <- tdf$num[which(ares==max(ares))[1]]
      X[i] <- tdf$x[which(ares==max(ares))[1]]
      p <- 1 - alpha/(2*(n-i+1))
      t <- qt(p,n-i-1)
      L[i] <- t*(n-i) / sqrt((n-i-1+t**2)*(n-i+1))
    }
  }
  D <- R - L
  if (is.infinite(max(which(D>0)))){
    return(NULL)
  }
  else if (!is.infinite(max(which(D>0)))){
    out.num <- max(which(D>0))
  }
  df <- data.frame(c(1:iteration),R,L)
  names(df) <- c("No. Outliers","Test Stat.","Critical Val.")
  out <- data.frame(I[1:out.num],X[1:out.num])
  met <- "Generalized ESD test"
  RVAL <- list(statistic = c(n = round(out.num), d = D[1:out.num]), 
               alternative = paste("there are up to", iteration, "outliers"),
               method = met, data.name = DNAME, 
               outlier.table = df,
               outliers = out)
  class(RVAL) <- "htest"
  return(RVAL)
}



library(RcisTarget)
library(umap)
data(motifAnnotations_mgi)
data(motifAnnotations_hgnc)
#mouse
motifRankings <- importRankings("~/cisTarget_databases/mm9-tss-centered-10kb-7species.mc9nr.feather")
motifRankings2 <- importRankings("~/cisTarget_databases/mm9-500bp-upstream-7species.mc9nr.feather")
#human
data_count<-read.csv()
motifRankingsh1 <- importRankings("/Users/mac/Desktop/cisTarget_databases/hg19-tss-centered-10kb-7species.mc9nr.feather")
motifRankingsh2 <- importRankings("/Users/mac/Desktop/cisTarget_databases/hg19-500bp-upstream-7species.mc9nr.feather")

data_count<-read.csv("~/data.csv",row.names=1) # the raw data
co1<-read.table("~/co_fc_ntf.txt")
co2<-read.table("~/co_fc_tf.txt")
coexpressed<-rbind(co1,co2)
gene_ntf<-read.csv("~/gene_ntf.csv")
gene_tf<-read.csv("~/gene_tf.csv")

# only retain the gene in the databases
gene <- colnames(motifRankings)
gene_ntf = as.character(gene_ntf[,2])
gene_tf = as.character(gene_tf[,2])
coexpressed = coexpressed[,-dim(coexpressed)[2]]
gene_total <- c(gene_ntf,gene_tf)
rownames(coexpressed) <-gene_total
colnames(coexpressed) <-gene_tf
coexpressed <- t(coexpressed)
gene_total<-gene_total[gene_total%in%gene]
coexpressed<-coexpressed[,gene_total]
gene_tf<-gene_tf[gene_tf%in%gene]
coexpressed<-coexpressed[gene_tf,]


# outlier screening, the alpha can be set
alpha = 0.05
O<-as.data.frame(matrix(0,nrow = dim(coexpressed)[1],ncol = dim(coexpressed)[2]))
rownames(O)<-rownames(coexpressed)
colnames(O)<-colnames(coexpressed)
#gESD test
O<-as.data.frame(matrix(0,nrow = dim(coexpressed)[1],ncol = dim(coexpressed)[2]))
rownames(O)<-rownames(coexpressed)
colnames(O)<-colnames(coexpressed)
iterations = floor(dim(coexpressed)[1]*0.3)
for (i in 1:dim(coexpressed)[2]) {
  res=gESD.test(as.numeric(coexpressed[,i]),alpha,iterations)
  print(i)
  if(!is.null(res)){
    outliers = res$outliers
    outliers = outliers[,1]
    O[outliers,i]=1
  }
}


copaired = Copaierd(O)
genesets = Genesets(copaired)
genesets_list <- levels(as.factor(copaired$V1))

# Regulatory modules screening
Regulons <- as.data.frame(matrix(0,dim(O)[1],dim(O)[2]))
rownames(Regulons)<-rownames(O)
colnames(Regulons)<-colnames(O)
ntf_gene = colnames(O)
tf_gene = rownames(O)

# mouse
num=0
for (q in 1:length(genesets_list)) {
  tf_name = rownames(genesets)[q]
  index = which(tf_gene%in%tf_name)
  genelist = unlist(strsplit(genesets[q,],split=","))
  motifEnrichmentTable_wGenes <- cisTarget(genelist, 
                                           motifRankings,
                                           motifAnnot=motifAnnotations_mgi)
  # motifAnnotations_mgi
  motifEnrichmentTable_wGenes2 <- cisTarget(genelist, 
                                            motifRankings2,
                                            motifAnnot=motifAnnotations_mgi)
  # 处理数据格式，将RcisTarget数据转化为data.frame
  if(!is.null(motifEnrichmentTable_wGenes$TF_highConf)){
    TF_hlists<-strsplit(motifEnrichmentTable_wGenes$TF_highConf,split = ";")
  }
  else{
    TF_hlists<-NULL
  }
  if(!is.null(motifEnrichmentTable_wGenes$TF_lowConf)){
    TF_llists<- strsplit(motifEnrichmentTable_wGenes$TF_lowConf,split = ";")
  }
  else{
    TF_llists<-NULL
  }
  if(!is.null(motifEnrichmentTable_wGenes$enrichedGenes)){
    regulated_gene<-strsplit(motifEnrichmentTable_wGenes$enrichedGenes,split = ";")
  }
  else{
    regulated_gene<-NULL
  }
  for (i in 1:length(TF_hlists)){
    str1 = TF_hlists[[i]]
    str2 = TF_llists[[i]]
    str3 = regulated_gene[[i]]
    str3 = unique(str3)
    if(length(str1)>0){
      str1 = gsub("\\(.*\\)","",str1)
      str1 = gsub("\\.","",str1)
      str1 = gsub(" ","",str1)
      str1 = tolower(str1)
      for(i in 1:length(str1)){
        str=unlist(strsplit(str1[i],""))
        str[1]=toupper(str[1])
        str=paste(str,collapse = "")
        str1[i]=str
      }
      if(tf_name %in% str1){
        for(k in 1:length(str3)){
          n = which(ntf_gene%in%str3[k])
          if(length(n)!=0){
            if(O[index,n]==1){
              Regulons[index,n]=1
              num = num + 1
            }
          }
        }
      }
    }
    if(length(str2)>0){
      str2 = gsub("\\(.*\\)","",str2)
      str2 = gsub("\\.","",str2)
      str2 = gsub(" ","",str2)
      str2 = tolower(str2)
      for(i in 1:length(str2)){
        str=unlist(strsplit(str2[i],""))
        str[1]=toupper(str[1])
        str=paste(str,collapse = "")
        str2[i]=str
      }
      if(tf_name %in% str2){
        for(k in 1:length(str3)){
          n = which(ntf_gene%in%str3[k])
          if(length(n)!=0){
            if(O[index,n]==1){
              Regulons[index,n]=1
              num = num + 1
            }
          }
        }
      }
    }
  }
  if(!is.null(motifEnrichmentTable_wGenes2$TF_highConf)){
    TF_hlists<-strsplit(motifEnrichmentTable_wGenes2$TF_highConf,split = ";")
  }
  else{
    TF_hlists<-NULL
  }
  if(!is.null(motifEnrichmentTable_wGenes2$TF_lowConf)){
    TF_llists<- strsplit(motifEnrichmentTable_wGenes2$TF_lowConf,split = ";")
  }
  else{
    TF_llists<-NULL
  }
  if(!is.null(motifEnrichmentTable_wGenes2$enrichedGenes)){
    regulated_gene<-strsplit(motifEnrichmentTable_wGenes2$enrichedGenes,split = ";")
  }
  else{
    regulated_gene<-NULL
  }
  for (i in 1:length(TF_hlists)){
    str1 = TF_hlists[[i]]
    str2 = TF_llists[[i]]
    str3 = regulated_gene[[i]]
    str3 = unique(str3)
    if(length(str1)>0){
      str1 = gsub("\\(.*\\)","",str1)
      str1 = gsub("\\.","",str1)
      str1 = gsub(" ","",str1)
      str1 = tolower(str1)
      for(i in 1:length(str1)){
        str=unlist(strsplit(str1[i],""))
        str[1]=toupper(str[1])
        str=paste(str,collapse = "")
        str1[i]=str
      }
      
      if(tf_name%in%str1){
        for(k in 1:length(str3)){
          n = which(ntf_gene%in%str3[k])
          if(length(n)!=0){
            if(O[index,n]==1){
              Regulons[index,n]=1
              num = num + 1
            }
          }
        }
      }
    }
    if(length(str2)>0){
      str2 = gsub("\\(.*\\)","",str2)
      str2 = gsub("\\.","",str2)
      str2 = gsub(" ","",str2)
      str2 = tolower(str2)
      for(i in 1:length(str2)){
        str=unlist(strsplit(str2[i],""))
        str[1]=toupper(str[1])
        str=paste(str,collapse = "")
        str2[i]=str
      }
      if(tf_name%in%str2){
        for(k in 1:length(str3)){
          n = which(ntf_gene%in%str3[k])
          if(length(n)!=0){
            if(O[index,n]==1){
              Regulons[index,n]=1
              num = num + 1
            }
          }
        }
      }
    }
  }
  
}

# human 
num = 0
for (q in 1:length(genesets_list)) {
  tf_name = rownames(genesets)[q]
  index = which(tf_gene%in%tf_name)
  print(tf_name)
  genelist = unlist(strsplit(genesets[q,],split=","))
  motifEnrichmentTable_wGenes <- cisTarget(genelist, 
                                           motifRankingsh1,
                                           motifAnnot=motifAnnotations_hgnc)
  # motifAnnotations_hgnc
  motifEnrichmentTable_wGenes2 <- cisTarget(genelist, 
                                            motifRankingsh2,
                                            motifAnnot=motifAnnotations_hgnc)
  # 处理数据格式，将RcisTarget数据转化为data.frame
  if(!is.null(motifEnrichmentTable_wGenes$TF_highConf)){
    TF_hlists<-strsplit(motifEnrichmentTable_wGenes$TF_highConf,split = ";")
  }
  else{
    TF_hlists<-NULL
  }
  if(!is.null(motifEnrichmentTable_wGenes$TF_lowConf)){
    TF_llists<- strsplit(motifEnrichmentTable_wGenes$TF_lowConf,split = ";")
  }
  else{
    TF_llists<-NULL
  }
  if(!is.null(motifEnrichmentTable_wGenes$enrichedGenes)){
    regulated_gene<-strsplit(motifEnrichmentTable_wGenes$enrichedGenes,split = ";")
  }
  else{
    regulated_gene<-NULL
  }
  for (i in 1:length(TF_hlists)){
    str1 = TF_hlists[[i]]
    str2 = TF_llists[[i]]
    str3 = regulated_gene[[i]]
    str3 = unique(str3)
    if(length(str1)>0){
      str1 = gsub("\\(.*\\)","",str1)
      str1 = gsub("\\.","",str1)
      str1 = gsub(" ","",str1)
      if(tf_name %in% str1){
        for(k in 1:length(str3)){
          n = which(ntf_gene%in%str3[k])
          if(length(n)!=0){
            if(O[index,n]==1){
              Regulons[index,n]=1
              num = num + 1
            }
          }
        }
      }
    }
    if(length(str2)>0){
      str2 = gsub("\\(.*\\)","",str2)
      str2 = gsub("\\.","",str2)
      str2 = gsub(" ","",str2)
      if(tf_name %in% str2){
        for(k in 1:length(str3)){
          n = which(ntf_gene%in%str3[k])
          if(length(n)!=0){
            if(O[index,n]==1){
              Regulons[index,n]=1
              num = num + 1
            }
          }
        }
      }
    }
  }
  if(!is.null(motifEnrichmentTable_wGenes2$TF_highConf)){
    TF_hlists<-strsplit(motifEnrichmentTable_wGenes2$TF_highConf,split = ";")
  }
  else{
    TF_hlists<-NULL
  }
  if(!is.null(motifEnrichmentTable_wGenes2$TF_lowConf)){
    TF_llists<- strsplit(motifEnrichmentTable_wGenes2$TF_lowConf,split = ";")
  }
  else{
    TF_llists<-NULL
  }
  if(!is.null(motifEnrichmentTable_wGenes2$enrichedGenes)){
    regulated_gene<-strsplit(motifEnrichmentTable_wGenes2$enrichedGenes,split = ";")
  }
  else{
    regulated_gene<-NULL
  }
  for (i in 1:length(TF_hlists)){
    str1 = TF_hlists[[i]]
    str2 = TF_llists[[i]]
    str3 = regulated_gene[[i]]
    str3 = unique(str3)
    if(length(str1)>0){
      str1 = gsub("\\(.*\\)","",str1)
      str1 = gsub("\\.","",str1)
      str1 = gsub(" ","",str1)
      if(tf_name%in%str1){
        for(k in 1:length(str3)){
          n = which(ntf_gene%in%str3[k])
          if(length(n)!=0){
            if(O[index,n]==1){
              Regulons[index,n]=1
              num = num + 1
            }
          }
        }
      }
    }
    if(length(str2)>0){
      str2 = gsub("\\(.*\\)","",str2)
      str2 = gsub("\\.","",str2)
      str2 = gsub(" ","",str2)
      if(tf_name%in%str2){
        for(k in 1:length(str3)){
          n = which(ntf_gene%in%str3[k])
          if(length(n)!=0){
            if(O[index,n]==1){
              Regulons[index,n]=1
              num = num + 1
            }
          }
        }
      }
    }
  }
}



num = 1 
copaired2 <- data.frame()
for (i in 1:dim(Regulons)[1]) {
  if(sum(Regulons[i,])>=3){
    for(j in 1:dim(Regulons)[2]){
      if(Regulons[i,j]!=0){
        copaired2[num,1] = rownames(Regulons)[i]
        copaired2[num,2] = colnames(Regulons)[j]
        num = num +1
      }
    }
  }
  
}
genesets2 = Genesets(copaired2)
genesets_list2 <- levels(as.factor(copaired2$V1))

data_auc <- as.data.frame(matrix(0,length(genesets_list2),dim(data_count)[2]))
library(AUCell)
cells_rankings <- AUCell_buildRankings(as.matrix(data_count))
for(i in 1:length(genesets_list2)){
  genes <- unlist(strsplit(genesets2[i,1],split = ","))
  geneSets <- list(geneSet1=genes)
  geneSets <- GSEABase::GeneSet(genes, setName=genesets_list2[i]) # alternative
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
  par(mfrow=c(1,1))
  set.seed(123)
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE))
data_auc[i,] <- as.numeric(getAUC(cells_AUC)) 

}
colnames(data_auc) = colnames(data_count)
rownames(data_auc) = genesets_list2

write.csv(data_auc,"~/data_auc.csv")




