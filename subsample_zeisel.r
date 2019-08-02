# module load hdf5/1.10.4
# wget https://storage.googleapis.com/linnarsson-lab-loom/l5_all.loom

# Data from: http://mousebrain.org/downloads.html
library(tidyverse)
library("rhdf5")
library("snow")
library(loomR)

path = "l5_all.loom"

lfile <- connect(filename = path, mode = "r+")
Lvl5=lfile$col.attrs$ClusterName[]
genes=lfile$row.attrs$Gene[]
cts = unique(Lvl5)

get_exp <- function(ct,lfile){
  library(tidyverse)
  library("rhdf5")
  library("snow")  
  library(loomR)
  library(Matrix)
  
  #lfile <- connect(filename = path, mode = "r+")
  Lvl5=lfile$col.attrs$ClusterName[]
  whichCT = Lvl5==ct
  exp = Matrix::Matrix(t(lfile[["matrix"]][whichCT,]))
  return(exp) 
}
allExp = lapply(as.list(cts),get_exp,lfile=lfile)

subsample_ct <- function(expIN,N = 30){
  library(tidyverse)
  library("rhdf5")
  library("snow")  
  library(loomR)
  
  print(dim(expIN))
  
  numGenes = dim(expIN)[1]
  numCells = dim(expIN)[2]
  #sampleGroups = split(1:numCells, sample(1:N, numCells, replace=T))
  newE = matrix(0,nrow=numGenes,ncol=30)
  for(i in 1:N){
    print(i)
    # Extract only cells from that subgroup
    #exp = expIN[,sampleGroups[[i]]]
    minSize = round(dim(expIN)[2]/N)
    if(minSize<3){minSize=3}
    whichI = sample(1:dim(expIN)[2],minSize)
    #exp = expIN[,sampleGroups[[i]]]
    exp = expIN[,whichI]
    
    #if(class(exp)=="numeric"){
    #  exp = t(Matrix::Matrix(rbind(exp,exp)))
    #}
    
    # Normalise the reads per cell
    exp2 = Matrix::t(Matrix::t(exp)/Matrix::colSums(exp))
    
    e = Matrix::rowSums(exp2)
    e = e/ dim(exp)[2]
    newE[,i] = e
  }
  #colnames(newE) = rep(ct,dim(newE)[2])
  return(newE)
}

library(clustermq)

options(
  clustermq.scheduler = "slurm",
  clustermq.template = "/nas/longleaf/home/nskene/.SLURMtemplate"
)

output = Q(subsample_ct, expIN=allExp,n_jobs=10, template=list(log_file="/nas/longleaf/home/nskene/ClusterMQ_logs/myLog.%a.log"),memory=2000)

for(i in 1:length(allExp)){colnames(output[[i]])=rep(cts[i],dim(output[[i]])[2])}

subsampledExp = output[[1]]
for(i in 2:length(allExp)){subsampledExp = cbind(subsampledExp,output[[i]])}
rownames(subsampledExp) = genes
save(subsampledExp,file="zeisel_subSampled.rda")

library(umap)
zeisel2018_pca = prcomp(t(subsampledExp))
plot(zeisel2018_pca$x[,1],zeisel2018_pca$x[,2])
aaa=umap(subsampledExp[1:20,])

