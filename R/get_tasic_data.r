get_tasic_data <- function(){
  # EXPRESSION DATA
  linURL = "https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/tasic/genes_rpkm.csv"
  path   = "genes_rpkm.csv"
  download.file(linURL, destfile=path)
  rpkm = read.csv(path,stringsAsFactors=FALSE,sep=",")
  rpkm = rpkm[!duplicated(rpkm$X),]
  rownames(rpkm)=rpkm$X
  rpkm = rpkm[,colnames(rpkm)!="X"]
  rpkm2 = Matrix::Matrix(data.matrix(rpkm))
  file.remove(path)
  
  # CELL META
  linURL = "https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/tasic/cell_metadata.csv"
  path   = "cell_metadata.csv"
  download.file(linURL, destfile=path)
  cell_meta = read.csv(path,stringsAsFactors=FALSE,sep=",")[,c("long_name","short_name")]
  file.remove(path)
  
  # CELL ANNOTATION
  linURL = "https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/tasic/cell_classification.csv"
  path   = "cell_classification.csv"
  download.file(linURL, destfile=path)
  cell_annot = read.csv(path,stringsAsFactors=FALSE,sep=",")[,c("X","primary")]
  colnames(cell_annot) = c("short_name","cluster_id")
  cell_annot2 = merge(cell_annot,cell_meta,by="short_name")
  file.remove(path)
  
  # CLUSTER ANNOTATION
  linURL = "https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/tasic/cluster_metadata.csv"
  path   = "cluster_metadata.csv"
  download.file(linURL, destfile=path)
  cluster_annot = read.csv(path,stringsAsFactors=FALSE,sep=",")
  all_annot = merge(cell_annot2,cluster_annot,by="cluster_id")
  file.remove(path)
  
  # ALIGN THE ANNOTATIONS AGAINST EXPRESSION
  rownames(all_annot) = all_annot$long_name
  all_annot=all_annot[colnames(rpkm),]
  
  output = list()
  output$exp = rpkm
  output$annot = all_annot
  return(output)
}