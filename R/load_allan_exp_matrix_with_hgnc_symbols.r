load_allan_exp_matrix_with_hgnc_symbols <- function(path="human_MTG_2018-06-14_intron-matrix.csv",prelim="human_MTG_2018-06-14"){
  # Rownames of the exp data are given as entrez_id... so first load the gene data
  geneAnnot = fread(sprintf("%s_genes-rows.csv",prelim))
  sampleAnnot = fread(sprintf("%s_samples-columns.csv",prelim))
  
  if(!"gene_symbol" %in% colnames(geneAnnot)){
    if("gene" %in% colnames(geneAnnot)){
      geneAnnot = geneAnnot %>% dplyr::rename(gene_symbol=gene)
    }else{
      stop("Cannot find gene_symbol column in gene annotation file")
    }
  }
  
  # Then load the expression data
  intronExp = fread(path)
  newDataset = list()
  newDataset$exp = Matrix::Matrix(data.matrix(intronExp))
  rownames(newDataset$exp) = geneAnnot$gene_symbol
  newDataset$exp = newDataset$exp[,-1]
  
  # Then prep the annotation matrix
  newDataset$annot = sampleAnnot %>% dplyr::rename(cellID=sample_name,celltype=cluster)  %>% dplyr::select(cellID,celltype)
  return(newDataset)
}