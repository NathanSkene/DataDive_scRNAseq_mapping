# Check all the Tasic celltypes are mapped
check_if_all_celltypes_mapped <- function(annotCTs,mappingCTs){
  annotCTs = unique(annotCTs)
  mappingCTs = unique(mappingCTs)
  missingFromAnnot = annotCTs[!annotCTs %in% mappingCTs]
  missingFromMapping = mappingCTs[!mappingCTs %in% annotCTs]
  if(length(missingFromAnnot)>0){
    print(sprintf("These cells are missing from the mapping file: %s",paste(missingFromAnnot,collapse=" \\\ ")))
  }
  if(length(missingFromMapping)>0){
    print(sprintf("These cells are missing from the annotation file: %s",paste(missingFromMapping,collapse=" \\\ ")))
  }
}