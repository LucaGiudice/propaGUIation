#' This function downloads the genes annotations for human and mouse
#'
#' @param organism ,string representing the orgnanism for which downloading the
#'  genes annotations
#' @param gene_names, vector of strings, gene names to search annotation for
#'
#' @return list with the gene id, chromosome, start and end positions.
#' @export
#'
#' @importFrom biomaRt useEnsembl getBM
#'

download_genes_annotation = function(organism = "human",gene_names){
  if(organism != "human" && organism != "mouse"){
    cat("Error: you are asking for an unsupported organism")
    return(NA)
  }
  if(length(gene_names) == 0){
    cat("Error: no gene names provided")
    return(NA)
  }
  dataset=""
  if(organism == "human"){
    dataset = "hsapiens_gene_ensembl"
    gene_id = "hgnc_symbol"
  }
  if(organism == "mouse"){
    dataset = "mmusculus_gene_ensembl"
    gene_id = "mgi_symbol"
  }

  ensembl <- useEnsembl(biomart = "ensembl", dataset = dataset)
  a = getBM(attributes=c(gene_id,"chromosome_name","start_position","end_position"),
    filters=c("chromosome_name",gene_id),values=list(c(1:19,"X","Y"),gene_names),mart=ensembl,useCache = FALSE)
  a = a[a[,1] != "",]
  a[,2] = paste("chr",a[,2],sep="")
  colnames(a) = c("ID","chr","start","end")
  differece_length = length(gene_names) - dim(a)[1]
  if(differece_length > 0 ){
    cat(differece_length, " gene(s) annotation(s) has/have not been found\n",sep="")
  }
  return(a)
}
