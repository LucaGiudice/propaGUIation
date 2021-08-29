#' @import igraph
#' 

#function to connect network nodes to metainformation
toolkitf <- function(net1, ann_net_b, frag_pattern = "F",ff_net = NULL) {
  toolnames=V(net1)$name
  if(is.null(ff_net)){
    symbol_position = 1:length(V(net1))
  }else {
    symbol_position = grep(frag_pattern,toolnames, invert=T)
  }
  symbol_names=toolnames[symbol_position]
  
  if(is.null(ff_net)){
    frag_position = 0
  }else {
    frag_position = grep(frag_pattern,toolnames, invert = F)
  }
  frag_names=data.frame(names=toolnames[frag_position])
  
  data_annotation=merge(frag_names, ann_net_b, by.x="names", by.y="ID", sort=F)

  #correct indices for annotation
  ann_frag_position=which(toolnames %in% data_annotation[,1]) 
  frag_link=paste("<a target='_blank' href=https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
    data_annotation$chr,"%3A",data_annotation$start,"%2D",data_annotation$end,
    "&hgsid=981993861_VmexClFuvk6xFB6uwwbxkArgAkC9> Annotation of ", 
    data_annotation$names,"</a>",sep="")
  vector_link=paste("<a target='_blank' href=https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
    sep = "",symbol_names, "> GeneCards of ", symbol_names,"</a>")
  
  nolink_vector=c(rep("No link available", length(toolnames)))
  toolkit_title=replace(nolink_vector,ann_frag_position,frag_link)
  toolkit_title=replace(toolkit_title,symbol_position,vector_link)
  
  V(net1)$title=toolkit_title
  
  return(net1)
}