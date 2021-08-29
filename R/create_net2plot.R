#' This function create an igraph network object for visualization puroposes
#'
#' @param g_net matrix, in each row there are starting and ending nodes of one edge 
#'  in the matrix. This matrix represent a gene-gene or gene-fragment network.
#' @param ff_net g_net matrix, in each row there are starting and ending nodes of one edge 
#'  in the matrix. This matrix represent a fragment-fragment network.
#' @param input_m vector of dimension nX1 representing the transcriptional profile befor propagation
#' @param gf_prop vector of dimension nX1 representing the transcriptional profile after 
#'  propagation using a gene-gene or gene-fragment network
#' @param ff_prop vector of dimension nX1 representing the transcriptional profile after
#'  propagation using a gene-fragment network
#' @param ann_net_b matrix, for each row presents the gene identifier, the chromosome in which the gene is,
#'  the starting end ending position in the sequence.
#' @param frag_pattern string, initial character of the fragments name
#'
#' @return igraph object
#' @export
#' 
#' @import igraph
#' 

create_net2plot = function(g_net, input_m, gf_prop, ann_net_b, frag_pattern="F",
                          ff_net = NULL,ff_prop = NULL){
  #Create the joint network
  if(is.null(ff_net)){
    data_net = g_net
  }
  else{
    data_net=rbind(g_net, ff_net)
    frag_pattern=paste("^",frag_pattern,sep="")
  }
  
  #Format and filter from duplicated rows the network edgelist 
  colnames(data_net)<-c("V1","V2")
  df=as.data.frame(data_net)
  df=df[!duplicated(t(apply(df,1,sort))),]
  net_m=as.matrix(df);rm(df)
  
  #Create igraph obj
  net1=igraph::graph_from_edgelist(net_m,directed = F)
  namesvector<-igraph::get.vertex.attribute(net1,"name")
  
  #Merge the two propagations
  if(is.null(ff_net)){
    merged_prop = gf_prop
  }
  else{
    check=!(rownames(gf_prop) %in% rownames(ff_prop))
    merged_prop=rbind(gf_prop[check,],ff_prop)
    rownames(merged_prop)[1:(sum(check))]=rownames(gf_prop)[check]
  }
  
  
  #Create propagation vector for network nodes: df with name | value : str and int
  prop=data.frame(name=rownames(merged_prop),value=merged_prop[,1],stringsAsFactors = F)
  rownames(prop)=seq(1,nrow(prop))
  miss_nodes=setdiff(as.vector(data_net),prop$name)
  if(length(miss_nodes)!=0){
    prop2=data.frame(name=miss_nodes,value=min(prop$value))
    prop=rbind(prop,prop2)
  }
  prop$value=round(prop$value,digits = 2)
  
  #Create starting vector for network genes: df with name | expr: str and others "-"
  starting_prop=data.frame(name=rownames(input_m),expr=input_m[,1])
  if(length(setdiff(as.vector(data_net),starting_prop$name)) > 0){
    starting_prop2=data.frame(name=setdiff(as.vector(data_net),starting_prop$name),expr="-")
    starting_prop=rbind(starting_prop,starting_prop2)
  }
  rownames(starting_prop)=seq(1,nrow(starting_prop))
  
  #Map of the edge colors
  edge_group=rep("orange",nrow(net_m))
  if(!is.null(ff_net)){
    gene_indxs=grep(frag_pattern,V(net1)$name,invert = T)
  }else {
     gene_indxs = 1:length(V(net1))
  }
  net1=igraph::set.vertex.attribute(net1,name="shape",value="circle",index = igraph::V(net1)[gene_indxs])
  net1=igraph::set.edge.attribute(net1, name="color", value=edge_group)
  if(!is.null(ff_net)){
    edge_group[grep(frag_pattern,net_m[,1])] = "purple" 
    edge_group[grep(frag_pattern,net_m[,2])] = "purple"
    net1=igraph::set.edge.attribute(net1, name="color", value=edge_group)
    #set genes as circles and fragments as squares and edges' colors
    frag_indxs=grep(frag_pattern,V(net1)$name,invert = F)
    net1=igraph::set.vertex.attribute(net1,name="shape",value="square",index = igraph::V(net1)[frag_indxs])
  }
  
  #assignment to network
  color_info=combine_netWprop(net1, prop, starting_prop)
  V(net1)$color=color_info[,3]
  V(net1)$propagation=color_info[,2]
  V(net1)$expr=color_info[,4]
  
  net1=toolkitf(net1,ann_net_b,frag_pattern,ff_net)
  
  return(net1)
}