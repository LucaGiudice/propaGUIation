#' @import igraph
#' 

combine_netWprop = function(net, propag, totaldataframe) { 
  totaldataframe[,2]=as.character(totaldataframe[,2])
  
  dtnames<-data.frame(names = names(V(net)))
  dtinfo=merge(dtnames,propag, by.x="names", by.y="name", sort=F)
  
  dtinfo[,2]=as.double(as.character(dtinfo[,2]))
  miniprop=dtinfo[order(dtinfo$value),]
  
  my_resolution = 5000
  my_palette    = colorRampPalette(c('white','orange1'))
  my_colors = my_palette(my_resolution)[as.numeric(cut(miniprop[,2], breaks=my_resolution))]
  
  miniprop[,3]=my_colors
  colnames(miniprop)<-c("Name","Propagation","Color")
  
  color_info=merge(dtnames,miniprop, by.x="names", by.y="Name", sort=F)
  
  color_info=merge(color_info, totaldataframe, by.x="names", by.y="name", sort=F )
  
  colnames(color_info)<-c("Name","Propagation","Color","Starting expression")
  
  return (color_info)
  
}