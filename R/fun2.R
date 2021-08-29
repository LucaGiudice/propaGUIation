#' function to assign the relative colors 
#' @param subnet igraph object
#' @importFrom grDevices colorRampPalette
#' 
fun2 <- function(subnet) {
  dtinfo=get.data.frame(subnet, what="vertices")
  miniprop=dtinfo[order(dtinfo$prop),]
  my_resolution = 5000
  my_palette = colorRampPalette(c('white','orange1'))
  my_colors = my_palette(my_resolution)[as.numeric(cut(miniprop[,4], breaks=my_resolution))]
  miniprop[,6]<-my_colors
  
  dtnames<-data.frame(names = names(V(subnet)))
  colordata=merge(dtnames, miniprop, by.x="names", by.y="name", sort=F)
  
  
  V(subnet)$color_rel<-colordata[,6]
  
  return(subnet)
}