# propaGUIation
Graphic user interface to visualize and explore the values of the patient's features (e.g. genes) after a network-based propagation

[![R-CMD-check](https://github.com/jokergoo/cola/workflows/R-CMD-check/badge.svg)](https://github.com/jokergoo/cola/actions)

## Features

1. It allows to explore a patient's/sample's profile (e.g. expression profile, somatic mutation profile, etc ...) after a network-based propagation
2. It allows to investigate the propagation score obtained by specific patient's features (e.g. genes) and their neighbourhood
3. It allows to download the propagated network and to import it in cytoscape

## Install

*propaGUIation* is available on R, you can install it by:

```r
devtools::install_github(
  repo="LucaGiudice/propaGUIation",
  ref = "main",
  dependencies = "Depends",
  upgrade = "always",
  quiet = F
)
```

## Vignettes

There are the following vignettes:

1. [A Quick Start of Using propaGUIation Package](https://github.com/LucaGiudice/supplementary-Simpati/blob/main/images/vignette.pdf)

## Interface

<img src="https://github.com/LucaGiudice/supplementary-Simpati/blob/main/images/gui_overall.png" />

  1. Select nodes by chromosomes: allows to filter the network and to keep only those nodes (e.g. genes) that belong to a specific chromosome
  2. Select nodes by genome region: allows to filter the network and to keep only those nodes that belong to a specific genome region
  3. Select or type node + Distance: allows to visualize the neighbourhood of one specific node of interest
  4. Select by propagation ranges: allows to filter the subnetwork generated with a "search". It allows to keep only the nodes with a value that falls inside a specific range
  5. Scale colours: allows to scale the colors of the visualized nodes as if the propagation would have been applied only on them
  6. More/Less: allow to increase or decrease the size of the neighbourhood around a node of interest based on the distance (e.g. first degree neighbourhood, second degree ...)
  7. Open in cytoscape: allows to open the network in cytoscape
  8. Download: allows to download the image of the network visualized and created with the GUI

### Usage

Few lines of code to run the GUI:

```r
library("propaGUIation")

#Load and set up the example data ----
data("example")

# get the original patient profile
input_m = example$input_gm
# get the propagated profile
gg_prop = example$gg_prop
# get the gene-gene interaction network
g_net = example$gg_net

# download the annotation data for the genes in the propagated profile
g_ann <- download_genes_annotation(organism = "human", gene_names = rownames(gg_prop))

#Create igraph object with all the information included
net=create_net2plot(g_net = g_net, input_m = input_m, gf_prop = gg_prop, ann_net_b = g_ann)

#Start GUI
start_GUI(net, g_ann, example$chr_len, example=T)
```

### Plots

Following plot shows you an example of how to interact with the GUI and its functionalities

<img src="https://github.com/LucaGiudice/supplementary-Simpati/blob/main/images/gui_scene.gif" />

## License

MIT @ Giudice Luca
