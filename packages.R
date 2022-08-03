#List
# BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

## Bioconductor packages
### Deseq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# Ordinary Orphans
## minpack
install.packages('minpack.lm', dependencies = T)
## biogrowth
install.packages('biogrowth', dependencies = T)
## fitdistrplus
install.packages('fitdistrplus', dependencies=T)
## logspline
install.packages('logspline', dependencies = T)
## mixtools
install.packages('mixtools', dependencies = T)
## data.tables
install.packages("data.table")
## install nlme
install.packages("nlme")
# ramify
install.packages('ramify')
#stringr
install.packages("stringr")
#dpylr
install.packages(c('dplyr', 'readxl'))
#inflection
install.packages('inflection')
#samSpectral
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SamSPECTRAL")

#growthcurver
install.packages("growthcurver")
#zoo
install.packages('zoo')
#slider
install.packages("slider")
#runner
install.packages('runner')
#akmedoids
install.packages('akmedoids')
#SciViews
install.packages('SciViews')
#baseline
install.packages('baseline')
#latex2exp
install.packages('latex2exp')
#ggplots
install.packages("ggplot2")
#ggridges
install.packages("ggridges")
#ggpubr
install.packages("ggpubr")
#forcats
install.packages("forcats")
#gmcm
install.packages("GMCM")
#complex heatmaps
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
#d3heatmap
install.packages('d3heatmap')
#text_repel
install.packages("ggrepel")
#gridExtra
install.packages("gridExtra")
#RColorBrewer
install.packages("RColorBrewer")
#cartocolor
install.packages("rcartocolor")
library(rcartocolor)  # package 
display_carto_all(colorblind_friendly = TRUE)  # show colours
#cartocolor
install.packages("factoextra")
#raincloudplots
install.packages('raincloudplots')
#staplr
install.packages("staplr")
#heatmaply
install.packages('heatmaply')
#pROC
install.packages('pROC')
#ROCR
install.packages('ROCR')
#wesanderson
install.packages("wesanderson")
#ggpmisc
install.packages('ggpmisc')
#shiny
install.packages("shiny")
#pracma
install.packages("pracma")
#proda
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("proDA")
#stabperf
BiocManager::install("StabPerf")
#ggbeeswarm
install.packages('ggbeeswarm')
#markdown
install.packages('markdown')
#DT
install.packages("DT")
#sendmailR
install.packages('sendmailR')
#beeswarm
install.packages("beeswarm")
#beeswarm
install.packages("tidyverse")
#modelsummary
install.packages('modelsummary')
#naijR
install.packages('naijR')
#rsconnect
install.packages('rsconnect')
#bioassay
install.packages('bioassays')
#phenoscreen
if (!require(devtools)) install.packages("devtools")
devtools::install_github('Swarchal/phenoScreen')
#ggplot2bdc
if(!require("devtools")) install.packages("devtools")
devtools::install_github("briandconnelly/ggplot2bdc")
#ggcorrplot
install.packages("ggcorrplot")
#venn diagram
install.packages('VennDiagram')
install.packages("ggVennDiagram")
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
#maps
install.packages("ggiraph")
install.packages('spData')
install.packages('tmap')
install.packages("leaflet")
