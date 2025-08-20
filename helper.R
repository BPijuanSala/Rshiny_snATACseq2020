########################################################################################
## Title: helper.R
## Author: Blanca Pijuan-Sala
## Description: Shiny app for snATAC-seq E8.25 embryos
## Date: 16 December 2019
########################################################################################
#setwd("/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/Experiments/PhD_BPS53/scripts/website_production/chromatin_early_organogenesis_split/")
# helper.R


#library(data.table)

ac = c(
  "1"="open",
  "0"="closed"
)

bluePal <- c("#BFBFBF","#6495ED","#000000")

all_colours = c(
  "Erythroid1" =  "#C72228",#[15] "Erythroid 1"                                 
  "Erythroid2" =  "#EF4E22",#[37] "Erythroid2"     
  "Erythroid3" =  "#f77b59",
  
  
  "Allantois" = "#532C8A",
  "Cardiomyocytes" =  "#B51D8D",  
  "Mixed mesoderm" =  "#ab80b7",#[26] E8.25 mixed mesoderm   
  "Endothelium" =  "#ff891c",                             
  "Erythroid" =  "#EF4E22",
  "ExE endoderm" = "#7F6874",                              
  "ExE mesoderm" =  "#8870ad",
  "Forebrain" = "#647a4f",
  "Mid/Hindbrain"="#8ca376",
  "Gut" =  "#EF5A9D",                                
  "Neural crest"= "#C3C388",
  "NMP" =  "#8EC792",                                        
  "Notochord" =  "#0F4A9C",                             
  "Paraxial mesoderm" =  "#8DB5CE",
  "Mesenchyme" = "#cc7818",
  "Somitic mesoderm" =  "#005579",                                   
  "Spinal cord" =  "#CDE088",                              
  "Surface ectoderm" = "#f7f79e",      
  "Pharyngeal mesoderm" =  "#C9EBFB"               
  
  
  
)
#all_colours = c(
 # "Allantois" = "#532C8A",
 # "Cardiomyocytes" =  "#B51D8D",  
#  "Mixed mesoderm" =  "#ab80b7",#[26] E8.25 mixed mesoderm   
#  "Endothelium" =  "#ff891c",                             
#  "Erythroid" =  "#EF4E22",
#  "ExE endoderm" = "#7F6874",                              
#  "ExE mesoderm" =  "#8870ad",
#  "Forebrain" = "#647a4f",
#  "Mid/Hindbrain"="#8ca376",
#  "Gut" =  "#EF5A9D",                                
#  "Neural crest"= "#C3C388",
#  "NMP" =  "#8EC792",                                        
#  "Notochord" =  "#0F4A9C",                             
#  "Paraxial mesoderm" =  "#8DB5CE",
#  "Mesenchyme" = "#cc7818",
#  "Somitic mesoderm" =  "#005579",                                   
#  "Spinal cord" =  "#CDE088",                              
#  "Surface ectoderm" = "#f7f79e",      
#  "Pharyngeal mesoderm" =  "#C9EBFB"               
  
  
  
#)


access_col =  c("closed"="#BFBFBF",
                "open"="#000000")
clust_endo_colours = c(
  "EC1"="#f9a602",#Endothelium
  "Ery1"="#8d021f",#Ery
  "Ery2"="#933a16",
  #"2"="#b43757",#Ery
  "Ery3"="#5e1914",#Ery
  "Al_EC"="#c06c84",#transition AL
  "Ery4"="#260805",#Ery
  "EC2"="#ef820d",#HE
  "Al1"="#81007f",#AL
  
  "Ery5"="#bf0a30",#Ery
  "Al2"="#7852a9",#AL
  "Al3"="#b200ed",#AL
  
  "Ery6"="#ea3c53",#Ery
  "Ery7"="#a45a52",#Ery
  "Al4"="#6f2da8",#AL
  "Ery8"= "#fa8072",#Ery
  "EC_Haem"="#ff0800",#transition blood
  "Ery9"="#ed2939"
  
)

#link_bin = HDF5Array("./www/data/snATAC_bin.hdf5", "snATAC_bin")
link_10X = HDF5Array("./www/data/counts10X.hdf5","counts10X")
link_chromVAR= HDF5Array("./www/data/chromVAR.hdf5", "chromVAR")
#link_10X = HDF5Array("./www/data/RNA.hdf5","RNA")
#load(file="./www/data/RNA_int.rda")
#loading data
#load("./www/data/dataTable_genomicRegions.rda")





load("www/data/genes10X.rda")
load("./www/data/meta_10X_celltype.rda")
#load(file="./www/data/RNA.rda")
#sort( sapply(ls(),function(x){object_size(get(x))})) 


#================
#rda
#================
#load("./www/data/meta.RData")

#load("./www/data/coords_snATAC_bin.rda")
#load("./www/data/TFs_chromVAR.rda")
#load("./www/data/TFs_chromVAR_abbrev_up.rda")
#load("./www/data/TFs_chromVAR_gene.rda")
#load("./www/data/genomic_regions.RData")
#load("./www/data/snATAC_coord_indices.rda")


#================
#feather
#================
library(feather)
meta <- as.data.frame(feather::read_feather("./www/data/meta.feather"))

chromVAR_meta <- as.data.frame(feather::read_feather("./www/data/chromVAR_meta.feather"))
TFs_chromVAR=chromVAR_meta[,1]
TFs_chromVAR_abbrev_up=chromVAR_meta[,2]
TFs_chromVAR_gene=chromVAR_meta[,3]

genomic_regions <- as.data.frame(feather::read_feather("./www/data/genomic_regions.feather"))

snATAC_coord_indices_df <- as.data.frame(feather::read_feather("./www/data/snATAC_coord_indices.feather"))
snATAC_coord_indices =snATAC_coord_indices_df[,2]
names(snATAC_coord_indices) =snATAC_coord_indices_df[,1]



plotLevels <- function(x, y, gene, cols=c("#BFBFBF","#6495ED","#000000"),
                       xlab="x",ylab="y",label_legend=expression('log'[10]*' normalized counts + 1'),titlePlot=gene,cexType=1,ybsub=0.1){
  
  redRamp <- colorRampPalette(cols[1:length(cols)])
  df <- data.frame(x = x, y = y, exp = gene)
  df <- df[order(df$exp,decreasing=F),]
  num_fragments = 20
  interval <- findInterval(df$exp, seq(min(df$exp), 
                                       max(df$exp), 
                                       (max(df$exp)-min(df$exp))/num_fragments))
  
  interval[interval==0]<-1
  colorsPlot <- redRamp((num_fragments+1))[interval]
  #bottom, left, top and right
  par(mar=c(8,2,4,2),xpd=NA)
  plot(df$x, df$y, col=colorsPlot, pch=20, cex=cexType,
       xlab="", ylab="", main=titlePlot, axes=F, cex.main=1.5)
  
  
  xl <- min(df$x)- (min(df$x)*(-5.57*10^(-4)))
  yb <- (min(df$y))-(-0.17*min(df$y))
  xr <- max(df$x)
  yt <- (min(df$y))-(min(df$y)*(-0.07))
  rect(xleft=head(seq(xl,xr,(xr-xl)/num_fragments),-1), ybottom = yb, 
       xright = tail(seq(xl,xr,(xr-xl)/num_fragments),-1), ytop = yt,
       col=redRamp(num_fragments), border=redRamp(num_fragments), lwd=0.5, xpd=NA)
  rect(xl,yb,xr,yt, xpd=NA)
  segments(seq(xl,xr,((xr-xl)/5)),yb,seq(xl,xr,((xr-xl)/5)),yb-(-yb*0.0381), xpd=NA)
  text(seq(xl,xr,((xr-xl)/5)), yb-(-yb*0.089), labels = round(seq(min(df$exp), max(df$exp), (max(df$exp)-min(df$exp))/5),3), cex=1.2, xpd=NA)
  text(stats::median(seq(xl,xr,((xr-xl)/5))), yb-(-yb*0.216), labels = label_legend,cex=1.4, xpd=NA)
  
  
  
}




plotTFLevels <- function(x, y, gene, cols=c("#BFBFBF",#mid
                                            "#de1404",#high
                                            "#0664cf"),#low
                         xlab="",ylab="",
                         label_legend=expression('log'[10]*' normalized counts + 1'),
                         titlePlot="",cexType=1,ybsub=0.1){
  
  redRamp <- colorRampPalette(cols[c(1,2)])
  blueRamp <- colorRampPalette(cols[c(1,3)])
  
  df <- data.frame(x = x, y = y, exp = gene)
  df <- df[order(df$exp,decreasing=F),]
  absolute = max(abs(df$exp))
  
  dfsub <- df[df$exp>=0,]
  num_fragments = 20
  interval <- findInterval(dfsub$exp, seq(0, 
                                          absolute, 
                                          (absolute-0)/num_fragments))
  
  interval[interval==0]<-1
  colorsPlot <- redRamp((num_fragments+1))[interval]
  df$cols = rep(0,nrow(df))
  df[df$exp>=0,"cols"] = as.character(colorsPlot)
  
  df <- df[order(df$exp,decreasing=T),]
  dfsub <- df[df$exp<=0,]
  dfsub$exp <- abs(dfsub$exp)
  num_fragments = 20
  interval <- findInterval(dfsub$exp, seq(0, 
                                          absolute, 
                                          (absolute-0)/num_fragments))
  
  interval[interval==0]<-1
  colorsPlot <- blueRamp((num_fragments+1))[interval]
  
  df[df$exp<=0,"cols"] = as.character(colorsPlot)
  df <- df[order(abs(df$exp),decreasing=F),]
  
  
  #bottom, left, top and right
  par(mar=c(8,2,2,2),xpd=NA)
  plot(df$x, df$y, col=as.character(df$cols), pch=20, cex=cexType,
       xlab="", ylab="", main=titlePlot, axes=F, cex.main=1.5)
  
  df <- df[order((df$exp),decreasing=F),]
  num_sec=length(unique(as.character(df$cols)))
  xl <- min(df$x)- (min(df$x)*(-5.57*10^(-4)))
  yb <- (min(df$y))-(-0.17*min(df$y))
  xr <- max(df$x)
  yt <- (min(df$y))-(min(df$y)*(-0.07))
  rect(xleft=head(seq(xl,xr,(xr-xl)/num_sec),-1), ybottom = yb, 
       xright = tail(seq(xl,xr,(xr-xl)/num_sec),-1), ytop = yt,
       col=unique(as.character(df$cols)), border=unique(as.character(df$cols)), lwd=0.5, xpd=NA)
  rect(xl,yb,xr,yt, xpd=NA)
  segments(seq(xl,xr,((xr-xl)/5)),yb,seq(xl,xr,((xr-xl)/5)),yb-(-yb*0.0381), xpd=NA)
  text(seq(xl,xr,((xr-xl)/5)), yb-(-yb*0.089), labels = round(seq(min(df$exp), max(df$exp), (max(df$exp)-min(df$exp))/5),3), cex=1, xpd=NA)
  text(stats::median(seq(xl,xr,((xr-xl)/5))), yb-(-yb*0.216), labels = label_legend,cex=1.2, xpd=NA)
  
  
  
}

