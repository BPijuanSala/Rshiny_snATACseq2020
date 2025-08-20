########################################################################################
## Title: server.R
## Author: Blanca Pijuan-Sala
## Description: Shiny app for snATAC-seq E8.25 embryos
## Date: 16 December 2019
########################################################################################



library(shinycssloaders)
library(shiny)
library(DT)
library(png)
library(Matrix)
library(data.table)
library(ggplot2)
library(HDF5Array)

source("helper.R",local=TRUE)


server <- function(input, output) {

  #########==================================
  # Accessibility tab
  #########==================================
  dataAccess <- reactive({
    
    idxGene <- unname(snATAC_coord_indices[input$OCR])
    coords_new = names(snATAC_coord_indices)[snATAC_coord_indices==idxGene]
    
    link_bin = HDF5Array::HDF5Array(paste0("./www/data/snATAC_bin_file",idxGene,".hdf5"), paste0("snATAC_bin_",idxGene))
    
    
    
    value = as.numeric(link_bin[match(as.character(input$OCR), 
                                                      as.character(coords_new)),])
    return(value)
  })
 
  plotAccess<-reactive({

   # idxGene <- unname(snATAC_coord_indices[input$OCR])
   # coords_new = names(snATAC_coord_indices)[snATAC_coord_indices==idxGene]
    
   # link_bin = HDF5Array(paste0("./www/data/snATAC_bin_file",idxGene,".hdf5"), paste0("snATAC_bin_",idxGene))
    
    
    
    df = data.frame(x = meta[,"umap_X"],
                    y=meta[,"umap_Y"],
                    value = dataAccess()
                    )

                      
    df$value2 = ac[as.character(df$value)]                 
    df = df[order(df$value,decreasing=F),]
    
    
    plot(df$x,df$y,axes=F,ylab="",xlab="",
             col=access_col[as.character(df$value2)],pch=20,
             main=input$OCR,cex.main=1.5)
    legend("bottomright", inset=.02, legend=c("open", "closed"),pch=16,
           col=c("#000000", "#BFBFBF"), cex=1.2)
   
  })
  
  
  barplotAccess = reactive({
    
   # idxGene <- unname(snATAC_coord_indices[input$OCR])
   # coords_new = names(snATAC_coord_indices)[snATAC_coord_indices==idxGene]
    
   # link_bin = HDF5Array(paste0("./www/data/snATAC_bin_file",idxGene,".hdf5"), paste0("snATAC_bin_",idxGene))
    
    
   # value = as.numeric(link_bin[match(as.character(input$OCR), 
   #                                   as.character(coords_new)),])
    value=dataAccess()
    tab0 = table(value,as.character(meta$ann))["1",]
    tab = tab0*100/table(meta$ann)
    df = data.frame(
      ann = as.vector(names(tab)),
      values = as.vector(unname(tab)),
      col = as.character(all_colours[as.character(names(tab))])
    )
    
    b = ggplot(data=df, aes(x=ann, y=values)) +
      geom_bar(stat="identity", width=0.5,fill= df$col) +
      theme_bw() +
      labs(y = "Nuclei with accessibility (%)") + 
      ggtitle(input$OCR) +
      theme(axis.title = element_text(face = "bold", size = 14),
            axis.text.y = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(size = 12, face = "bold", 
                                       angle = 45, hjust = 1),
            legend.position = "none",
            axis.title.x = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank()
            #panel.grid.minor = element_blank()
            ) +
      annotate("text", 
               x = sort(names(tab)), 
               y = rep_len(c(max(tab)*1.1, max(tab) * 1.2), 
                           length.out = length(tab)), 
               label = tab0[sort(names(tab))]) 
    
    return(b)
    
  })
  
  
  output$Access_plot <- renderPlot({
    
    validate(
      need(input$OCR %in% names(snATAC_coord_indices),
           "The coordinate is not valid. Please check the 'Genomic regions' tab to explore valid coordinates and enter one in the form of e.g. chr1_3008841_3009342." )
    )
    plotAccess()
    
  })
  
  output$Access_barplot <- renderPlot({
    
    validate(
      need(input$OCR %in% names(snATAC_coord_indices),
           "The coordinate is not valid. Please check the 'Genomic regions' tab to explore valid coordinates and enter one in the form of e.g. chr1_3008841_3009342." )
    )
    barplotAccess()
    
  })
  

output$downloadAccess <- downloadHandler(
    filename = function() { paste0(input$OCR, "_UMAP.pdf") },
    content = function(file) {
      pdf(file, width=7,height=8)
      
      df = data.frame(x = meta[,"umap_X"],
                      y=meta[,"umap_Y"],
                      value = dataAccess()
      )
      
      
      df$value2 = ac[as.character(df$value)]                 
      df = df[order(df$value,decreasing=F),]
      
      
      plot(df$x,df$y,axes=F,ylab="",xlab="",
           col=access_col[as.character(df$value2)],pch=20,
           main=input$OCR,cex.main=1.5)
      legend("bottomright", inset=.02, legend=c("open", "closed"),pch=16,
             col=c("#000000", "#BFBFBF"), cex=1.2)
      dev.off()
        }
  )

  
  output$download_barplotAccess <- downloadHandler(
    filename = function() { paste0(input$OCR, "_barplot.pdf") },
    content = function(file) {
      pdf(file, width=10,height=6)
      print(barplotAccess())
      dev.off()
    }
  )
  
  
  
  #########==================================
  # Genomic regions tab
  #########==================================

  
  
  output$OCR_DT <- renderDataTable({
   # dt
    DT::datatable(genomic_regions,
                 filter = list(position = 'top', clear = FALSE),
                
               options = list(
                deferRender = TRUE,
               search = list(regex = TRUE, caseInsensitive = TRUE, search = ''),
              pageLength = 10
           ),rownames= FALSE)    
    
  })

  
  # Downloadable csv of selected dataset ----
  output$downloadOCR <- downloadHandler(
    filename = function() {"OCR_table.csv"},
    content = function(file) {
      write.csv(genomic_regions[input$OCR_DT_rows_all,], file, row.names = FALSE)
    }
  )
  
  
  #########==================================
  # Topics tab
  #########==================================
  
  plotTopic<-reactive({
    df = data.frame(x = meta[,"umap_X"],
                    y=meta[,"umap_Y"],
                    value=meta[,gsub(" ","_",input$Topic)])
    df = df[order(df$value,decreasing=F),]
    
    
    
    plotLevels(df$x, df$y, df$value, cols=c("#BFBFBF","#6495ED","#000000"),
               xlab="",ylab="",titlePlot=input$Topic,label_legend="Value",
               cexType=1,ybsub=0.1)
    
    
  })
  
  output$Topic_plot <- renderPlot({
    plotTopic()
    
  })
  

  output$topics_motif1 <-renderText({
    
    validate(
      need(gsub(" ","",input$Topic) %in% paste0("Topic",c(2:14,16:100)),
           "This topic did not present regions uniquely contributing to it." )
    )
    known=paste0("./links/Topics_motifs/",gsub(" ","",input$Topic),".html")
    html_pattern_known = paste0(input$Topic," selected: <a href=\"",known,"\" target=\"_blank\">Click here for enriched motifs</a>")
    
    return(html_pattern_known)
  })
  
  
  output$downloadTopic <- downloadHandler(
    filename = function() { paste0(input$Topic, ".pdf") },
    content = function(file) {
      pdf(file, width=7,height=9)
      df = data.frame(x = meta[,"umap_X"],
                      y=meta[,"umap_Y"],
                      value=meta[,gsub(" ","_",input$Topic)])
      df = df[order(df$value,decreasing=F),]
      
      
      
      plotLevels(df$x, df$y, df$value, cols=c("#BFBFBF","#6495ED","#000000"),
                 xlab="",ylab="",titlePlot=input$Topic,label_legend="Value",
                 cexType=1,ybsub=0.1)
      dev.off()
    }
  )
  
  
  
  
  
  
  #########==================================
  # Enriched TFs tab
  #########==================================
  data_TF <-reactive({
    return(as.numeric(link_chromVAR[match(as.character(input$TF), 
                                   as.character((TFs_chromVAR_abbrev_up))),]))
  })
  plotTF <-reactive({
    

    df = data.frame(x = meta$umap_X,
                    y=meta$umap_Y,
                    value=data_TF())
    
                      
   # df = df[order(abs(df$value),decreasing=F),]
    plotTFLevels(df$x, df$y, df$value,
                 xlab="",ylab="",titlePlot=input$TF,
                 cexType=1,ybsub=0.1,label_legend="Z-score")
    
   
    
  })
  
  
  output$chromVAR_plot <- renderPlot({
    
    validate(
      need(input$TF %in% (TFs_chromVAR_abbrev_up),
           "This transcription factor does not exist. Please check the drop-down list." )
    )
    plotTF()
    
  })
  
  output$download_chromVAR <- downloadHandler(
    filename = function() { paste0(input$TF, ".pdf") },
    content = function(file) {
      pdf(file, width=7,height=9)
      
      df = data.frame(x = meta$umap_X,
                      y=meta$umap_Y,
                      value=data_TF())
      
      
      df = df[order(abs(df$value),decreasing=F),]
      
      plotTFLevels(df$x, df$y, df$value,
                   xlab="",ylab="",titlePlot=input$TF,
                   cexType=1,ybsub=0.1,label_legend="Z-score")
      dev.off()
    }
  )
  
  output$selected_chromVAR <- renderText({
    
    validate(
      need(input$TF %in% (TFs_chromVAR_abbrev_up),
           " " )
    )
    d = TFs_chromVAR[match(as.character(input$TF), 
                           as.character((TFs_chromVAR_abbrev_up)))]
    return(paste("<b>ID in motif collection 'mouse_pwms_v2' from 'chromVARmotifs' v0.2.0 R package:</b>",d))
      })
  
  
  #Set images
  output$motif <- renderImage({
    
    validate(
      need(input$TF %in% (TFs_chromVAR_abbrev_up),
           " " )
    )
    d = TFs_chromVAR[match(as.character(input$TF), 
                           as.character((TFs_chromVAR_abbrev_up)))]
    outfile <- tempfile(fileext='.png')
    p <- png::readPNG(paste0("./www/images/logos/",d,".png"))
    png::writePNG(p, target=outfile)
    list(src = outfile, contentType = 'image/png', width = 450, height = 200)
  }, deleteFile = TRUE)
  
  
    
    
    output$download_motif <- downloadHandler(
      
      filename <- function() {
        d = TFs_chromVAR[match(as.character(input$TF), 
                               as.character((TFs_chromVAR_abbrev_up)))]
        return(paste0(d,".png"))
      },
      
      content <- function(file) {
        d = TFs_chromVAR[match(as.character(input$TF), 
                               as.character((TFs_chromVAR_abbrev_up)))]
        p <- png::readPNG(paste0("./www/images/logos/",d,".png"))
        png::writePNG(p, target=file)
        list(src = file, contentType = 'image/png', width = 450, height = 200)
        
        },
      
      contentType = 'image/png'
    )
  
  
    dataTF_RNA <- reactive({
      
      d = TFs_chromVAR_gene[match(as.character(input$TF), 
                                  as.character((TFs_chromVAR_abbrev_up)))]
      
      
      TF=data_TF()
      
      
      y <- aggregate(TF,by=list(meta$ann),FUN=mean)
      TF_mean <- y$x[order(as.character(y$Group.1))]
      names(TF_mean) <- as.character(y$Group.1)[order(as.character(y$Group.1))]
      names(TF_mean) <- as.character(y$Group.1)[order(as.character(y$Group.1))]
      
      RNA= as.numeric(link_10X[match(as.character(d),as.character(genes10X)),])
      
      y <- aggregate(RNA,by=list(meta_10X_celltype),FUN=mean)
      RNA_mean <- y$x[order(as.character(y$Group.1))]
      names(RNA_mean) <- as.character(y$Group.1)[order(as.character(y$Group.1))]
      names(RNA_mean) <- as.character(y$Group.1)[order(as.character(y$Group.1))]                
      
      
      
      
      TF_mean["Erythroid1"] =TF_mean[names(TF_mean)=="Erythroid"]
      TF_mean["Erythroid2"] =TF_mean[names(TF_mean)=="Erythroid"]
      TF_mean["Erythroid3"] =TF_mean[names(TF_mean)=="Erythroid"]
      TF_mean["Intermediate mesoderm"] =TF_mean[names(TF_mean)=="Mixed mesoderm"]
      TF_mean["ExE mesoderm"] =TF_mean[names(TF_mean)=="Mixed mesoderm"]
      
      RNA_mean["Forebrain"] =RNA_mean[names(RNA_mean)=="Forebrain/Midbrain/Hindbrain"]
      RNA_mean["Mid/Hindbrain"] =RNA_mean[names(RNA_mean)=="Forebrain/Midbrain/Hindbrain"]
      
      TF_mean = TF_mean[(names(TF_mean)%in%c("Intermediate mesoderm","ExE mesoderm","Forebrain/Midbrain/Hindbrain"))==F]
      RNA_mean = RNA_mean[(names(RNA_mean)%in%c("Intermediate mesoderm","ExE mesoderm","Forebrain/Midbrain/Hindbrain"))==F]
      
      celltypesOv = intersect(names(TF_mean),names(RNA_mean))
      TF_mean = TF_mean[celltypesOv]
      RNA_mean = RNA_mean[celltypesOv]
      df = data.frame(
        TF = TF_mean,
        RNA = RNA_mean,
        ann = names(TF_mean),
        col = all_colours[celltypesOv]
      )
      #bottom, left, top and right
      df <- df[order(as.character(df$ann)),]
      return(df)
    })
    
    plotTF_RNA <- reactive({
      df = dataTF_RNA()
      par(mar=c(4,5,2,28),xpd=NA)
      
    plot(df$RNA,df$TF,
           pch=20,cex=3,cex.lab=1.4,cex.axis=1.2,cex.main=1.5,#main=gene,
           col=as.character(df$col),
           ylab="Chromatin Accessibility Z-score",
           xlab=expression('log'[2]*'(RNA expression)'),
         main=input$TF
          )
 
      legend("topright", inset=c(-1.85,0),bty="n",
                    legend=unique(df$ann),
                    col=as.character(unique(df$col)),
                    pch=16,ncol=2,cex=1.3,pt.cex=1.5)
      par(xpd = FALSE) #Only because we made TRUE above 
      
      abline(h=0,col="red")
    # p = ggplot(data=df, aes(x=RNA, y=TF)) +
     #   geom_point(aes(col= df$ann),size = 5, 
    #               alpha = 0.9) +
    #   scale_color_manual(values = as.character(df$col)[order(as.character(df$ann))], 
    #                      labels = as.character(df$ann)[order(as.character(df$ann))], 
     #                     drop = FALSE, name = "") +
     #  theme_bw() +
     #  theme (legend.text=element_text(size=rel(1.2)),
     #         axis.title = element_text(face = "bold", size = 14),
     #         axis.text.y = element_text(size = 12, face = "bold"),
      #        axis.text.x = element_text(size = 12, face = "bold"),
      #        plot.title = element_text(size=16)) +
       
       
     #  geom_hline(yintercept = 0,colour="red") +
     #
     #   labs(y = "Chromatin Accessibility Z-score", 
     #        x= "log2(RNA expression)") + 
    #    ggtitle(input$TF) +
     # guides(colour=guide_legend(nrow=10))
      
       
      
     # return(p)
      
      
    })
    
    
    output$TF_RNA <- renderPlot({
      validate(
        need(input$TF %in% (TFs_chromVAR_abbrev_up),
             " " ),
        need(TFs_chromVAR_gene[match(as.character(input$TF), 
                                     as.character((TFs_chromVAR_abbrev_up)))] %in% as.character(genes10X),
             "The corresponding gene cannot be found." )
        )
      
      plotTF_RNA()
      
    })
    
    
    output$download_TF_RNA <- downloadHandler(
      filename = function() { paste0(input$TF, "_TF_RNA.pdf") },
      content = function(file) {
        pdf(file, width=10,height=4.5)
        df = dataTF_RNA()
        par(mar=c(4,5,2,27),xpd=NA)
        
        plot(df$RNA,df$TF,
             pch=20,cex=3,cex.lab=1.4,cex.axis=1.2,cex.main=1.5,#main=gene,
             col=as.character(df$col),
             ylab="Chromatin Accessibility Z-score",
             xlab=expression('log'[2]*'(RNA expression)'),
             main=input$TF
        )
        
        legend("topright", inset=c(-1.5,0),bty="n",
               legend=unique(df$ann),
               col=as.character(unique(df$col)),
               pch=16,ncol=2,cex=1.3,pt.cex=1.5)
        par(xpd = FALSE) #Only because we made TRUE above 
        
        abline(h=0,col="red")
        dev.off()
      }
    )
  
  
  
    
    
    
    
    
    
    #########==================================
    # Accessibility SUBSET tab
    #########==================================
    
    #Set images
   # output$UMAP_subset_celltype <- renderImage({
   #   outfile <- tempfile(fileext='.png')
   #   p <- readPNG("./www/images/UMAP_AL_EC_HAEM_celltypes_website.png")
   #   writePNG(p, target=outfile)
   #   list(src = outfile, contentType = 'image/png', width = 400, height = 350)
   # }, deleteFile = TRUE)
    
    #output$UMAP_subset_celltype_subcolor <- renderImage({
   #   outfile <- tempfile(fileext='.png')
    #  p <- readPNG("./www/images/UMAP_AL_EC_HAEM_website.png")
    #  writePNG(p, target=outfile)
   #   list(src = outfile, contentType = 'image/png', width = 400, height = 350)
    #}, deleteFile = TRUE)
    dataAccess_subset<-reactive({
      
      idxGene <- unname(snATAC_coord_indices[input$OCR_subset])
      coords_new = names(snATAC_coord_indices)[snATAC_coord_indices==idxGene]
      
      link_bin = HDF5Array::HDF5Array(paste0("./www/data/snATAC_bin_subset_file",idxGene,".hdf5"), paste0("snATAC_bin_subset_",idxGene))
      
      return(as.numeric(link_bin[match(as.character(input$OCR_subset), 
                                       as.character(coords_new)),])
      )
      
      
    })
    
    plotAccess_subset<-reactive({
      ac = c(
        "1"="open",
        "0"="closed"
      )
      
      meta_endo=meta[is.na(meta$al_haem_endo_clusters)==F,]
      
      
      
      df = data.frame(x = meta_endo[,"umap_X"],
                      y=meta_endo[,"umap_Y"],
                      value = dataAccess_subset())
      
      
      df$value2 = ac[as.character(df$value)]                 
      df = df[order(df$value,decreasing=F),]
      
      par(mar=c(1,2,2,2),xpd=NA)
      
      p = plot(df$x,df$y,axes=F,ylab="",xlab="",
               col=access_col[as.character(df$value2)],pch=20,
               main=input$OCR,cex.main=1.5,cex=0.9)
      legend("topright",inset=0.05, legend=c("open", "closed"),pch=16,
             col=c("#000000", "#BFBFBF"), cex=1.2)
      
      
    })
    
    
    barplotAccess_subset = reactive({
      meta_endo=meta[is.na(meta$al_haem_endo_clusters)==F,]
      
      
      value = dataAccess_subset()
      
      tab0 = table(value,as.character(meta_endo$al_haem_endo_clusters))["1",]
      tab = tab0*100/table(meta$al_haem_endo_clusters)
      df = data.frame(
        ann = as.vector(names(tab)),
        values = as.vector(unname(tab)),
        col = as.character(clust_endo_colours[as.character(names(tab))])
      )
      
      b = ggplot(data=df, aes(x=ann, y=values)) +
        geom_bar(stat="identity", width=0.5,fill= df$col) +
        theme_bw() +
        labs(y = "Nuclei with accessibility (%)") + 
        ggtitle(input$OCR_subset) +
        theme(axis.title = element_text(face = "bold", size = 14),
              axis.text.y = element_text(size = 12, face = "bold"),
              axis.text.x = element_text(size = 12, face = "bold", 
                                         angle = 45, hjust = 1),
              legend.position = "none",
              axis.title.x = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_blank()
              #panel.grid.minor = element_blank()
        ) +
        annotate("text", 
                 x = sort(names(tab)), 
                 y = rep_len(c(max(tab)*1.1, max(tab) * 1.2), 
                             length.out = length(tab)), 
                 label = tab0[sort(names(tab))]) 
      
      return(b)
      
    })
    
    
    output$Access_plot_subset <- renderPlot({
      
      validate(
        need(input$OCR_subset %in% names(snATAC_coord_indices),
             "The coordinate is not valid. Please check the 'Genomic regions' tab to explore valid coordinates and enter one in the form of e.g. chr1_3008841_3009342." )
      )
      plotAccess_subset()
      
    })
    
    output$Access_barplot_subset <- renderPlot({
      
      validate(
        need(input$OCR_subset %in% names(snATAC_coord_indices),
             "The coordinate is not valid. Please check the 'Genomic regions' tab to explore valid coordinates and enter one in the form of e.g. chr1_3008841_3009342." )
      )
      barplotAccess_subset()
      
    })
    
    output$downloadAccess_subset <- downloadHandler(
      filename = function() { paste0(input$OCR_subset, "_AL_EC_HAEM_UMAP.pdf") },
      content = function(file) {
        pdf(file, width=7,height=7.5)
        ac = c(
          "1"="open",
          "0"="closed"
        )
        
        meta_endo=meta[is.na(meta$al_haem_endo_clusters)==F,]
        
        
        
        df = data.frame(x = meta_endo[,"umap_X"],
                        y=meta_endo[,"umap_Y"],
                        value = dataAccess_subset()
        )
        
        
        df$value2 = ac[as.character(df$value)]                 
        df = df[order(df$value,decreasing=F),]
        #bottom, left, top and right
        
        par(mar=c(1,2,2,2),xpd=NA)
        
        p = plot(df$x,df$y,axes=F,ylab="",xlab="",
                 col=access_col[as.character(df$value2)],pch=20,
                 main=input$OCR,cex.main=1.5)
        legend("topright",inset=0.05, legend=c("open", "closed"),pch=16,
               col=c("#000000", "#BFBFBF"), cex=1.2)
        dev.off()
      }
    )
    
    
    output$download_barplotAccess_subset <- downloadHandler(
      filename = function() { paste0(input$OCR_subset, "_AL_EC_HAEM_barplot.pdf") },
      content = function(file) {
        pdf(file, width=10,height=6)
        print(barplotAccess_subset())
        dev.off()
      }
    )
    
 
    
    
    
    #########==================================
    # Patterns
    #########==================================
    
    #Set images
   # output$patterns <- renderImage({
    #  outfile <- tempfile(fileext='.png')
    #  p <- readPNG("./www/images/patterns.png")
    #  writePNG(p, target=outfile)
    #  list(src = outfile, contentType = 'image/png', width = 550, height = 450)
    #}, deleteFile = TRUE)
    
    output$selected_pattern1 <-renderText({
      
      validate(
        need(input$pattern_input %in% seq(1:12),
             "Pattern number not valid. Select one from the drop-down menu." )
      )
      html_pattern=paste0("<b> Pattern ",input$pattern_input," selected:</b>")
      return(html_pattern)
    })
    
    output$selected_pattern2 <-renderText({
      
      validate(
        need(input$pattern_input %in% seq(1:12),
             "Pattern number not valid. Select one from the drop-down menu." )
      )
      known=paste0("./links/HOMER_patterns/cluster_",input$pattern_input,"/knownResults.html")
      html_pattern_known = paste0("<a href=\"",known,"\" target=\"_blank\">Click here for known motifs</a>")
      
      return(html_pattern_known)
    })
    
    output$selected_pattern3 <-renderText({
      
      validate(
        need(input$pattern_input %in% seq(1:12),
             "Pattern number not valid. Select one from the drop-down menu." )
      )
      unknown=paste0("./links/HOMER_patterns/cluster_",input$pattern_input,"/homerResults.html")
      html_pattern_unknown = paste0("<a href=\"",unknown,"\" target=\"_blank\">Click here for <i>de novo</i> motifs</a>")
      return(html_pattern_unknown)
    })
    
    
  
  }#close all







