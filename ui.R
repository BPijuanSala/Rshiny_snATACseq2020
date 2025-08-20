########################################################################################
## Title: ui.R
## Author: Blanca Pijuan-Sala
## Description: Shiny app for snATAC-seq E8.25 embryos
## Date: 16 December 2019
########################################################################################


library(shiny)
library(Matrix)
library(data.table)
load("./www/data/TFs_chromVAR_abbrev_up.rda")

ui <- fluidPage(

  
  
  tags$head(
    tags$style(HTML("
                    
                    h1 {
                    padding-left: 20px;
                    padding-top: 40px;
                    padding-bottom: 40px;
                    
                    background: #de8c12; /* For browsers that do not support gradients */
                    background: -webkit-linear-gradient(#de8c12, white); /* For Safari 5.1 to 6.0 */
                    background: -o-linear-gradient(#de8c12, white); /* For Opera 11.1 to 12.0 */
                    background: -moz-linear-gradient(#de8c12, white); /* For Firefox 3.6 to 15 */
                    background: linear-gradient(#de8c12, white); /* Standard syntax */
                    background-color: #de8c12;
                    color: black;
                    }
                    .explain {
                    
                    padding-left: 20px;
                    padding-top: 20px;
                    padding-bottom: 20px;
                    }
                    
                    .footer{ 
                    text-align: right;    
                    bottom: 0px; 
                    color: gray;
                    margin-right: 20px
                    margin-left: 20px
                    margin-top: 200px
                    margin-bottom: 20px
                    
                    
                    }
                    .downloadBut{ 
                    text-align: center; 
                    
                    
                    }
                    .downloadHeader{
                    align:center;
                    background-color: #f0f0f5;
                    margin-right: 20px;
                    margin-left: 20px;
                    margin-top: 20px;
                    margin-bottom: 20px;
                    padding-left: 20px;
                    padding-top: 10px;
                    padding-bottom: 10px;
                    
                    }
                    
                    "))
                    ),
#  tags$head(tags$link(rel="shortcut icon", href="./images/favicon.ico")),
  
  
 # titlePanel("Single-nucleus ATAC-seq of E8.25 mouse embryos"),
  
  headerPanel(
    "Single-nucleus ATAC-seq of E8.25 mouse embryos"),
  
   tags$div(
     class="explain",
     tags$p("This is the website accompanying the manuscript",tags$a(href="",target="_blank",
                                                                     "Pijuan-Sala, B., et. al., Nat. Cell Biol., 2020."),
            " A custom UCSC Genome Browser session is available at",tags$a(href="https://tinyurl.com/snATACseq-GSE133244-UCSC",target="_blank",
                                                                          "https://tinyurl.com/snATACseq-GSE133244-UCSC."))#,
     
    # tags$p(tags$span(style="color:red","Temporary note: External links may not work until publication. All data has been made available to reviewers."))
     
     ),
  tabsetPanel(
    tabPanel("Main",


             fluidRow(
               column(5, offset=0.5,
                      HTML('<img src="./images/UMAP_website.png", 
                    width="600px" height="500px" style="float:right"/>',
                           '<p style="color:black"></p>')
                     # imageOutput("snATACseq")
                      
               ),#close column for images
               column(4,offset=1,
                      tags$br(),
                      tags$br(),
                      
                      tags$p("This website contains the following tabs:"),
                      tags$br(),
                      
                      tags$p(tags$b("Experimental design:"),"Here, you can find a brief outline of the experimental design of this dataset."),
                      tags$br(), 

                      tags$p(tags$b("Topics:"),"Here, you can explore the contribution of each topic to each nucleus and the transcription factor enrichment results for each topic."),
                      tags$br(),
                      
                      tags$p(tags$b("Single-cell TF enrichment:"),"Here, you can visualise the enrichment scores for transcription factors at the single-cell level, computed using chromVAR",
                             tags$a(href="https://www.nature.com/articles/nmeth.4401",target="_blank",
                                    "(Schep et al., Nature Methods, 2017)"),"."),
                      tags$br(),
                      
                      tags$p(tags$b("Genomic regions:"),"Here, you can browse the list of genomic regions and their associated details."),
                      tags$br(),                     
                      
                      tags$p(tags$b("Accessibility of OCRs:"),"Here, you can check the accessibility for specific genomic coordinates in the UMAP visualisaiton."),
                      tags$br(),
                      
                      tags$p(tags$b("Allantoic-haemato-endothelial landscape:"),"Here, you can check the accessibility for specific genomic coordinates in the UMAP visualisaiton for the allantoic-haemato-endothelial subset."),
                      
                      tags$br(),
                      
                      tags$p(tags$b("Dynamic patterns towards endothelium:"),"Here, you will visualise the dynamic patterns of accessibility and access transcription factor enrichment results."),
                      
                      tags$br(),
                      
                      tags$p(tags$b("External links:"),"Here, we provide links to the paper, the data and the code used for the analysis.")
                      
                )#close column
             ) #close fluidrow
    ),#close tab panel
    
    
    
    
    tabPanel("Experimental design",
             fluidRow(
               column(7,offset=1,
               tags$br(),
               tags$p("We dissected mouse embryos at days 8.25 post-fertilisation and snap-froze them. Samples were then processed using",
                      tags$a(href="https://www.nature.com/articles/s41593-018-0079-3",target="_blank","single-nucleus ATAC-seq"),
                      "with slight modifications (please see the",tags$a(href="",target="_blank","publication"),"for further details)."),
               tags$br(),
               tags$br(),
               HTML('<img src="./images/diagram.png", 
                    width="1000px" style="float:right"/>',
                    '<p style="color:black"></p>')
               #imageOutput("diagram")
               )#close column
               )#Explanation
             ),#close tabPanel experimental design
   
    tabPanel( "Topics",
              tags$br(),

              fluidRow(
                
                column(3, offset=1,
                       fluidRow(
                         selectInput("Topic", label = h3("Select topic"),
                                     
                                     choices = paste("Topic",seq(1,100)))
                       ),#close fluidrow for input topics
                       fluidRow(
                         shinycssloaders::withSpinner(plotOutput("Topic_plot",
                                                                 height = "550px",
                                                                 width = "370px"
                             
                         ),color="#0dc5c1")#close spinner
                       ),#close fluidRow for plot topic
                       fluidRow(
                         downloadButton("downloadTopic", label = "Download plot")
                         
                       )#close fluidRow for download

                       
                       
                       ),#close entire first column
                    
                column(7,offset=0.5,
                       tags$br(),
                       tags$br(),
                       tags$h4("TF enrichment analysis on regions uniquely contributing to the topic"),
                       tags$br(),
                       #DT::dataTableOutput("TF_topics_DT")

                       htmlOutput("topics_motif1"),
                       tags$br(),
                       tags$br(),
                       
                       tags$p(tags$span(style="color:red","Check the regions contributing to each topic in the ",tags$b('Genomic regions')," tab!"))
                       
                )#close column motifs
             
               
                
                       
                )#close fluidRow
      
    ),#Close tabPanel topics
    
    
    tabPanel("Single-cell TF enrichment",
             tags$br(),
             tags$br(),
             
             fluidRow(
               
               column(3, offset=1,
                      fluidRow(
                        selectizeInput("TF", label = h3("Select transcription factor"),
                                       choices=sort(TFs_chromVAR_abbrev_up), selected="GATA1")
                      ),#close fluidRow input
                      fluidRow(
                        shinycssloaders::withSpinner(plotOutput("chromVAR_plot",
                                                                height = "480px",
                                                                width = "350px"),#close plot output
                                  
                        color="#0dc5c1")#close spinner
                      ),#close fluidrow for chromVAR plot
                      
                      tags$br(),

                      fluidRow(
                        downloadButton("download_chromVAR", label = "Download UMAP")
                        
                      ),#close fluidrow for download button
                      tags$br(),

                      fluidRow(
                        downloadButton("download_TF_RNA", label = "Download RNA-TF plot")
                        
                      ),#close fluidrow for download button
                      tags$br(),
                      
                      fluidRow(
                        downloadButton("download_motif", label = "Download motif logo")
                        
                      )#close fluidrow for download button
                      
                      
                      ),#close first column
               column(5,offset=1,
                      tags$br(),
                      tags$br(),
                      
                      fluidRow(
                        htmlOutput("selected_chromVAR")
                      ),
                      tags$br(),
                      tags$br(),
                      
                      fluidRow(
                        shinycssloaders::withSpinner(plotOutput("TF_RNA",
                                                                height = "310px",
                                                                width = "700px"
                        ),color="#0dc5c1")
                      ),
                      tags$br(),
                      
                      tags$br(),
                      
                      fluidRow(
                        
                        imageOutput("motif")
                       )
                      
                      
                      
                      )#close second column
               
               
               
               
               
               ) #close first fluidRow
           
             
    ),#Close tab panel TF enrichment
    
    tabPanel("Genomic regions",
             fluidRow(
               column(10,offset=1,
                      tags$br(),
                      tags$p("You can browse the list of genomic coordinates below. By pressing the download button, you will be able to export your desired list of coordinates with their details into a .csv file. 
                             Please note that if a genomic region is assigned to multiple genes, the entry for this region will be repeated as many times as the number of associated genes.")
                      )),#close column and fluidrow
             fluidRow(
               column(5,offset=4,
                      tags$br(),
                      downloadButton("downloadOCR", "Download custom table")
               )#close column
             ),#close fluidrow download
             fluidRow(
               column(11.5,offset=1,
                      tags$br(),
                      shinycssloaders::withSpinner(
                        DT::dataTableOutput("OCR_DT"),color="#0dc5c1")
               )#close column
             )#close fluidrow for input genomic regions
               ),#close tab panel genomic regions
    
    
    tabPanel("Accessibility of OCRs",
             fluidRow(
               tags$br(),
               
               column(4, offset=1,
                      fluidRow(
                        textInput("OCR", label = h3("Select a genomic region"),
                                  placeholder="Enter a PeakID",
                                  value = "chr5_147657804_147658307")
                      ),#close fluidrow for input topics
                      fluidRow(
                        tags$br(),
                        shinycssloaders::withSpinner(plotOutput("Access_plot",
                                                                height = "500px",
                                                                width = "400px"
                        ),color="#0dc5c1")
                      ),#close fluidRow for plot OCR
                      tags$br(),
                      tags$br(),
                      
                      fluidRow(
                        column(3, offset=1.5,
                               downloadButton("downloadAccess", label = "Download plot")
                        )
                      )#close fluidRow for download
                      
               ),#close entire first column
               
               column(6,offset=0.9,
                      fluidRow(
                        tags$br(),
                        tags$br(),
                        tags$br(),
                        
                        tags$p("The plot below shows the percentage of nuclei with the genomic region accessible per cell type. On top of the bars, you will find the actual number of nuclei with accessibility at the specified locus."),
                        tags$br(),
                        
                        shinycssloaders::withSpinner(plotOutput("Access_barplot",
                                                                height = "400px",
                                                                width = "700px"
                        ),color="#0dc5c1"),
                        tags$br()
                      ),#close fluid row
                      fluidRow(
                        column(7,offset=3,
                               downloadButton("download_barplotAccess", label = "Download barplot")
                        )
                      )#close fluidRow for download
                      
               )#close column barplots
               
             )#close fluidRow
             
    ),#close tab panel accessibility
    
    
    
    tabPanel("Allantoic-haemato-endothelial landscape",
             tags$br(),
             fluidRow(
               column(5,offset=1,
                      
                 textInput("OCR_subset", label = h3("Select a genomic region"),
                           placeholder="Enter a PeakID",
                           value = "chr5_147657804_147658307")
               )#close column
               ),#close fluidrow for input topics
           
             
             fluidRow(
               column(3,
                      tags$br(),

                      HTML('<img src="./images/UMAP_AL_EC_HAEM_celltypes_website.png", 
                    width="400px" height="350px" style="float:right"/>',
                           '<p style="color:black"></p>')
                      #imageOutput("UMAP_subset_celltype")
                      
                      ),#close column for celltype plot
               column(3,
                      tags$br(),

                      
                      HTML('<img src="./images/UMAP_AL_EC_HAEM_website.png", 
                    width="400px" height="350px" style="float:right"/>',
                           '<p style="color:black"></p>')
                     # imageOutput("UMAP_subset_celltype_subcolor")
   
               ),#close column for cell subtype plot
               column(4,
                     
                      fluidRow(
                        shinycssloaders::withSpinner(plotOutput("Access_plot_subset",
                                                                height = "390px",
                                                                width = "370px"
                        ),color="#0dc5c1")
                      ),#close fluidRow for plot OCRsubset

                      fluidRow(
                        column(3, offset=1,
                               downloadButton("downloadAccess_subset", label = "Download plot")
                        )
                      )#close fluidRow for download
                      
               )#close column plotting access
             ),#close fluidRow for plots celltypes and access
             fluidRow(
               column(10,offset=1,
                      fluidRow(
                        tags$br(),
                        tags$br(),
                        
                        
                        tags$p("The plot below shows the percentage of nuclei with the genomic region accessible per cell type. On top of the bars, you will find the actual number of nuclei with accessibility at the specified locus."),
                        tags$br(),
                        
                        tags$br()
                        
                      ))),#close fluidRow for sentence
             fluidRow(
               column(10,offset=3,
                      
                        shinycssloaders::withSpinner(plotOutput("Access_barplot_subset",
                                                                height = "400px",
                                                                width = "700px"
                        ),color="#0dc5c1"),
                        tags$br()
               )#close column for barplot
                      ),#close fluid row for barplot
                      fluidRow(
                        column(9,offset=5,
                               downloadButton("download_barplotAccess_subset", label = "Download barplot")
                        )#close column for download
                      )#close fluidRow for download
                      
               
 ),#close tab panel allantoic haematoendothelial landscape
   tabPanel("Dynamic patterns towards endothelium",
            fluidRow(
              column(5,offset=1,
                     tags$br(),
                     HTML('<img src="./images/patterns.png", 
                    width="600px" style="float:right"/>',
                          '<p style="color:black"></p>')
                    # imageOutput("patterns")
                     
                     
                     ),#close column image output
              
              
              
              column(4,offset=1,
                     fluidRow(
                     selectInput("pattern_input", label = h3("Check out TF enrichment results for pattern"),
                                 
                                 choices = seq(1:12))
                     ),#close fluidRow for input
                     fluidRow(
                       htmlOutput("selected_pattern1"),
                       tags$br(),
                       shinycssloaders::withSpinner(htmlOutput("selected_pattern2"),
                       color="#0dc5c1"),#close spinner
                       htmlOutput("selected_pattern3")
                     
                       
                     )#close fluidRow for links
                     
                     )
            )#close first fluidRow
     
   ), #close tabPanel patterns endothelium
    tabPanel("External links",
             fluidRow(
               column(10,offset=1,
                      tags$br(),
                    
                      tags$br(),
                      
                      
                      tags$p(tags$b("Paper:"),"This analysis is part of:",tags$a(href="",target="_blank",
                                                                                      "Pijuan-Sala, B., et. al., Nat. Cell Biol., 2020")),
                      tags$p(tags$b("Data:"),"Data has been deposited at GEO:",tags$a(href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133244",target="_blank",
                                                                                      "GSE133244")),
                      tags$p(tags$b("UCSC session:"),"Explore the Bigwig tracks for each cell type in ",tags$a(href="https://tinyurl.com/snATACseq-GSE133244-UCSC",target="_blank",
                                                                                                               "our UCSC session.")), 
                      
               
                      tags$p(tags$b("Code:"),"Code used to analyse this dataset is available at",tags$a(href="https://github.com/BPijuanSala/MouseOrganogenesis_snATACseq_2020",target="_blank",
                                                                                                        "https://github.com/BPijuanSala/MouseOrganogenesis_snATACseq_2020."),
                             "It contains a README.md file explaining each step of the code. "),
                      
                      tags$br(),
                      tags$p(tags$b("Citing this website:"),"If you use this website, we would be grateful if you could cite: ",tags$a(href="",target="_blank",
                                                                                                        "Pijuan-Sala, B., et. al., Nat. Cell Biol., 2020"))
                            
                      
                      )
             )
             
             
             )
    
    
    

  ),#tabset close
  tags$br(),
  tags$br(),
  tags$br(),
  tags$br(),
  
  tags$div(class="footer",
           tags$hr(),
           tags$p("Göttgens Laboratory 2020"))
                    )

