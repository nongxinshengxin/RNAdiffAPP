#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinythemes
#' @import colourpicker
#' @noRd
app_ui <- function(request) {
  # Define UI for application
  navbarPage(title=div(a("RNAdiff v1.0.0")),
             theme =shinythemes::shinytheme("flatly"),
             tabPanel("Home", sidebarLayout(
               sidebarPanel(h3("RNAdiff App",style = "font-family: 'times'"),
                            p("Make your RNA-Seq downstream analysis easy!",style = "font-family: 'times'; font-size:12pt;color:grey"),
                            br(),
                            strong("Welcome to follow our Wechat Official Account.", style = "font-family: 'times'; font-size:16pt"),
                            hr(),
                            p("Our Wechat Official Account is dedicated to sharing knowledge, sharing tools and sharing experiences.",style = "font-family: 'times'; font-size:14pt"),
                            br(),
                            p("For any feedback and tool suggestions, please feel free to leave a message at", a("Github issues.", href="https://github.com/nongxinshengxin/RNAdiffAPP/issues",style = "font-family: 'times'; font-size:14pt"), style = "font-family: 'times'; font-size:14pt"),
                            hr(),
                            hr(),
                            p(em("Email: nongxinshengxin@163.com",style = "font-family: 'times'; font-size:12pt;color:grey"))
               ),
               mainPanel(
                 h3("What is RNAdiff App ?",style = "font-family: 'times'"),
                 p("This app is used to perform",span("RNA-Seq downstream analysis",style="font-weight:bold"), "and is available to the user through an interactive interface.", style = "font-family: 'times'; font-size:14pt"),
                 h3("How is this RNAdiff App built ?",style = "font-family: 'times'"),
                 p("The app is built on the R language and shiny packages, calling DESeq2, edgeR, ggplot2 and other R packages.", style = "font-family: 'times'; font-size:14pt"),
                 h3("What can the RNAdiff App do ?",style = "font-family: 'times'"),
                 p("The app is divided into five main sections.", style = "font-family: 'times'; font-size:16pt"),
                 p(em("Function 1: Analysis of differentially expressed genes, which can be performed using both DESeq2 and edgeR methods."), style = "font-family: 'times'; font-size:14pt"),
                 p(em("Function 2: Plotting volcano maps, based on the results of differentially expressed gene analysis."), style = "font-family: 'times'; font-size:14pt"),
                 p(em("Function 3: Calculating TPM and plotting heatmap, or only plotting heatmap."), style = "font-family: 'times'; font-size:14pt"),
                 p(em("Function 4: GO or KEGG enrichment analysis. Enrichment analysis based on the clusterProfiler package."), style = "font-family: 'times'; font-size:14pt"),
                 p(em("Function 5: Plotting bubble maps, based on the results of GO or KEGG enrichment analysis."), style = "font-family: 'times'; font-size:14pt"),
                 h3("How to use the RNAdiff App ?",style = "font-family: 'times'"),
                 p("You can view the help documentation on our Wechat Official Accounts --",span("nongxinshengxin",style="font-weight:bold;color:blue;font-family: 'SimSun'"), style = "font-family: 'times'; font-size:14pt"),
                 includeMarkdown(system.file("app/www/wx.rmd", package = "RNAdiffAPP")),
                 hr()
               )
             )
             ),
             tabPanel("DEGanalysis",
                      sidebarLayout(
                        sidebarPanel(
                          div(
                            fileInput("matFile", "Choose Reads Matrix Data File"
                            )
                          ),
                          div(checkboxInput('header', 'Header', TRUE)),
                          radioButtons('sep','Sep',c(Tab='\t',Comma=',',Semicolon=';'),selected = '\t',inline=T),
                          div(fileInput("conditionFile","Choose Condition File",accept = c(".txt"))),
                          fluidRow(column(6,textInput("control","Control Group",value = "control")),
                                   column(6,textInput("case","Test Group",value = "case"))),
                          selectInput("dataset","pick ALL or Up or Down",choices = c("ALL","Up","Down")),
                          selectInput("tools","Choose DESeq2 or edgeR ?",choices = c("DESeq2","edgeR")),
                          div(actionButton("action","Start")),
                          hr(),
                          div(downloadButton("downlodData","Download"))
                          #div(downloadButton("downlodData","Download"))
                        ),
                        mainPanel(actionLink("win1","Click here show you FeatureCounts Data Format"),
                                  br(),
                                  actionLink("win2","Click here show you Condition Data Format"),
                                  div(helpText("Show the results in this view."),
                                      tags$hr(),
                                      dataTableOutput("Result"))))),
             tabPanel("Volcanoplot",sidebarLayout(
               sidebarPanel(
                 fileInput("volcanoFile", "Choose DESeq2 or edgeR Result File"),
                 helpText("The output data of DEGanalysis is the input data."),
                 checkboxInput('header2', 'Header', TRUE),
                 radioButtons('sep2','Sep',c(Tab='\t',Comma=',',Semicolon=';'),selected = '\t',inline=T),
                 fluidRow(column(4,textInput("name1","column name of geneID",value = "geneid")),
                          column(4,textInput("name2","column name of FDR",value = "FDR")),
                          column(4,textInput("name3","column name of logFC",value = "logFC"))),
                 textInput("fdrvalue","Choose your FDR threshold",value = 0.05),
                 textInput("FCvalue","Choose your logFC threshold",value = 1),
                 fileInput("geneIDlist", "Choose geneID list"),
                 actionButton("plotstart","Start to Plot")

               ),
               mainPanel(splitLayout(div(helpText("Show the results in this view."),
                                         plotOutput("plot1"),helpText("Show preview of input data"),tableOutput("text1")),div(colourpicker::colourInput("col1","Select down genes color","red"),
                                                                                                                              colourpicker::colourInput("col2","Select nosig genes color","#00B2FF"),
                                                                                                                              colourpicker::colourInput("col3","Select up genes color","orange"),
                                                                                                                              sliderInput("alpha","Choose points alpha",min = 0.1,max = 1,value = 1),
                                                                                                                              sliderInput("size","Choose points size",min = 0.5,max = 6,value = 2),
                                                                                                                              colourpicker::colourInput("textcol","Select text color","black"),
                                                                                                                              sliderInput("textsize","Choose text size",min = 1,max = 7,value = 3),

                                                                                                                              textInput("volcanotiltle","Set your plot tilte",value = "Volcano Plot"),
                                                                                                                              fluidRow(column(6,textInput("plot1width","Plot width",value = 7)),
                                                                                                                                       column(6,textInput("plot1height","Plot height",value = 5))),
                                                                                                                              radioButtons('extPlot', 'Plot output format',choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                                                                                                                              downloadButton("plot1downloadData","Download Plot")
                                         ),
                                     cellWidths = c("70%","30%"))
               ))),
             tabPanel("TPM&heatmap",sidebarLayout(
               sidebarPanel(
                 fileInput("matrixFile", "Choose Matrix File"),
                 helpText("The output data of FeatureCounts or Expression Matrix is the input data."),
                 checkboxInput('headerT', 'Header', TRUE),
                 radioButtons('sepT','Sep',c(Tab='\t',Comma=',',Semicolon=';'),selected = '\t',inline=T),
                 selectInput("TP","Calculate and Plot or Only Plot?",choices = c("only","both")),
                 actionButton("TPplotStart","Start")

               ),
               mainPanel(
                 splitLayout(
                   div(actionLink("winht","Click here show you input Data Format"),helpText("Show the results in this view."),plotOutput("plot_ht"),helpText("Show preview of input data"),tableOutput("table_tpm")),
                   div(colourpicker::colourInput("col_tp1","Select higher value color","red"),
                       colourpicker::colourInput("col_tp2","Select lower value color","blue"),
                       colourpicker::colourInput("col_tp3","Select middle value color","white"),
                       fluidRow(column(6,checkboxInput('row_c', 'Col cluster?', TRUE)),
                                column(6,checkboxInput('col_c', 'Row cluster?', TRUE))),
                       fluidRow(column(6,textInput("htwidth","Plot width",value = 7)),
                                column(6,textInput("htheight","Plot height",value = 5))),
                       radioButtons('ext2ht', 'Plot output format',choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                       helpText("Download Heatmap"),
                       downloadButton("download_ht","Download"),
                       helpText("Download TPM"),
                       downloadButton("download_tpm","Download")),
                   cellWidths = c("70%","30%")
                 )
               ))),
             tabPanel("ClusterProfiler",sidebarLayout(
               sidebarPanel(
                 div(
                   fileInput("genelist", "Select the gene list file you want to enrich.")
                 ),
                 div(checkboxInput('header4', 'Header', FALSE)),
                 div(fileInput("term2gene","Select the list of term2gene for a particular species.",accept = c(".txt",
                                                                                                               ".csv"))),
                 checkboxInput('header5', 'Header', TRUE),
                 radioButtons('sep4','Sep',c(Tab='\t',Comma=',',Semicolon=';'),selected = '\t',inline=T),
                 fileInput("term2name","Select the list of term2name for a particular species.",accept = c(".txt",
                                                                                                           ".csv")),
                 checkboxInput('header6', 'Header', TRUE),
                 radioButtons('sep5','Sep',c(Tab='\t',Comma=',',Semicolon=';'),selected = '\t',inline=T),
                 selectInput("ea_tools","Choose GOea or KEGGea ?",choices = c("goea","keggea")),
                 div(actionButton("ea_action","Start")),
                 hr(),
                 textInput("prefix","Add a file prefix",value = "ClusterProfiler"),
                 div(downloadButton("ea_downloadData","Download"))
               ),
               mainPanel(
                 actionLink("cwin1","Click here show you genelist Data Format"),
                 br(),
                 actionLink("cwin2","Click here show you term2gene Data Format"),
                 br(),
                 actionLink("cwin3","Click here show you term2name Data Format"),
                 div(helpText("Show the results in this view."),
                     tags$hr(),
                     dataTableOutput("ea_Result"))))),
             tabPanel("BubblePlot",sidebarLayout(
               sidebarPanel(
                 fileInput("bubbleFile", "Choose GO/KEGG Enrichment Result File"),
                 helpText("The output data of ClusterProfiler is the input data."),
                 checkboxInput('header3', 'Header', TRUE),
                 radioButtons('sep3','Sep',c(Tab='\t',Comma=',',Semicolon=';'),selected = '\t',inline=T),
                 fluidRow(column(4,textInput("b_name1","column name of GO/KEGG term",value = "Description")),
                          column(4,textInput("b_name2","column name of FDR/p.adjust",value = "p.adjust")),
                          column(4,textInput("b_name3","column name of Ontology (only in GOEA)",value = "Ontology"))),
                 fluidRow(column(6,textInput("b1","column name of Study.term",value = "Study.term")),
                          column(6,textInput("b2","column name of Study.total",value = "Study.total"))),
                 fluidRow(column(6,textInput("b3","column name of Pop.term",value = "Pop.term")),
                          column(6,textInput("b4","column name of Pop.total",value = "Pop.total"))),
                 textInput("fdr2value","Choose your FDR threshold",value = 0.05),
                 actionButton("plot2start","Start to Plot")

               ),
               mainPanel(splitLayout(div(helpText("Show the results in this view."),
                                         plotOutput("plot2"),helpText("Show preview of input data"),tableOutput("text2")),div(colourpicker::colourInput("collow","Select lower value color","red"),
                                                                                                                              colourpicker::colourInput("colhigh","Select higher value color","#00B2FF"),
                                                                                                                              sliderInput("sizerange","Choose points size range",min = 0.5,max = 6,value = c(0.5,6)),
                                                                                                                              textInput("bubbletiltle","Set your plot tilte",value = "Bubble Plot"),
                                                                                                                              fluidRow(column(6,textInput("plot2width","Plot width",value = 7)),
                                                                                                                                       column(6,textInput("plot2height","Plot height",value = 5))),
                                                                                                                              radioButtons('ext2Plot', 'Plot output format',choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                                                                                                                              downloadButton("plot2downloadData","Download Plot")
                                         ),
                                     cellWidths = c("80%","20%"))
               ))),
             tags$footer(p("Contact: nongxinshengxin@163.com"), align="center")
  )

}


#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "RNAdiffAPP"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
