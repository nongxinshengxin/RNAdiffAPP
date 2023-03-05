#' The application server-side
#'
#' @param input,output Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinythemes
#' @import colourpicker
#' @import DESeq2
#' @import edgeR
#' @import ggplot2
#' @import ggrepel
#' @import clusterProfiler
#' @import tidyr
#' @import aplot
#' @import shinyalert
#' @import pheatmap
#' @import corrplot
#' @import Hmisc
#' @noRd
app_server <- function(input, output) {
  options(shiny.maxRequestSize = 50 * 1024^2)
  #############模块TPM&heatmap响应
  
  ##示例文件
  
  rawMat<-matrix(sample(0:100,size = 18),3,6)
  rownames(rawMat)<-paste("AT1G0",c(101:103),"0",sep = "")
  colnames(rawMat)<-c(paste("Control",c(1:3),sep = ""),paste("Treat",c(1:3),sep = ""))
  rawFC<-data.frame(chr=c(rep("Chr1",3)),start=c(1250,25966,35689),end=c(1350,27976,38689),strand=c(rep("+",3)),Length=c(1520,569,6541))
  rawFC<-cbind(rawFC,as.data.frame(rawMat))
  
  #通过点击按钮，可以弹出界面，这个界面可以显示ui的输出界面，因此要做tagList中添加ui输出控件
  #
  observeEvent(input$winht, {
    shinyalert(
      html = TRUE,
      text = tagList(
        "Raw FeatureCounts Data Format",
        tableOutput("raw_example"),
        "Matrix Data Format",
        tableOutput("mat_example")
      ),
      size = "l")
  })
  
  ##关联示例控件的函数
  output$raw_example <- renderTable({
    
    head(rawFC)
  },rownames = T)
  
  output$mat_example <- renderTable({
    
    head(rawMat)
  },rownames = T)
  
  ###矩阵数据处理、cor计算、TPM计算
  mat<-eventReactive(input$TPplotStart,{
    if (input$TP=="heatmap"){
      rawmat<-input$matrixFile
      if (is.null(rawmat))
        return(NULL)
      htmat <- read.delim(file = rawmat$datapath, header = input$headerT, row.names = 1,as.is = TRUE, sep = input$sepT,encoding='UTF-8')
      htmat
    }else if(input$TP=="cor"){
      rawmat<-input$matrixFile
      if (is.null(rawmat))
        return(NULL)
      countdata <- read.delim(file = rawmat$datapath, header = input$headerT, row.names = 1,as.is = TRUE, sep = input$sepT,encoding='UTF-8')
      #countdata<-as.matrix(countdata)
      countdata<-as.matrix(countdata)
      #cors<-cor(countdata, method = c("pearson"))
      
      if (input$matNum=="one"){
        cors<-Hmisc::rcorr(countdata,type = input$cor_m)
        cors$P[is.na(cors$P)]<-1
        cors
      }else{
        rawmat2<-input$secondMatFile
        if (is.null(rawmat2))
          return(NULL)
        countdata2 <- read.delim(file = rawmat2$datapath, header = input$headerT, row.names = 1,as.is = TRUE, sep = input$sepT,encoding='UTF-8')
        countdata2<-as.matrix(countdata2)
        cors<-Hmisc::rcorr(countdata,countdata2, type = input$cor_m)
        cors$P[is.na(cors$P)]<-1
        cors
      }
      
    }
    
  })
  
  ###计算TPM CPM FPKM
  calculatemat<-eventReactive(input$TPcalculateStart,{
    rawmat<-input$matrixFile
    if (is.null(rawmat))
      return(NULL)
    countdata <- read.delim(file = rawmat$datapath, header = input$headerT, row.names = 1,as.is = TRUE, sep = input$sepT,encoding='UTF-8')
    metadata <- countdata[,1:5]#提取基因信息count数据前的几列
    countdata <- countdata[,6:ncol(countdata)]#提取counts数，counts数据主题部分
    cpm <- t(t(countdata)/colSums(countdata) * 1000000)#参考cpm定义
    #avg_cpm <- data.frame(avg_cpm=rowMeans(cpm))
    #-----TPM Calculation------
    kb <- metadata[,5] / 1000
    rpk <- countdata / kb
    tpm <- t(t(rpk)/colSums(rpk) * 1000000)
    fpkm <- t(t(rpk)/colSums(countdata) * 10^6)
    if (input$CP=="TPM")
      return(tpm)
    if (input$CP=="CPM")
      return(cpm)
    if (input$CP=="FPKM")
      return(fpkm)
  })
  
  
  #列分组信息输入
  colgroup<-eventReactive(input$TPplotStart,{
    colg<-input$colgroupFile
    if (is.null(colg))
      return(NULL)
    cg <- read.delim(file = colg$datapath, header = T, row.names = 1,as.is = TRUE, sep = "\t",encoding='UTF-8')
    cg
  })
  
  
  #行分组信息输入
  rowgroup<-eventReactive(input$TPplotStart,{
    colg<-input$rowgroupFile
    if (is.null(colg))
      return(NULL)
    cg <- read.delim(file = colg$datapath, header = T, row.names = 1,as.is = TRUE, sep = "\t",encoding='UTF-8')
    cg
  })
  
  output$table_tpm <- renderTable({
    
    if (input$plotOrcal=="plot"){
      if (is.null(mat()))
        return(NULL)
      head(mat())
    }else{
      if (is.null(calculatemat()))
        return(NULL)
      head(calculatemat())
    }
  },rownames = T)
  
  
  drawht<-reactive({
    ###热图绘制
    scale_test <- apply(mat(), 2, function(x){log2(x+1)})
    
    pheatmap(mat =scale_test,
             cluster_cols = input$row_c,
             cluster_rows = input$col_c,
             angle_col = "45",
             cellwidth=15,
             colorRampPalette(colors = c(input$col_tp2,input$col_tp3,input$col_tp1))(100),
             annotation_col = colgroup(),
             annotation_row = rowgroup(),
             annotation_names_col = input$col_a,
             annotation_names_row = input$row_a
             # annotation_names_row = F
    )#适用于平均值
  })
  
  
  drawcor<-reactive({
    corrplot::corrplot(mat()$r, type = input$corType, method =input$corMethod,order = input$corOrder,tl.col = "black", tl.srt = 45,
             col = colorRampPalette(colors = c(input$col_tp2,input$col_tp3,input$col_tp1))(100),
             is.corr = input$corlim,p.mat = mat()$P, insig = "label_sig",
             sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = input$sig_col)
    
  })
  
  
  output$plot_ht<-renderPlot({
    if (input$plotOrcal=="calculate")
      return(NULL)
    
    if (input$TP=="cor"){
      if (is.null(drawcor()))
        return(NULL)
      drawcor()
    }else if(input$TP=="heatmap"){ 
      if (is.null(drawht()))
        return(NULL)
      drawht()}
  })
  
  
  #下载图片
  output$download_ht<-downloadHandler(
    filename = function() {
      if (input$TP=="cor"){
        paste('corHeatmap', Sys.Date(), '.',input$ext2ht, sep='')
      }else if (input$TP=="heatmap"){
        paste('Heatmap', Sys.Date(), '.',input$ext2ht, sep='')
      }
      
    },
    content=function(file){
      if (input$TP=="cor"){
        if (input$ext2ht=="png"){
          png(file)
        }else if(input$ext2ht=="pdf"){
          pdf(file,width =as.numeric(input$htwidth), height = as.numeric(input$htheight) )
        }else{
          jpeg(file)
        }
        corrplot::corrplot(mat()$r, type = input$corType, method =input$corMethod,order = input$corOrder,tl.col = "black", tl.srt = 45,
                 col = colorRampPalette(colors = c(input$col_tp2,input$col_tp3,input$col_tp1))(100),
                 is.corr = input$corlim,p.mat = mat()$P, insig = "label_sig",
                 sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = input$sig_col)
        dev.off()
        #ggsave(file,plot =replayPlot(p1()), width = as.numeric(input$htwidth), height = as.numeric(input$htheight))
      }else if(input$TP=="heatmap"){
        ggsave(file,plot = drawht(), width = as.numeric(input$htwidth), height = as.numeric(input$htheight))
      }
      
    })
  
  
  ##下载tpm
  output$download_tpm<-downloadHandler(filename = function() {paste("TPM", Sys.Date(), ".csv", sep="")},content=function(file) {
    write.csv(calculatemat(), file)
  })
  
  #DESeq2流程
  #eventReactive 隔离流程
  result<-eventReactive(input$action,{
    if (input$tools=="DESeq2"){
      ### DESeq2流程
      rawmat<-input$matFile
      if (is.null(rawmat))
        return(NULL)
      mycounts <- read.delim(file = rawmat$datapath, header = input$header, row.names = 1,as.is = TRUE, sep = input$sep,encoding='UTF-8')
      rawfile<-input$conditionFile
      if (is.null(rawfile))
        return(NULL)
      conditiondata<-read.table(file = rawfile$datapath,header = F,encoding='UTF-8')
      condition<-c(conditiondata[,1])
      condition <- factor(condition)
      colData <- data.frame(row.names = colnames(mycounts), condition)
      colData
      dds <- DESeqDataSetFromMatrix(mycounts, colData, design = ~condition)
      dds <- DESeq(dds)
      res <- results(dds, contrast=c("condition", input$case, input$control))
      res = res[order(res$pvalue),]
      res<-na.omit(res)
      up<-subset(res, padj < 0.05 & log2FoldChange > 1)
      #uplist<-rownames(up)
      down<-subset(res, padj < 0.05 & log2FoldChange < -1)
      #downlist<-rownames(down)
      if (input$dataset=="ALL"){
        datasetOutput<-res
      }
      else if (input$dataset=="Up"){
        datasetOutput<-up
      }else {
        datasetOutput<-down
      }
      datasetOutput}
    
    else {
      ### edgeR流程
      rawmat<-input$matFile
      if (is.null(rawmat))
        return(NULL)
      x<-read.delim(file = rawmat$datapath, header=input$header,row.names=1, stringsAsFactors=FALSE, as.is = TRUE, sep = input$sep,encoding='UTF-8')
      rawfile<-input$conditionFile
      if (is.null(rawfile))
        return(NULL)
      group<-read.table(file = rawfile$datapath,header = F,encoding='UTF-8')
      
      y <- DGEList(counts=x, group=group[,1])
      keep <- rowSums(cpm(y)>1) >= 2
      y <- y[keep,]
      y$samples$lib.size <- colSums(y$counts)
      y <- calcNormFactors(y)
      
      group<-factor(group[,1])
      design<-model.matrix(~0+group)
      colnames(design)<-levels(group)
      rownames(design)<-colnames(x)
      
      y <- estimateGLMCommonDisp(y,design)
      y <- estimateGLMTrendedDisp(y,design)
      y <- estimateGLMTagwiseDisp(y,design)
      fit<-glmFit(y,design)
      
      design_df<-as.data.frame(design)
      control_num<-which(names(design_df)==input$control)
      case_num<-which(names(design_df)==input$case)
      if (control_num <case_num){
        case_num=-1
        
      }
      if (control_num >case_num){
        control_num=-1
      }
      lrt<-glmLRT(fit,contrast = c(case_num,control_num))
      
      
      # Diff table
      diff <- as.data.frame(topTags(lrt, n = nrow(y)))
      
      # Summary table
      diff<-na.omit(diff)
      up<-subset(diff, FDR < 0.05 & logFC > 1)
      #uplist<-rownames(up)
      down<-subset(diff, FDR < 0.05 & logFC < -1)
      #downlist<-rownames(down)
      if (input$dataset=="ALL"){
        datasetOutput<-diff
      }
      else if (input$dataset=="Up"){
        datasetOutput<-up
      }
      else {
        datasetOutput<-down
      }
      datasetOutput
    }
  })
  
  
  
  output$Result <- renderDataTable({
    
    if (is.null(result()))
      return(NULL)
    result()
  })
  
  ###设置示例数据
  test1<-matrix(sample(0:100,size = 60),10,6)
  rownames(test1)<-paste("AT1G0",c(101:110),"0",sep = "")
  colnames(test1)<-c(paste("Control",c(1:3),sep = ""),paste("Treat",c(1:3),sep = ""))
  
  test2<-c(rep("Control",3),rep("Treat",3))
  
  #通过点击按钮，可以弹出界面，这个界面可以显示ui的输出界面，因此要做tagList中添加ui输出控件
  #弹窗一
  observeEvent(input$win1, {
    shinyalert(
      html = TRUE,
      text = tagList(
        "FeatureCounts Data Format",
        tableOutput("example1"),
      ),
      size = "l")
  })
  
  ##关联示例控件的函数
  output$example1 <- renderTable({
    
    head(test1)
  },rownames = T)
  
  #弹窗2及关联函数
  observeEvent(input$win2, {
    shinyalert(
      html = TRUE,
      text = tagList(
        "Condition Data Format",
        tableOutput("example2")
      ),
      size = "l")
  })
  
  output$example2 <- renderTable({
    
    head(test2)
  },rownames = F,colnames = F)
  
  
  output$downlodData<-downloadHandler(filename = function() {paste(input$dataset, input$tools, ".csv", sep="")},content=function(file) {
    write.csv(result(), file)
  })
  
  ########################
  #火山图
  volcanodf<-eventReactive(input$plotstart,{
    rawdf<-input$volcanoFile
    if (is.null(rawdf))
      return(NULL)
    dataset<-read.delim(file = rawdf$datapath, header=input$header2, stringsAsFactors=FALSE, as.is = TRUE, sep = input$sep2,encoding='UTF-8')
    cut_off_fdr =as.numeric(input$fdrvalue)
    cut_off_logFC= as.numeric(input$FCvalue)
    # padj<-input$name2
    # log2FoldChange<-input$name3
    idpos<-which(names(dataset)==input$name1)
    fdrpos<-which(names(dataset)==input$name2)
    fcpos<-which(names(dataset)==input$name3)
    
    dataset$change = ifelse(dataset[,fdrpos] < cut_off_fdr & abs(dataset[,fcpos]) >= cut_off_logFC,
                            ifelse(dataset[,fcpos]> cut_off_logFC ,'Up','Down'),
                            'Stable')
    dataset<-na.omit(dataset)
    dataset<-dataset[,c(idpos,fdrpos,fcpos,ncol(dataset))]
    dataset
  })
  
  
  iddf<-eventReactive(input$plotstart,{
    rawdf<-input$geneIDlist
    if (is.null(rawdf))
      return(NULL)
    IDdata<-read.delim(file = rawdf$datapath, header=F, stringsAsFactors=FALSE, as.is = TRUE,sep = "\t",encoding='UTF-8')
    iddf<-merge(IDdata,volcanodf(),by.x=names(IDdata),by.y=names(volcanodf())[1],all.x)
    iddf
  })
  
  
  output$text1 <- renderTable({
    
    if (is.null(input$volcanoFile))
      return(NULL)
    raw<-input$volcanoFile
    df<-read.delim(file = raw$datapath, header=input$header2, stringsAsFactors=FALSE, as.is = TRUE, sep = input$sep2,encoding='UTF-8')
    head(df)
  })
  
  
  plot1_obj<-reactive({
    if (is.null(volcanodf()))
      return(NULL)
    if (is.null(iddf())){
      ggplot(
        # draw plot
        volcanodf(), aes(x = volcanodf()[,3], y = -log10(volcanodf()[,2]), colour=change)) +
        geom_point(alpha=input$alpha, size=input$size) +
        scale_color_manual(values=c(input$col1, input$col2,input$col3))+
        # draw line
        geom_vline(xintercept=c(-as.numeric(input$FCvalue),as.numeric(input$FCvalue)),lty=4,col="#990000",lwd=0.8) +
        geom_hline(yintercept = -log10(as.numeric(input$fdrvalue)),lty=4,col="#990000",lwd=0.8) +
        # change labs
        labs(x="log2(fold change)",
             y="-log10 (FDR)",
             title = input$volcanotiltle)+
        theme_bw()+
        # set theme
        theme(plot.title = element_text(hjust = 0.5),
              legend.position="right",
              legend.title = element_blank())
    }
    else {
      ggplot(
        # draw plot
        volcanodf(), aes(x = volcanodf()[,3], y = -log10(volcanodf()[,2]), colour=change)) +
        geom_point(alpha=input$alpha, size=input$size) +
        scale_color_manual(values=c(input$col1, input$col2,input$col3))+
        # draw line
        geom_vline(xintercept=c(-as.numeric(input$FCvalue),as.numeric(input$FCvalue)),lty=4,col="#990000",lwd=0.8) +
        geom_hline(yintercept = -log10(as.numeric(input$fdrvalue)),lty=4,col="#990000",lwd=0.8) +
        #add gene text
        geom_text_repel(
          data = iddf(),
          aes(x = iddf()[,3], y = -log10(iddf()[,2]), label = iddf()[,1]),
          colour=input$textcol,
          size = input$textsize,
          force= 20,box.padding = 1, point.padding = 1,hjust = 0.5,
          min.segment.length = 0,
          arrow = arrow(length = unit( 0.01, "npc"), type = "open", ends = "last"),
          segment.color= "grey20",segment.size= 0.5,segment.alpha= 0.8,nudge_y= 1)+
        # change labs
        labs(x="log2(fold change)",
             y="-log10 (FDR)",
             title = input$volcanotiltle)+
        theme_bw()+
        # set theme
        theme(plot.title = element_text(hjust = 0.5),
              legend.position="right",
              legend.title = element_blank())
    }
  })
  
  
  output$plot1<-renderPlot({
    if (is.null(plot1_obj()))
      return(NULL)
    plot1_obj()
    
  })
  
  output$plot1downloadData<-downloadHandler(
    filename = function() {
      paste('volcanoPlot', Sys.Date(), '.',input$extPlot, sep='')
    },
    content=function(file){
      ggsave(filename = file,plot = plot1_obj(), width = as.numeric(input$plot1width), height = as.numeric(input$plot1height), dpi = 300)
    })
  
  ######Bubble plot
  #eventReactive延迟反应式保证和按钮绑定
  bubbledf<-eventReactive(input$plot2start,{
    rawdf<-input$bubbleFile
    if (is.null(rawdf))
      return(NULL)
    go<-read.delim(file = rawdf$datapath, header=input$header3, stringsAsFactors=FALSE, as.is = TRUE, sep = input$sep3,encoding='UTF-8')
    cut_off_fdr =as.numeric(input$fdr2value)
    GOtermpos<-which(names(go)==input$b_name1)
    fdrpos<-which(names(go)==input$b_name2)
    ontpos<-which(names(go)==input$b_name3)
    
    go <- go[go[,fdrpos]<cut_off_fdr,]
    go<-go[order(go[,fdrpos],decreasing = T),]
    if (length(ontpos)!=0){
      go<-go[order(go[,ontpos]),]
      go
    }else{
      go
    }
    
  })
  
  #将输入的表格显示出来
  output$text2 <- renderTable({
    
    if (is.null(input$bubbleFile))
      return(NULL)
    raw<-input$bubbleFile
    df<-read.delim(file = raw$datapath, header=input$header2, stringsAsFactors=FALSE, as.is = TRUE, sep = input$sep2,encoding='UTF-8')
    head(df)
  })
  
  plot2_obj<-reactive({
    if (is.null(bubbledf()))
      return(NULL)
    pos1<-which(names(bubbledf())==input$b1)
    pos2<-which(names(bubbledf())==input$b2)
    pos3<-which(names(bubbledf())==input$b3)
    pos4<-which(names(bubbledf())==input$b4)
    
    
    gofold<-(as.numeric(bubbledf()[,pos1])/as.numeric(bubbledf()[,pos2]))/(as.numeric(bubbledf()[,pos3])/as.numeric(bubbledf()[,pos4]))
    
    fdrpos<-which(names(bubbledf())==input$b_name2)
    ontpos<-which(names(bubbledf())==input$b_name3)
    
    gof<--log10(bubbledf()[,fdrpos])
    
    #Category <- go$name.space
    Number <- bubbledf()[,pos1]
    
    GOtermpos<-which(names(bubbledf())==input$b_name1)
    y<-factor(bubbledf()[,GOtermpos],levels = bubbledf()[,GOtermpos])
    
    p1<-ggplot(bubbledf(),aes(gofold,y))+
      geom_point(aes(size=gof,color=Number))+
      scale_y_discrete(position = "right")+
      scale_color_gradient(low = input$collow, high = input$colhigh)+ 
      labs(color='Counts',alpha="Number of enriched genes",size="−log10(FDR)",x="Enrichment factor (fold)",y="",title="")+
      theme_bw()+theme(axis.text = element_text(color = "black",size = 14),legend.text = element_text(size = 14),legend.title=element_text(size=14),axis.title.x = element_text(size = 14))+
      scale_size_continuous(range=c(input$sizerange[1],input$sizerange[2]))
    
    if (length(ontpos)!=0){
      p2<-ggplot(bubbledf(),aes(x="",y=y,fill=bubbledf()[,ontpos]))+
        geom_tile()+
        labs(fill="Ontology")+
        theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),panel.background = element_blank(),legend.text = element_text(size = 14),legend.title=element_text(size=14))
      
      insert_left(p1,p2,width = 0.05)
    }else{
      ggplot(bubbledf(),aes(gofold,y))+
        geom_point(aes(size=gof,color=Number))+
        scale_y_discrete(position = "right")+
        scale_color_gradient(low = input$collow, high = input$colhigh)+ 
        labs(color='Counts',alpha="Number of enriched genes",size="−log10(FDR)",x="Enrichment factor (fold)",y="",title="")+
        theme_bw()+theme(axis.text = element_text(color = "black",size = 14),legend.text = element_text(size = 14),legend.title=element_text(size=14),axis.title.x = element_text(size = 14))+
        scale_size_continuous(range=c(input$sizerange[1],input$sizerange[2]))
    }
    
    
  })
  
  output$plot2<-renderPlot({
    if (is.null(plot2_obj()))
      return(NULL)
    plot2_obj()
    
  })
  
  output$plot2downloadData<-downloadHandler(
    filename = function() {
      paste('BubblePlot', Sys.Date(), '.',input$ext2Plot, sep='')
    },
    content=function(file){
      ggsave(filename = file,plot = plot2_obj(), width = as.numeric(input$plot2width), height = as.numeric(input$plot2height))
    })
  
  #ClusterProfiler
  
  cp<-eventReactive(input$ea_action,{
    if (input$ea_tools=="goea"){
      ### goea流程
      raw1<-input$genelist
      if (is.null(raw1))
        return(NULL)
      gene1 <- read.delim(file = raw1$datapath, header = input$header4,as.is = TRUE,stringsAsFactors=FALSE,encoding='UTF-8')
      raw2<-input$term2gene
      if (is.null(raw2))
        return(NULL)
      term2gene<-read.delim(file = raw2$datapath, header = input$header5,as.is = TRUE,stringsAsFactors=FALSE,sep = input$sep4,encoding='UTF-8')
      raw3<-input$term2name
      if (is.null(raw2))
        return(NULL)
      term2name<-read.delim(file = raw3$datapath, header = input$header6,as.is = TRUE,stringsAsFactors=FALSE,sep = input$sep5,encoding='UTF-8')
      gene1 <- gene1$V1[1:nrow(gene1)]
      df <- enricher(gene = gene1, TERM2GENE = term2gene, TERM2NAME = term2name, pvalueCutoff = 1, qvalueCutoff = 1)
      input<-df@result
      input<-separate(input,GeneRatio,c('Study.term','Study.total'),sep = '\\/')
      input<-separate(input,BgRatio,c('Pop.term','Pop.total'),sep = '\\/')
      go2ont <- go2ont(term2gene$GO)
      mergedf<-merge(input,go2ont,by.x = 'ID',by.y = 'go_id',all.x = T)%>%arrange(p.adjust)
      mergedf
    }
    
    else {
      ###
      raw1<-input$genelist
      if (is.null(raw1))
        return(NULL)
      gene1 <- read.delim(file = raw1$datapath, header = input$header4,as.is = TRUE,stringsAsFactors=FALSE,encoding='UTF-8')
      raw2<-input$term2gene
      if (is.null(raw2))
        return(NULL)
      term2gene<-read.delim(file = raw2$datapath, header = input$header5,as.is = TRUE,stringsAsFactors=FALSE,sep = input$sep4,encoding='UTF-8')
      raw3<-input$term2name
      if (is.null(raw2))
        return(NULL)
      term2name<-read.delim(file = raw3$datapath, header = input$header6,as.is = TRUE,stringsAsFactors=FALSE,sep = input$sep5,encoding='UTF-8')
      gene1 <- gene1$V1[1:nrow(gene1)]
      df <- enricher(gene = gene1, TERM2GENE = term2gene, TERM2NAME = term2name, pvalueCutoff = 1, qvalueCutoff = 1)
      input<-df@result
      input<-separate(input,GeneRatio,c('Study.term','Study.total'),sep = '\\/')
      input<-separate(input,BgRatio,c('Pop.term','Pop.total'),sep = '\\/')%>%arrange(p.adjust)
      input
      
    }
  })
  
  output$ea_Result <- renderDataTable({
    
    if (is.null(cp()))
      return(NULL)
    cp()
  })
  
  output$ea_downloadData<-downloadHandler(filename = function() {paste(input$prefix, input$ea_tools, ".txt", sep="")},content=function(file) {
    write.table(cp(), file,quote = F,row.names = F,sep = '\t')
  })
  
  ###设置示例数据
  genelisttest<-c(paste("AT1G0",c(101:105),"0",sep = ""))
  term2genetest<-cbind(c(paste("GO:0000",c(411:415),"0",sep = "")),genelisttest)
  term2nametest<-cbind(c(paste("GO:0000",c(411:413),"0",sep = "")),c("tRNA binding","exocyst","Golgi membrane"))
  
  
  #通过点击按钮，可以弹出界面，这个界面可以显示ui的输出界面，因此要做tagList中添加ui输出控件
  #弹窗一
  observeEvent(input$cwin1, {
    shinyalert(
      html = TRUE,
      text = tagList(
        "Genelist Data Format",
        tableOutput("example3"),
      ),
      size = "l")
  })
  
  ##关联示例控件的函数
  output$example3 <- renderTable({
    
    head(genelisttest)
  },rownames = F,colnames = F)
  
  #弹窗2及关联函数
  observeEvent(input$cwin2, {
    shinyalert(
      html = TRUE,
      text = tagList(
        "Term2gene Data Format",
        tableOutput("example4")
      ),
      size = "l")
  })
  
  output$example4 <- renderTable({
    
    head(term2genetest)
  },rownames = F,colnames = F)
  
  #弹窗3及关联函数
  observeEvent(input$cwin3, {
    shinyalert(
      html = TRUE,
      text = tagList(
        "Term2gene Data Format",
        tableOutput("example5")
      ),
      size = "l")
  })
  
  output$example5 <- renderTable({
    
    head(term2nametest)
  },rownames = F,colnames = F)
}
