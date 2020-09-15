library("shinythemes")
library("ggplot2")
library("edgeR")
library("clusterProfiler")
library("scales")
library("org.Hs.eg.db")
library("pheatmap")
library("biomaRt")
library("shiny")
library("png")
library("grid")

extract_tissue <<- function (tissue) {
  tissueSample <- Sample[Sample$SMTSD == tissue,] # Extract the whole blood tissue sample ID
  tissueSampleID <- tissueSample$SAMPID
  tissueSampleID <- gsub("-", ".", tissueSampleID)
  allSampleID <- as.character(colnames(gctdf))
  position2 <- c() # Extract blood sample position in all the expression matrix
  for (i in 1:length(tissueSampleID)) {
    for (j in 1:length(allSampleID)) {
      if (tissueSampleID[i] == allSampleID[j]) {
        position2 [j] = j
      }
    }
  }
  position2 <- position2[!is.na(position2)]
  return(position2)
}
## Function to extract DE genes
findCOVID_DE_Genes <<- function(DE.genes) {
  DE.genes <- DE.genes$table
  rownames(DE.genes) <-  gsub("\\.[0-9]*$", "", rownames(DE.genes))
  DE.genes_COVID.19 <- DE.genes[rownames(DE.genes) %in% COVID_19_ID,]
  DE.genes_COVID.19$padj <- p.adjust(DE.genes_COVID.19$PValue, method = 'fdr', n = length(COVID_19$hgnc_symbol))
  return(DE.genes_COVID.19)
}



server <- function(input, output, session){
 
  observeEvent(input$run, {
    
    withProgress(message = 'Start analysis:', value = 0, { 
      n <- 7
      
      incProgress(1/n, detail = c("Loading data..."))
      
      COVID_19 <<- read.csv(file = "data/COVID-19 related genes.csv")
      
      mart <<- useMart("ENSEMBL_MART_ENSEMBL")
      mart <<- useDataset("hsapiens_gene_ensembl", mart)
      COVID_19_name <<- as.character(COVID_19$hgnc_symbol)
      COVID_19_ID <<- getBM(attributes='ensembl_gene_id',
                            filters = 'hgnc_symbol',
                            values = COVID_19_name,
                            mart = mart)
      COVID_19_ID <<- as.data.frame(COVID_19_ID)
      COVID_19_ID <<- unique(as.character(COVID_19_ID$ensembl_gene_id))
      
      Sample <<- read.delim("data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt") # The ID here is the sample ID
      annotationFIle <<- read.delim("data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt") ## The ID here is the subject ID
     
      #tissue <- "Whole Blood"
      tissue <- input$Tissue # need to change to any tissue you like 
      
      if (tissue == "Whole Blood") {
        gctdf <<- read.delim(file = "data/gctFileBlood.txt")
      } else if (tissue == "Stomach") {
        gctdf <<- read.delim(file = "data/gctFileStomach.txt")
      } else if (tissue == "Kidney - Cortex") {
        gctdf <<- read.delim(file = "data/gctFileKidney.txt")
      } else if (tissue == "Lung") {
        gctdf <<- read.delim(file = "data/gctFileLung.txt")
      } else if (tissue == "Heart - Left Ventricle") {
        gctdf <<- read.delim(file = "data/gctFileHeartLeftVentricle.txt")
      } else if (tissue == "Heart - Atrial Appendage") {
        gctdf <<- read.delim(file = "data/gctFileHeartAtrialAppendage.txt")
      } else if (tissue == "Brain - Cortex") {
        gctdf <<- read.delim (file = "data/gctFileBrainCortex.txt") 
      } else if (tissue == "Brain - Hippocampus") {
        gctdf <<- read.delim(file = "data/gctFileBrainHippocampus.txt")
      }
      
      #Pvalue <- 0.05
      #ageSelect <- "age>60"
      
      ageSelect <<- input$ageSelect
      PValue <<- as.numeric(input$PValue)
      LogFC <<- as.numeric(input$LogFC)
      if (input$factor == "Gender") {
        factorValue = TRUE
      } else {
        factorValue = FALSE
      }
    
      
      gctFIle <- gctdf
      

      
      #gctFile= gctFIle

      analysisGct <<- function(gctFile, gctFileTPM) {
        
        gctFile <- gctFile[rowSums(gctFile) !=0,] ## Remove the 0 reads row
        #gctFileTPM <- gctFileTPM[rowSums(gctFileTPM) !=0,] ## Remove the 0 reads row
        
        ## Extract the ID from the annotation file
        AllID <- annotationFIle$SUBJID
        sex <- annotationFIle$SEX
        age <- annotationFIle$AGE
        
        #gctFile <- gctdf[,position2]
        gctFileID <- colnames(gctFile)
        
        # change sample ID to subject ID
        SampleNewID <- gsub("\\.", "-", colnames(gctFile)) # change sample ID to subject ID
        
        SampleNewID <- gsub("-[0-9][0-9][0-9][0-9]-.*", "", SampleNewID)
        
        colnames(gctFile) <- SampleNewID
        ## Extract the annotation file
        annotationTissue <- annotationFIle[as.character(AllID) %in% as.character(SampleNewID),]
        
        ## Remove the data without the annotation
        removeData <<- which(is.na(annotationTissue$DTHHRDY)) # remove the NA data
        
        ## If NA data is existed, remove it
        if (length(removeData) != 0){
          annotationTissue <- annotationTissue[-c(removeData),]
          gctFile <- gctFile[,-c(removeData)]
          #gctFileTPM <- gctFileTPM[,-c(removeData)]
        }
        
        incProgress(1/n, detail = c("Normalizing data..."))
        
        
        
        gender <- as.character(annotationTissue$SEX)
        gender[gender == 1] <- "Male"
        gender[gender == 2] <- "Female"
        dieReason <- as.character(annotationTissue$DTHHRDY)
        ageGroup <- as.character(annotationTissue$AGE)

        if (ageSelect == "age>30") {
          
          ageGroup[ageGroup != "30-39" & ageGroup != "40-49" & ageGroup != "50-59" & ageGroup != "60-69" & ageGroup != "70-79"] <- "Young"
          ageGroup[ageGroup == "30-39" | ageGroup == "40-49"| ageGroup == "50-59" | ageGroup == "60-69" | ageGroup == "70-79"] <- "Old"
          
        } else if (ageSelect == "age>40") {
          
          ageGroup[ageGroup != "40-49" & ageGroup != "50-59" & ageGroup != "60-69" & ageGroup != "70-79"] <- "Young"
          ageGroup[ageGroup == "40-49"| ageGroup == "50-59" | ageGroup == "60-69" | ageGroup == "70-79"] <- "Old"
          
        } else if (ageSelect == "age>50") {
          
          ageGroup[ageGroup != "50-59" & ageGroup != "60-69" & ageGroup != "70-79"] <- "Young"
          ageGroup[ageGroup == "50-59" | ageGroup == "60-69" | ageGroup == "70-79"] <- "Old"
          
        } else if (ageSelect == "age>60") {
          
          ageGroup[ageGroup != "60-69" & ageGroup != "70-79"] <- "Young"
          ageGroup[ageGroup == "60-69" | ageGroup == "70-79"] <- "Old"
          
        } else if (ageSelect == "age>70") {
          ageGroup[ageGroup != "70-79"] <- "Young"
          ageGroup[ageGroup == "70-79"] <- "Old"
        }
        
        dieReason[dieReason == "0" | dieReason == "N"] <- "A"
        
        dieReason[dieReason != "A"] <- "N"
        
        ## Make it global 
        gender <<- gender
        dieReason <<- dieReason
        ageGroup <<- ageGroup
        
        gctFileNormal <- gctFile[,dieReason == "N"]
            
        rownames(gctFileNormal) <-  gsub("\\.[0-9]*$", "", rownames(gctFileNormal))
            
        genderNormal <- as.factor(gender[dieReason != "A"])
            
        ageGroupNormal <- as.factor(ageGroup[dieReason != "A"])
         
        # Experiment design   
        designNormal <- model.matrix(~ genderNormal + ageGroupNormal)
            
        GroupNormal <<- factor(paste(genderNormal, ageGroupNormal,sep="."))
            
        DE.genes.Normal <- DGEList(gctFileNormal)
            
        incProgress(1/n, detail = c("calculating normal factors..."))
            
        DE.genes.Normal <- calcNormFactors(DE.genes.Normal)
            
        incProgress(1/n, detail = c("estimating binomial dispersions..."))
            
        DE.genes.Normal <- estimateDisp(DE.genes.Normal, designNormal)
            
        incProgress(1/n, detail = c("Fit linear model"))
            
        DE.genes.Normal <- glmQLFit(DE.genes.Normal, designNormal)
            
        incProgress(1/n, detail = c("Conducting genewise statistical tests..."))
        
        if (factorValue) {
          DE.genes.normal.list <<- glmQLFTest(DE.genes.Normal, coef = 2)
        } else {
          DE.genes.normal.list <<- glmQLFTest(DE.genes.Normal, coef = 3)
        }
            
            
        DE.genes.normal.list_p <- findCOVID_DE_Genes(DE.genes.normal.list)
            
        DE.genes.normal.list_p <- DE.genes.normal.list_p[DE.genes.normal.list_p$PValue < PValue,]
            
        fittedValuesRaw <<- DE.genes.Normal$fitted.values
            
        DE.genes.set <- DE.genes.normal.list_p[abs(DE.genes.normal.list_p$logFC) > LogFC,]
        
        if (nrow(DE.genes.set) == 0) {
          return(FALSE)
        } else {
          fittedValues <- DE.genes.Normal$fitted.values[rownames(DE.genes.set),]
          
          medians.normal <<- t(apply(fittedValues, 1, function (x) tapply(x, GroupNormal, median)))
          
          name <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),
                        filters = 'ensembl_gene_id',
                        values = rownames(DE.genes.set),
                        mart = mart,uniqueRows = FALSE)
          
          temp <- rownames(DE.genes.set)
          colnames(name) <- c("SYMBOL", "ID")
          GroupName.normal <- name[match(temp,name$ID),]
          
          rownames(medians.normal) <- GroupName.normal$SYMBOL
          HM.data.normal <- t(scale(t(medians.normal)))
          rownames(DE.genes.set) <- GroupName.normal$SYMBOL
          
          DE.genes.set$SYMBOL <- GroupName.normal$SYMBOL
          DE.genes.set$logFC <- round(DE.genes.set$logFC,2)
          DE.genes.set$PValue <- scientific(DE.genes.set$PValue, digits = 2)
          DE.genes.set$padj <- scientific(DE.genes.set$padj, digits = 2)
          DE.genes.set$logCPM <- scientific(DE.genes.set$logCPM, digits = 2)
          DE.genes.set$F <- round(DE.genes.set$F, 2)
          
          col_order <- c("SYMBOL", "logFC","logCPM", "F", "PValue", "padj")
          DE.genes.set <- DE.genes.set[,col_order]
          resultAll <- list(HM.data.normal, DE.genes.set)
          return(resultAll)
        }
        
            
        
      }
     
      
      result <<- analysisGct(gctFile = gctFIle, gctFileTPM = gctFIleTPM)
      if (is.logical(result)) {
        img1 <- readPNG("www/WrongGene.PNG")
        output$myHeatmap <- renderPlot({grid.raster(img1)})
        error <- c("We can't find DE genes. Please higher the threshold!")
        output$geneName <- renderTable({grid.raster(error)})
        thedata <<- reactive(geneName)
        
        
      } else {
        output$myHeatmap  <- renderPlot(
          {
            pheatmap(result[[1]], angle_col = c("0"), fontsize = 15, display_numbers = TRUE, number_color = "blue", cluster_cols = FALSE)
          }
        )
        geneName <<- as.data.frame(result[[2]])
        output$geneName <<- renderTable ({geneName})
        thedata <<- reactive(geneName)
      }
      
      
      })
    
    
    
    showNotification("Done. Please check the result.", type= "message")
    
  })
  

  
  observeEvent(input$showBoxPlot, {
    
    genderSubGroup <- gender[dieReason != "A"]
    ageSubGroup <- ageGroup[dieReason != "A"]
    subGroup <- paste0(genderSubGroup,".",ageSubGroup)
    
    DE.genes_COVID.19 <- fittedValuesRaw[rownames(fittedValuesRaw) %in% COVID_19_ID,]
    #COVID_Name <- select(EnsDb.Hsapiens.v79,
    #                                keys = rownames(DE.genes_COVID.19),
    #                                keytype = "GENEID", columns = c("SYMBOL","GENEID"))
    
    name <<- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),
                   filters = 'ensembl_gene_id',
                   values = rownames(DE.genes_COVID.19),
                   mart = mart,uniqueRows = FALSE)
    
    
    temp <- rownames(DE.genes_COVID.19)
    colnames(name) <- c("SYMBOL", "ID")
    COVID_Name <- name[match(temp,name$ID),]
    
    rownames(DE.genes_COVID.19) <- COVID_Name$SYMBOL
    
    ## The position of the gene you like to return
    geneSymbol <- input$GeneSymbol
    i <- which(rownames(DE.genes_COVID.19) %in% geneSymbol)

    tempGroup <- as.data.frame(t(as.data.frame(t(DE.genes_COVID.19[i,]))))
    tempGroup$group <- subGroup
    tempGroup$ageGroup <- ageSubGroup
    tempGroup$Gender <- genderSubGroup
    
    output$myBoxPlot  <- renderPlot(
      {
        ggplot(tempGroup, aes(x = group , y = log(tempGroup[,1]), fill = Gender)) +
          geom_boxplot() + 
          ylab(rownames(DE.genes_COVID.19)[i]) +
          scale_fill_brewer(palette = "Accent") +
          theme(axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text = element_text(size = 14),
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 16))
        
        }
    )
    


  })# close show boxplot 
  
  observeEvent(input$EnrichmentAnalysis,{
    withProgress(message = 'Generating Plot...(be patient)', value = 0,{
      n = 5
      
      # P value and sample size
      PValue2<<-as.numeric(input$PValue2)
      if (input$factor2 == "All DE genes") {
        target_gene_id <- as.vector(rownames(result[[2]]))
      } else if (input$factor2 == "DE genes expressed more in males") {
        target_gene_id <- as.vector(rownames(result[[2]]))[geneName$logFC > 0]
      } else {
        target_gene_id <- as.vector(rownames(result[[2]]))[geneName$logFC < 0]
      }

      
      target_gene_id <- as.character(mapIds(org.Hs.eg.db, target_gene_id, 'ENTREZID', 'SYMBOL'))
      
      incProgress(1/n, detail = c("GO enrichment for molecular function"))
      
      COVID_19_Background <- as.character(mapIds(org.Hs.eg.db, COVID_19_name,'ENTREZID', 'SYMBOL'))
      ## GO enrichment with clusterProfiler
      ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",
                         universe = COVID_19_Background,
                         gene = target_gene_id,
                         pvalueCutoff = PValue2,
                         ont = "MF",
                         readable=TRUE)
      if (nrow(as.data.frame(ego_MF)) > 5) {
        MF = 5
      } else {
        MF = nrow(as.data.frame(ego_MF))
      }
      ego_result_MF <- as.data.frame(ego_MF)[1:MF, ]
      # ego_result_MF <- ego_result_MF[order(ego_result_MF$Count),]
      
      incProgress(1/n, detail = c("GO enrichment for cellular component"))

      ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                         gene = target_gene_id,
                         universe = COVID_19_Background,
                         pvalueCutoff = PValue2,
                         ont = "CC",
                         readable=TRUE)
      if (nrow(as.data.frame(ego_CC)) > 5) {
        CC = 5
      } else {
        CC = nrow(as.data.frame(ego_CC))
      }
      ego_result_CC <- as.data.frame(ego_CC)[1:CC, ]
      # ego_result_CC <- ego_result_CC[order(ego_result_CC$Count),]
      incProgress(1/n, detail = c("GO enrichment for biological process"))

            
      ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                         gene = target_gene_id,
                         universe = COVID_19_Background,
                         pvalueCutoff = PValue2,
                         ont = "BP",
                         readable=TRUE)
      if (nrow(as.data.frame(ego_BP)) > 5) {
        BP = 5
      } else {
        BP = nrow(as.data.frame(ego_BP))
      }

            
      ego_result_BP <- na.omit(as.data.frame(ego_BP)[1:BP, ])
      # ego_result_BP <- ego_result_BP[order(ego_result_BP$Count),]
      incProgress(1/n, detail = c("Drawing plot"))
      if (CC == 0 & MF == 0 & BP == 0) {
        cat("No Go term is matched")
        img <- readPNG("www/wrong.PNG")
        output$myGO  <- renderPlot({grid.raster(img)})

      } else {
        if (CC == 0 & MF == 0) { 
        go_enrich_df <- data.frame(ID=c(ego_result_BP$ID),
                                   Description=c(ego_result_BP$Description),
                                   GeneNumber=c(ego_result_BP$Count),
                                   type=factor(c(rep("biological process", BP)), levels=c("biological process")))
        CPCOLS <- c("#66C3A5")#BP
        
      } else if (CC == 0 & BP == 0) { 
        go_enrich_df <- data.frame(ID=c(ego_result_MF$ID),
                                   Description=c(ego_result_MF$Description),
                                   GeneNumber=c(ego_result_MF$Count),
                                   type=factor(c(rep("molecular function", MF)), levels=c("molecular function")))
        CPCOLS <- c("#FD8D62")#MF
        
      } else if (MF == 0 & BP == 0) { 
        go_enrich_df <- data.frame(ID=c(ego_result_CC$ID),
                                   Description=c(ego_result_CC$Description),
                                   GeneNumber=c(ego_result_CC$Count),
                                   type=factor(c(rep("cellular component", CC)), levels=c("cellular component")))
        CPCOLS <- c("#8DA1CB")# CC
      } else if (CC == 0) {
        go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_MF$ID),
                                   Description=c(ego_result_BP$Description, ego_result_MF$Description),
                                   GeneNumber=c(ego_result_BP$Count, ego_result_MF$Count),
                                   type=factor(c(rep("biological process", BP), 
                                                 rep("molecular function", MF)), levels=c("biological process","molecular function")))
        CPCOLS <- c("#66C3A5","#FD8D62")
        
      } else if (MF == 0) {
        go_enrich_df <-  data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID),
                                    Description=c(ego_result_BP$Description, ego_result_CC$Description),
                                    GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count),
                                    type=factor(c(rep("biological process", BP),
                                    rep("cellular component", CC)), levels=c("biological process", "cellular component")))
        CPCOLS <- c("#66C3A5", "#8DA1CB")
        
      } else if (BP == 0) {
        go_enrich_df <- data.frame(ID=c(ego_result_CC$ID, ego_result_MF$ID),
                                   Description=c(ego_result_CC$Description, ego_result_MF$Description),
                                   GeneNumber=c(ego_result_CC$Count, ego_result_MF$Count),
                                   type=factor(c(rep("molecular function", MF),
                                   rep("cellular component", CC)), levels=c("molecular function", "cellular component")))
        CPCOLS <- c("#8DA1CB", "#FD8D62")
      } else {
        go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                                   Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                                   GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                                   type=factor(c(rep("biological process", BP), rep("cellular component", CC),
                                           rep("molecular function", MF)), levels=c("biological process","cellular component", "molecular function")))
        CPCOLS <- c( "#66C3A5", "#8DA1CB", "#FD8D62")
      }
        

                
      ## numbers as data on x axis
      go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
      ## shorten the names of GO terms
      shorten_names <- function(x, n_word=4, n_char=40){
        if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
        {
          if (nchar(x) > 40) x <- substr(x, 1, 40)
          x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                           collapse=" "), "...", sep="")
          return(x)
        } 
        else
        {
          return(x)
        }
      }
      
      labels=(sapply(
        levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
        shorten_names))
      names(labels) = rev(1:nrow(go_enrich_df))
      
      ## colors for bar // green, blue, orange
      library(ggplot2)
      myGO <<- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
        geom_bar(stat="identity", width=0.8) + coord_flip() + 
        scale_fill_manual(values = CPCOLS) + theme_bw() + 
        scale_x_discrete(labels=labels) +
        xlab("GO term") + 
        theme(axis.text=element_text(face = "bold", color="gray50", size = 15),
              legend.title = element_text(size = 15),
              legend.text  = element_text(size = 15),
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15)) +
        labs(title = "The Most Enriched GO Terms - COVID-19 Disease Genes")
      output$myGO  <- renderPlot(
        {myGO}
      )
        }
      
    })
    showNotification("Done. Already Finish the GO plot", type= "message")
    
   
    
  })
  output$downloadData1 <- downloadHandler(
    filename = function() {"GeneTable.csv"},
    content = function(file) {
      write.csv(thedata(), file, row.names = FALSE)
    }
  )

}

ui <- navbarPage("DeCovid App",
                 theme = shinytheme("cosmo"),
                 tabPanel(
                   HTML(
                     '
                      <div align="justify">
                      <h3>Tutorial Video </h3>
                      </div>
                     '
                   ),
                   HTML('<iframe width="1120" height="630" src="https://www.youtube.com/embed/vwtjcw7cJG4" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>'),
                   HTML('<div align="justify">
                          <h3>Description </h3>
This app provides information on gene expression differences between man and women and old versus young individuals for genes in the COVID19 disease map database. Gene expression data was obtained from the GTEx project, that collects expression data for multiple human tissues.
<ol style="list-style-type: decimal">
<li><p>You can select tissue of interest, companions type (man vs female or young vs old) and a p.value of differential expression to identify regulated COVID-19 genes.
<li><p>Also, you can simply introduce a gene name to see the expressions data across these population groups.
                          </div>
                        '
                   ),
                   HTML(
                     '
                      <div align="justify">
                      <h3>Pipeline </h3>
                      </div>
                     '
                   ),
                   title = "Introduction",
                   
                   mainPanel(img(src="pipeline.PNG", align = "center"))
                 ),
                 tabPanel(title='Data Analysis',
                          sidebarLayout(
                            
                            sidebarPanel(
                              h3("Step 1: Differential expression analysis"),
                              h5("*In this step, you can choose a tissue you are interested in and perform differential expression analysis to find COVID-19 disease genes with a gender difference or age difference!"),
                              selectInput("Tissue", "select the tissue you want to analysis", choices = c("Whole Blood", "Kidney - Cortex", "Stomach",
                                                                                                          "Lung", "Heart - Left Ventricle", "Heart - Atrial Appendage",
                                                                                                          "Brain - Cortex", "Brain - Hippocampus" 
                              )),
                              textInput("PValue", "P-value"),
                              textInput("LogFC", "Absolute Log FoldChange"),
                              radioButtons("factor", "select the factor you are interested in", choices = c("Gender", "Age")),
                              selectInput("ageSelect", "select the age threshold to determine old group", choices = c("age>30", "age>40", "age>50", "age>60", "age>70"),selected = c("age>60")),
                            
                              #radioButtons("dataSize", "select the data size you want to analysis", choices = c("Uniform data(Recommend)", "All data")),
                              actionButton("run", "1-RUN!"),
                              downloadButton("downloadData1", "Download Table"),
                              
                              hr(),
                              h3("Step 2: See the boxplot for COVID-19 disease genes"),
                              h5("*In this step, you can check the expression of the COVID-19 disease genes you are interested in!"),
                              selectInput("GeneSymbol", "Enter a COVID-19 disease gene to see the boxplot", choices = c(
                                "ISG15","DVL1","MMP23B","CDK11B","CDK11A","PRKCZ","TNFRSF14","H6PD","PIK3CD","CORT","MTOR","TNFRSF1B","CASP9","EPHA2","PLA2G2E","PLA2G2A","PLA2G5","PLA2G2D","PLA2G2F","PLA2G2C","EIF4G3","NBPF3","CDC42","WNT4","E2F2","IFNLR1","TRIM63","RPS6KA1","NR0B2","WDTC1","MAP3K6","RPA2","PTAFR","ATPIF1","RAB42","EPB41","EIF3I","LCK","HDAC1","AZIN2","SFPQ","ZC3H12A","BMP8A","PPIE","BMP8B","EXO5","NFYC","FOXO6","EDN2","PPIH","SLC2A1","CDC20","EIF2B3","IPP","PIK3R3","AGBL4","RAB3B","ZFYVE9","NDC1","PCSK9","PRKAA2","JUN","JAK1","PDE4B","GADD45A","PRKACB","SH3GLB1","GBP3","VAV3","GNAI3","WNT2B","NRAS","SIKE1","NGF","ATP1A1","NOTCH2","PPIAL4B","PIAS3","TXNIP","PRKAB2","PPIAL4G","PPIAL4C","MCL1","CTSS","CTSK","ARNT","PSMD4","RFX5","CREB3L4","RAB13","ADAR","PMVK","SHC1","RAB25","NES","IFI16","AIM2","ATP1A2","ATP1A4","USP21","APOA2","ATF6","RXRG","ATP1B1","KIFAP3","FASLG","TNN","TNR","RNASEL","DHX9","NCF2","EDEM3","IVNS1ABP","TPR","PDC","PLA2G4A","UCHL5","ELF3","REN","CNTN2","CDK18","RAB29","RAB7B","IKBKE","RASSF5","IL10","IL19","CD46","IRF6","TRAF5","ATF3","BATF3","TGFB2","DUSP10","TLR5","WNT9A","WNT3A","ARF1","RAB4A","ACTA1","LYST","ACTN2","AKT3","NLRP3","RSAD2","ID2","ADAM17","ROCK2","E2F6","GEN1","GDF7","POMC","RAB10","AGBL5","CAD","EIF2B4","FOSL2","NLRC4","EIF2AK2","PRKD3","SOS1","CDKL4","MAP4K3","PRKCE","EPAS1","SOCS5","REL","XPO1","RAB1A","SPRED2","PPP3R1","BMP10","EGR4","ACTG2","CTNNA2","CD8A","CD8B","EIF2AK3","MAL","ZAP70","EIF5B","RFX8","MAP4K4","IL1R1","BCL2L11","ANAPC1","IL1A","IL1B","IL1RN","RABL2A","PROC","RAB6C","TUBA3E","TUBA3D","CXCR4","ZEB2","NMI","PRPF40A","TANK","DPP4","IFIH1","SSB","PDK1","pk","ATF2","PDE11A","ITGA4","PDE1A","GLS","STAT1","STAT4","HSPE1","PPIL3","CASP10","CASP8","CDK15","FZD7","ICOS","EEF1B2","GPR1","CPO","CREB1","FZD5","RPE","FN1","XRCC5","PLCD4","TTLL4","PRKAG3","WNT6","WNT10A","TUBA4A","IRS1","CCL20","SP100","CAB39","INPP5D","ATG16L1","SAG","RAB17","HDAC4","AGXT","PDCD1","EDEM1","ATG7","RAF1","HDAC11","WNT7A","NR2C2","RAB5A","THRB","EOMES","PDCD6IP","PLCD1","MYD88","CTNNB1","XCR1","CCR1","CDC25A","ATRIP","PRKAR2A","USP4","RHOA","AMT","DAG1","RNF123","GNAT1","GNAI2","CISH","TLR9","PBRM1","PRKCD","WNT5A","PTPRG","TOMM70A","SENP7","ATG3","CD80","GSK3B","RABL3","KPNA1","RAB7A","RAB43","PIK3R4","TOPBP1","RAB6B","PIK3CB","ATP1B3","XRN1","ATR","CPB1","CP","SIAH2","IL12A","PRKCI","CLDN11","TNFSF10","PIK3CA","EIF2B5","DVL3","EIF4G1","MAP3K13","SENP2","ETV5","EIF4A2","CLDN1","FGF12","LRRC15","GP5","RNF168","DLG1","PDE6B","FGFRL1","FGFR3","LRPAP1","RAB28","PPARGC1A","TLR10","TLR1","TLR6","UBE2K","RHOH","UCHL1","TXK","TEC","PDGFRA","KDR","GC","ALB","CXCL8","CXCL6","CXCL1","PF4","PPBP","CXCL5","CXCL3","CXCL2","CDKL2","G3BP2","CXCL9","CXCL10","CXCL11","CXCL13","FGF5","BMP3","MAPK10","PTPN13","SPP1","HERC5","MTTP","PPP3CA","NFKB1","UBE2D3","AIMP1","SGMS2","LEF1","CASP6","PLA2G12A","EGF","CAMK2D","IL2","IL21","FGF2","PLK4","ELF2","RAB33B","GYPB","FBXW7","TLR2","FGB","FGA","FGG","PPID","KLHL2","MSMO1","CPE","DDX60","VEGFC","IRF2","CASP3","TLR3","SLC9A3","NKD2","DROSHA","RICTOR","PRKAA1","C6","PAIP1","FGF10","DHX29","MAP3K1","RAB3C","PDE4D","CWC27","PPWD1","PIK3R1","CDK7","OCLN","HMGCR","PDE8B","SERINC5","RASGRF2","ATG10","XRCC4","EDIL3","RASA1","PAM","FER","TICAM2","ATG12","PPIC","IL3","CSF2","IRF1","IL4","CDKL3","CXCL14","WNT8A","EGR1","ETF1","CTNNA1","NRG2","EIF4EBP3","CD14","IK","HDAC3","RNF14","FGF1","SPINK1","PDE6A","PDGFRB","CAMK2A","CD74","TNIP1","G3BP1","ITK","IL12B","FGF18","DUSP1","CLTB","FGFR4","RAB24","DDX41","ADAMTS2","CANX","MAPK9","FLT4","MGAT1","TRIM41","IRF4","FOXC1","RIPK1","TUBB2A","TUBB2B","DSP","BMP6","EEF1E1","PAK1IP1","MAK","SIRT5","GMPR","NUP153","E2F3","TRIM38","TRIM27","UBD","TRIM40","TRIM26","MDC1","TUBB","LTA","TNF","BAG6","CSNK2B","HSPA1B","C2","C4B","TNXB","NOTCH4","TAP2","PSMB8","TAP1","RXRB","DAXX","BAK1","FANCE","SRPK1","MAPK14","MAPK13","RAB44","PPIL1","PIM1","RNF8","NFYA","TREM2","MEA1","SRF","VEGFA","HSP90AB1","PLA2G7","IL17A","IL17F","BMP5","BAG2","RAB23","MB21D1","EEF1A1","UBE2J1","ATG5","AIM1","FOXO3","PPIL6","WASF1","CDK19","FYN","TUBE1","HDAC2","RFX6","PDE7B","MAP3K5","IFNGR1","TNFAIP3","UTRN","RAB32","TAB2","PPIL4","ESR1","LPA","PLG","MAP3K4","PARK2","PDE10A","MPC1","RPS6KA2","MLLT4","PRKAR1B","GET4","EIF3B","GNA12","CARD11","ACTB","EIF2AK1","RAC1","RPA3","AHR","HDAC9","IL6","CREB5","NOD1","PDE1C","RALA","CDK13","GCK","CAMK2B","NPC1L1","PPIA","UPP1","DDC","EGFR","TRIM50","FZD9","MLXIPL","EIF4H","LAT2","NCF1","HSPB1","YWHAG","GNAI1","CD36","GNAT3","CROT","CDK14","FZD1","KRIT1","CDK6","PDK4","ASNS","TRRAP","TRIM4","EPO","ACHE","TRIM56","SERPINE1","IFT22","SH2B2","SLC26A5","KMT2E","SRPK2","PIK3CG","PRKAR2B","CAV1","WNT2","CFTR","WNT16","GPR37","IMPDH1","IRF5","SMO","CREB3L2","TRIM24","PARP12","RAB19","BRAF","MGAM","CASP2","TPK1","NOS3","CDK5","FASTK","RHEB","PRKAG2","SHH","NCAPG2","CTSB","FGF20","ATP6V1B2","FGF17","BMP1","PPP3CC","PDLIM2","EGR3","TNFRSF10B","TNFRSF10C","TNFRSF10D","TNFRSF10A","NEFL","GNRH1","PTK2B","CLU","FZD3","DUSP4","NRG1","RNF122","EIF4EBP1","FGFR1","PLAT","IKBKB","CEBPD","PRKDC","SNAI2","RB1CC1","TCEA1","LYN","MOS","PENK","TOX","RAB2A","PDE7A","LY96","E2F5","CA2","WWP1","MMP16","RIPK2","GEM","GDF6","RNF19A","PABPC1","YWHAZ","FZD6","EIF3E","EIF3H","SNTB1","DERL1","MYC","TG","AGO2","PTK2","EEF1D","SCRIB","GPT","RFX3","JAK2","CD274","IFNB1","IFNW1","IFNA21","IFNA4","IFNA7","IFNA10","IFNA16","IFNA17","IFNA14","IFNA5","IFNA6","IFNA13","IFNA2","IFNA8","IFNA1","IFNE","CDKN2A","CDKN2B","TEK","IFNK","DDX58","DNAJA1","CREB3","CLTA","PRKACG","HNRNPK","DAPK1","CDK20","SHC3","SYK","PTCH1","TRIM14","TGFBR1","NR4A3","BAAT","PPP3R2","TNC","TLR4","TRAF1","C5","RAB14","GSN","HSPA5","MAPKAP1","CDK9","ENG","CRAT","ASS1","TSC1","RALGDS","CACFD1","VAV2","RXRA","NOTCH1","TRAF2","NPDC1","TUBB4B","TUBAL3","PRKCQ","GATA3","OPTN","DCLRE1C","VIM","STAM","ARL5B","RAB18","CREM","FZD8","CXCL12","GDF2","MAPK8","SGMS1","DKK1","MBL2","CDK1","EGR2","CTNNA3","SIRT1","TYSND1","EIF4EBP2","PRF1","DDIT4","PLA2G12B","PPP3CB","PLAU","VCL","PPIF","SFTPA2","SFTPA1","NRG3","GHITM","PTEN","ACTA2","FAS","IFIT2","IFIT3","IFIT1","PDE6C","PLCE1","BLNK","PI4K2A","CHUK","WNT8B","BTRC","FGF8","GBF1","NFKB2","TRIM8","DUSP5","CASP7","EIF3A","BAG3","FGFR2","HMX2","MMP21","SIRT3","IFITM2","IFITM1","IFITM3","HRAS","IRF7","PIDD1","TOLLIP","DUSP8","INS","TH","CD81","ART1","NUP98","RHOG","TRIM6","EIF3F","ARNTL","PDE3B","TSG101","CSRP3","E2F8","BDNF","PRRG4","CAT","APIP","CD44","TRAF6","CREB3L1","AMBRA1","ATG13","ARHGAP1","RTN4RL2","CTNND1","CD6","DDB1","DAK","EEF1G","COX8A","VEGFB","PLCB3","BAD","TRMT112","RPS6KA4","GPHA2","MRPL49","MAP3K11","RELA","KAT5","CFL1","FOSL1","RAB1B","ACTN3","PC","PPP1CA","RPS6KB2","LRP5","PPP6R3","GAL","TPCN2","FGF19","FGF4","FGF3","FADD","INPPL1","PDE2A","RAB6A","ARRB1","UVRAG","WNT11","RAB30","FZD4","RAB38","CCDC67","MRE11A","BIRC3","MMP7","MMP20","MMP27","MMP8","MMP10","MMP1","MMP3","MMP12","MMP13","CASP12","CASP4","CASP5","CASP1","CARD16","RAB39A","ATM","IL18","PTS","APOA5","APOA4","APOA1","PCSK7","ARCN1","DDX6","HYOU1","HMBS","CBL","CHEK1","TIRAP","ETS1","JAM3","RAD52","WNT5B","FGF23","FGF6","CD9","TNFRSF1A","LAG3","CD4","USP5","PTPN6","PHB2","GDF3","C3AR1","CLEC7A","ETV6","LRP6","DUSP16","CDKN1B","PLCZ1","PDE3A","ETNK1","KRAS","IRAK4","RAPGEF3","HDAC7","WNT10B","WNT1","PRKAG1","TUBA1B","TUBA1A","TUBA1C","TMBIM6","FAIM2","ATF1","KRT8","EIF4B","SP1","MAP3K12","ATF7","NFE2","PDE1B","MMP19","CDK2","RAB5B","IL23A","STAT2","STAT6","LRP1","GLI1","DDIT3","OS9","CDK4","USP15","TBK1","IFNG","IL22","MDM2","RAB21","KRR1","E2F7","KITLG","SOCS2","CRADD","CDK17","APAF1","CHPT1","DRAM1","IGF1","HSP90B1","NFYB","RFX4","SART3","DAO","ACACB","ATP2A2","PPP1CC","MAPKAPK5","ERP29","PTPN11","OAS1","OAS3","OAS2","TPCN1","NOS1","PRKAB1","CIT","RAB35","PXN","SIRT4","PLA2G1B","DYNLL1","OASL","CAMKK2","EIF2B1","SCARB1","UBC","FZD10","MMP17","ULK1","TUBA3C","FGF9","SGCG","CDK8","FLT1","HMGB1","BRCA2","FOXO1","TNFSF11","RB1","TRIM13","UCHL3","DNAJC3","STK24","FGF14","TNFSF13B","RAB20","TFDP1","ATP4B","RAB2B","MMP14","CEBPE","BCL2L2","IRF9","NFATC4","GZMB","CFL2","NFKBIA","SOS2","CDKL1","BMP4","ATG14","TBPL2","GPR135","PPM1A","PRKCH","HIF1A","RAB15","MAX","EIF2S1","ACTN1","SLC10A1","ZFYVE1","PSEN1","PGF","EIF2B2","FOS","BATF","C14orf1","TGFB3","PTPN21","RPS6KA5","SERPINA1","DICER1","HSP90AA1","TRAF3","EIF5","BAG5","C14orf2","AKT1","ACTC1","SPRED1","RASGRP1","THBS1","EIF2AK4","PLCB2","ITPKA","PLA2G4D","TP53BP1","PPIP5K1","PDIA3","MFAP1","EIF3J","B2M","SHC4","FGF7","USP8","RAB27A","NEDD4","RFX7","TPM1","RAB8B","USP3","PPIB","RAB11A","MAP2K1","SMAD3","MAP2K5","PIAS1","KIF23","RPLP1","PML","MPI","SIN3A","NRG4","RASGRF1","ARNT2","PDE8A","ISG20","FURIN","ARRDC4","IGF1R","AXIN1","RAB40C","STUB1","MAPK8IP3","ZNF598","SLC9A3R2","TSC2","RAB26","TRAF7","MLST8","MMP25","IL32","TRAP1","PPL","USP7","SOCS1","PLA2G10","EEF2K","PLK1","PRKCB","EIF3C","MVP","MAPK3","ITGAL","PYCARD","GPT2","SIAH1","NKD1","NOD2","CYLD","MMP2","GNAO1","AMFR","MT2A","NLRC5","MMP15","CSNK2A2","NAE1","TRADD","E2F4","CTRL","NFATC3","CDH1","VPS4A","NFAT5","NQO1","TERF2IP","MAF","PLCG2","MBTPS1","IRF8","CDK10","TUBB3","CRK","SERPINF2","RPA1","PAFAH1B1","CLUH","ARRB2","C1QBP","MED31","SLC13A5","DVL2","FGF11","EIF4A1","SHBG","ATP1B2","TP53","PIK3R6","PIK3R5","MAP2K4","ULK2","SPECC1","MAP2K3","KCNJ12","NOS2","NLK","SARM1","RAB34","TRAF4","RNF135","CCL2","SLFN11","MMP28","CCL5","CCL16","CCL4","ACACA","CDK12","THRA","CCR7","EIF1","GAST","DHX58","RAB5C","STAT5B","STAT5A","STAT3","TUBG1","TUBG2","IFI35","BRCA1","HDAC5","FZD2","PLCD3","MAPT","WNT3","WNT9B","KPNB1","SP2","TTLL6","CALCOCO2","PHB","PDK2","HLF","TRIM25","SCPEP1","CLTC","TUBD1","RPS6KB1","KCNH6","MAP3K3","STRADA","CSHL1","ICAM2","ERN1","GNA13","AXIN2","PRKCA","PRKAR1A","MAP2K6","RAB37","SLC9A3R1","SUMO2","GRB2","UNK","TRIM65","CDK3","FOXJ1","SOCS3","RPTOR","ACTG1","SIRT7","RAC3","FASN","RAB40B","USP14","YES1","RAB12","RAB31","VAPA","GNAL","TUBB6","PTPN2","ROCK1","RIOK3","RNF125","PIK3C3","SMAD2","SMAD7","SMAD4","RAB27B","TCF4","PHLPP1","BCL2","SOCS6","NFATC1","SHC2","CDC34","BSG","FGF22","ELANE","SBNO2","STK11","TCF3","SGTA","EEF2","MAP2K2","CREB3L3","SIRT6","RFX2","TUBB4A","C3","SH2D3A","VAV1","INSR","CD209","MAP2K7","RAB11B","EIF3G","ICAM1","ICAM5","ICAM3","TYK2","PDE4A","KEAP1","CDKN2D","LDLR","RAB3D","ECSIT","JUNB","CALR","TRMT1","RFX1","PRKACA","GIPC1","TECR","CASP14","NOTCH3","RAB8A","BST2","JAK3","IFI30","RAB3A","PDE4C","JUND","LRRC25","ELL","CEBPA","CEBPG","USF2","HCST","TYROBP","SPRED3","ACTN4","SIRT2","IFNL3","IFNL4","IFNL2","AKT2","RAB4B","TGFB1","CD79A","ATP1A3","GSK3A","ERF","LIPE","APOE","RELB","FOSB","VASP","HIF3A","NPAS1","SAE1","C5AR1","PLA2G4C","RPL18","DBP","FGF21","PPP1R15A","BAX","IRF3","NUP62","PRKCG","TFPT","LILRA4","KIR3DL2","NLRP2","TRIM28","ZCCHC3","TRIB3","RBCK1","CSNK2A1","PSMF1","MAVS","BMP2","PLCB1","PLCB4","MKKS","PCSK2","DSTN","BCL2L1","SNTA1","E2F1","EIF2S2","MAP1LC3A","EDEM2","MMP24","EIF6","GDF5","SAMHD1","SRC","LBP","PLCG1","FITM2","HNF4A","SERINC3","MMP9","CD40","UBE2V1","CEBPB","PTPN1","NFATC2","BMP7","ZBP1","RAB22A","GNAS","NELFCD","TUBB1","ATP5E","ADRM1","EEF1A2","JAM2","APP","IFNAR2","IL10RB","IFNAR1","IFNGR2","ERG","ETS2","MX2","MX1","UBASH3A","PDE9A","ICOSLG","SUMO3","ITGB2","BID","TUBA8","USP18","DGCR8","ZNF74","PPIL2","MAPK1","GNAZ","RAB36","BCR","MMP11","CHEK2","XBP1","NF2","PLA2G3","RNF185","HMOX1","EIF3D","NCF4","RAC2","PDXP","PLA2G6","TAB1","ATF4","ST13","EP300","XRCC6","NHP2L1","SREBF2","WNT7B","PPARA","HDAC10","MAPK12","MAPK11","RABL2B","CD99","STS","TLR7","TLR8","RAB9A","FIGF","ACE2","CDKL5","PHKA2","EIF1AX","RPS6KA3","MBTPS2","EIF2S3","PDK3","TAB3","DMD","CYBB","USP9X","DDX3X","CASK","UBA1","CDK16","ARAF","TIMP1","WAS","SUV39H1","HDAC6","BMP15","GSPT2","RAB41","FOXO4","MED12","HDAC8","FGF16","ATRX","RPS6KA6","NOX1","BTK","RAB40AL","RAB40A","PLP1","RAB9B","ATP1B4","ELF4","AIFM1","RAB33A","FGF13","AVPR2","IRAK1","OPN1LW","FLNA","UBL4A","G6PD","IKBKG","RAB39B","EIF1AY","MT-ATP8","MT-ATP6"
                              ),multiple = FALSE),
                              
                              actionButton("showBoxPlot", "2-Show Boxplot"),
                              hr(),
                              h3("Step 3: GO Enrichment Analysis"),
                              h5("*In this step, you can perform a GO enrichment analysis for DE genes found in step1 under the COVID-19 diseases genes background!"),
                              textInput("PValue2", "P-value"),
                              radioButtons("factor2", "select the gene size you want to analysis", choices = c("All DE genes", "DE genes expressed more in males", "DE genes expressed more in females")),
                              actionButton("EnrichmentAnalysis", "3-GO enrichment analysis")
                              
                            ),
                            mainPanel(
                              tabsetPanel(type = "tab",
                                          tabPanel("Gene name",tableOutput("geneName")),
                                          tabPanel("Heatmap", plotOutput("myHeatmap",width = "800px", height = "1200px")),
                                          tabPanel("Boxplot", plotOutput("myBoxPlot", width = "800px", height = "800px")),
                                          tabPanel("Enrichment Analysis", plotOutput("myGO", width = "800px", height = "800px"))
                              )
                              
                              
                            )
                            
                          )))






shinyApp(ui = ui, server = server)









