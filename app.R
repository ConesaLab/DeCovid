library("shinythemes")
library("ggplot2")
library("edgeR")
library("clusterProfiler")
library("scales")
library("org.Hs.eg.db")
library("pheatmap")
library("shiny")
library("shinydashboard")
library("DT")
library("png")
library("grid")
library("AnnotationDbi")

COVID_19_1 <<- read.csv(file = "data/COVID-19 related genes.csv")
COVID_19_2 <<- read.csv(file = "data/COVID-19 Disease Map genes.csv")
COVID_19_1 <<- COVID_19_1$hgnc_symbol
COVID_19_2 <<- COVID_19_2$hgnc_symbol
COVID_19_3 <<- union(COVID_19_1, COVID_19_2)
COVID_19_name <<- c()

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
  DE.genes_COVID.19$padj <- p.adjust(DE.genes_COVID.19$PValue, method = 'fdr', n = length(COVID_19))
  return(DE.genes_COVID.19)
}



server <- function(input, output, session){
  observeEvent(input$tabs == "introduction", {
    if (req(input$tabs) == "DataAnalysis") {
      
      values <- reactiveValues()
      values$successStep1 <- FALSE
      observeEvent(input$run, {
        
        withProgress(message = 'Start analysis:', value = 0, { 
          n <- 7
          
          incProgress(1/n, detail = c("Loading data..."))
          
          backgroundGene <<- input$background
          if (backgroundGene == "All") {
            COVID_19 <<- COVID_19_3
          } else if (backgroundGene == "COVID-19 Disease Map genes") {
            COVID_19 <<- COVID_19_2
          } else {
            COVID_19 <<- COVID_19_1
          }
          
          COVID_19_ID <- as.character(mapIds(org.Hs.eg.db,
                                             keys=as.character(COVID_19),
                                             column="ENSEMBL",
                                             keytype="SYMBOL",
                                             multiVals="first"))
          
          COVID_19_ID <<- as.character(COVID_19_ID)[!is.na(COVID_19_ID)]
          #mart <<- useMart("ENSEMBL_MART_ENSEMBL")
          #mart <<- useDataset("hsapiens_gene_ensembl", mart)
          COVID_19_name <<- as.character(COVID_19)
          updateSelectInput(session, "GeneSymbol", choices = COVID_19_name)
          
          #COVID_19_ID <<- getBM(attributes='ensembl_gene_id',
          #                      filters = 'hgnc_symbol',
          #                      values = COVID_19_name,
          #                      mart = mart)
          #COVID_19_ID <<- as.data.frame(COVID_19_ID)
          #COVID_19_ID <<- unique(as.character(COVID_19_ID$ensembl_gene_id))
          
          Sample <<- read.delim("data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt") # The ID here is the sample ID
          annotationFIle <<- read.delim("data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt") ## The ID here is the subject ID
          
          #tissue <- "Whole Blood"
          tissue <<- input$Tissue # need to change to any tissue you like 
          
          if (tissue == "Whole Blood") {
            gctdf <<- read.delim(file = "data/gctFileBlood.txt")
          } else if (tissue == "Spleen") {
            gctdf <<- read.delim(file = "data/Spleen.txt")
          } else if (tissue == "lymphocytes") {
            gctdf <<- read.delim(file = "data/lymphocytes.txt")
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
            if (tissue == "Whole Blood")  {
              gctFileNormal <- gctFile[,dieReason == "N"]
              rownames(gctFileNormal) <-  gsub("\\.[0-9]*$", "", rownames(gctFileNormal))
              genderNormal <- as.factor(gender[dieReason != "A"])
              ageGroupNormal <- as.factor(ageGroup[dieReason != "A"])
            } else {
              #gctFileNormal <- gctFile[,dieReason == "N"]
              gctFileNormal <- gctFile
              rownames(gctFileNormal) <-  gsub("\\.[0-9]*$", "", rownames(gctFileNormal))
              #genderNormal <- as.factor(gender[dieReason != "A"])
              genderNormal <- as.factor(gender)
              #ageGroupNormal <- as.factor(ageGroup[dieReason != "A"])
              ageGroupNormal <- as.factor(ageGroup)
            }
              
            
            
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
            
            DE.genes.set <<- DE.genes.normal.list_p[abs(DE.genes.normal.list_p$logFC) > LogFC,]
            
            if (nrow(DE.genes.set) < 2) {
              return(FALSE)
            } else {
              fittedValues <- DE.genes.Normal$fitted.values[rownames(DE.genes.set),]
              
              medians.normal <<- t(apply(fittedValues, 1, function (x) tapply(x, GroupNormal, median)))
              
              GroupName.normal <- as.character(mapIds(org.Hs.eg.db,
                                                      keys=rownames(DE.genes.set),
                                                      column="SYMBOL",
                                                      keytype="ENSEMBL",
                                                      multiVals="first"))
              
              #name <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),
              #              filters = 'ensembl_gene_id',
              #              values = rownames(DE.genes.set),
              #              mart = mart,uniqueRows = FALSE)
              
              #temp <- rownames(DE.genes.set)
              #colnames(name) <- c("SYMBOL", "ID")
              #GroupName.normal <- name[match(temp,name$ID),]
              
              rownames(medians.normal) <- GroupName.normal
              HM.data.normal <- t(scale(t(medians.normal)))
              rownames(DE.genes.set) <- GroupName.normal
              
              DE.genes.set$SYMBOL <- GroupName.normal
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
            output$geneName <- renderDataTable({grid.raster(error)})
            thedata <<- reactive(geneName)
            
            
          } else {
            output$myHeatmap  <- renderPlot(
              {
                pheatmap(result[[1]], angle_col = c("0"), fontsize_col = 14 ,fontsize = 9,cluster_cols = FALSE)
              }
            )
            geneName <<- as.data.frame(result[[2]])
            rownames(geneName) <- NULL
            output$geneName <<- renderDataTable({geneName})
            thedata <<- reactive(geneName)
          }
          
          
        })
        
        
        values$successStep1<-TRUE
        showNotification("Done. Please check the result.", type= "message")
        
      })
      
      
      
      observeEvent(input$showBoxPlot, {
        shiny::validate(
          need( values$successStep1==TRUE, message = ('Missing succesful Step 1')),
          errorClass =  showNotification("Missing succesful Step 1", type= "error")
        )
        
        if (tissue == "Whole Blood") {
          genderSubGroup <- gender[dieReason != "A"]
          ageSubGroup <- ageGroup[dieReason != "A"]
        } else {
          genderSubGroup <- gender
          ageSubGroup <- ageGroup
        }
       
        subGroup <- paste0(genderSubGroup,".",ageSubGroup)
        
        DE.genes_COVID.19 <- fittedValuesRaw[rownames(fittedValuesRaw) %in% COVID_19_ID,]
        #COVID_Name <- select(EnsDb.Hsapiens.v79,
        #                                keys = rownames(DE.genes_COVID.19),
        #                                keytype = "GENEID", columns = c("SYMBOL","GENEID"))
        
        #name <<- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),
        #               filters = 'ensembl_gene_id',
        #               values = rownames(DE.genes_COVID.19),
        #               mart = mart,uniqueRows = FALSE)
        #
        
        #temp <- rownames(DE.genes_COVID.19)
        #colnames(name) <- c("SYMBOL", "ID")
        #COVID_Name <- name[match(temp,name$ID),]
        
        COVID_Name <- as.character(mapIds(org.Hs.eg.db,
                                          keys=rownames(DE.genes_COVID.19),
                                          column="SYMBOL",
                                          keytype="ENSEMBL",
                                          multiVals="first"))
        
        
        rownames(DE.genes_COVID.19) <- COVID_Name
        
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
        shiny::validate(
          need( values$successStep1==TRUE, message = ('Missing succesful Step 1')),
          errorClass =  showNotification("Missing succesful Step 1", type= "error")
        )
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
            go_enrich_df$Description <- as.factor(go_enrich_df$Description)
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
              labs(title = "The Most Enriched GO Terms - COVID-19 Disease Map Genes")
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
    else if (req(input$tabs) == "AdvancedAnalysis") {
      
      values_A <- reactiveValues()
      values_A$successStep1 <- FALSE
      values_A$successStep2 <- FALSE
      
      observeEvent(input$addGroupInf, {
        
        controlGroup <<- input$controlGroup
        experimentalGroup <<- input$experimentalGroup
        shiny::validate(
          need(all(!controlGroup%in%experimentalGroup) == TRUE, message = ('A Group is existed in both control and experimental group!')),
          errorClass = showNotification("Group existed in both control and experimental group!", type= "error")
        )
        
        if (input$factor_A == "Gender") {
          factorValue <<- TRUE
        } else {
          factorValue <<- FALSE
        }
        if (factorValue == TRUE) {
          tempControl <- gsub("\ .*", "", controlGroup)
          tempExp <- gsub("\ .*", "", experimentalGroup)
          shiny::validate(
            need(all(!tempControl%in%tempExp) == TRUE, message = ('Group Wrong: Please check if experimental/control group have same gender')),
            errorClass =  showNotification("Group Wrong: Please check if experimental/control group have same gender!", type= "error")
          )
        } else {
          tempControl <<-gsub(".* ", "", controlGroup)
          tempExp <<- gsub(".* ", "", experimentalGroup)
          shiny::validate(
            need(all(!tempControl%in%tempExp) == TRUE, message = ('Group Wrong: Please check if experimental/control group have same age')),
            errorClass = showNotification('Group Wrong: Please check if experimental/control group have same age', type = "error")
          )
        } 
        ageIndexC <- paste(gsub(".* ", "", controlGroup), collapse = " ")
        matches <- regmatches(ageIndexC, gregexpr("[[:digit:]]+", ageIndexC))
        ageIndexC <- max(as.numeric(unlist(matches)))
        
        ageIndexE <-  gsub(".* ", "", experimentalGroup)
        matches <- regmatches(ageIndexE, gregexpr("[[:digit:]]+", ageIndexE))
        ageIndexE <-  max(as.numeric(unlist(matches)))
        
        if (ageIndexC > ageIndexE) {
          ageIndexC <- paste(gsub(".* ", "", controlGroup), collapse = " ")
          matches <- regmatches(ageIndexC, gregexpr("[[:digit:]]+", ageIndexC))
          ageIndexC <<- min(as.numeric(unlist(matches)))
        } else {
          ageIndexE <-  gsub(".* ", "", experimentalGroup)
          matches <- regmatches(ageIndexE, gregexpr("[[:digit:]]+", ageIndexE))
          ageIndexE <<-  min(as.numeric(unlist(matches)))
        }

        ageSelect <<- input$ageSelect_A
        
        ageIndexS <- as.numeric(unlist(regmatches(ageSelect, gregexpr("[[:digit:]]+", ageSelect))))
        
        PValue <<- as.numeric(input$PValue_A)
        LogFC <<- as.numeric(input$LogFC_A)

        showNotification("Success!", type = "message" )
        values_A$successStep1 <- TRUE
        
        
      })
      observeEvent(input$run_A, {
        
        shiny::validate(
          need(values_A$successStep1, message = ('Missing experimental/control Group')),
          errorClass = showNotification('Missing experimental/control Group', type = "error")
        )
        
        
        
        if (factorValue == FALSE) {
          shiny::validate(
            need(ageIndexC >= ageIndexS && ageIndexS >= ageIndexE || ageIndexC <= ageIndexS && ageIndexS <= ageIndexE,
                 message = ('Old threshold should between the expermintal group and control group')),
            errorClass = showNotification('Please select a correct old people threshold', type = "error")
          )
        }
        
        withProgress(message = 'Start analysis:', value = 0, { 
          n <- 7
          backgroundGene <<- input$background_A
          if (backgroundGene == "All") {
            COVID_19 <<- COVID_19_3
          } else if (backgroundGene == "COVID-19 Disease Map genes") {
            COVID_19 <<- COVID_19_2
          } else {
            COVID_19 <<- COVID_19_1
          }

          incProgress(1/n, detail = c("Loading data..."))
          COVID_19_ID <- as.character(mapIds(org.Hs.eg.db,
                                             keys=as.character(COVID_19),
                                             column="ENSEMBL",
                                             keytype="SYMBOL",
                                             multiVals="first"))
          
          COVID_19_ID <<- as.character(COVID_19_ID)[!is.na(COVID_19_ID)]
          #mart <<- useMart("ENSEMBL_MART_ENSEMBL")
          #mart <<- useDataset("hsapiens_gene_ensembl", mart)
          COVID_19_name <<- as.character(COVID_19)
          updateSelectInput(session, "GeneSymbol_A", choices = COVID_19_name)
          ?updateSelectInput
          #COVID_19_ID <<- getBM(attributes='ensembl_gene_id',
          #                      filters = 'hgnc_symbol',
          #                      values = COVID_19_name,
          #                      mart = mart)
          #COVID_19_ID <<- as.data.frame(COVID_19_ID)
          #COVID_19_ID <<- unique(as.character(COVID_19_ID$ensembl_gene_id))
          
          Sample <<- read.delim("data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt") # The ID here is the sample ID
          annotationFIle <<- read.delim("data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt") ## The ID here is the subject ID
          
          #tissue <- "Whole Blood"
          tissue <- input$Tissue_A # need to change to any tissue you like 
          
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
          } else if (tissue == "Spleen") {
            gctdf <<- read.delim(file = "data/Spleen.txt")
          } else if (tissue == "lymphocytes") {
            gctdf <<- read.delim(file = "data/lymphocytes.txt")
          } 
          
          #Pvalue <- 0.05
          #ageSelect <- "age>60"
          gctFIle <- gctdf
          
          #gctFile <- gctFIle
          analysisGct_A <<- function(gctFile, gctFileTPM) {
            
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
            
            if (tissue == "Whole Blood")  {
              gctFileNormal <- gctFile[,dieReason == "N"]
              rownames(gctFileNormal) <-  gsub("\\.[0-9]*$", "", rownames(gctFileNormal))
              genderNormal <- as.factor(gender[dieReason != "A"])
              ageGroupNormal <- as.factor(ageGroup[dieReason != "A"])
            } else {
              #gctFileNormal <- gctFile[,dieReason == "N"]
              gctFileNormal <- gctFile
              rownames(gctFileNormal) <-  gsub("\\.[0-9]*$", "", rownames(gctFileNormal))
              #genderNormal <- as.factor(gender[dieReason != "A"])
              genderNormal <- as.factor(gender)
              #ageGroupNormal <- as.factor(ageGroup[dieReason != "A"])
              ageGroupNormal <- as.factor(ageGroup)
            }
            
            ## Extract User Select Group
            expFactor1 <- gsub(" .*","",gsub("sex:","",experimentalGroup))
            if (expFactor1 == "male") {
              expFactor1 = "Male"
            } else {
              expFactor1 = "Female"
            }
            expFactor2 <- gsub(".*:", "",experimentalGroup)
            
            contFactor1 <- gsub(" .*","",gsub("sex:","",controlGroup))
            if (contFactor1 == "male"){
              contFactor1 = "Male"
            } else {
              contFactor1 = "Female"
            }
            contFactor2 <- gsub(".*:", "",controlGroup)
            if (as.numeric(gsub("..-","",expFactor2)) - as.numeric(gsub("...>","",ageSelect)) <= 0) {
              expFactor2 = c("Young")
            } else {
              expFactor2 = c("Old")
            }
            
            if (as.numeric(gsub("..-","",contFactor2)) - as.numeric(gsub("...>","",ageSelect)) <= 0) {
              contFactor2 = c("Young")
            } else {
              contFactor2 = c("Old")
            }
            
            NormalAll <- paste0(as.character(genderNormal), ".",as.character(ageGroupNormal))
            ExperAll <- paste0(as.character(expFactor1), ".", as.character(expFactor2))
            contrAll <- paste0(as.character(contFactor1), ".", as.character(contFactor2))
            SelectAll <- c(ExperAll, contrAll)
            
            genderSelect <<- as.factor(as.character(genderNormal[NormalAll %in% SelectAll]))
            ageGroupSelect <<- as.factor(as.character(ageGroupNormal[NormalAll %in% SelectAll]))
            gctFileSelect <<- gctFileNormal[,NormalAll %in% SelectAll]
            
            # Experiment design 
            
            ## Determine how many contrast I need to put in the experiment

            if (factorValue == TRUE) {
              if (length(levels(ageGroupSelect)) == 1 ){
                designSelect <- model.matrix(~genderSelect)
              } else {
                designSelect <- model.matrix(~genderSelect + ageGroupSelect)
              }
              
            } else {
              if (length(levels(genderSelect)) == 1) {
                designSelect <- model.matrix(~ageGroupSelect)
              } else {
                designSelect <- model.matrix(~genderSelect + ageGroupSelect)
              } 
            }
            
            
            GroupSelect <<- factor(paste(genderSelect, ageGroupSelect,sep="."))
            
            DE.genes.Select <- DGEList(gctFileSelect)
            
            incProgress(1/n, detail = c("calculating normal factors..."))
            
            
            DE.genes.Select <- calcNormFactors(DE.genes.Select)
            
            incProgress(1/n, detail = c("estimating binomial dispersions..."))
            
            DE.genes.Select <- estimateDisp(DE.genes.Select, designSelect)
            
            incProgress(1/n, detail = c("Fit linear model"))
            
            DE.genes.Select <- glmQLFit(DE.genes.Select, designSelect)
            
            incProgress(1/n, detail = c("Conducting genewise statistical tests..."))
            
            ## Test position
            if(factorValue) {
              if (length(levels(ageGroupSelect)) == 1){
                DE.genes.Select.list <<- glmQLFTest(DE.genes.Select, coef = 2)
              } else {
                DE.genes.Select.list <<- glmQLFTest(DE.genes.Select, coef = 2)
              }
            } else {
              if (length(levels(genderSelect)) == 1){
                DE.genes.Select.list <<- glmQLFTest(DE.genes.Select, coef = 2)
              } else {
                DE.genes.Select.list <<- glmQLFTest(DE.genes.Select, coef = 3)
              }
            }
            
            
            
            DE.genes.Select.list_p <- findCOVID_DE_Genes(DE.genes.Select.list)
            
            DE.genes.Select.list_p <- DE.genes.Select.list_p[DE.genes.Select.list_p$PValue < PValue,]
            
            fittedValuesRaw <<- DE.genes.Select$fitted.values

            DE.genes.set <<- DE.genes.Select.list_p[abs(DE.genes.Select.list_p$logFC) > LogFC,]
            
            if (nrow(DE.genes.set) < 2) {
              return(FALSE)
            } else {
              fittedValues <- DE.genes.Select$fitted.values[rownames(DE.genes.set),]
              
              medians.Select <<- t(apply(fittedValues, 1, function (x) tapply(x, GroupSelect, median)))
              
              GroupName.Select <- as.character(mapIds(org.Hs.eg.db,
                                                      keys=rownames(DE.genes.set),
                                                      column="SYMBOL",
                                                      keytype="ENSEMBL",
                                                      multiVals="first"))
              
              #name <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),
              #              filters = 'ensembl_gene_id',
              #              values = rownames(DE.genes.set),
              #              mart = mart,uniqueRows = FALSE)
              
              #temp <- rownames(DE.genes.set)
              #colnames(name) <- c("SYMBOL", "ID")
              #GroupName.Select <- name[match(temp,name$ID),]
              
              rownames(medians.Select) <- GroupName.Select
              HM.data.Select <- t(scale(t(medians.Select), center = FALSE))
              rownames(DE.genes.set) <- GroupName.Select
              
              DE.genes.set$SYMBOL <- GroupName.Select
              DE.genes.set$logFC <- round(DE.genes.set$logFC,2)
              DE.genes.set$PValue <- scientific(DE.genes.set$PValue, digits = 2)
              DE.genes.set$padj <- scientific(DE.genes.set$padj, digits = 2)
              DE.genes.set$logCPM <- scientific(DE.genes.set$logCPM, digits = 2)
              DE.genes.set$F <- round(DE.genes.set$F, 2)
              
              col_order <- c("SYMBOL", "logFC","logCPM", "F", "PValue", "padj")
              DE.genes.set <- DE.genes.set[,col_order]
              resultAll <- list(HM.data.Select, DE.genes.set)
              return(resultAll)
            }
            
            
            
          }
          
          
          result <<- analysisGct_A(gctFile = gctFIle, gctFileTPM = gctFIleTPM)
          if (is.logical(result)) {
            
            img1 <- readPNG("www/WrongGene.PNG")
            output$myHeatmap_A <- renderPlot({grid.raster(img1)})
            error <- c("We can't find DE genes. Please higher the threshold!")
            output$geneName_A <- renderDataTable({grid.raster(error)})
            thedata <<- reactive(geneName_A)
            
          } else {
            
            output$myHeatmap_A  <- renderPlot(
              {
                pheatmap(result[[1]], angle_col = c("0"), fontsize_col = 14 ,fontsize = 9,cluster_cols = FALSE)
              }
            )
            geneName_A <<- as.data.frame(result[[2]])
            rownames(geneName_A) <- NULL
            output$geneName_A <<- renderDataTable({geneName_A})
            thedata <<- reactive(geneName_A)
          }
          
          
        })
        
        
        values_A$successStep2 <- TRUE
        showNotification("Done. Please check the result.", type= "message")
        
      })
      
      observeEvent(input$showBoxPlot_A, {
        shiny::validate(
          need( values_A$successStep2 == TRUE, message = ('Missing succesful Step 1')),
          errorClass =  showNotification("Missing succesful Step 1", type= "error")
        )
        
        genderSubGroup <- as.character(genderSelect)
        ageSubGroup <- as.character(ageGroupSelect)
        subGroup <- paste0(genderSubGroup,".",ageSubGroup)
        
        DE.genes_COVID.19 <- fittedValuesRaw[rownames(fittedValuesRaw) %in% COVID_19_ID,]
        #COVID_Name <- select(EnsDb.Hsapiens.v79,
        #                                keys = rownames(DE.genes_COVID.19),
        #                                keytype = "GENEID", columns = c("SYMBOL","GENEID"))
        
        #name <<- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),
        #               filters = 'ensembl_gene_id',
        #               values = rownames(DE.genes_COVID.19),
        #               mart = mart,uniqueRows = FALSE)
        #
        
        #temp <- rownames(DE.genes_COVID.19)
        #colnames(name) <- c("SYMBOL", "ID")
        #COVID_Name <- name[match(temp,name$ID),]
        
        COVID_Name <- as.character(mapIds(org.Hs.eg.db,
                                          keys=rownames(DE.genes_COVID.19),
                                          column="SYMBOL",
                                          keytype="ENSEMBL",
                                          multiVals="first"))
        
        
        rownames(DE.genes_COVID.19) <- COVID_Name
        
        ## The position of the gene you like to return
        geneSymbol <- input$GeneSymbol_A
        i <- which(rownames(DE.genes_COVID.19) %in% geneSymbol)
        
        tempGroup <- as.data.frame(t(as.data.frame(t(DE.genes_COVID.19[i,]))))
        tempGroup$group <- subGroup
        tempGroup$ageGroup <- ageSubGroup
        tempGroup$Gender <- genderSubGroup
        
        output$myBoxPlot_A  <- renderPlot(
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
      
      
      
      
      output$downloadData1_A <- downloadHandler(
        filename = function() {"GeneTable.csv"},
        content = function(file) {
          write.csv(thedata(), file, row.names = FALSE)
        }
      )
      
      
    }
  }) 


  
}

ui <- dashboardPage(skin = "black",
                    dashboardHeader(title = "DeCovid App"),
                    dashboardSidebar(
                      sidebarMenu(
                        id = "tabs",
                        menuItem("Introduction", tabName = "Introduction", icon = icon("dashboard"), selected = TRUE),
                        menuItem("Data Analysis", tabName = "DataAnalysis", icon = icon("th")),
                        menuItem("Advanced Analysis", tabName = "AdvancedAnalysis", icon = icon("line-chart"))
                      )
                    ),
                    dashboardBody(
                      tabItems(
                        # First tab content
                        tabItem(tabName = "Introduction",
                                tabPanel(
                                  title = "Introduction",
                                  HTML(
                                    '
                      <div align="justify">
                      <h3>Tutorial Video </h3>
                      </div>
                     '
                                  ),
                                  HTML('<iframe width="1120" height="630" src="https://www.youtube.com/embed/aBwrSgVLSqQ" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>'),
                                  HTML('<div align="justify">
                        <h3>Description </h3>
This app provides information on gene expression differences between man and women and old versus young individuals for genes in the COVID19 Disease Map database. Gene expression data was obtained from the GTEx project, that collects expression data for multiple human tissues.
<ol style="list-style-type: decimal">
<li><p>You can select tissue of interest, companions type (man vs female or young vs old) and a p.value of differential expression to identify regulated COVID-19 Disease Map genes.
<li><p>Also, you can simply introduce a gene name to see the expressions data across these population groups.
                          </div>
                        '
                                  ),
                                  HTML('<img style="width: 70%; display: block; margin-left: left; margin-right: auto;" 
                   src="Pipeline.JPG"/>')
                                  
                                )
                        ),
                        
                        # Second tab content
                        tabItem(tabName = "DataAnalysis",
                                fluidRow(
                                  column(width = 3,
                                         box(title = 'Step 1: Differential expression analysis',width = 13,
                                             status = 'primary',solidHeader = T,
                                             h5("*In this step, you can choose a tissue you are interested in and perform differential expression analysis to find COVID-19 Disease Map genes with a gender difference or age difference!"),
                                             br(),
                                             selectInput("Tissue", "select the tissue you want to analysis", choices = c("Whole Blood", "Spleen",
                                                                                                                         "lymphocytes",
                                                                                                                         "Stomach",
                                                                                                                         "Kidney - Cortex",
                                                                                                                         "Lung",
                                                                                                                         "Heart - Left Ventricle",
                                                                                                                         "Heart - Atrial Appendage",
                                                                                                                         "Brain - Cortex",
                                                                                                                         "Brain - Hippocampus")),
                                             textInput("PValue", "P-value", value = "0.05"),
                                             textInput("LogFC", "Absolute Log FoldChange", value = "0.5"),
                                             selectInput("ageSelect", "select the age threshold to determine old group", choices = c("age>30", "age>40", "age>50", "age>60", "age>70"),selected = c("age>60")),
                                             radioButtons("factor", "select the factor you are interested in", choices = c("Gender", "Age")),
                                             radioButtons("background", "Please select background gene", choices = c("All", "COVID-19 Disease Map genes", "COVID-19 Related Pathway genes")),
                                             actionButton('run','1-RUN!',icon("send outline icon"),class="btn btn-primary",
                                                          style="color: #ffffff; background-color: #35424d; border-color: #35424d")
                                         )
                                  ),
                                  column(width = 9,
                                         tabBox(title = 'DE Genes',width = 13,
                                                tabPanel(title = 'DE Genes Table',
                                                         dataTableOutput('geneName')
                                                ),
                                                tabPanel(title = 'Heatmap',
                                                         plotOutput('myHeatmap', height = 1000)
                                                )
                                         ),
                                         column(width = 12,
                                                downloadButton("downloadData1", "Download Table")
                                         )
                                  )
                                ), 
                                fluidRow(
                                  column(width = 3,
                                         box(title = 'Step 2:View the expression of COVID-19 Disease Map genes',width = 13,
                                             status = 'primary',solidHeader = T,
                                             h5("*In this step, you can check the expression of the COVID-19 Disease Map genes in a boxplot!"),
                                             br(),
                                             selectInput("GeneSymbol", "Enter a COVID-19 Disease Map gene to see the boxplot", choices = COVID_19_name,multiple = FALSE),
                                             
                                             actionButton("showBoxPlot", "2-Show Boxplot", icon("send outline icon"),class="btn btn-primary",
                                                          style="color: #ffffff; background-color: #35424d; border-color: #35424d"),
                                         )
                                  ),
                                  column(width = 6,
                                         tabBox(title = 'BoxPlot',width = 16,
                                                tabPanel(title = 'BoxPlot',
                                                         plotOutput("myBoxPlot")                              
                                                )
                                         )
                                  )
                                ),
                                fluidRow(
                                  column(width = 3,
                                         box(title = 'Step 3: GO term Enrichment Analysis',width = 13,
                                             status = 'primary',solidHeader = T,
                                             h5("*In this step, you can perform a GO enrichment analysis for DE genes found in step1 under the COVID-19 Diseases Map genes background!"),
                                             br(),
                                             textInput("PValue2", "P-value",value = "0.1"),
                                             radioButtons("factor2", "select the gene size you want to analysis", choices = c("All DE genes", "DE genes expressed more in males", "DE genes expressed more in females")),
                                             actionButton("EnrichmentAnalysis", "3-GO enrichment analysis", icon("send outline icon"), class="btn btn-primary",
                                                          style="color: #ffffff; background-color: #35424d; border-color: #35424d")
                                         )
                                  ),
                                  column(width = 9,
                                         tabBox(title = 'GO plot',width = 13,
                                                tabPanel(title = 'GO plot',
                                                         plotOutput("myGO")
                                                )
                                         )
                                  )
                                  
                                )
                        ),
                        tabItem(tabName = "AdvancedAnalysis",
                                fluidRow(
                                  column(width = 3,
                                         box(title = 'Step 1: Differential expression analysis',width = 13,
                                             status = 'primary',solidHeader = T,
                                             h5("*In this step, you can choose a tissue you are interested in and perform differential expression analysis to find COVID-19 Disease Map genes with a gender difference or age difference!"),
                                             br(),
                                             selectInput("experimentalGroup", "select the experimental group", choices = c("sex:male age:20-30", "sex:male age:30-40", "sex:male age:50-60", "sex:male age:60-70", "sex:male age:70-80",
                                                                                                                           "sex:female age:20-30", "sex:female age:30-40", "sex:female age:50-60", "sex:female age:60-70", "sex:female age:70-80"), multiple = TRUE),
                                             selectInput("controlGroup", "select the control group", choices = c("sex:male age:20-30", "sex:male age:30-40", "sex:male age:50-60", "sex:male age:60-70", "sex:male age:70-80",
                                                                                                                 "sex:female age:20-30", "sex:female age:30-40", "sex:female age:50-60", "sex:female age:60-70", "sex:female age:70-80"), multiple = TRUE),
                                             radioButtons("factor_A", "select the factor you are interested in", choices = c("Gender", "Age")),
                                             actionButton("addGroupInf", "add group information", icon("table"), class="btn btn-primary",
                                                          style="color: #ffffff; background-color: #35424d; border-color: #35424d"),
                                             br(),
                                             selectInput("Tissue_A", "select the tissue you want to analysis", choices = c("Whole Blood", "Spleen",
                                                                                                                           "lymphocytes",
                                                                                                                           "Stomach",
                                                                                                                           "Kidney - Cortex",
                                                                                                                           "Lung",
                                                                                                                           "Heart - Left Ventricle",
                                                                                                                           "Heart - Atrial Appendage",
                                                                                                                           "Brain - Cortex",
                                                                                                                           "Brain - Hippocampus"
                                             )),
                                             textInput("PValue_A", "P-value", value = "0.05"),
                                             textInput("LogFC_A", "Absolute Log FoldChange", value = "0.5"),
                                             selectInput("ageSelect_A", "select the age threshold to determine old group", choices = c("age>30", "age>40", "age>50", "age>60", "age>70"),selected = c("age>60")),
                                             radioButtons("background_A", "Please select background gene", choices = c("All", "COVID-19 Disease Map genes", "COVID-19 Pathway genes")),
                                             
                                             actionButton('run_A','1-RUN!',icon("send outline icon"),class="btn btn-primary",
                                                          style="color: #ffffff; background-color: #35424d; border-color: #35424d")
                                         )
                                  ),
                                  column(width = 9,
                                         tabBox(title = 'DE Genes',width = 13,
                                                tabPanel(title = 'DE Genes Table',
                                                         dataTableOutput('geneName_A')
                                                ),
                                                tabPanel(title = 'Heatmap',
                                                         plotOutput('myHeatmap_A', height = 1000)
                                                )
                                         ),
                                         column(width = 12,
                                                downloadButton("downloadData1_A", "Download Table")
                                         )
                                  )
                                  ),
                                fluidRow(
                                  column(width = 3,
                                         box(title = 'Step 2:View the expression of COVID-19 Disease Map genes',width = 13,
                                             status = 'primary',solidHeader = T,
                                             h5("*In this step, you can check the expression of the COVID-19 Disease Map genes in a boxplot!"),
                                             br(),
                                             selectInput("GeneSymbol_A", "Enter a COVID-19 Disease Map gene to see the boxplot", choices = COVID_19_name, multiple = FALSE),
                                             
                                             actionButton("showBoxPlot_A", "2-Show Boxplot", icon("send outline icon"),class="btn btn-primary",
                                                          style="color: #ffffff; background-color: #35424d; border-color: #35424d"),
                                         )
                                  ),
                                  column(width = 4,
                                         tabBox(title = 'BoxPlot',width = 13,
                                                tabPanel(title = 'BoxPlot',
                                                         plotOutput("myBoxPlot_A")                              
                                                )
                                         )
                                  )
                                )
                                
                      ))))






shinyApp(ui = ui, server = server)
