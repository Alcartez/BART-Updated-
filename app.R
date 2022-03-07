## Only run in interactive R sessions
if (interactive()) {
  
# Automatic install if necessary
  
  if(!require("shiny")) {
    install.packages('shiny', dependencies = TRUE)
    library(shiny)
  }
  if(!require("DT")) {
    install.packages('DT', dependencies = TRUE)
    library(DT)
  }
  if(!require("GEOquery")) {
    BiocManager::install("GEOquery")
    library(GEOquery)
  }
  if(!require("oligo")) {
    BiocManager::install("oligo")
    library(oligo)
  }
  if(!require("affy")) {
    BiocManager::install("affy")
    suppressPackageStartupMessages(library(affy))
  }
  if(!require("limma")) {
    BiocManager::install("limma")
    suppressPackageStartupMessages(library(limma))
  }
  if(!require("rhandsontable")) {
    install.packages("rhandsontable", dependencies = TRUE)
    library(rhandsontable)
  }
  if(!require("gplots")) {
    install.packages('gplots', dependencies = TRUE)
    library(gplots)
  }
  if(!require("ggfortify")) {
    install.packages('ggfortify', dependencies = TRUE)
    library(ggfortify)
  }
  if(!require("ggplot2")) {
    install.packages('ggplot2', dependencies = TRUE)
    library(ggplot2)
  }
  if(!require("plotly")) {
    install.packages('plotly', dependencies = TRUE)
    library(plotly)
  }
  if(!require("RColorBrewer")) {
    install.packages('RColorBrewer', dependencies = TRUE)
    library(RColorBrewer)
  }
  if(!require("gage")) {
    BiocManager::install("gage")
    library(gage)
  }
  if(!require("gageData")) {
    BiocManager::install("gageData")
    library(gageData)
  }
  if(!require("WebGestaltR")) {
    install.packages('WebGestaltR', dependencies = TRUE)
    library(WebGestaltR)
  }
  if(!require("calibrate")) {
    install.packages('calibrate', dependencies = TRUE)
    library(calibrate)
  }
  
  options(shiny.reactlog = TRUE)
  options(shiny.maxRequestSize = 100*1024^2)
  options(shiny.error = browser)
  
  system("rm -r GSE*")
  system("rm -r CEL*")
  system("rm -r Project*")
  system("rm -r ./www/*.png")
  system("rm ./phenodata.txt")
  system("rm ./phenodata.txt")
  system("rm *.html")
  
  
  TestFile <-
    read.table("normalized_expression.txt",
               sep = '\t',
               header = T)  
  
  

  ui <- fluidPage(
    titlePanel(div("BART - Bioinformatics Array Research Tool", img(src="bart.jpg", height = 100, width = 80)), windowTitle = "BART - Bioinformatics Array R Tool"),
    
    sidebarLayout(
      sidebarPanel(
        
        textInput("geoID","Enter GEO accession number","GSE20986"),
        actionButton("goButton", "Go"),
        fileInput('file1', 'Or choose file to upload',multiple = TRUE,
                  accept = c(
                    'text/csv',
                    '.txt',
                    '.csv',
                    '.gz',
                    '.zip'
                  )
        ),
        
        uiOutput("annotation"),
        uiOutput("annSubmit"),
        uiOutput("choosePlatform"),
        uiOutput("groups"),
        uiOutput("useCEL"),
        uiOutput("groupsSubmit"),
        uiOutput("comparison")
        
      ),
      
      mainPanel(
        #textOutput("text1"),
        tableOutput("contents"),
        tabsetPanel(
          tabPanel('Description',
                   # h3("BART (Bioinformatics Array R Tool) is a tool developed by the Bioinformatics Core at the Salk Institute."),
                   # h4("BART automates the download and analysis process across diverse microarray platforms. To begin, enter a GSE accession ID or upload an expression table."),
                   # h5("BART automatically parses any experimental data available on GEO and suggests groupings to use for differential expression testing."),
                   # h5("RMA normalization is performed using the Affy (Gautier et al., 2004) or Oligo (Carvalho et al., 2010) package if CEL files are available. When an expression table is used, a log 2 transformation is automatically applied if a skew is detected."),
                   # h5("For visualization, Hierarchical Clustering is performed using the hclust R function with the 1000 most highly expressed genes. BART also generates a Principle Component Analysis (PCA) plot using the top two principal components, which show the most sample variance."),
                   # h5("The limma package (Ritchie et al., 2015) is leveraged for differential expression testing. BART annotates the results with gene names and symbols and displays the full differentially expressed gene lists for each pairwise comparison of the groups specified. Volcano Plots are available to visualize the results."),
                   cat(file=stderr(), paste(system("which R"))),
                   # textOutput("selected_var"),
                   #  h4(sessionInfo()),
                   h3("BART is an R shiny web application for microarray analysis"),
                   downloadButton("downloadTestFile","Download a test input file for analysis"),
                   downloadButton("downloadTutorial","Download a tutorial for BART"),
                   helpText(   a("Click here to learn more about BART",     href="https://bitbucket.org/Luisa_amaral/bart", target="_blank")),
                   h4("To use BART with automatic download from GEO:"),
                   h5("1. Enter a GSE Accession ID for a microarray experiment from GEO and press 'Go'"),
                   h5("2. Wait for BART to download the files to begin analysis. When download is done, you will be prompted to visit the 'Sample Grouping' tab."),
                   h5("3. Choose the platform you wish to use if there is more than one available."),
                   h5("4. Browse the phenodata BART has found in the 'Sample Grouping' tab. You can choose to use manual entry if your desired grouping is not displayed by BART. Check which grouping you would like to use for differential expression testing and press 'Submit'." ),
                   h5("5. Wait for BART to perform normalization (if necessary) and differential expression testing."),
                   h5("6. When analysis is complete, you can navigate to any tab to see your results. Some data may take a few seconds to load."),
                   #h4("6. Use the drop-down menu entitled 'View differential expression results' to see the different comparisons"),
                   
                   h3(""),
                   h4("To use BART with your own expression table:"),
                   h5("1. Click the 'Browse' button and select the expression table you wish to use. This table should have a column for IDs and columns with expression values for each sample. A sample expression table is available for download below."),
                   h5("2. Once the table is loaded, you will be prompted to navigate to the 'Sample Grouping' tab and enter a grouping for differential expression via manual entry. Once you have entered a grouping, press 'Submit'."),
                   h5("3. Wait for BART to perform a log2 transformation (if necessary) and differential expression testing."),
                   h5("4. You can now go to any tab and view your results."),
                   
                   # downloadButton("downloadTestFile","Download a test input file for analysis"),
                   # downloadButton("downloadTutorial","Download a tutorial for BART"),
                   
                   helpText("Please refresh BART in between analyses")
                   
                   
          ),
          tabPanel('Sample grouping',
                   helpText("Please enter a grouping for differential expression testing and click 'Submit'"),
                   rHandsontableOutput("hot"),
                   uiOutput("fileGroupsSubmit")
          ),
          tabPanel('Normalized Expression',
                   h4("Click on a row to see a bar graph of the data below"),
                   DT::dataTableOutput("exprs"),
                   uiOutput("logmes"),
                   plotlyOutput('x2'),
                   DT::dataTableOutput("file"),
                   helpText("It may take a few seconds to load table"),
                   uiOutput("downloadNormExprs")),
          # downloadButton('downloadExprs', 'Download Table')),
          tabPanel('Clustering',
                   uiOutput("heatmap"),
                   uiOutput("downloadHeat")),
          #downloadButton('downloadHM', 'Download Heatmap')),
          tabPanel('PCA',
                   plotly::plotlyOutput("plotPCA", width = "150%", height = "600px"),
                   uiOutput("downloadPCA")),
          tabPanel('Differential Expression',
                   uiOutput("diffTitle"),
                   textInput("pcutoff", "Enter an adj. P value cutoff", "0.05"),
                   DT::dataTableOutput("diffExprs"),
                   helpText("It may take a few seconds to load table"),
                   uiOutput("downloadDiff")),
          # downloadButton('downloadDiffExprs', 'Download Table')),
          tabPanel('Volcano Plot',
                   #uiOutput("volcano_fast"),
                   plotlyOutput("volcano"),
                   #downloadButton('downloadVolcano', 'Download Volcano Plot')),
                   uiOutput("downloadVol")),
          
          tabPanel('Pathway Analysis (GAGE)',
                   uiOutput("diff_Title"),
                   #uiOutput("species"),
                   #uiOutput("gageSubmit"),
                   h4("Up regulated KEGG pathways:"),
                   DT::dataTableOutput("KEGGpathwaysUp"),
                   uiOutput("downloadGageUp"),
                   #downloadButton('downloadkegg_greater', 'Download Full Results'),
                   h4("Down regulated KEGG pathways:"),
                   DT::dataTableOutput("KEGGpathwaysDown"),
                   uiOutput("downloadGageDown")),
          #downloadButton('downloadkegg_less', 'Download Full Results')),
          tabPanel("Pathway Analysis (WebGestalt)",
                   uiOutput("species2"),
                   uiOutput("IDtype"),
                   uiOutput("wegSubmit"),
                   h4("Significantly Over-Represented KEGG Pathways"),
                   DT::dataTableOutput("enrich_results"), 
                   uiOutput("GOSlim"),
                   uiOutput("downloadWeb")
                   #downloadButton('downloadWebgestalt', 'Download Full Analysis')
          )
          
        )
      )
    )
  )
  
  server <- function(input, output, session) {
    v <- reactiveValues(data = NULL)
    v$SysInfo <- "Initializing"
    options(shiny.maxRequestSize = 9 * 1024 ^ 2)
    s <- sessionInfo()
    output$selected_var <- renderText({ 
      s$otherPkgs$GEOquery$Version
    })
    # s <- sessionInfo()
    
    #Functions for analyses
    
    #downloads series matrix from GSE using GEOquery and determines if it is a cross-platform study
    GSEdownload <- reactive({
      if (is.null(input$geoID)) {
        return(NULL)
      }
      geoID <- input$geoID
      if (is.null(geoID))
        return(NULL)
      else {
        just_celfiles <- FALSE
        gse <- NULL
        #try to download 100 times
        for (i in 1:100) {
          results = tryCatch({
            gse <- getGEO(geoID, GSEMatrix = TRUE, AnnotGPL = TRUE)
          },
          error = function(e) {
          },
          finally = {
          })
          if (!is.null(gse)) {
            break
          }
        }
        if (is.null(gse)) {
          cat(file = stderr(),
              "GEO download failed")
          showModal(
            modalDialog(
              title = "GEO download failed",
              "Please ensure that your GSE accession number is correct and that it has a series matrix file which is downloadable on GEO. If so, refresh the page and try again in a few minutes. (The GEO server is likely down)"
            )
          )
          return(NULL)
        }
        v$plat_choices <- c()
        v$choose_platform <- FALSE
        v$gse <- gse
        v$tableidx <- 1
        if (length(gse) > 0) {
          v$choose_platform <- TRUE
          for (i in 1:length(gse)) {
            v$plat_choices <- c(v$plat_choices, annotation(gse[[i]]))
          }
        }
      }
    })
    
    #parses phenodata to display groupings
    #grabs the CEL files if available
    #Uses input$geoID
    #Makes v$pdata, v$filenames, v$exitDir, v$gse, v$pdataALL
    pdata <- reactive({
      if (is.null(input$choosePlatform)) {
        return(NULL)
      }
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Performing ", value = 0)
      progress$inc(0.70, detail = "Downloading and parsing supplemental data")
      geoID <- input$geoID
      just_celfiles <- FALSE
      gse <- v$gse
      cross <- FALSE
      if (length(gse) > 1) {
        cross <- TRUE
        message <- paste(geoID, "Cross Platform")
        cat(message, "\n")
        cat(input$choosePlatform, v$plat_choices)
        idx <- match(input$choosePlatform, v$plat_choices)
        cat(idx)
      }
      if (length(gse) == 1) {
        cat("Not cross platform")
        idx <- 1
      }
      
      phenodata <- pData(gse[[idx]])
      
      org_idx <- grep('organism', colnames(phenodata))[1]
      v$organism <- phenodata[,org_idx][1]
      cat(rownames(phenodata), nrow(phenodata), "\n")
      
      if (nrow(phenodata) <= 1) {
        stop("Not enough samples in dataset", call. = FALSE)
      }
      fvarLabels(gse[[idx]]) <- make.names(fvarLabels(gse[[idx]]))
      v$geoFile <-
        paste("./", geoID, sep = "") #folder with TAR file and filelist
      geoFile <- v$geoFile
      results <- ""
      exitDir <-
        paste("CEL_", geoID, sep = "") #folder with phenodata & untarred files
      supp = NULL
      supp = tryCatch({
        (getGEOSuppFiles(geoID))
      },
      error = function(e) {
      },
      warning = function(w) {
      },
      finally = {
      })
      
      progress$inc(0.25, detail = "Downloading and parsing supplemental data")
      v$no_supp <- FALSE
      if (is.null(supp)) {
        v$no_supp <- TRUE
        cat("SUPPLEMENTARY FILES NOT OBTAINED")
        filenames <- sampleNames(gse[[idx]])
      } else {
        if (!file.exists(exitDir)) {
          string <- paste(geoID, "/", geoID, "_RAW.tar", sep = "")
          untar(string, exdir = exitDir)
          cels <- list.files(exitDir, pattern = "[CEL.gz]")
          sap <- sapply(paste(exitDir, cels, sep = "/"), gunzip)
        }
        setwd(exitDir)
        setwd("..")
      }
      
      if (cross) {
        rows <- rownames(phenodata)
        for (i in 1:length(rows)) {
          cat("row", i)
          filenam = paste(rows[i], ".CEL", sep = "")  # maybe change this later to supplementary_file name
          if (i == 1) {
            p <-
              paste("ls ", exitDir, "/", filenam,  "> filenames.txt", sep = "")
          } else {
            p <-
              paste("ls ", exitDir, "/", filenam,  ">> filenames.txt", sep = "")
          }
          system(p)
        }
        filenames = tryCatch({
          read.table("filenames.txt")
        },
        error = function(e) {
          message <-
            paste("Supplemental files not found, using GSE matrix instead")
          cat(message)
        },
        warning = function(w) {
          message <-
            paste("Supplemental files not found, using GSE matrix instead")
          cat(message)
        },
        finally = {
        })
      } else {
        p <- paste("ls ", exitDir, "/*.CEL > filenames.txt", sep = "")
        system(p)
        filenames = tryCatch({
          read.table("filenames.txt")
        },
        error = function(e) {
          message <-
            paste("Supplemental files not found, using GSE matrix instead")
          cat(message)
        },
        warning = function(w) {
          message <-
            paste("Supplemental files not found, using GSE matrix instead")
          cat(message)
        },
        finally = {
        })
      }
      if (is.null(filenames)) {
        v$no_supp <- TRUE
        filenames <- sampleNames(gse[[idx]])
      }
      groupFile <- paste(exitDir, "/phenodata.txt", sep = "")
      filenames <- as.data.frame(filenames)
      for (i in 1:ncol(phenodata)) {
        name <- paste(colnames(phenodata[i]))
        if (name == "source_name_ch1") {
          colnames(phenodata)[i] <- "source"
          name <- "source"
        }
        col <- phenodata[i]
        grp <- grep("characteristics", name)
        if (i == 1) {
          if (length(grp) == 0) {
            col <- sub("^", paste(name, ": ", sep = ""), col[, c(name)])
          }
          pdata <- as.data.frame(col)
          colnames(pdata)[i] <- name
        }
        else {
          if (length(grp) == 0) {
            col <- sub("^", paste(name, ": ", sep = ""), col[, c(name)])
          }
          pdata <- cbind(pdata, col)
          colnames(pdata)[i] <- name
        }
      }
      pdata <- as.data.frame(pdata)
      pdataRel <- pdata #relevant P data to be formed
      nameVec = c()
      cat("************************************************* \n")
      show(head(pdata))
      for (i in names(pdata)) {
        col <- pdata[, c(i)]
        title <- strsplit(paste(col[1]), ":")[[1]][1]
        numGroups <- length(unique(col))
        if ((numGroups < length(col) && numGroups > 1 && !(title %in% nameVec)) || title == "title") {
          nameVec <- c(nameVec, title)
          cat(
            "Column name: ",
            title,
            "\n",
            "Number of different groups: ",
            numGroups,
            "\n",
            sep = ""
          )
          if (numGroups < 21) {
            cat("Here is a list of groups and their frequencies: ",
                "\n",
                "\n")
          }
          if (numGroups >= 21) {
            cat("Too many groups to display. Here is a sample of the groups: ",
                "\n")
          }
          freq_table <- as.data.frame(table(col))
          freq_table <- freq_table[order(-freq_table$Freq), c(1, 2)]
          colnames(freq_table) <- c(title, "Freq")
          for (j in 1:21) {
            if (j > numGroups) {
              break
            }
            group <- paste(freq_table[j, 1])
            cat(group, ": ", paste(freq_table[j, 2]), "\n")
          }
          cat("===============================================",
              "\n")
        }
        else {
          pdataRel[, c(i)] <- NULL
        }
      }
      colnames(pdataRel) <- nameVec
      pdata <- as.data.frame(pdataRel, stringsAsFactors = FALSE)
      newpdata <- sapply(pdata, as.character)
      show(head(pdata))
      if (ncol(pdata) == 0) {
        pdata =  as.data.frame(matrix(
          data = "",
          nrow = nrow(pdata),
          ncol = 1
        ))
        newpdata = pdata
        colnames(pdata[1]) = "No characteristics found"
        colnames(newpdata[1]) = "No characteristics found"
      } else {
        for (i in 1:(ncol(pdata))) {
          for (j in 1:(nrow(pdata))) {
            # newpdata[j, i] <-
            #    as.character(strsplit(paste(pdata[j, i]), ":")[[1]][2])
            cont <- as.character(strsplit(paste(pdata[j, i]), ":")[[1]])
            newpdata[j, i] <- cont[length(cont)]
          }
        }
      }
      
      rownames(newpdata) <- rownames(phenodata)
      Manual_Entry <-
        as.data.frame(matrix(
          data = "edit here",
          nrow = nrow(pdata),
          ncol = 1
        ))
      colnames(Manual_Entry) <- "manual entry"
      newpdata <- cbind(newpdata, Manual_Entry)
      Batch_Effect <-
        as.data.frame(matrix(
          data = "edit here",
          nrow = nrow(pdata),
          ncol = 1
        ))
      colnames(Batch_Effect) <- "batch effect"
      newpdata <- cbind(newpdata, Batch_Effect)
      v$pdata <- newpdata
      v$filenames <- filenames
      v$gse <- gse
      v$exitDir <- exitDir
      colnames(newpdata) = make.unique(colnames(newpdata))
      v$pdataALL <- newpdata
      v$dataALL <- newpdata
      v$idx <- idx
      cat("\nCURRENT PDATA COLS", colnames(v$pdata))
    })
    
    
    #performs normalization
    #Uses v$exitDir, v$pdata, v$filenames
    #Makes phenodata.txt, v$eSet
    normalizedExpression <- reactive({
      if (is.null(v$exitDir) ||
          is.null(v$pdata) || is.null(v$filenames)) {
        return(NULL)
      }
      cat("\n", v$exitDir)
      sink("./phenodata.txt")
      cat("Name", "FileName", "Target", sep = "\t")
      cat("\n", sep = "")
      for (i in 1:(nrow(v$pdata))) {
        c1 <- paste(v$filenames[[i, 1]])
        c2 <- paste(v$filenames[[i, 1]])
        c3 <- ""
        for (j in 1:ncol(v$pdata)) {
          g <- v$pdata[i, j]
          if (j == 1) {
            g <- gsub(" ", "", g)
            c3 <- paste(c3, g, sep = "")
          }
          else {
            g <- gsub(" ", "", g)
            c3 <- paste(c3, g, sep = "_")
          }
        }
        c3 <- make.names(c3)
        cat(c1, c2, c3, sep = "\t")
        cat("\n", sep = "")
      }
      sink()
      v$useMatrix <- FALSE
      
      if (input$useCEL == "Expres_table") {
        v$useMatrix <- TRUE
        celfiles <- NULL
        celFiles <- NULL
      }
      desc <- "phenodata.txt"
      if (v$useMatrix == FALSE) {
        celfiles = tryCatch({
          read.affy(covdesc = desc, path = ".")
        },
        warning = function(w) {
          cat(file = stderr(), "Trying Oligo package \n")
        },
        error = function(e) {
          cat(file = stderr(), "Trying Oligo package \n")
        },
        finally = {
          
        })
      }
      
      useOligo <- is.null(celfiles)
      if (useOligo) {
        if (v$useMatrix == FALSE) {
          celFiles = tryCatch({
            list.celfiles(v$exitDir, full.names = TRUE)
          },
          warning = function(w) {
            
          },
          error = function(e) {
            
          },
          finally = {
            
          })
        }
        if (!is.null(celFiles) && !isEmpty(celFiles)) {
          celfiles <- read.celfiles(celFiles)
          eSet <- oligo::rma(celfiles)
        } else {
          #no supplementary CEL files found, so using user supplied expression values
          cat(file = stderr(), "Using GSE matrix 1\n")
          if(input$useCEL == "CEL") {
            showModal(modalDialog(title = "CEL file analysis failed", "Supplementary CEL files were not found or not downloaded correctly. Using expression table for analysis"))
          }
          v$useMatrix <- TRUE
          eSet <- v$gse[[v$idx]]
          ex <- exprs(v$gse[[v$idx]])
          #Detecting if Log transformation needed
          qx <-
            as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
          LogC <- (qx[5] > 100) ||
            (qx[6] - qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
          if (LogC) {
            ex[which(ex <= 0)] <- NaN
            exprs(v$gse[[v$idx]]) <- log2(ex)
            # exprs(v$gse[[v$idx]]) <- normalizeBetweenArrays(exprs(v$gse[[v$idx]]))
            cat("Log 2 transformation applied to data \n")
            output$logmes <- renderUI({
              h5("Log 2 transformation was applied to data")
            })
          } else {
            output$logmes <- renderUI({
              h5("Using expression table from GEO. No log transformation was needed.")
            })
          }
          v$no_supp <- TRUE
        }
      }
      if (!useOligo && !v$useMatrix) {
        cat("\n Using affy package \n")
        eSet = tryCatch({
          affy::rma(celfiles)
        },
        error = function(e) {
          cat(file = stderr(),
              "Problems encountered in affy rma normalization")
        },
        finally = {
          
        })
        if (is.null(eSet)) {
          cat(file = stderr(), "Using GSE matrix 2\n")
          v$no_supp <- TRUE
          v$useMatrix <- TRUE
          eSet <- v$gse[[v$idx]]
          ex <- exprs(v$gse[[v$idx]])
          if (nrow(ex) == 0) {
            showModal(
              modalDialog(
                title = "Oops",
                "Normalization failed and no supplementary expression matrix found. BART may be incompatible with this experiment"
              )
            )
            return(NULL)
          }
          #Detecting if Log transformation needed
          qx <-
            as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
          LogC <- (qx[5] > 100) ||
            (qx[6] - qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
          if (LogC) {
            ex[which(ex <= 0)] <- NaN
            exprs(v$gse[[v$idx]]) <- log2(ex)
            cat("Log 2 transformation applied to data \n")
            output$logmes <- renderUI({
              h5("Log 2 transformation was applied to data")
            })
          } else {
            output$logmes <- renderUI({
              h5(
                "No supplementary CEL files found, using user supplied expression values from GEO. Log transformation was not needed."
              )
              
            })
          }
          v$no_supp <- TRUE
        }
      }
      v$eSet <- eSet
      v$exprs <- exprs(eSet)
      v$ptable <- read.table("./phenodata.txt", header = T)
      v$samples <- make.names(v$ptable$Target)
      v$samples <- as.factor(v$samples)
      batch <- factor(v$batch)
      if (!is.null(v$batch)) {
        v$design <- model.matrix( ~ 0 + v$samples + batch)
      } else {
        v$design <- model.matrix( ~ 0 + v$samples)
      }
      colnames(v$design) = make.unique(make.names(colnames(v$design)))
      
    })
    
    #Performs differential expression and makes volcano plots
    #uses v$gse, v$eSet, v$samples, v$design
    #makes v$diffList
    getDiffExpr <- reactive ({
      if (is.null(v$gse) || is.null(v$eSet) || is.null(v$design)) {
        return(NULL)
      }
      for (i in 1:length(unique(v$samples))) {
        colnames(v$design)[i] <- paste(unique(sort(v$samples))[i])
      }
      ret <- list()
      fit <- lmFit(v$gse[[v$idx]], v$design)
      genelist <- fit$genes
      show(head(genelist))
      entrez = genelist$Gene.ID
      geneinfo = data.frame(matrix(nrow = nrow(genelist), ncol = 1))
      if ('ID' %in% colnames(genelist)) {
        geneinfo[, ncol(geneinfo) + 1] <- genelist$ID
        colnames(geneinfo)[ncol(geneinfo)] = 'ID'
      }
      if ('Gene.symbol' %in% colnames(genelist)) {
        geneinfo[, ncol(geneinfo) + 1] <- genelist$`Gene.symbol`
        colnames(geneinfo)[ncol(geneinfo)] = 'Gene.symbol'
      } else if ('Gene.Symbol' %in% colnames(genelist)) {
        geneinfo[, ncol(geneinfo) + 1] <- genelist$`Gene.Symbol`
        colnames(geneinfo)[ncol(geneinfo)] = 'Gene.symbol'
      } else if ('GeneSymbol' %in% colnames(genelist)) {
        geneinfo[, ncol(geneinfo) + 1] <- genelist$`GeneSymbol`
        colnames(geneinfo)[ncol(geneinfo)] = 'Gene.symbol'
      } else if ('GeneName' %in% colnames(genelist)) {
        geneinfo[, ncol(geneinfo) + 1] <- genelist$`GeneName`
        colnames(geneinfo)[ncol(geneinfo)] = 'Gene.symbol'
      }
      if ('Gene.title' %in% colnames(genelist)) {
        geneinfo[, ncol(geneinfo) + 1] <- genelist$`Gene.title`
        colnames(geneinfo)[ncol(geneinfo)] = 'Gene.title'
      } else if ('Gene.Title' %in% colnames(genelist)) {
        geneinfo[, ncol(geneinfo) + 1] <- genelist$`Gene.Title`
        colnames(geneinfo)[ncol(geneinfo)] = 'Gene.title'
      }
      if ('Description' %in% colnames(genelist)) {
        geneinfo[, ncol(geneinfo) + 1] <- genelist$`Description`
        colnames(geneinfo)[ncol(geneinfo)] = 'Description'
      }
      if ('Gene.ID' %in% colnames(genelist)) {
        geneinfo[, ncol(geneinfo) + 1] <- genelist$Gene.ID
        colnames(geneinfo)[ncol(geneinfo)] = 'Entrez.ID'
      } else if ('ENTREZ_GENE_ID' %in% colnames(genelist)) {
        geneinfo[, ncol(geneinfo) + 1] <- genelist$ENTREZ_GENE_ID
        colnames(geneinfo)[ncol(geneinfo)] = 'Entrez.ID'
      } else if ('Entrez.Gene' %in% colnames(genelist)) {
        geneinfo[, ncol(geneinfo) + 1] <- genelist$Entrez.Gene
        colnames(geneinfo)[ncol(geneinfo)] = 'Entrez.ID'
      }
      show(head(geneinfo))
      rownames(geneinfo) = geneinfo$ID
      geneinfo = as.data.frame(geneinfo)
      geneinfo = geneinfo[, -c(1), drop = F]
      if (!v$no_supp) {
        cat("Using RMA normalized CEL files \n")
        output$logmes <- renderUI({
          h5("Used RMA normalized CEL files")
        })
        fit <- lmFit(v$exprs, v$design)
      }
      #Making list of args for makeContrasts
      myargs <- list()
      pair <- 1
      if (is.null(v$batch)) {
        for (i in 1:(ncol(v$design) - 1)) {
          num <- i + 1
          for (j in (num):(ncol(v$design))) {
            myargs[[pair]] <-
              paste(colnames(v$design)[i], "-", colnames(v$design)[j], sep = "")
            pair <- pair + 1
          }
        }
      } else {
        for (i in 1:(ncol(v$design) - v$num_batch)) {
          num <- i + 1
          for (j in (num):(ncol(v$design) - (v$num_batch - 1))) {
            myargs[[pair]] <-
              paste(colnames(v$design)[i], "-", colnames(v$design)[j], sep = "")
            pair <- pair + 1
          }
        }
      }
      show(myargs)
      show(v$design)
      show(list(levels = v$design))
      
      myargs <- append(myargs, list(levels = v$design))
      contrast.matrix <- do.call(makeContrasts, myargs)
      
      fits <- contrasts.fit(fit, contrast.matrix)
      ebFit <- eBayes(fits, 0.01)
      no_symbol <- FALSE
      
      cat("Output filenames: \n")
      tT <- c()
      normExp <- v$exprs
      for (i in 1:(length(myargs) - 1)) {
        tT <-
          topTable(
            ebFit,
            coef = i,
            adjust = "fdr",
            sort.by = "p",
            number = Inf
          )
        rem = c(which(!(
          rownames(tT) %in% rownames(geneinfo)
        )))
        if (!isEmpty(rem)) {
          tT = tT[-rem, ]
        }
        rem2 = c(which(!(
          rownames(geneinfo) %in% rownames(tT)
        )))
        if (!isEmpty(rem2)) {
          geneinfo = geneinfo[-rem2, ]
        }
        st = tT[order(rownames(tT)), ]
        st <- subset (st,
                      select = c("adj.P.Val",
                                 "P.Value",
                                 "logFC",
                                 "AveExpr",
                                 "t",
                                 "B"))
        tT =  cbind(geneinfo, st)
        tT = tT[order(tT$adj.P.Val, decreasing = F), ]
        show(head(tT))
        f <- gsub("-", "-vs-", myargs[[i]])
        outfile <- paste(f, ".txt", sep = "")
        cat(paste(outfile), "\n")
        head(tT, row.names = FALSE, n = 10)
        ret <- c(ret, list(tT, paste(f)))
        #for use later
        #topgenes = results[results[, "adj.P.Val"] < 0.05, ]
        #topups = topgenes[topgenes[, "logFC"] > 1, ]
        #topdowns = topgenes[topgenes[, "logFC"] < -1, ]
      }
      if ('Gene.symbol' %in% colnames(geneinfo)) {
        normExp <- normExp[, order(v$samples)]
        currNames <- colnames(normExp)
        for (i in 1:length(colnames(normExp))) {
          colnames(normExp)[i] <-
            paste(currNames[i], ' ', v$samples[order(v$samples)[i]])
        }
        rem3 = c(which(!(
          rownames(normExp) %in% rownames(geneinfo)
        )))
        if (!isEmpty(rem3)) {
          normExp = normExp[-rem3, ]
        }
        SYMBOL <- as.data.frame(paste(geneinfo$Gene.symbol))
        normExp <- cbind(SYMBOL, normExp)
        colnames(normExp)[1] = "Symbol"
      }
      v$norm_exprs <- normExp
      v$diffList <- ret
      # system(paste("rm -r ", v$exitDir))
      # system(paste("rm -r ", v$geoFile))
    })
    
    # Runs WebGestaltR
    RunWebGestaltAnalysis <- reactive({
      if (!is.null(input$org)) {
        cat("IN WEBGESTALT")
        curr_organism <- input$org
      }
      if (is.null(v$organism) || is.null(input$comparison)) {
        return(NULL)
      }
      if (is.null(v$currDiffList) || nrow(v$currDiffList) <= 1 ) {
        showModal(modalDialog(title = "Not enough significantly differentially expressed genes for analysis"))
        return(NULL)
      }
      if (is.null(input$pcutoff)) {
        pcutoff = 0.05
      } else {
        pcutoff = as.double(input$pcutoff)
      }
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Performing WebGestalt analysis", value = 0)
      progress$inc(0.10)
      cat(pcutoff)
      if (is.null(input$org)) {
        full_org = tolower(v$organism)
        curr_organism = paste(substring(full_org,1,1), strsplit(full_org, " ")[[1]][2], sep = '')
        IDtype = "genesymbol"
        sig_genes = v$currDiffList$Gene.symbol[1:(which(v$currDiffList$adj.P.Val>pcutoff)[1]-1)]
      } else {
        IDtype = input$IDtype
        sig_genes = rownames(v$currDiffList)[1:(which(v$currDiffList$adj.P.Val>pcutoff)[1]-1)]
        show(head(v$currDiffList))
        cat("sig_gene ",nrow(sig_genes) )
        
      }
      write.table(sig_genes, file = "webgest1.txt", quote = F, row.names = F, col.names = F)
      sig_genes[which((sig_genes == ''))] = paste0(as.character('NA'))
      sig_genes = sapply(strsplit(as.character(sig_genes), "///"), "[[", 1)
      write.table(sig_genes, file = "webgest.txt", quote = F, row.names = F, col.names = F)
      
      if(curr_organism %in% listOrganism()) {
        org = curr_organism
        paste(as.vector(sig_genes))
        v$enrichResult<-WebGestaltR(enrichMethod="ORA",organism=org,
                                    enrichDatabase="pathway_KEGG",interestGene = as.vector(sig_genes),
                                    interestGeneType=IDtype, referenceSet = "genome_protein-coding", outputDirectory = getwd(), projectName = input$comparison)
        if(length(v$enrichResult) <= 1) {
          showModal(modalDialog(title = "KEGG pathway analysis failed", "Likely not enough significantly differentially expressed genes for this comparison, or the gene symbols are incorrect. Try adjusting the adj p value cutoff or try a different differential expression list."))
        }
        system(paste("mv ./Project_",input$comparison,"/goslim_summary_",input$comparison,".png", " ./www/goslim_summary_",input$comparison,".png",sep = ""))
      }
      else {
        cat("Organism not available")
        showModal(modalDialog(title = "WebGestalt analysis is not available for this organism"))
        return(NULL)
      }
    }
    )
    
    #########
    #Run pathway analysis
    # only works for Homo sapiens, Mus Muculus, Rattus norvegicus, Saccharomyces cerevisiae
    RunPathwayAnalysis <- reactive({
      if (is.null(v$organism)) {
        return(NULL)
      }
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Performing KEGG analysis", value = 0)
      progress$inc(0.10)
      data(bods)
      bods <- as.data.frame(bods)
      tT <- v$currDiffList
      foldchanges = as.matrix(tT$logFC)
      foldchanges <- mapply(foldchanges, FUN = as.numeric)
      if ("Gene.symbol" %in% colnames(tT) &&
          !("Entrez.ID" %in% colnames(tT))) {
        if (v$organism == "Homo Sapiens" || v$organism == "Homo sapiens") {
          library(org.Hs.eg.db)
          tT$Entrez.ID = mapIds(
            org.Hs.eg.db,
            keys = tT$Gene.symbol,
            column = "ENTREZID",
            keytype = "SYMBOL",
            multiVals = "first"
          )
        }
        if (v$organism == "Mus Muculus") {
          library(org.Mm.eg.db)
          tT$Entrez.ID = mapIds(
            org.Mm.eg.db,
            keys = tT$Gene.symbol,
            column = "ENTREZID",
            keytype = "SYMBOL",
            multiVals = "first"
          )
        }
        if (v$organism == "Rattus norvegicus") {
          library(org.Rn.eg.db)
          tT$Entrez.ID = mapIds(
            org.Rn.eg.db,
            keys = tT$Gene.symbol,
            column = "ENTREZID",
            keytype = "SYMBOL",
            multiVals = "first"
          )
        }
        if (v$organism == "Saccharomyces cerevisiae") {
          library(org.Sc.sgd.db)
          tT$Entrez.ID = mapIds(
            org.Sc.sgd.db,
            keys = tT$Gene.symbol,
            column = "ENTREZID",
            keytype = "SYMBOL",
            multiVals = "first"
          )
        }
      }
      if ("Entrez.ID" %in% colnames(tT)) {
        names(foldchanges) = tT$Entrez.ID
        cat(v$organism)
        #To extract the Entrez ID (for pathway analysis)
        progress$inc(0.10)
        if (v$organism == "Homo sapiens") {
          cat(file = stderr(), "This is human!!!", "\n")
          data(kegg.sets.hs)
          data(sigmet.idx.hs)
          kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
          v$keggres = gage(
            foldchanges,
            gsets = kegg.sets.hs,
            same.dir = TRUE,
            use.stouffer = TRUE
          )
          progress$inc(0.10)
        }
        else if (v$organism == "Mus musculus") {
          cat(file = stderr(), "This is mouse!!!", "\n")
          data(kegg.sets.mm)
          data(sigmet.idx.mm)
          kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
          v$keggres = gage(
            foldchanges,
            gsets = kegg.sets.mm,
            same.dir = TRUE,
            use.stouffer = TRUE
          )
          progress$inc(0.10)
        }
        
        else if (v$organism == "Rattus norvegicus") {
          cat(file = stderr(), "This is rat!!!", "\n")
          data(kegg.sets.rn)
          data(sigmet.idx.rn)
          kegg.sets.rn = kegg.sets.rn[sigmet.idx.rn]
          v$keggres = gage(
            foldchanges,
            gsets = kegg.sets.rn,
            same.dir = TRUE,
            use.stouffer = TRUE
          )
          progress$inc(0.10)
        }
        
        else if (v$organism == "Saccharomyces cerevisiae") {
          cat(file = stderr(), "This is yeast!!!", "\n")
          data(kegg.sets.sc)
          data(sigmet.idx.sc)
          kegg.sets.sc = kegg.sets.sc[sigmet.idx.sc]
          v$keggres = gage(
            foldchanges,
            gsets = kegg.sets.sc,
            same.dir = TRUE,
            use.stouffer = TRUE
          )
          progress$inc(0.10)
        }
        else {
          cat("Not available for this species")
          showModal(modalDialog(title = "GAGE analysis is not available for this species"))
        }
      } else {
        cat("Kegg analysis not available for this GPL")
        showModal(modalDialog(title = "GAGE analysis is not available for this GPL"))
        return(NULL)
      }
    })
    #Start running analysis here
    
    #if a file is uploaded (no GEO)
    observeEvent(input$file1, {
      if (is.null(input$file1)) {
        return(NULL)
      }
      if (!is.null(v$gse) || !is.null(v$inexprs)) {
        showModal(
          modalDialog(
            title = "Hey!",
            "Please refresh the page if you want to perform a different analysis"
          )
        )
        Sys.sleep(20)
      }
      # Progress Bar
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Performing ", value = 0)
      progress$inc(0.2, detail = "initialization")
      infile <- input$file1
      if (length(input$file1[, 1]) == 1) {
        v$inexprs <- read.csv(infile$datapath, sep = "\t", header = T,row.names = NULL)
      }
      ex <- v$inexprs
      show(head(ex))
      # if ("SYMBOL" %in% colnames(ex)) {
      #   rownames(ex) = make.names(ex$SYMBOL, unique = T)
      #   ex = ex[, -c(which(colnames(ex) == "SYMBOL"))]
      #   v$inexprs <- ex
      # }
      # else if ("Symbol" %in% colnames(ex)) {
      #   rownames(ex) = make.names(ex$Symbol, unique = T)
      #   ex = ex[, -c(which(colnames(ex) == "Symbol"))]
      #   v$inexprs <- ex
      # }
      # else if ("ID" %in% colnames(ex)) {
      #   rownames(ex) = ex$ID
      #   ex = ex[, -c(which(colnames(ex) == "ID"))]
      #   v$inexprs <- ex
      # }
      # else if ("PROBE_ID" %in% colnames(ex)) {
      #   rownames(ex) = ex$PROBE_ID
      #   ex = ex[, -c(which(colnames(ex) == "PROBE_ID"))]
      #   v$inexprs <- ex
      # }
      # else if ("X" %in% colnames(ex)) {
      #   rownames(ex) = ex$X
      #   ex = ex[, -c(which(colnames(ex) == "X"))]
      #   v$inexprs <- ex
      # } 
      if (length(unique(ex[,1])) == nrow(ex)) {
        show(head(ex))
        rownames(ex) = ex[,1]
        ex = ex[, -c(1)]
        v$inexprs <- ex
      } else {
        showModal(
          modalDialog(
            title = "Input file not accepted",
            "Please refer to the BART tutorial available on this page for more information about formatting the expression table. The first column must contain unique identifiers. Refresh the page and try again."
          )
        )
        return(NULL)
      }
      if (typeof(ex[1,1]) == "character" || (as.numeric(ex[3,1]) == 3)) {
        showModal(
          modalDialog(
            title = "Input file not accepted",
            "Please refer to the BART tutorial available on this page for more information about formatting the expression table. All but the first column must contain numbers. Refresh the page and try again."
          )
        )
        return(NULL)
      }
      
      #show(head(ex))
      progress$inc(0.05, detail = "determining if normalization needed")
      qx <-
        as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
      LogC <- (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
      if (LogC) {
        ex[which(ex <= 0)] <- NaN
        v$inexprs <- log2(ex)
        progress$inc(0.65, detail = "normalization")
        cat("Log 2 transformation applied to data \n")
        output$logmes <- renderUI({
          h5("Log 2 transformation was applied to data")
        })
      } else {
        output$logmes <- renderUI({
          h5("Log 2 transformation was not needed")
        })
      }
      
      sampNot <-
        showNotification(
          "Please go to the 'Sample grouping' tab to proceed with analysis",
          type = "message",
          duration = 35
        )
      
      
      output$hot <- renderRHandsontable({
        mat <-
          data.frame(
            character = rep("edit here", ncol(v$inexprs)),
            character = rep("edit here", ncol(v$inexprs))
          )
        rownames(mat) <- colnames(v$inexprs)
        colnames(mat) <- c("Enter grouping", "Batches (if necessary)")
        width = nchar(rownames(mat)[1]) * 10
        cat(width)
        rhandsontable(
          mat,
          useTypes = FALSE ,
          readOnly = FALSE,
          rowHeaderWidth = width,
          columnSorting = TRUE
        ) %>%
          hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
      })
      output$fileGroupsSubmit <- renderUI ({
        actionButton("fileGroupsSubmit", "Submit")
      })
      observeEvent(input$fileGroupsSubmit, {
        if(is.null(input$hot)) {
          showNotification("Please enter a grouping before submitting")
          return(NULL)
        }
        pdata <-
          data.frame(hot_to_r(input$hot)[, 1], row.names = rownames(hot_to_r(input$hot)))
        batch <- factor(hot_to_r(input$hot)[, 2])
        if (1 %in% table(pdata) || nrow(unique(pdata)) <= 1) {
          cat("grouping not accepted")
          showModal(
            modalDialog(
              title = "Grouping not accepted",
              "At least two different groups with at least two replicates per group are needed for testing, please enter a different grouping"
            )
          )
          
        } else {
          removeNotification(sampNot)
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "Performing ", value = 0)
          cat("grouping accepted")
          samples <- make.names(pdata[, 1])
          output$exprs <- DT::renderDataTable({
            if (is.null(v$inexprs)) {
              return(NULL)
            }
            v$displayExprs <- v$inexprs[, order(samples)]
            exprs <- v$displayExprs
            for (i in 1:ncol(v$displayExprs)) {
              colnames(exprs)[i] = paste(colnames(v$displayExprs)[i], samples[order(samples)][i])
            }
            exprs
          }, selection = list(mode = "single", selected = 1))
          
          output$downloadNormExprs <- renderUI ({
            if (!is.null(v$inexprs)) {
              downloadButton('downloadExprs', 'Download Table')
            }
          })
          
          output$downloadExprs <- downloadHandler(
            filename = function() {
              paste("normalized_expression.txt", sep = "")
            },
            content = function(file) {
              write.table(
                v$inexprs,
                file,
                sep = "\t",
                row.names = T,
                quote = FALSE
              )
            }
          )
          
          
          output$x2 <- renderPlotly({
            if (is.null(v$displayExprs)) {
              return(NULL)
            }
            s = input$exprs_rows_selected
            if (!length(s)) {
              return(NULL)
            } else {
              m = v$displayExprs[s[length(s)], ]
            }
            plot_ly(
              x = factor(colnames(v$inexprs)[order(samples)], levels = colnames(v$inexprs)[order(samples)]),
              y = as.double(m),
              type = "bar",
              color = factor(samples[order(samples)])
            ) %>%
              layout(title = rownames(v$displayExprs)[s],
                     margin = list(b = 100))
          })
          
          progress$inc(0.40, detail = "differential expression testing")
          if (length(unique(batch)) == 1) {
            design <- model.matrix( ~ 0 + samples)
          } else {
            design <- model.matrix( ~ 0 + samples  + batch)
          }
          for (i in 1:length(unique(samples))) {
            colnames(design)[i] <- paste(unique(sort(samples))[i])
          }
          fit <- lmFit(v$inexprs, design)
          cat("fit")
          myargs <- list()
          pair <- 1
          num_batch = length(unique(batch))
          if (num_batch <= 1) {
            for (i in 1:(ncol(design) - 1)) {
              num <- i + 1
              for (j in (num):(ncol(design))) {
                myargs[[pair]] <-
                  paste(colnames(design)[i], "-", colnames(design)[j], sep = "")
                pair <- pair + 1
              }
            }
          } else {
            for (i in 1:(ncol(design) - num_batch)) {
              num <- i + 1
              for (j in (num):(ncol(design) - (num_batch - 1))) {
                myargs[[pair]] <-
                  paste(colnames(design)[i], "-", colnames(design)[j], sep = "")
                pair <- pair + 1
              }
            }
          }
          cat("made args\n")
          myargs <- append(myargs, list(levels = design))
          contrast.matrix <- do.call(makeContrasts, myargs)
          fits <- contrasts.fit(fit, contrast.matrix)
          ebFit <- eBayes(fits, 0.01)
          diffList <- c()
          cat("Output filenames: \n")
          for (i in 1:(length(myargs) - 1)) {
            tT <-
              topTable(
                ebFit,
                coef = i,
                adjust = "fdr",
                sort.by = "p",
                number = Inf
              )
            f <- gsub("-", "-vs-", myargs[[i]])
            progress$inc(0.10, detail = f)
            # name = paste("www/", f, "_volcanoplot.png", sep = "")
            # png(name,
            #     res = 100,
            #     width = 800,
            #     height = 800)
            # with(
            #   tT,
            #   plot(
            #     logFC,
            #     -log10(adj.P.Val),
            #     pch = 16,
            #     cex = 0.45,
            #     main = "Volcano plot",
            #     xlab = "Log Fold Change",
            #     ylab = "-log10(adj.P.Val)"
            #   )
            # )
            # if (nrow(subset(tT, adj.P.Val < .05 &
            #                 abs(logFC) > 1)) < 100) {
            #   with(
            #     subset(tT, adj.P.Val < .05 &
            #              abs(logFC) > 1),
            #     points(
            #       logFC,
            #       -log10(adj.P.Val),
            #       pch = 16,
            #       cex = 0.45,
            #       col = "red"
            #     )
            #   )
            #   with(
            #     subset(tT, adj.P.Val < .05 &
            #              abs(logFC) > 1),
            #     textxy(
            #       logFC,
            #       -log10(adj.P.Val),
            #       labs = rownames(tT),
            #       cex = .6
            #     )
            #   )
            # } else {
            #   with(
            #     subset(tT, adj.P.Val < .05 &
            #              abs(logFC) > 1),
            #     points(
            #       logFC,
            #       -log10(adj.P.Val),
            #       pch = 16,
            #       cex = 0.45,
            #       col = "red"
            #     )
            #   )
            #   #with(
            #   #  subset(tT, adj.P.Val < .05 &
            #   #           abs(logFC) > 1)[1:20, ],
            #   #  textxy(
            #   #    logFC,
            #   #    -log10(adj.P.Val),
            #   #    labs = rownames(tT),
            #   #    cex = .6
            #   #  )
            #   #)
            # }
            # dev.off()
            outfile <- paste(f, ".txt", sep = "")
            cat(paste(outfile), "\n")
            diffList <- c(diffList, list(tT, paste(f)))
          }
          v$diffList <- diffList
          output$comparison <- renderUI({
            if (is.null(v$diffList)) {
              return(NULL)
            }
            selectInput(
              "comparison" ,
              "View differential expression results:",
              choices = v$diffList[c(FALSE, TRUE)],
              selected = v$diffList[[2]]
            )
          })
          showNotification("Analysis done. You can now navigate to any tab to see your results",
                           duration = 0)
          
          #output$pcutoff <- renderUI({
          #  if (is.null(v$diffList)) {
          #    return(NULL)
          #  }
          #  textInput("pcutoff", "Enter an adj. P value cutoff for the differential expression table", ".05")
          #})
          
          output$diffExprs <- DT::renderDataTable({
            if (is.null(input$comparison) || is.null(input$pcutoff)) {
              return(NULL)
            } else {
              indx <- grep(input$comparison, v$diffList[c(FALSE,TRUE)])
              indx = indx*2
            }
            #v$currDiffList <- v$diffList[[indx-1]]
            out <-
              data.frame(rownames(v$diffList[[indx - 1]]), v$diffList[[indx - 1]])
            colnames(out)[1] = "ID"
            v$currDiffList <- out
            cat("HERE1")
            out <-
              out[which(out$adj.P.Val < as.double(input$pcutoff)), ]
            show(which(out$adj.P.Val < as.double(input$pcutoff)))
            if ("Entrez.ID" %in% colnames(out)) {
              out$Entrez.ID <-
                sapply(strsplit(as.character(out$Entrez.ID), "G"), "[[", 1)
              #show$Entrez.ID[1:2000] <- paste0(a(show$Entrez.ID[1:2000],href=paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=", show$Entrez.ID[1:2000], sep = ""), target="_blank"))
              out$Entrez.ID <-
                paste0(a(
                  out$Entrez.ID,
                  href = paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",out$Entrez.ID,sep = ""),
                  target = "_blank"
                ))
            }
            out
          }, rownames = FALSE, escape = FALSE)
          
          output$diffTitle <- renderUI({
            if (is.null(input$comparison)) {
              return(NULL)
            } else {
              indx <- grep(input$comparison, v$diffList[c(FALSE,TRUE)])
              indx = indx*2
            }
            h4(v$diffList[[indx]])
          })
          
          output$downloadDiff <- renderUI ({
            downloadButton('downloadDiffExprs', 'Download Table')
          })
          
          output$downloadDiffExprs <- downloadHandler(
            filename = function() {
              paste(input$comparison,
                    "_differential_expression.txt",
                    sep = "")
            },
            content = function(file) {
              show(head(v$currDiffList))
              cat("IN DOWNLOAD")
              v$diffdownloaddone <- TRUE
              write.table(v$currDiffList ,
                          file,
                          sep = "\t",
                          row.names = F, quote = FALSE)
            }
          )
          
          
          output$diff_Title <- renderUI({
            helpText("Pathway analysis is not currently available for user uploaded data tables")
          })
          
          # FUTURE GAGE CODE
          # output$species <- renderUI({
          #   #helpText("Pathway analysis is not currently available for user uploaded data tables")
          #   selectInput("org2", "Organism: ", c("Homo sapiens", "Mus muculus", "Rattus norvegicus", "Saccharomyces cerevisiae"),selected = character(0))
          # })
          # 
          # output$gageSubmit <- renderUI ({
          #   if (is.null(input$org2)) {
          #     return(NULL)
          #   }
          #   actionButton("gage_go", "Submit")
          # }) 
          # 
          # observeEvent(input$gage_go, {
          #   if (!is.null(input$org2)) {
          #     v$organism <- input$org2
          #   }
          #   output$KEGGpathwaysUp <- DT::renderDataTable({
          #     if (is.null(input$comparison)) {
          #       return(NULL)
          #     }
          #     indx <- grep(input$comparison, v$diffList[c(FALSE,TRUE)])
          #     indx = indx*2
          #     v$currDiffList <- v$diffList[[indx - 1]]
          #     tT <- v$currDiffList
          #     RunPathwayAnalysis()
          #     temp.data <- as.data.frame(v$keggres$greater)
          #     v$kegg_greater_all <- temp.data
          #     show(head(temp.data))
          #     filtered.data <-
          #       temp.data[temp.data$q.val <= .5 &
          #                   !(is.na(temp.data$q.val)),]
          #     datatable(filtered.data, rownames = TRUE)
          #   })
          #   output$KEGGpathwaysDown <- DT::renderDataTable({
          #     temp.data <- as.data.frame(v$keggres$less)
          #     v$kegg_lesser_all <- temp.data
          #     filtered.data <-
          #       temp.data[temp.data$q.val <= .5 &
          #                   !(is.na(temp.data$q.val)),]
          #     datatable(filtered.data, rownames = TRUE)
          #   })
          # })
          
          output$species2 <- renderUI({
            #helpText("Pathway analysis is not currently available for user uploaded data tables")
            selectInput("org", "Organism: ", listOrganism(),selected = character(0))
          })
          
          
          output$IDtype <- renderUI({
            selectInput("IDtype", "ID Type in First Column of uploaded table: ", listIDType(), selected = "genesymbol")
          })
          
          output$wegSubmit <- renderUI ({
            if (is.null(input$IDtype)) {
              return(NULL)
            }
            actionButton("webg_go", "Submit")
          }) 
          
          observeEvent(input$webg_go, {
            cat(input$org)
            output$enrich_results <- renderDataTable({
              if (is.null(input$org) || length(input$org)<=0) {
                cat("NULL")
                cat(length(input$org))
                return(NULL)
              }
              v$organism <- input$org
              cat("v$org ", v$organism)
              if (is.null(input$comparison)) {
                cat("WHTA")
                return(NULL)
              }
              cat("IN ENRICH RESULT")
              #if (is.null(v$currDiffList)) {return(NULL)}
              indx <- grep(input$comparison, v$diffList[c(FALSE,TRUE)])
              indx = indx*2
              v$currDiffList <- v$diffList[[indx - 1]]
              RunWebGestaltAnalysis()
              if(length(v$enrichResult) <= 1) {return(NULL)}
              table_path <- paste("./Project_",input$comparison,"/enrichment_results_",input$comparison,".txt", sep = "") 
              dt <- read.csv(table_path, sep = "\t")
              
              dt$geneset <-
                paste0(a(
                  dt$geneset,
                  href = paste(dt$link),
                  target = "_blank"
                ))
              dt <- dt[,c(-3,-10)]
              datatable(dt, rownames = F, escape = F)
            })
          })
          
          output$GOSlim <- renderUI({
            if (is.null(v$currDiffList)) {return(NULL)}
            if(length(v$enrichResult) <= 1) {return(NULL)}
            source <- paste("goslim_summary_",input$comparison,".png", sep = "")
            img(src = source, height = 696, width = 1500)
          })
          
          output$downloadWeb <- renderUI ({
            if(!is.null(v$enrichResult))
              downloadButton('downloadWebgestalt', 'Download Full Analysis')
          })
          
          output$downloadWebgestalt <- downloadHandler(
            filename = "bart_webgestalt_results.html", 
            content = function(file) {
              file.copy(paste("./Project_",input$comparison,"/Report_",input$comparison,".html", sep = ""), file) 
            }
          )
          # RunWebGestaltAnalysis()
          
          #if (!is.null(input$org)) {
          #  v$organism <- input$org
          #}
          
          
          
          
          
          output$volcano <- renderPlotly ({
            if (is.null(input$comparison)) {
              return(NULL)
            }
            
            progress <- shiny::Progress$new()
            on.exit(progress$close())
            progress$set(message = "Generating Volcano Plot ", value = 0)
            progress$inc(0.3)
            indx <- grep(input$comparison, v$diffList[c(FALSE,TRUE)])
            indx = indx*2
            v$currDiffList <- v$diffList[[indx - 1]]
            tT <- v$currDiffList
            cat("NROW tT", nrow(tT), "\n")
            color = matrix(data = "black",
                           nrow = nrow(tT),
                           ncol = 1)
            for (i in 1:nrow(tT)) {
              if (tT$adj.P.Val[i] < .05 && abs(tT$logFC[i]) > 1) {
                color[i,] = 'red'
              }
            }
            text = ''
            if ("Gene.symbol" %in% colnames(tT)) {
              text = tT$Gene.symbol
            } else {
              text = rownames(tT)
            }
            progress$inc(0.3)
            g = ggplot(tT) + geom_point(aes(
              x = logFC,
              y = -log10(adj.P.Val),
              text = text
            ), color = color)
            v$gg <- ggplotly(g) %>% layout(title = input$comparison,
                                           margin = list(t = 40))
            htmlwidgets::saveWidget(v$gg, file = paste(input$comparison, '_volcanoplot.html', sep = ""))
            showNotification(
              "It may take an additional ten seconds to display plot",
              type = "message",
              duration = 10
            )
            v$gg
          })
          
          output$downloadVol <- renderUI({
            if (!is.null(v$gg)) {
              downloadButton('downloadVolcano', 'Download Volcano Plot')
            }
          })
          
          output$downloadVolcano <- downloadHandler(
            filename = function() {
              paste(input$comparison, '_volcanoplot.html', sep = "")
            },
            content = function(file) {
              cat(paste(input$comparison, "_volcanoplot.html", sep = ""))
              source = paste(input$comparison, "_volcanoplot.html", sep = "")
              file.copy(source, file, overwrite = T)
            }
          )
          
          
          #Faster volcano plot, but not labelled 
          #  output$volcano <- renderUI({
          #    if (is.null(input$comparison)) {
          #      return(NULL)
          #    }
          #    source <-
          #      paste(input$comparison, "_volcanoplot.png", sep = "")
          #    img(src = source)
          #  })
          #  output$downloadVolcano <- downloadHandler(
          #    filename = function() {
          #      paste(input$comparison, '_volcanoplot.png', sep = "")
          #    },
          #    content = function(file) {
          #      cat(paste(input$comparison, "_volcanoplot.png", sep = ""))
          #      source = paste("www/", input$comparison, "_volcanoplot.png", sep = "")
          #      file.copy(source, file, overwrite = T)
          #    }
          # )
          
          output$plotPCA <- renderPlotly ({
            if (is.null(v$inexprs)) {
              return(NULL)
            }
            cat("\nPCA\n")
            exprs <- v$inexprs
            data.matrix <- as.matrix(exprs)
            
            data.matrix = as.matrix(data.matrix[order(rowMeans(data.matrix), decreasing =
                                                        TRUE)[1:1000], ])
            t = t(data.matrix)
            t = data.frame(t, "labels" = colnames(data.matrix))
            grouping = samples
            p <-
              ggplot(prcomp(t[1:(ncol(t) - 1)])) + geom_point(aes(
                x = PC1,
                y = PC2,
                colour = grouping,
                text = t$labels
              ))
            v$p <- ggplotly(p) # %>%
            htmlwidgets::saveWidget(v$p, file = paste('PCA.html', sep = ""))
            v$p
          })
          
          output$downloadPCA <- renderUI ({
            if(!is.null(v$p)) {
              downloadButton('downloadP', 'Download PCA Plot')
            }
          })
          
          output$downloadP <- downloadHandler(
            filename = function() {
              paste('PCA.html', sep = "")
            },
            content = function(file) {
              cat(paste("PCA.html", sep = ""))
              source = paste("PCA.html", sep = "")
              file.copy(source, file, overwrite = T)
            }
          )
          
          output$heatmap <- renderUI ({
            if (is.null(v$inexprs)) {
              return(NULL)
            }
            cat("Clustering\n")
            progress <- shiny::Progress$new()
            on.exit(progress$close())
            progress$set(message = "Clustering", value = 0)
            progress$inc(0.10, detail = "this may take a few seconds")
            data <- as.matrix(v$inexprs)
            currNames <- colnames(data)
            for (i in 1:length(colnames(data))) {
              colnames(data)[i] <- paste(samples[i], " ", currNames[i])
            }
            data2 = data[apply(data, 1, function(x)
              var(x) > 0.005),] + 1
            data2 = as.matrix(data2[order(rowMeans(data2), decreasing =
                                            TRUE)[1:1000], ])
            progress$inc(0.20, detail = "")
            my_palette <-
              colorRampPalette(c("green", "black", "red"))(n = 299)
            HMfile <- paste("www/heatmap.png", sep = "")
            png(HMfile,
                res = 200,
                height = 1200,
                width = 1200)
            heatmap.2(
              data2,
              margins = c(10, 10) ,
              trace = "none",
              col = my_palette,
              cexRow = 0.1,
              cexCol = 0.5,
              distfun = function(c)
                as.dist(1 - cor(t(c))),
              scale = "row",
              labRow = c(""),
              key.title = "",
              hclustfun = function(x)
                hclust(x, method = "ward.D2")
            )
            dev.off()
            v$clust_done <- TRUE
            cat("Cluster DONE\n")
            source <- paste("heatmap.png", sep = "")
            img(src = source,
                height = 900,
                width = 900)
          })
          
          
          output$downloadHeat <- renderUI ({
            if (!is.null(v$clust_done)) {
              downloadButton('downloadHM', 'Download Heatmap')
            }
          })
          
          output$downloadHM <- downloadHandler(
            filename = function() {
              paste('Heatmap.png', sep = '')
            },
            content = function(file) {
              source = paste("www/heatmap.png", sep = "")
              file.copy(source, file, overwrite = T)
            }
          )
          
        }
      })
    })
    
    #if a GEO is going to be used
    observeEvent(input$goButton, {
      if (!is.null(v$gse) || !is.null(v$inexprs)) {
        showModal(
          modalDialog(
            title = "Hey!",
            "Please refresh the page if you want to perform a different analysis"
          )
        )
        Sys.sleep(20)
        return(NULL)
        # session$reload()
      } else {
        cat("ok")
      }
      # Progress Bar
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Performing ", value = 0)
      progress$inc(0.2, detail = "initialization")
      progress$inc(0.05, detail = "Downloading data from GEO")
      GSEdownload()
      if (is.null(v$plat_choices)) {
        showModal(
          modalDialog(
            title = "GEO download failed",
            "Please ensure that your GSE accession number is correct and that it has a series matrix file which is downloadable on GEO. If so, refresh the page and try again in a few minutes. (The GEO server is likely down)"
          )
        )
        return(NULL)
      }
      output$choosePlatform <- renderUI({
        if (is.null(v$plat_choices)) {
          return(NULL)
        }
        selectInput(
          "choosePlatform" ,
          "Platform:",
          choices = v$plat_choices,
          selected = v$plat_choices[1]
        )
      })
      
      
      observeEvent(input$choosePlatform, {
        progress$inc(0.40, detail = "Parsing data from GEO")
        pdata()
        #makes checkbox to choose phenodata to group
        output$groups <- renderUI({
          cat("\nIN MAKING GROUPS")
          v$currgroups <- colnames(v$pdataALL)
          choice = c()
          curcol = c()
          for( i in 1:ncol(v$pdataALL) ) {
            cat(colnames((v$pdataALL)))
            cat("I", i)
            curcol = v$pdataALL[,i]
            if (length(unique(curcol)) > 1 && !any(table(curcol)<2) ) {
              choice = c(choice, colnames(v$pdataALL)[i])
            } 
            if (!is.na(colnames(v$pdataALL)[i]) && (colnames(v$pdataALL)[i] == "manual entry" || colnames(v$pdataALL)[i] == "batch effect")) {
              choice = c(choice, colnames(v$pdataALL)[i])
            }
          }
          sampNot <-
            showNotification(
              "Please go to the 'Sample grouping' tab to proceed with analysis",
              type = "message",
              duration = 20
            )
          
          checkboxGroupInput(
            "groups",
            "Group based on:",
            choices = choice,
            selected = 1
          )
        })
      }, priority = 2)
      if (!is.null(v$pdata)) {
        sampNot <-
          showNotification(
            "Please go to the 'Sample grouping' tab to proceed with analysis",
            type = "message",
            duration = 0
          )
      } else {
        sampNot = NULL
      }
      observe({
        
        v$idx
        if (is.null(v$tableidx) || is.null(v$pdata)) {
          return(NULL)
        }
        if (v$idx != v$tableidx) {
          cat("\nUPDATING\n")
          
          updateCheckboxGroupInput(
            session,
            "groups",
            "Group based on:",
            choices = colnames(v$pdataALL),
            selected = 1
          )
          cat(input$groups)
          v$tableidx = v$idx
        }
        
      }, priority = 100)
      
      output$hot <- renderRHandsontable({
        # if (is.null(v$dataALL)) {
        #  return(NULL)
        #}
        if ("manual entry" %in% colnames(v$dataALL)) {
          show(v$dataALL)
          cat("HERE 1")
          rhandsontable(
            v$dataALL,
            useTypes = F,
            readOnly = TRUE,
            columnSorting = TRUE,
            rowHeaderWidth = 100
          ) %>%
            hot_col(col = c("manual entry"), readOnly = FALSE) %>%
            hot_col(col = c("batch effect"), readOnly = FALSE) %>%
            hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
        } else {
          cat("HERE !")
          rhandsontable(
            v$dataALL,
            useTypes = F,
            readOnly = TRUE,
            columnSorting = TRUE,
            rowHeaderWidth = 100
          ) %>%
            hot_col(col = c("batch effect"), readOnly = FALSE) %>%
            hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
        }
      })
      
      output$groupsSubmit <- renderUI ({
        # if (is.null(input$groups) || is.null(input$hot)) {
        if (is.null(input$groups)) {
          return(NULL)
        }
        actionButton("groupsSubmit", "Submit")
      })
      
      output$useCEL <- renderUI({
        if (is.null(input$groups)) {
          return(NULL)
        }
        radioButtons(
          "useCEL",
          "Select which starting material to use for analysis:",
          c(
            "CEL files (if available)" = "CEL",
            "Expression table" = "Expres_table"
          ),
          selected = "CEL"
        )
      })
      
      observeEvent(input$groupsSubmit, {
        group_list = c()
        if (is.null(input$hot)) {
          cat("INPUT$hot is NULL")
          v$pdata <- v$dataALL[,input$groups, drop = FALSE]
        } else {
          v$pdata <- hot_to_r(input$hot)[, input$groups, drop = FALSE]
        }
        v$batch <- NULL
        v$num_batch <- NULL
        if ('batch effect' %in% colnames(v$pdata)) {
          v$batch <- v$pdata[, c("batch effect")]
          v$num_batch <- length(unique(v$batch))
          v$pdata <- v$pdata[,-c(ncol(v$pdata)), drop = FALSE]
        } 
        show(v$pdata)
        cat(nrow(v$pdata))
        for (i in 1:(nrow(v$pdata))) {
          group <- ""
          for (j in 1:ncol(v$pdata)) {
            g <- v$pdata[i, j]
            if (j == 1) {
              g <- gsub(" ", "", g)
              group <- paste(group, g, sep = "")
            }
            else {
              g <- gsub(" ", "", g)
              group <- paste(group, g, sep = "_")
            }
          }
          group <- make.names(group)
          group_list = c(group_list, group)
        }
        if (length(unique(group_list)) <= 1) {
          showModal(
            modalDialog(
              title = "Grouping not accepted",
              "At least two different groups are needed for testing, please enter a different grouping"
            )
          )
          cat("Less than two groups, please enter a different grouping")
          return(NULL)
        }
        
        if (any(as.double(table(group_list)) < 2)) {
          show(table(group_list))
          showModal(
            modalDialog(
              title = "Grouping not accepted",
              "Each group must have at least two replicates for testing, please enter a different grouping"
            )
          )
          cat("Each group must have at least two replicates for testing")
          return(NULL)
        }
        if ('manual entry' %in% colnames(v$pdata)) {
          if (length(unique(v$pdata[,'manual entry'])) <= 1) {
            v$pdata <- v$pdata[,-which(colnames(v$pdata) == "manual entry"), drop = FALSE]
            show(v$pdata)
            showNotification("manual entry column is being ignored.")
          }
        }
        if (!is.null(v$num_batch) &&
            any(as.double(table(group_list)) < v$num_batch)) {
          cat(as.double(table(group_list)))
          #  showModal(
          #    modalDialog(title = "Grouping not accepted",
          #                "There cannot be more batches than replicates")
          #  )
          # return(NULL)
        }
        if (!is.null(v$num_batch) && v$num_batch <= 1) {
          showModal(
            modalDialog(title = "Grouping not accepted",
                        "There must be more than one batch. If you do not wish to use batch effect correction, do not check that option.")
          )
          return(NULL)
        }
        if (!is.null(sampNot)) {
          removeNotification(sampNot)
        }
        cat(file = stderr(), "Grouping based on:", paste(input$groups), "\n")
        showNotification(paste("Grouping based on:", paste(input$groups)))
        cat(file = stderr(), "Running Analysis\n")
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Performing ", value = 0)
        progress$inc(0.3, detail = "Normalization")
        normalizedExpression()
        cat(file = stderr(), "Finished Normalizing\n")
        
        output$heatmap <- renderUI ({
          if (is.null(v$exprs)) {
            return(NULL)
          }
          cat("Clustering\n")
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "Clustering", value = 0)
          progress$inc(0.10, detail = "this may take a few seconds")
          data <- as.matrix(v$exprs)
          currNames <- colnames(data)
          for (i in 1:length(colnames(data))) {
            colnames(data)[i] <- paste(v$samples[i], currNames[i])
          }
          data2 = data[apply(data, 1, function(x)
            var(x) > 0.005),] + 1
          data2 = as.matrix(data2[order(rowMeans(data2), decreasing =
                                          TRUE)[1:1000], ])
          progress$inc(0.20, detail = "this may take a few seconds")
          my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)
          # my_palette<- colorRampPalette(brewer.pal(9, 'YlOrRd'))(299)
          
          HMfile <-
            paste("www/", input$geoID, "_heatmap.png", sep = "")
          png(HMfile,
              res = 200,
              height = 1200,
              width = 1200)
          heatmap.2(
            data2,
            #main = "Clustering",
            margins = c(10, 10) ,
            trace = "none",
            col = my_palette,
            cexRow = 0.1,
            cexCol = 0.5,
            distfun = function(c)
              as.dist(1 - cor(t(c))),
            scale = "row",
            labRow = c(""),
            key.title = "",
            hclustfun = function(x)
              hclust(x, method = "ward.D2")
          )
          dev.off()
          cat("Cluster DONE\n")
          v$clust_done <- TRUE
          source <- paste(input$geoID, "_heatmap.png", sep = "")
          img(src = source,
              height = 900,
              width = 900)
        })
        
        output$downloadHeat <- renderUI ({
          if (!is.null(v$clust_done)) {
            downloadButton('downloadHM', 'Download Heatmap')
          }
        })
        
        
        output$downloadHM <- downloadHandler(
          filename = function() {
            paste(input$geoID, '_Heatmap.png', sep = '')
          },
          content = function(file) {
            source = paste("www/", input$geoID, "_heatmap.png", sep = "")
            file.copy(source, file, overwrite = T)
          }
        )
        
        output$plotPCA <- renderPlotly ({
          if (is.null(v$exprs)) {
            return(NULL)
          }
          cat("\nPCA\n")
          exprs <- v$exprs
          data <- as.matrix(exprs)
          data.matrix = as.matrix(data[order(rowMeans(data), decreasing =
                                               TRUE)[1:1000], ])
          t = t(data.matrix)
          t = data.frame(t, 'labels' = rownames(t))
          show(head(data.matrix))
          grouping = v$samples
          p <-
            ggplot(prcomp(t[1:(ncol(t) - 1)])) + geom_point(aes(
              x = PC1,
              y = PC2,
              colour = grouping,
              text = t$labels
            ))
          v$p <- ggplotly(p)
          htmlwidgets::saveWidget(v$p, file = paste('PCA.html', sep = ""))
          v$p
        })
        
        output$downloadPCA <- renderUI ({
          if(!is.null(v$p)) {
            downloadButton('downloadP', 'Download PCA Plot')
          }
        })
        
        output$downloadP <- downloadHandler(
          filename = function() {
            paste('PCA.html', sep = "")
          },
          content = function(file) {
            cat(paste("PCA.html", sep = ""))
            source = paste("PCA.html", sep = "")
            file.copy(source, file, overwrite = T)
          }
        )
        
        # output$PCAplot <- renderUI ({
        #   if (is.null(v$exprs) || is.null(v$ptable)) {
        #     return(NULL)
        #   }
        #   cat("\nPCA\n")
        #   exprs <- v$exprs
        #   nGroups <- length(unique(v$ptable$Target))
        #   legVec <- vector()
        #   for (i in 1:nGroups) {
        #     legVec <- c(legVec, c(1))
        #   }
        #   
        #   legSize <- 1
        #   if (max(nchar(paste(v$ptable$Target))) > 40) {
        #     legSize <- .5
        #   }
        #   data.matrix <- exprs
        #   data.PC = prcomp(t(data.matrix), scale. = TRUE)
        #   PCAfile <-
        #     paste("www/", input$geoID, "PCAplot.png", sep = "")
        #   png(PCAfile,
        #       res = 160,
        #       height = 1000,
        #       width = 1000)
        #   plot(data.PC$x[, 1:2],
        #        col = factor(v$samples),
        #        pch = 1)
        #   legend(
        #     "bottomright",
        #     legend = unique(v$ptable$Target),
        #     pch = legVec,
        #     col = factor(unique(v$samples)),
        #     cex = legSize,
        #     horiz = FALSE
        #   )
        #   dev.off()
        #   cat("PCA done")
        #   source = paste(input$geoID, "PCAplot.png", sep = "")
        #   img(src = source,
        #       height = 600,
        #       width = 600)
        #   
        # })
        
        if (!is.null(v$pdata)) {
          if (input$useCEL ==  "CEL" && v$useMatrix == TRUE) {
            cat("CEL file analysis failed, using GSE expression table instead")
            s = showNotification(
              "CEL file analysis failed, used GSE expression table instead",
              type = "warning",
              duration = 35
            )
          }
          progress$inc(0.29, detail = "Differential Expression")
          
          getDiffExpr()
          showNotification("Analysis done. You can now navigate to any tab to see your results",
                           duration = 0)
        }
        output$comparison <- renderUI({
          if (is.null(v$diffList)) {
            return(NULL)
          }
          selectInput(
            "comparison" ,
            "View differential expression results:",
            choices = v$diffList[c(FALSE, TRUE)],
            selected = v$diffList[[2]]
          )
        })
        output$diffExprs <- DT::renderDataTable({
          if (is.null(input$comparison) ||
              is.null(input$pcutoff) || input$pcutoff <= 0) {
            return(NULL)
          } else {
            indx <- grep(input$comparison, v$diffList[c(FALSE,TRUE)])
            indx = indx*2
            #indx <-
            #  (which(
            #    sapply(
            #      v$diffList[c(FALSE,TRUE)],
            #      FUN = function(X)
            #        input$comparison %in% X
            #    )
            #  ))
          }
          v$currDiffList <- v$diffList[[indx - 1]]
          show <- v$diffList[[indx - 1]]
          show <-
            show[which(show$adj.P.Val < as.double(input$pcutoff)), ]
          if ("Entrez.ID" %in% colnames(show) && nrow(show) > 0) {
            show$Entrez.ID[which((show$Entrez.ID == ''))] = paste0(as.character('NA'))
            show$Entrez.ID <-
              sapply(strsplit(as.character(show$Entrez.ID), "///"), "[[", 1)
            show$Entrez.ID = paste0(
              '<a ',
              'href=',
              paste(
                "http://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                show$Entrez.ID,
                sep = ''
              ),
              '>',
              show$Entrez.ID,
              '</a>'
            )
            #show$Entrez.ID[1:2000] <- paste0(a(show$Entrez.ID[1:2000],href=paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=", show$Entrez.ID[1:2000], sep = ""), target="_blank"))
          }
          show
        }, rownames = FALSE, escape = FALSE)
        
        output$diff_Title <- renderUI({
          if (is.null(input$comparison)) {
            return(NULL)
          } else {
            indx <- grep(input$comparison, v$diffList[c(FALSE,TRUE)])
            indx = indx*2
          }
          h4(v$diffList[[indx]])
        })
        
        output$downloadDiff <- renderUI ({
          downloadButton('downloadDiffExprs', 'Download Table')
        })
        
        output$downloadDiffExprs <- downloadHandler(
          filename = function() {
            paste(input$comparison,
                  "_differential_expression.txt",
                  sep = "")
          },
          content = function(file) {
            show(head(v$currDiffList))
            v$diffdownloaddone <- TRUE
            write.table(v$currDiffList ,
                        file,
                        sep = "\t",
                        row.names = F, quote = FALSE)
          }
        )
        
        output$exprs <- DT::renderDataTable({
          if (is.null(v$norm_exprs)) {
            return(NULL)
          }
          norm_table <- v$norm_exprs
          v$norm_table <- as.data.frame(norm_table, drop = F)
          show(head(norm_table))
          norm_table
        }, selection = list(mode = "single", selected = 1))
        
        output$x2 <- renderPlotly({
          if (is.null(v$norm_table)) {
            return(NULL)
          }
          s = input$exprs_rows_selected
          if (!length(s)) {
            return(NULL)
          }
          title = ''
          if ("Symbol" %in% colnames(v$norm_table)) {
            title = v$norm_table$Symbol[s]
            normPlot = v$norm_table[, -c(1)]
            m = v$norm_table[s[length(s)], -c(1)]
          } else {
            title = rownames(v$norm_table)[s]
            normPlot = v$norm_table
            m = v$norm_table[s[length(s)], ]
          }
          plot_ly(
            x = factor(colnames(normPlot), levels = colnames(normPlot)),
            y = as.double(m),
            type = "bar",
            color = factor(v$samples[order(v$samples)])
          ) %>%
            layout(title = title,
                   margin = list(b = 100))
          
        })
        
        output$downloadNormExprs <- renderUI ({
          if (!is.null(v$norm_exprs)) {
            downloadButton('downloadExprs', 'Download Table')
          }
        })
        
        output$downloadExprs <- downloadHandler(
          filename = function() {
            paste(input$geoID, "normalized_expression.txt", sep = "")
          },
          content = function(file) {
            write.table(
              v$norm_exprs,
              file,
              sep = "\t",
              row.names = T,
              quote = FALSE
            )
          }
        )
        
        
        output$volcano <- renderPlotly ({
          if (is.null(input$comparison)) {
            return(NULL)
          }
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "Generating Volcano Plot ", value = 0)
          progress$inc(0.3)
          indx <- grep(input$comparison, v$diffList[c(FALSE,TRUE)])
          indx = indx*2
          v$currDiffList <- v$diffList[[indx - 1]]
          tT = v$currDiffList
          color = matrix(data = "black",
                         nrow = nrow(tT),
                         ncol = 1)
          progress$inc(0.1)
          for (i in 1:nrow(tT)) {
            if (tT$adj.P.Val[i] < .05 && abs(tT$logFC[i]) > 1) {
              color[i,] = 'red'
            }
          }
          text = ''
          if ("Gene.symbol" %in% colnames(tT)) {
            text = tT$Gene.symbol
          } else {
            text = tT$ID
          }
          progress$inc(0.3)
          g = ggplot(tT) + geom_point(aes(
            x = logFC,
            y = -log10(adj.P.Val),
            text = text
          ), color = color)
          
          v$gg <- ggplotly(g) %>% layout(title = input$comparison,
                                         margin = list(t = 40))
          
          htmlwidgets::saveWidget(v$gg, file = paste(input$comparison, '_volcanoplot.html', sep = ""))
          progress$inc(0.2)
          showNotification(
            "It may take an additional ten seconds to display plot",
            type = "message",
            duration = 10
          )
          v$gg
        })
        
        
        output$downloadVol <- renderUI({
          if (!is.null(v$gg)) {
            downloadButton('downloadVolcano', 'Download Volcano Plot')
          }
        })
        
        output$downloadVolcano <- downloadHandler(
          filename = function() {
            paste(input$comparison, '_volcanoplot.html', sep = "")
          },
          content = function(file) {
            cat(paste(input$comparison, "_volcanoplot.html", sep = ""))
            source = paste(input$comparison, "_volcanoplot.html", sep = "")
            file.copy(source, file, overwrite = T)
          }
        )
        
        
        output$KEGGpathwaysUp <- DT::renderDataTable({
          if (is.null(input$comparison)) {
            return(NULL)
          }
          indx <- grep(input$comparison, v$diffList[c(FALSE,TRUE)])
          indx = indx*2
          #indx <-
          #  (which(
          #    sapply(
          #      v$diffList,
          #      FUN = function(X)
          #        input$comparison %in% X
          #    )
          #  ))
          v$currDiffList <- v$diffList[[indx - 1]]
          tT <- v$currDiffList
          RunPathwayAnalysis()
          temp.data <- as.data.frame(v$keggres$greater)
          v$kegg_greater_all <- temp.data
          show(head(temp.data))
          filtered.data <-
            temp.data[temp.data$q.val <= .5 &
                        !(is.na(temp.data$q.val)),]
          datatable(filtered.data, rownames = TRUE)
        })
        output$KEGGpathwaysDown <- DT::renderDataTable({
          temp.data <- as.data.frame(v$keggres$less)
          v$kegg_lesser_all <- temp.data
          filtered.data <-
            temp.data[temp.data$q.val <= .5 &
                        !(is.na(temp.data$q.val)),]
          datatable(filtered.data, rownames = TRUE)
        })
        
        
        output$enrich_results <- renderDataTable({
          if (is.null(input$comparison)) {
            return(NULL)
          }
          #if (is.null(v$currDiffList)) {return(NULL)}
          indx <- grep(input$comparison, v$diffList[c(FALSE,TRUE)])
          indx = indx*2
          #indx <-
          #  (which(
          #    sapply(
          #     v$diffList,
          #      FUN = function(X)
          #        input$comparison %in% X
          #    )
          #  ))
          v$currDiffList <- v$diffList[[indx - 1]]
          RunWebGestaltAnalysis()
          if(length(v$enrichResult) <= 1) {return(NULL)}
          table_path <- paste("./Project_",input$comparison,"/enrichment_results_",input$comparison,".txt", sep = "") 
          dt <- read.csv(table_path, sep = "\t")
          
          dt$geneset <-
            paste0(a(
              dt$geneset,
              href = paste(dt$link),
              target = "_blank"
            ))
          dt <- dt[,c(-3,-10)]
          datatable(dt, rownames = F, escape = F)
        })
        
        output$GOSlim <- renderUI({
          if (is.null(v$currDiffList)) {return(NULL)}
          if(length(v$enrichResult) <= 1) {return(NULL)}
          source <- paste("goslim_summary_",input$comparison,".png", sep = "")
          img(src = source, height = 696, width = 1500)
        })
        
      })
      
      output$downloadGageUp <- renderUI ({
        if(!is.null(v$kegg_greater_all))
          downloadButton('downloadkegg_greater', 'Download Full Results')
      })
      
      output$downloadkegg_greater <- downloadHandler(
        filename = function() {
          paste(input$comparison,
                "_pathway_analysis_up.txt",
                sep = "")
        },
        content = function(file) {
          write.table(as.data.frame(v$kegg_greater_all),
                      file,
                      sep = "\t",
                      row.names = T, quote = FALSE)
        }
      )
      
      
      output$downloadGageDown <- renderUI ({
        if(!is.null(v$kegg_lesser_all)) {
          downloadButton('downloadkeggres', 'Download Full Results')
        }
      })
      
      
      output$downloadkeggres <- downloadHandler(
        filename = function() {
          paste(input$comparison,
                "_pathway_analysis_down.txt",
                sep = "")
        },
        content = function(file) {
          write.table(as.data.frame(v$kegg_lesser_all),
                      file,
                      sep = "\t",
                      row.names = T, quote = FALSE)
        }
      )
      
      output$downloadWeb <- renderUI ({
        if(!is.null(v$enrichResult))
          downloadButton('downloadWebgestalt', 'Download Full Analysis')
      })
      
      output$downloadWebgestalt <- downloadHandler(
        filename = "bart_webgestalt_results.html", 
        content = function(file) {
          file.copy(paste("./Project_",input$comparison,"/Report_",input$comparison,".html", sep = ""), file) 
        }
      )
      
    })
    output$downloadTestFile <- downloadHandler(
      filename = "GSE20986normalized_expression.txt",
      content = function(con) {
        write.table(
          TestFile,
          con,
          row.names = F,
          sep = "\t",
          quote = F
        )
      }
    )
    
    output$downloadTutorial <- downloadHandler(
      filename = "BART_Tutorial.pdf",
      content = function(con) {
        file.copy("./BART_Tutorial.pdf", con)
      }
    )
  }
  
  shinyApp(ui, server)
}
