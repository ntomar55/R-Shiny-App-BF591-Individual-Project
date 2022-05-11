## Author: Taylor Falk
## tfalk@bu.edu
## BU BF591
## Assignment 7

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(dplyr)
library(tidyverse)
library(ggrepel)
library(genefilter)
library(DT)
library(tidyverse)

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("BF591 Final Project"),
    h3("Nikita Tomar"),
    # HTML("<p>To use this application, download the CSV <code>deseq_res.csv</code> from the data directory of this app's repository.</p>"),
    verbatimTextOutput("debug"),
    tabsetPanel(
        # title = "NAV TITLE",
        tabPanel(
            "Sample Info",
            sidebarLayout(
                sidebarPanel(
                    fileInput("c1_fileinput", "Input Sample Information data",
                              accept = ".csv", placeholder = "sample_info.csv"),
                ),
                mainPanel(
                    # "Main Panel 1",
                    tabsetPanel(
                        tabPanel(
                            "Summary",
                            tableOutput("c1_table1")
                        ),
                        tabPanel(
                            "Table",
                            DTOutput("c1_table2")
                        ),
                        tabPanel(
                            "Plots",
                            sidebarLayout(
                                sidebarPanel(
                                    radioButtons("c1_radio1", "Choose the column name to plot",
                                                 choices = ""),
                                    submitButton("Plot", width = '100%')
                                ), # sidebarPanel
                                
                                mainPanel(
                                    plotOutput("c1_hist_plot")
                                ) # end of mainPanel
                            )
                        )
                    )
                )
            ),
            
        ),
        tabPanel(
            "Count Matrix",
            sidebarLayout(
                sidebarPanel(
                    # "Input Controls",
                    fileInput("c2_fileinput", "Load normalized count matrix results",
                              accept = ".csv", placeholder = "norm_count_matrix.csv"),
                    sliderInput(inputId = "c2_slider_ptile", min = 0, max = 100,
                                label = "Select the percentile of variance to filter genes.",
                                value = 50, step = 1),
                    sliderInput(inputId = "c2_slider_nz_sampl", min = 1, max = 2,
                                label = "Select X so as to include genes with at least X samples that are non-zero",
                                value = 2, step = 1),
                    submitButton("Plot", width = '100%')
                ),
                mainPanel(
                    tabsetPanel(
                        tabPanel("Summary", 
                                 textOutput("c2_sample_text"),
                                 tableOutput("c2_table")),
                        tabPanel("Scatter Plot",
                                 mainPanel(
                                     plotOutput("c2_scatter1"),
                                     plotOutput("c2_scatter2")
                                 )),
                        tabPanel(
                            "Heatmap",
                            sidebarLayout(
                                sidebarPanel(
                                    radioButtons("c2_radio_log", "Log scale options for heatmap",
                                                 choices = c("Enable", "Disable"),
                                                 selected = "Enable"),
                                ),
                                mainPanel(
                                    plotOutput("c2_heatmap")
                                )
                            ),
                        ),
                        tabPanel(
                            "PCA",
                            sidebarLayout(
                                sidebarPanel(
                                    radioButtons("c2_radio_pc1", "Select Principal Component (X)",
                                                 choices = c(1,2,3,4,5,6,7),
                                                 selected = 1),
                                    radioButtons("c2_radio_pc2", "Select Principal Component (Y)",
                                                 choices = c(1,2,3,4,5,6,7),
                                                 selected = 2),
                                    sliderInput(inputId = "c2_slider_pca", min = 1, max = 7,
                                                label = "Select number of top components for beeswarm plot",
                                                value = 5, step = 1),
                                ),
                                mainPanel(
                                    plotOutput("c2_scatter_pca"),
                                    # plotOutput("c2_beeswarm_pca")
                                )
                            ),
                        )
                    )
                )
            ),
            
        ),
        tabPanel(
            "Differential Expression",
            sidebarLayout(
                sidebarPanel(
                    fileInput("c3_fileinput", "Load differential expression results",
                              accept = ".csv", placeholder = "diff_exp.csv"),
                ),
                mainPanel(
                    tabsetPanel(
                        tabPanel(
                            "Table",
                            DTOutput("c3_table_de")
                        ),
                        tabPanel(
                            "Analyses",
                            sidebarLayout(
                                sidebarPanel(
                                    HTML(paste(rep("<p>A volcano plot can be generated with <b>'log<sub>2</sub> fold-change'</b> on the x-axis and <b>'p-adjusted'</b> on the y-axis.</p>"),
                                               collapse = "")),
                                    
                                    radioButtons("c3_radio_1", "Choose the column for the x-axis",
                                                 choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"),
                                                 selected = "log2FoldChange"),
                                    radioButtons("c3_radio_2", "Choose the column for the y-axis",
                                                 choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"),
                                                 selected = "padj"),
                                    
                                    colourInput("c3_col_1", label = "Base point color" ,"#22577A"),
                                    colourInput("c3_col_2", label = "Highlight point color" ,"#FFCF56"),
                                    
                                    sliderInput(inputId = "c3_slider_padj", min = -150, max = 0,
                                                label = "Select the magnitude of the p adjusted coloring:", 
                                                value = -75, step = 1),
                                    
                                    submitButton("Plot", width = '100%')
                                ), # sidebarPanel
                                
                                mainPanel(
                                    tabsetPanel(
                                        tabPanel("Plot", plotOutput("c3_volcano_plot")),
                                        tabPanel("Table", tableOutput("c3_data_table"))
                                    )
                                ) # end of mainPanel
                                
                            ) # end of sidebarLayout
                        )
                    )
                )
            ),
        ),
        tabPanel(
            "Visualization of Gene Expression",
            sidebarLayout(
                sidebarPanel(
                    fileInput("c4_fileinput1", "Load normalized count matrix data",
                              accept = ".csv", placeholder = "norm_count_matrix.csv"),
                    fileInput("c4_fileinput2", "Load normalized count matrix data",
                              accept = ".csv", placeholder = "sample_info.csv"),
                ),
                mainPanel(
                    sidebarLayout(
                        sidebarPanel(
                            radioButtons("c4_radio1", "Choose the categorical column name to plot",
                                         choices = ""),
                            # textInput("c4_search", "", placeholder = "Input gene name"),
                            # htmlOutput("c4_search_text"),
                            selectizeInput(
                                inputId = "c4_search_gene", 
                                label = "Input gene name",
                                multiple = FALSE,
                                choices = "",
                                #c("Search Bar" = "", paste0(LETTERS,sample(LETTERS, 26)))
                                options = list(
                                    create = FALSE,
                                    placeholder = "Search gene",
                                    maxItems = '1',
                                    onDropdownOpen = I("function($dropdown) {if (!this.lastQuery.length) {this.close(); this.settings.openOnFocus = false;}}"),
                                    onType = I("function (str) {if (str === \"\") {this.close();}}")
                                )
                            ),
                            radioButtons("c4_radio2", "Choose the type of plot",
                                         choices = c("bar plot", "boxplot", "violin plot", "beeswarm plot")),
                            submitButton("Plot", width = '100%')
                        ), # sidebarPanel
                        
                        mainPanel(
                            textOutput("c4_error"),
                            plotOutput("c4_plot")
                        ) # end of mainPanel
                    )
                )
            )
            
        ),
    ),
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    shinyOptions(cache = cachem::cache_disk(file.path(dirname(tempdir()), "myapp-cache")))
    options(shiny.maxRequestSize=30*1024^2)
    # output$debug <- renderPrint({
    #     sessionInfo()
    # })
    #' load_Data
    #' #'
    #' #' @details Okay this one is a little weird but bear with me here. This is 
    #' #' still a "function", but it will take no arguments. The `reactive({})` bit 
    #' #' says "if any of my inputs (as in, input$...) are changed, run me again". 
    #' #' This is useful when a user clicks a new button or loads a new file. In 
    #' #' our case, look for the uploaded file's datapath argument and load it with 
    #' #' read.csv. Return this data frame in the normal return() style.
    
    c1.load_data <- reactive({
        req(input$c1_fileinput)
        file <- input$c1_fileinput
        if (is.null(file)) {
            return(NULL)
        }
        else {
            data <- read.csv(file$datapath, header= TRUE, sep=",")
        }
        observe({ #for histogram
            dataf <- data %>%
                select_if(function(x) (is.integer(x) || is.numeric(x)) && (length(unique(x)) > 2))
            choices <- colnames(dataf) 
            updateRadioButtons(session, "c1_radio1", choices = choices)
        })
        
        # observe({ #buttons for DE
        #     dechoices <- get_de_names(data)
        #     updateRadioButtons(session, "de_var", choices = dechoices)
        # })
        return(data)
    })
    c1.radio1 <- reactive({
        # if (is.null(input$c1_radio1)) return(2)
        input$c1_radio1
    })
    output$c1_table1 <- renderTable({
        datam <- c1.load_data()
        summarize <- function(x){
            desc = ""
            col_type = class(x)
            if(col_type == "integer" || col_type == "numeric") {
                x <- na.omit(x)
                m <- round(mean(x),2)
                avg <- as.character(round(mean(x),2))
                sdev <- as.character(round(sd(x), 2))
                desc <- paste(c(avg, " (+/- ", sdev, ")"), collapse="")
            }
            else if((length(unique(x)) != length(x)) && (length(unique(x)) > 1) && (length(unique(x)) < 5)) {
                col_type = "factor"
                desc <- paste(levels(as.factor(unique(x))),collapse=",")
            }
            else if(col_type == "character") {
                if(length(unique(x)) == length(x)) {
                    desc <- "ID/Unique"
                } else {
                    desc <- x[1]
                }
            }
            return(list(col_type, desc))
        }
        col1 <- colnames(datam)
        col1
        length(col1)
        col2 <- lapply(datam, function(x){summarize(x)[[1]]})
        col3 <- lapply(datam, function(x){summarize(x)[[2]]})
        length(col2)
        length(col3)
        tbl <- as.data.frame(cbind(col1, col2, col3))
        colnames(tbl) <- c("Column Name", "Type", "Mean (sd) or Distinct Values")
        return(tbl)
    }, bordered = T)
    output$c1_table2 <- renderDT({
        datam <- c1.load_data()
        return(datam)
    }, options = list(rownames = TRUE))
    output$c1_hist_plot <- renderPlot({
        datam <- c1.load_data()
        plot <- ggplot(datam, aes(x=!!sym(input$c1_radio1))) + 
            geom_histogram(fill = "blue", color = "black", alpha=0.6) +
            theme(legend.position="bottom") +
            theme_minimal() 
        return(plot)
    })
    c2.load_data <- reactive({
        req(input$c2_fileinput)
        file <- input$c2_fileinput
        if (is.null(file)) {
            return(NULL)
        }
        else {
            data <- read.csv(file$datapath, header= TRUE, sep=",")
        }
        data <- data %>% rename(Gene = X)
        observe({ #for component 4
            updateSliderInput(session, "c2_slider_nz_sampl", max = dim(data)[2], value = as.integer(dim(data)[2]/2))
        })
        
        return(data)
    })
    c2.slider_ptile <- reactive({
        if (is.null(input$c2_slider_ptile)) return(90)
        input$c2_slider_ptile
    })
    c2.slider_nz_sampl <- reactive({
        if (is.null(input$c2_slider_nz_sampl)) return(2)
        input$c2_slider_nz_sampl
    })
    c2.radio_log <- reactive({
        if (is.null(input$c2_radio_log)) return("Enable")
        input$c2_radio_log
    })
    c2.radio_pc1 <- reactive({
        if (is.null(input$c2_radio_pc1)) return(1)
        input$c2_radio_pc1
    })
    c2.radio_pc2 <- reactive({
        if (is.null(input$c2_radio_pc2)) return(2)
        input$c2_radio_pc2
    })
    c2.slider_pca <- reactive({
        if (is.null(input$c2_slider_pca)) return(5)
        input$c2_slider_pca
    })
    c2_modify_data <- function(datam, ptile, nsampl) {
        m <- datam[c(-1)]
        rv <- rowVars(m)
        cutoff_v <- quantile(rowVars(m), ptile)
        idx = which(rv > cutoff_v & rowSums(m != 0) >= nsampl)
        datam$Pass = 'fail'
        datam$Pass[idx] = 'pass'
        datam$NumZero = rowSums(m == 0)
        datam$Median = apply(m, 1, median)
        datam$RowVar = rv
        datam = relocate(datam, c(Pass, NumZero, Median, RowVar), .after = 1)
        return(datam)
    }
    output$c2_table <- renderTable({
        datam = c2.load_data()
        m = datam[c(-1)]
        ngene = dim(m)[1]
        nsampl = dim(m)[2]
        datam <- c2_modify_data(datam, c2.slider_ptile()/100, c2.slider_nz_sampl())
        datam.filt = datam[datam$Pass=='pass',]
        m2 <- datam.filt[-c(1,2,3,4,5)]
        npass = dim(m2)[1]
        nfail = ngene - npass
        tbl = data.frame(X=c('Total', 'Pass', 'Fail'), Genes=c(ngene, npass, nfail),'Percentage (%)'=c(100, round(npass/ngene,2), round(nfail/ngene,2)))
        output$c2_sample_text <- renderText({ 
            paste(paste(c("Total samples in dataset ", dim(m)[2]),collapse=""),
                  paste(c("Filter criteria 1: Genes with at least ", c2.slider_ptile()/100, " percentile variance"),collapse=""), 
                  paste(c("Filter criteria 2: Genes with at least ", c2.slider_nz_sampl(), " non-zero samples"), collapse=""), collapse='\n')
        })
        return(tbl)
    })
    output$c2_scatter1 <- renderPlot({
        datam = c2.load_data()
        datam = c2_modify_data(datam, c2.slider_ptile()/100, c2.slider_nz_sampl())
        plot <- ggplot(data = datam, aes(x = log(Median+1), y=NumZero)) +
            geom_point(aes(color = Pass)) +
            labs( color = c('Filtered')) +
            theme(legend.position = "bottom") +
            scale_color_manual(values = c("black", "red")) +
            theme_minimal()
        return(plot)
    })
    output$c2_scatter2 <- renderPlot({
        datam = c2.load_data()
        datam = c2_modify_data(datam, c2.slider_ptile()/100, c2.slider_nz_sampl())
        plot <- ggplot(data = datam, aes(x = log(Median+1), y=log(RowVar))) +
            geom_point(aes(color = Pass)) +
            labs( color = c('Filtered')) +
            theme(legend.position = "bottom") +
            scale_color_manual(values = c("black", "red")) +
            theme_minimal()
        return(plot)
    })
    output$c2_heatmap <- renderPlot({
        datam <- c2.load_data()
        datam <- c2_modify_data(datam, c2.slider_ptile()/100, c2.slider_nz_sampl())
        filt_data <- datam[datam$Pass=='pass',]
        filt_data <- filt_data[-c(1,2,3,4,5)]
        label_name = "counts"
        if(c2.radio_log() == 'Enable') {
            filt_data <- log(filt_data+1)
            label_name <- "log(counts)"
        }
        plot <- Heatmap(filt_data, cluster_rows = F, name=label_name, show_row_names = F, heatmap_legend_param = list())
        return(plot)
    })
    output$c2_scatter_pca <- renderPlot({
        datam <- c2.load_data()
        datam <- c2_modify_data(datam, c2.slider_ptile()/100, c2.slider_nz_sampl())
        filt_data <- datam[datam$Pass=='pass',]
        filt_data <- filt_data[-c(1,2,3,4,5)]
        print(c2.slider_ptile())
        print(c2.slider_nz_sampl())
        print(dim(filt_data))
        print('filt data')
        print(head(filt_data))
        if(dim(filt_data)[1]==0 || dim(filt_data)[2]==0) {
            return(NULL)
        }
        pca <- prcomp(filt_data)
        pca_df <- as.data.frame(pca$x)
        pca_var <- sapply(pca$sdev, function(x) x^2 )
        pcnt_var <- sapply(pca_var, function(x) x/sum(pca_var))
        
        PC1 <- as.integer(c2.radio_pc1())
        PC2 <- as.integer(c2.radio_pc2())
        x_lab <- paste0("PC", PC1, ", variance: ", as.character(round(pcnt_var[PC1]*100,2)), "%")
        y_lab <- paste0("PC", PC2, ", variance: ", as.character(round(pcnt_var[PC2]*100,2)), "%")
        
        plot <- ggplot(pca_df, aes(x=!!sym(str_c("PC", PC1)), y=!!sym(str_c("PC", PC2)))) +
            geom_point(stat="identity") +
            xlab(x_lab) +
            ylab(y_lab) +
            ggtitle("PCA Scatter Plot") +
            theme(legend.position="bottom") +
            theme_minimal()
        PCN <- as.integer(c2.slider_pca())
        # output$c2_beeswarm_pca <- renderPlot({
        #     bw_plot <- beeswarm(log(pca_df[,1:PCN]+1), pch = 20, col = 2:PCN+1, ylab = "log(PC+1)")
        #     return(bw_plot)
        # })
        return(plot)
    })
    c3.load_data <- reactive({
        req(input$c3_fileinput)
        file <- input$c3_fileinput
        if (is.null(file)) {
            return(NULL)
        }
        else {
            data <- read.csv(file$datapath, header= TRUE, sep=",")
        }
        data <- data %>% rename(Gene = X)
        return(data)
    })
    c3.radio_1 <- reactive({
        if (is.null(input$c3_radio_1)) return("log2FoldChange")
        input$c3_radio_1
    })
    c3.radio_2 <- reactive({
        if (is.null(input$c3_radio_2)) return("padj")
        input$c3_radio_2
    })
    c3.slider_padj <- reactive({
        if (is.null(input$c3_slider_padj)) return(-210)
        input$c3_slider_padj
    })
    c3.col_1 <- reactive({
        if (is.null(input$c3_col_1)) return("blue")
        input$c3_col_1
    })
    c3.col_2 <- reactive({
        if (is.null(input$c3_col_2)) return("green")
        input$c3_col_2
    })
    volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
        plot <- ggplot(data = dataf, aes(x = !!sym(x_name), y=-log10(!!sym(y_name)))) +
            geom_point(aes(color = padj< 1*10^(slider))) +
            labs( color = str_glue('{y_name} 1 x 10^ {slider}')) +
            theme(legend.position = "bottom") +
            theme_minimal()
        return(plot)
    }
    draw_table <- function(dataf, slider) {
        dataf <- arrange(dataf, pvalue)
        dataf <- filter(dataf, between(padj, 0, 10 ^slider ))
        dataf <- mutate(dataf, pvalue = formatC(dataf$pvalue, digits = 2, format = "E" ),
                        padj = formatC(dataf$padj, digits = 2, format = "E"))
        return(dataf)
    }
    output$c3_table_de <- renderDT({
        datam <- c3.load_data()
        return(datam)
    }, options = list(rownames = TRUE))
    output$c3_volcano_plot <- renderPlot({
        volcano_plot(dataf = c3.load_data(), slider = c3.slider_padj(), 
                     x_name=c3.radio_1(), y_name=c3.radio_2(), 
                     color1=c3.col_1(), color2=c3.col_2())
    }) # replace this NULL
    # Same here, just return the table as you want to see it in the web page
    output$c3_data_table <- renderTable({
        draw_table(dataf = c3.load_data(), slider = c3.slider_padj())
    }) # replace this NULL
    genetext = " "
    c4.load_data1 <- reactive({
        req(input$c4_fileinput1)
        file <- input$c4_fileinput1
        if (is.null(file)) {
            return(NULL)
        }
        else {
            data <- read.csv(file$datapath, header= TRUE, sep=",")
        }
        data <- data %>% rename(Gene = X)
        observe({ #for component 4
            # genetext <- paste(data[,1], collapse = ' ')
            updateSelectizeInput(inputId="c4_search_gene", choices = c("Search Bar" = "", sample(data[,1],5000)))
        })
        return(data)
    })
    c4.load_data2 <- reactive({
        req(input$c4_fileinput2)
        file <- input$c4_fileinput2
        if (is.null(file)) {
            return(NULL)
        }
        else {
            data <- read.csv(file$datapath, header= TRUE, sep=",")
        }
        observe({ #for histogram
            dataf <- data_info %>%
                select_if(function(x) ((length(unique(x)) != length(x)) && (length(unique(x)) > 1) && (length(unique(x)) < 5)
                ))
            choices <- colnames(dataf) 
            updateRadioButtons(session, "c4_radio1", choices = choices)
        })
        return(data)
    })
    c4.radio1 <- reactive({
        # if (is.null(input$c4_radio1)) return("padj")
        input$c4_radio1
    })
    c4.radio2 <- reactive({
        if (is.null(input$c4_radio2)) return("bar plot")
        input$c4_radio2
    })
    # insert_mark_tag <- function(s, loc_index, all_locs) {
    #     str_sub(s, all_locs[loc_index, 2] + 1, all_locs[loc_index, 2]) <- "</mark>"
    #     str_sub(s, all_locs[loc_index, 1], all_locs[loc_index, 1] - 1) <- "<mark>"
    #     s
    # }
    # output$c4_search_text <- renderText(HTML(
    #     if (nchar(input$c4_search_text)>0)
    #         str_replace_all(genetext, sprintf("(%s)", input$c4_search_text), "<mark>\\1</mark>") else
    #             genetext
    # ))
    # Show Selected Value in Console
    observe({
        print(input$searchme)
    })
    
    get_plot <- function(dat, type, colname) {
        plot <- NULL
        print(head(dat))
        if(type=="bar plot") {
            plot <- ggplot(dat, aes_string(x=colname, fill=colname )) + 
                geom_bar( ) +
                theme_minimal()
        } else if(type=="boxplot") {
            plot <- ggplot(dat, aes_string(x=colname, y="GE_Values", fill=colname)) + 
                geom_boxplot() +
                theme_minimal()
        } else if(type=="violin plot") {
            plot <- ggplot(dat, aes_string(x = colname, y = "GE_Values", fill = colname)) +
                geom_violin(trim = FALSE) +
                geom_boxplot(width = 0.07) +
                theme_minimal()
        } else if(type=="beeswarm plot") {
            plot <- ggplot(dat, aes_string(x = colname, y = "GE_Values", color = colname)) +
                geom_beeswarm(cex = 3) +
                theme_minimal()
        }
        return(plot+theme(legend.position="bottom"))
    }
    output$c4_plot <- renderPlot({
        data_cnt <- c4.load_data1()
        data_info <- c4.load_data2()
        gene <- input$c4_search_gene
        genelist <- data_cnt[,1]
        if(gene=="" || gene %in% genelist == FALSE) {
            output$c4_error <- renderText({
                paste(c("The input gene value ",
                        input$c4_search_gene, 
                        " is either invalid or not present in the count matrix. Please try another gene name"),
                      collapse = "")})
            return(NULL)
        }
        output$c4_error <- renderText({""})
        idx <- which(genelist == gene)
        # observeEvent(input$c4_radio1,{ 
        #     colname <- input$c4_radio1
        # })
        colname <- c4.radio1()
        
        plot_data <- data_cnt[idx,-c(1)] %>% t() %>% as.data.frame() %>%
            mutate(as.factor(data_info[,colname]))
        
        colnames(plot_data) <- c("GE_Values", colname)
        plot <- get_plot(plot_data, c4.radio2(), colname)
        return(plot)
    })
}


# Run the application
shinyApp(ui = ui, server = server)

