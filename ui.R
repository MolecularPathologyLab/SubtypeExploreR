
### Options for Spinner
options(spinner.color="#00a65a", spinner.color.background="#ffffff", spinner.size=0.5)

# build ui.R ----
## 1. header ----
header <- 
    dashboardHeader(
        title = "SubtypeExploreR"
    )


## 2. sidebar ----
sidebar <-
    dashboardSidebar(
        
        sidebarMenu(
            menuItem("Home", tabName = "home", icon = icon("home")),
            menuItem("Expression Boxplots", tabName = "boxplots", icon = icon("chart-column")),
            menuItem("ssGSEA", tabName = "ssgsea", icon = icon("chart-line")),
            
            helpText("Developed by",
                     br(),
                     a("Ryan Byrne", href="https://www.researchgate.net/profile/Ryan-Byrne-14", target="_blank"),
                     br(),
                     "&",
                     br(),
                     a("Sudhir Malla", href="https://www.researchgate.net/profile/Sudhir-Malla-2", target="_blank"),
                     br(),
                     br(),
                     a("Molecular Pathology", href="https://dunne-lab.com/", target="_blank"),
                     br(),
                     a("Research Group", href="https://dunne-lab.com/", target="_blank"),
                     br(),
                     "Queen's University Belfast",
                     style="padding-left:2em; 
                padding-right:2em;
                position:absolute; 
                bottom:1em;
                text-align:center")
        )
    )


## 3. body ----
body <- 
    dashboardBody(
        tabItems(
            
            ### HOME TAB ----
            tabItem(tabName = "home",
                    fluidRow(
                        tags$iframe(src = './introduction.html', 
                                    width = '100%', height = '800px', 
                                    frameborder = 0, scrolling = 'auto'
                        )
                    )
            ),
            
            
            ### BOXPLOTS TAB ----
            tabItem(tabName = "boxplots",
                    
                    fluidRow(
                        column(2,
                               box(width = NULL, status = "success", solidHeader = TRUE,
                                   title = "Boxplot Parameters",
                                   div(style="display: inline-block;vertical-align:top; width: 100%;",
                                       selectInput("boxplot_gene", "Select Gene:", 
                                                   choices = NULL, 
                                                   selected = NULL, multiple = FALSE)),
                                   ### Select variable-of-interest
                                   div(style="display: inline-block;vertical-align:top; width: 100%",
                                       selectInput("boxplot_var", "Select Subtype:",
                                                   choices = c("PDS", "CMS", "iCMS"), 
                                                   selected = "PDS", multiple = FALSE)),
                                   # Button user clicks to view gene expression boxplots after selecting an
                                   # existing gene list/signature or entering their own gene list
                                   div(style="display: inline-block;vertical-align:top;",
                                       actionButton("boxplot_submit", "View Expression"))
                               )),
                        ### FOCUS Boxplot Output
                        column(5,
                               box(width = NULL, status = "success", solidHeader = TRUE,
                                   title = "FOCUS",
                                   withSpinner(
                                       plotOutput("focus_gene_exp_boxplot"),
                                       type = 1, 
                                       size = 0.5
                                   )
                               )),
                        
                        ### GSE39582 Boxplot Output
                        column(5,
                               box(width = NULL, status = "success", solidHeader = TRUE,
                                   title = "GSE39582",
                                   withSpinner(
                                       plotOutput("gse39582_gene_exp_boxplot"),
                                       type = 1, 
                                       size = 0.5
                                   )
                               ))
                    )
            ),
            
            
            ### ssGSEA TAB ----
            tabItem(tabName = "ssgsea",
                    
                    fluidRow(
                        ### Example or User dataset selectiion
                        column(3,
                               box(width = NULL, status = "success", solidHeader = TRUE,
                                   title = "ssGSEA Parameters",
                                   ### Select gene-of-interest
                                   div(style="display: inline-block;vertical-align:top; width: 100%;",
                                       selectInput("ssgsea_collection", "Select Gene Set Collection:", 
                                                   choices = c("Custom",
                                                               "Hallmark",
                                                               "KEGG",
                                                               "BioCarta",
                                                               "Reactome",
                                                               "PID"),
                                                   selected = "Hallmark", # Hallmark selected on initialisation
                                                   multiple = FALSE)),
                                   
                                   conditionalPanel(
                                       condition = "input.ssgsea_collection != 'Custom'",  # condition so panel only displays if user wishes to use an existing gene set
                                       # Input to select the gene set form a drop down menu
                                       # Options in drop down menu change based on the
                                       # ontology/gene set collection select in selectInput above
                                       div(style="display: inline-block;vertical-align:top; width: 100%;",
                                           selectizeInput("gene_set", 
                                                          label = "Select Gene Set:",
                                                          choices = NULL))
                                   ),
                                   # Panel that displays if user wishes to use an existing
                                   # signature/gene list
                                   conditionalPanel(
                                       condition = "input.ssgsea_collection == 'Custom'",  # condition so panel only displays if user wishes to use their own gene list
                                       
                                       div(style="display: inline-block;vertical-align:top; width: 160px;",
                                           # Input to enter a list of gene symbols to use for GSEA
                                           # Each gene symbol should be on a new line
                                           textAreaInput("ssgsea_gene_list",
                                                         "Enter a list of gene symbols (as shown)",
                                                         value = "", # value on initialisation
                                                         width = NULL,
                                                         placeholder = "TP53\nCCND2\nCDKN2A\nMDM2", # Example gene list which user will see
                                                         rows = 10 # Number of rows to display at once. If more genes (rows) entered a scroll bar appears
                                           )),
                                       br(),
                                       div(style="display: inline-block;vertical-align:top; width: 160px;",
                                           # Input to enter a name for the gene list entered
                                           # This is optional and if entered the name is used as
                                           # title of GSEA plots
                                           textInput(inputId = "custom_gene_set_name",
                                                     label = "Enter a name for your gene set (optional):",
                                                     placeholder = "Enter a name")
                                       )
                                   ),
                                   
                                   
                                   ### Select variable-of-interest
                                   div(style="display: inline-block;vertical-align:top; width: 100%;",
                                       selectInput("ssgsea_var_1", "Select Subtype:", 
                                                   choices = c("PDS", "CMS", "iCMS"), 
                                                   selected = "PDS", 
                                                   multiple = FALSE)),
                                   # Button user clicks to perform ssGSEA after selecting an
                                   # existing gene list/signature or entering their own gene
                                   # list
                                   br(),
                                   div(style="display: inline-block;vertical-align:top;",
                                       actionButton("ssgsea_submit", "Perform ssGSEA"))
                               )),
                        ### Boxplot Output
                        column(4,
                               box(width = NULL, status = "success", solidHeader = TRUE,
                                   title = "FOCUS",
                                   withSpinner(
                                       plotOutput("focus_ssgsea_boxplot"),
                                       type = 1, 
                                       size = 0.5
                                   )
                               )),
                        
                        ### Boxplot Output
                        column(4,
                               box(width = NULL, status = "success", solidHeader = TRUE,
                                   title = "GSE39582",
                                   withSpinner(
                                       plotOutput("gse39582_ssgsea_boxplot"),
                                       type = 1, 
                                       size = 0.5
                                   )
                               ))
                    ),
            )
            
        )
    )


## 4. put UI together ----
ui <-
    dashboardPage( 
        header, sidebar, body,
        ### aesthetics
        skin = "green"
    )

