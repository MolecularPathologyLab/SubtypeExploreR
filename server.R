
## Define server ----
server <- function(session, input, output) {
    
    ### FOCUS expression data ----
    focus_exp <- reactive({
        as.data.frame(focus_express)
    })
    
    ### FOCUS clinical/molecular data ----
    focus_clin <- reactive({
        as.data.frame(focus_clino)
    })
    
    ### stored GSE39582 expression data ----
    gse39582_exp <- reactive({
        as.data.frame(gse39582_express)
    })
    
    ### stored GSE39582 clinical/molecular data ----
    gse39582_clin <- reactive({
        as.data.frame(gse39582_clino)
    })

    
    #### BOXPLOTS TAB ####
    
    # Create a vector with all the genes
    available_genes <- reactive({
        available_genes <- c(focus_exp()$Gene.Symbol, gse39582_exp()$Gene.Symbol)
        available_genes <- unique(available_genes)
        available_genes <- available_genes[!grepl("///", available_genes)]
        available_genes <- sort(available_genes)
        available_genes
    })
    
    # Update the choice of genes for the user
    observe({
        
        ### Select gene for expression boxplots
        updateSelectizeInput(session, 'boxplot_gene', choices = available_genes(),
                             options = list(maxOptions = length(available_genes())),
                             selected = "A1BG",
                             server = TRUE)
    })
    
    
    # Check if at least one gene symbol in the user's custom gene set is
    # present in the FOCUS expression data
    focus_boxplot_gene_present <- reactive({
        input$boxplot_gene %in% focus_exp()$Gene.Symbol
    })
    
    # Get expression values for the chosen gene
    focus_boxplot_gene_exp <- eventReactive(input$boxplot_submit, {
        
        # Ensure at gene symbol has been selected
        validate(need(input$boxplot_gene != "",
                 "Please select a gene"))
        
        # Ensure at least gene symbol is
        # present in the expression data and if not inform the user
        validate(
            need(focus_boxplot_gene_present(),
                 "The gene you have selected is not present in the gene expression data for this dataset"))
        
        exp_values <- as.numeric(focus_exp()[focus_exp()$Gene.Symbol == input$boxplot_gene, colnames(focus_exp())[colnames(focus_exp()) != "Gene.Symbol"]])
        exp_df <- data.frame("Sample_ID" = colnames(focus_exp())[colnames(focus_exp()) != "Gene.Symbol"],
                             "Gene_Exp" = exp_values)
        
        exp_df
        
    })
    
    focus_boxplot_groups <- eventReactive(input$boxplot_submit, {
        if (input$boxplot_var == "CMS"){
            boxplot_groups <- focus_clin()[,c("Sample_ID", "CMS_RF_0_4")]
            colnames(boxplot_groups)[2] <- "CMS"
            boxplot_groups
        }
        else if (input$boxplot_var == "PDS"){
            boxplot_groups <- focus_clin()[,c("Sample_ID", "PDS_call")]
            colnames(boxplot_groups)[2] <- "PDS"
            boxplot_groups
        }
        else if (input$boxplot_var == "iCMS"){
            boxplot_groups <- focus_clin()[,c("Sample_ID", "iCMS_final_prediction")]
            colnames(boxplot_groups)[2] <- "iCMS"
            boxplot_groups
        }
    })
    
    focus_gene_exp_boxplot_data <- eventReactive(input$boxplot_submit, {
        merge(focus_boxplot_gene_exp(), focus_boxplot_groups(), by = "Sample_ID")
    })
    
    focus_gene_exp_boxplot_obj <- eventReactive(input$boxplot_submit,{
        
        
        req(focus_gene_exp_boxplot_data())
        
        # Set colours to use
        if (input$boxplot_var == "CMS"){
            cols <- cms_cols
        }
        else if (input$boxplot_var == "PDS"){
            cols <- pds_cols
        }
        else if (input$boxplot_var == "iCMS"){
            cols <- icms_cols
        }
        
        # Set comparisons for statistical testing
        if (input$boxplot_var == "CMS"){
            stat_comparisons <- list(c("CMS2", "CMS1"), 
                                     c("CMS2", "CMS3"),
                                     c("CMS2", "CMS4"))
        }
        else if (input$boxplot_var == "PDS"){
            stat_comparisons <- list(c("PDS2", "PDS3"), 
                                     c("PDS1", "PDS3"),
                                     c("PDS1", "PDS2"))
        }
        else if (input$boxplot_var == "iCMS"){
            stat_comparisons <- list(c("iCMS2", "iCMS3"))
        }
        
        
        ggplot(focus_gene_exp_boxplot_data(),
               aes_string(x = input$boxplot_var, y = "Gene_Exp", 
                          col = input$boxplot_var, fill = input$boxplot_var)) +
            geom_violin(width=0.8, alpha = 0.3, colour = NA) +
            geom_beeswarm(alpha = 0.5, size = 3) +
            geom_boxplot(width=0.1, color="white", alpha=0.2) +
            scale_fill_manual(name=input$boxplot_var, values = cols, guide = 'none') +
            scale_colour_manual(name=input$boxplot_var, values = cols, guide = 'none') +
            stat_compare_means(comparisons = stat_comparisons,
                               method = 'wilcox.test',
                               ## Add pairwise comparisons p-value
                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001,
                                                                0.01, 0.05, 1),
                                                  symbols = c("****", "***", "**",
                                                              "*", "ns"))) +
            xlab(NULL) +
            ylab('Gene Expression') +
            ggtitle(str_trunc(paste(stri_wrap(input$boxplot_gene, width = 36, whitespace_only = TRUE), collapse = "\n"), width = 85, side = "right", ellipsis = "...")) +
            box_theme.1
        
    })
    
    output$focus_gene_exp_boxplot <- renderPlot({
        focus_gene_exp_boxplot_obj()
    })
    
    
    # Check if at least one gene symbol in the user's custom gene set is
    # present in the GSE39582 expression data
    gse39582_boxplot_gene_present <- reactive({
        input$boxplot_gene %in% gse39582_exp()$Gene.Symbol
    })
    
    # Get expression values for the chosen gene
    gse39582_boxplot_gene_exp <- eventReactive(input$boxplot_submit, {
        
        # Ensure at gene symbol has been selected
        validate(need(input$boxplot_gene != "",
                      "Please select a gene"))
        
        # Ensure at least gene symbol is
        # present in the expression data and if not inform the user
        validate(
            need(gse39582_boxplot_gene_present(),
                 "The gene you have selected is not present in the gene expression data for this dataset"))
        
        exp_values <- as.numeric(gse39582_exp()[gse39582_exp()$Gene.Symbol == input$boxplot_gene, colnames(gse39582_exp())[colnames(gse39582_exp()) != "Gene.Symbol"]])
        exp_df <- data.frame("Sample_ID" = colnames(gse39582_exp())[colnames(gse39582_exp()) != "Gene.Symbol"],
                             "Gene_Exp" = exp_values)
        exp_df
        
    })
    
    gse39582_boxplot_groups <- eventReactive(input$boxplot_submit, {
        if (input$boxplot_var == "CMS"){
            boxplot_groups <- gse39582_clin()[,c("Sample_ID", "CMS")]
            colnames(boxplot_groups)[2] <- "CMS"
            boxplot_groups
        }
        else if (input$boxplot_var == "PDS"){
            boxplot_groups <- gse39582_clin()[,c("Sample_ID", "PDS_call")]
            colnames(boxplot_groups)[2] <- "PDS"
            boxplot_groups
        }
        else if (input$boxplot_var == "iCMS"){
            boxplot_groups <- gse39582_clin()[,c("Sample_ID", "iCMS_final_prediction")]
            colnames(boxplot_groups)[2] <- "iCMS"
            boxplot_groups
        }
    })
    
    gse39582_gene_exp_boxplot_data <- eventReactive(input$boxplot_submit, {
        merge(gse39582_boxplot_gene_exp(), gse39582_boxplot_groups(), by = "Sample_ID")
    })
    
    gse39582_gene_exp_boxplot_obj <- eventReactive(input$boxplot_submit,{
        
        
        req(gse39582_gene_exp_boxplot_data())
        
        # Set colours to use
        if (input$boxplot_var == "CMS"){
            cols <- cms_cols
        }
        else if (input$boxplot_var == "PDS"){
            cols <- pds_cols
        }
        else if (input$boxplot_var == "iCMS"){
            cols <- icms_cols
        }
        
        # Set comparisons for statistical testing
        if (input$boxplot_var == "CMS"){
            stat_comparisons <- list(c("CMS2", "CMS1"), 
                                     c("CMS2", "CMS3"),
                                     c("CMS2", "CMS4"))
        }
        else if (input$boxplot_var == "PDS"){
            stat_comparisons <- list(c("PDS2", "PDS3"), 
                                     c("PDS1", "PDS3"),
                                     c("PDS1", "PDS2"))
        }
        else if (input$boxplot_var == "iCMS"){
            stat_comparisons <- list(c("iCMS2", "iCMS3"))
        }
        
        
        ggplot(gse39582_gene_exp_boxplot_data(),
               aes_string(x = input$boxplot_var, y = "Gene_Exp", 
                          col = input$boxplot_var, fill = input$boxplot_var)) +
            geom_violin(width=0.8, alpha = 0.3, colour = NA) +
            geom_beeswarm(alpha = 0.5, size = 3) +
            geom_boxplot(width=0.1, color="white", alpha=0.2) +
            scale_fill_manual(name=input$boxplot_var, values = cols, guide = 'none') +
            scale_colour_manual(name=input$boxplot_var, values = cols, guide = 'none') +
            stat_compare_means(comparisons = stat_comparisons,
                               method = 'wilcox.test',
                               ## Add pairwise comparisons p-value
                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001,
                                                                0.01, 0.05, 1),
                                                  symbols = c("****", "***", "**",
                                                              "*", "ns"))) +
            xlab(NULL) +
            ylab('Gene Expression') +
            ggtitle(str_trunc(paste(stri_wrap(input$boxplot_gene, width = 36, whitespace_only = TRUE), collapse = "\n"), width = 85, side = "right", ellipsis = "...")) +
            box_theme.1
        
    })
    
    output$gse39582_gene_exp_boxplot <- renderPlot({
        gse39582_gene_exp_boxplot_obj()
    })
    
    
    #### ssGSEA TAB ####
    
    # When user selects a gene set collection set the corresponding string
    # which appears at the start of that gene set collection's gene sets
    collection_identifier <- eventReactive(input$ssgsea_collection,{
        if (input$ssgsea_collection == "Hallmark"){
            "HALLMARK"
        }
        else if (input$ssgsea_collection == "BioCarta"){
            "BIOCARTA"
        }
        else if (input$ssgsea_collection == "KEGG"){
            "KEGG"
        }
        else if (input$ssgsea_collection == "Reactome"){
            "REACTOME"
        }
        else if (input$ssgsea_collection == "PID"){
            "PID"
        }
        else if (input$ssgsea_collection == "WikiPathways"){
            "WP"
        }
    })
    
    # When user selects an gene set collection create a vector with the name of
    # all the gene sets in that collection
    collection_gene_sets <- eventReactive(input$ssgsea_collection,{
        collection_gene_set_pattern <- paste0("^", collection_identifier(), "_")
        collection_inds <- grep(collection_gene_set_pattern, names(gene_sets_list))
        collection_gene_sets <- names(gene_sets_list)[collection_inds]
        collection_gene_sets
    })
    
    # Update the gene set selectizeInput to reflect the available gene sets
    # based on the chosen ontology
    observeEvent(input$ssgsea_collection, {
        collection_gene_set_pattern <- paste0("^", collection_identifier(), "_")
        collection_gene_sets_available <- gsub(collection_gene_set_pattern, "", collection_gene_sets())
        collection_gene_sets_available <- gsub("_", " ", collection_gene_sets_available)
        updateSelectizeInput(session = session,
                             "gene_set",
                             label = "Select Gene Set:",
                             choices = collection_gene_sets_available,
                             options = list(maxOptions = length(collection_gene_sets_available)),
                             server = TRUE)
    })
    
    ssgsea_custom_genes <- eventReactive(input$ssgsea_submit,{
        if (input$ssgsea_collection != "Custom"){
            NULL
        }
        else if (input$ssgsea_collection == "Custom"){
            # Split gene list input on new line character
            genes <- strsplit(input$ssgsea_gene_list, "\n")[[1]]
            # Keep only unique gene symbols
            genes <- unique(genes)
            # Remove any horizontal whitespace (space, tab) at start or end
            # i.e. before or after gene symbol
            genes <- trimws(genes, 
                            which = "both", # leading and trailing whitespace
                            whitespace = "[ \t]" # space and tab
            )
            # Remove any empty strings
            genes <- genes[genes != ""]
            # Convert gene symbols to upper case in case the user entered them
            # partially or completely in lower case 
            genes <- toupper(genes)
            
            # If no input then make genes an empty string
            if (length(genes) == 0){
                ""
            }
            # Else return the unique gene symbols (with any empty strings
            # removed and all in upper case
            else {
                genes
            }
        }
    })
    
    # Check if at least one gene symbol in the user's custom gene set is
    # present in the FOCUS expression data
    focus_custom_genes_present <- reactive({
        any(ssgsea_custom_genes() %in% focus_exp()$Gene.Symbol)
    })
    
    # Get ssGSEA scores for the chosen/custom gene set
    focus_ssgsea_scores <- eventReactive(input$ssgsea_submit,{
        if (input$ssgsea_collection == "Hallmark"){
            validate(need(input$gene_set != "", "Please select a gene set"))
            full_geneset_name <- gsub(" ", "_", input$gene_set)
            full_geneset_name <- paste0(collection_identifier(), "_", full_geneset_name)
            focus_clin()[,c("Sample_ID", full_geneset_name)]
        }
        else {
            if (input$ssgsea_collection != "Custom"){
                validate(need(input$gene_set != "", "Please select a gene set"))
                full_geneset_name <- gsub(" ", "_", input$gene_set)
                full_geneset_name <- paste0(collection_identifier(), "_", full_geneset_name)
                genes_list <- gene_sets_list[full_geneset_name]
            }
            else if (input$ssgsea_collection == "Custom"){
                # Ensure at least one gene symbol in the user's custom gene set is
                # present in the expression data and if not inform the
                # user (as ssGSEA can not be be performed if no gene symbols
                # match those in the expression data)
                validate(
                    need(focus_custom_genes_present(),
                         "None of the gene symbols in your custom gene set are present in the gene expression data for this dataset"))
                
                # Put the genes into a list
                genes_list <- list(ssgsea_custom_genes())
                names(genes_list) <- "Custom_gene_set"
            }
            
            # Make Gene Symbols the rownames of expression data 
            # rather than first column
            express <- focus_exp()
            rownames(express) <- express$Gene.Symbol
            express$Gene.Symbol <- NULL
            
            # Perform ssGSEA
            ssgsea_scores <- gsva(expr = as.matrix(express),
                                  gset.idx.list = genes_list,
                                  method = "ssgsea",
                                  kcdf = "Gaussian",
                                  min.sz = 1,
                                  max.sz = Inf,
                                  ssgsea.norm = TRUE,
                                  verbose = FALSE)
            
            # Transpose and format the ssGSEA result
            ssgsea_scores <- as.data.frame(t(ssgsea_scores))
            ssgsea_scores <- data.frame(Sample_ID = rownames(ssgsea_scores),
                                        ssgsea_scores)
            rownames(ssgsea_scores) <- NULL
            
            # Return the ssGSEA scores
            ssgsea_scores
        }
        
    })
    
    focus_ssgsea_groups <- eventReactive(input$ssgsea_submit, {
        if (input$ssgsea_var_1 == "CMS"){
            ssgsea_groups <- focus_clin()[,c("Sample_ID", "CMS_RF_0_4")]
            colnames(ssgsea_groups)[2] <- "CMS"
            ssgsea_groups
        }
        else if (input$ssgsea_var_1 == "PDS"){
            ssgsea_groups <- focus_clin()[,c("Sample_ID", "PDS_call")]
            colnames(ssgsea_groups)[2] <- "PDS"
            ssgsea_groups
        }
        else if (input$ssgsea_var_1 == "iCMS"){
            ssgsea_groups <- focus_clin()[,c("Sample_ID", "iCMS_final_prediction")]
            colnames(ssgsea_groups)[2] <- "iCMS"
            ssgsea_groups
        }
    })
    
    focus_ssgsea_boxplot_data <- eventReactive(input$ssgsea_submit, {
        merge(focus_ssgsea_scores(), focus_ssgsea_groups(), by = "Sample_ID")
    })
    
    focus_ssgsea_boxplot_obj <- eventReactive(input$ssgsea_submit,{
        
        req(focus_ssgsea_boxplot_data())
        
        if (input$ssgsea_collection != "Custom"){
            full_geneset_name <- gsub(" ", "_", input$gene_set)
            full_geneset_name <- paste0(collection_identifier(), "_", full_geneset_name)
            plot_title <- str_trunc(paste(stri_wrap(input$gene_set, width = 30, whitespace_only = TRUE), collapse = "\n"), width = 85, side = "right", ellipsis = "...")
        }
        else {
            full_geneset_name <- colnames(focus_ssgsea_boxplot_data())[2]
            if (is.null(input$custom_gene_set_name) | input$custom_gene_set_name == ""){
                plot_title <- NULL
            }
            else{
                plot_title <- str_trunc(paste(stri_wrap(input$custom_gene_set_name, width = 30, whitespace_only = TRUE), collapse = "\n"), width = 85, side = "right", ellipsis = "...")
            }
            
        }
        
        # Set colours to use
        if (input$ssgsea_var_1 == "CMS"){
            cols <- cms_cols
        }
        else if (input$ssgsea_var_1 == "PDS"){
            cols <- pds_cols
        }
        else if (input$ssgsea_var_1 == "iCMS"){
            cols <- icms_cols
        }
        
        
        # Set comparisons for statistical testing
        if (input$ssgsea_var_1 == "CMS"){
            stat_comparisons <- list(c("CMS2", "CMS1"), 
                                     c("CMS2", "CMS3"),
                                     c("CMS2", "CMS4"))
        }
        else if (input$ssgsea_var_1 == "PDS"){
            stat_comparisons <- list(c("PDS2", "PDS3"), 
                                     c("PDS1", "PDS3"),
                                     c("PDS1", "PDS2"))
        }
        else if (input$ssgsea_var_1 == "iCMS"){
            stat_comparisons <- list(c("iCMS2", "iCMS3"))
        }
        
        
        ggplot(focus_ssgsea_boxplot_data(),
               aes_string(x = input$ssgsea_var_1, y = full_geneset_name, 
                          col = input$ssgsea_var_1, fill = input$ssgsea_var_1)) +
            geom_violin(width=0.8, alpha = 0.3, colour = NA) +
            geom_beeswarm(alpha = 0.5, size = 3) +
            geom_boxplot(width=0.1, color="white", alpha=0.2) +
            scale_fill_manual(name=input$ssgsea_var_1, values = cols, guide = 'none') +
            scale_colour_manual(name=input$ssgsea_var_1, values = cols, guide = 'none') +
            stat_compare_means(comparisons = stat_comparisons,
                               method = 'wilcox.test',
                               ## Add pairwise comparisons p-value
                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001,
                                                                0.01, 0.05, 1),
                                                  symbols = c("****", "***", "**",
                                                              "*", "ns"))) +
            xlab(NULL) +
            ylab('ssGSEA score') +
            ggtitle(plot_title) +
            box_theme.1
        
    })
    
    output$focus_ssgsea_boxplot <- renderPlot({
        focus_ssgsea_boxplot_obj()
    })
    
    
    
    ###### GSE39582 #####
    
    # Check if at least one gene symbol in the user's custom gene set is
    # present in the GSE39582 expression data
    gse39582_custom_genes_present <- reactive({
        any(ssgsea_custom_genes() %in% gse39582_exp()$Gene.Symbol)
    })
    
    # Get ssGSEA scores for the chosen/custom gene set
    gse39582_ssgsea_scores <- eventReactive(input$ssgsea_submit,{
        if (input$ssgsea_collection == "Hallmark"){
            validate(need(input$gene_set != "", "Please select a gene set"))
            full_geneset_name <- gsub(" ", "_", input$gene_set)
            full_geneset_name <- paste0(collection_identifier(), "_", full_geneset_name)
            gse39582_clin()[,c("Sample_ID", full_geneset_name)]
        }
        else {
            if (input$ssgsea_collection != "Custom"){
                validate(need(input$gene_set != "", "Please select a gene set"))
                full_geneset_name <- gsub(" ", "_", input$gene_set)
                full_geneset_name <- paste0(collection_identifier(), "_", full_geneset_name)
                genes_list <- gene_sets_list[full_geneset_name]
            }
            else if (input$ssgsea_collection == "Custom"){
                # Ensure at least one gene symbol in the user's custom gene set is
                # present in the expression data and if not inform the
                # user (as ssGSEA can not be be performed if no gene symbols
                # match those in the expression data)
                validate(
                    need(gse39582_custom_genes_present(),
                         "None of the gene symbols in your custom gene set are present in the gene expression data for this dataset"))
                
                # Put the genes into a list
                genes_list <- list(ssgsea_custom_genes())
                names(genes_list) <- "Custom_gene_set"
            }
            
            # Make Gene Symbols the rownames of expression data 
            # rather than first column
            express <- gse39582_exp()
            rownames(express) <- express$Gene.Symbol
            express$Gene.Symbol <- NULL
            
            # Perform ssGSEA
            ssgsea_scores <- gsva(expr = as.matrix(express),
                                  gset.idx.list = genes_list,
                                  method = "ssgsea",
                                  kcdf = "Gaussian",
                                  min.sz = 1,
                                  max.sz = Inf,
                                  ssgsea.norm = TRUE,
                                  verbose = FALSE)
            
            # Transpose and format the ssGSEA result
            ssgsea_scores <- as.data.frame(t(ssgsea_scores))
            ssgsea_scores <- data.frame(Sample_ID = rownames(ssgsea_scores),
                                        ssgsea_scores)
            rownames(ssgsea_scores) <- NULL
            
            # Return the ssGSEA scores
            ssgsea_scores
        }
        
    })
    
    gse39582_ssgsea_groups <- eventReactive(input$ssgsea_submit, {
        if (input$ssgsea_var_1 == "CMS"){
            ssgsea_groups <- gse39582_clin()[,c("Sample_ID", "CMS")]
            colnames(ssgsea_groups)[2] <- "CMS"
            ssgsea_groups
        }
        else if (input$ssgsea_var_1 == "PDS"){
            ssgsea_groups <- gse39582_clin()[,c("Sample_ID", "PDS_call")]
            colnames(ssgsea_groups)[2] <- "PDS"
            ssgsea_groups
        }
        else if (input$ssgsea_var_1 == "iCMS"){
            ssgsea_groups <- gse39582_clin()[,c("Sample_ID", "iCMS_final_prediction")]
            colnames(ssgsea_groups)[2] <- "iCMS"
            ssgsea_groups
        }
    })
    
    gse39582_ssgsea_boxplot_data <- eventReactive(input$ssgsea_submit, {
        merge(gse39582_ssgsea_scores(), gse39582_ssgsea_groups(), by = "Sample_ID")
    })
    
    gse39582_ssgsea_boxplot_obj <- eventReactive(input$ssgsea_submit,{
        
        
        req(gse39582_ssgsea_boxplot_data())
        
        if (input$ssgsea_collection != "Custom"){
            full_geneset_name <- gsub(" ", "_", input$gene_set)
            full_geneset_name <- paste0(collection_identifier(), "_", full_geneset_name)
            plot_title <- str_trunc(paste(stri_wrap(input$gene_set, width = 30, whitespace_only = TRUE), collapse = "\n"), width = 85, side = "right", ellipsis = "...")
        }
        else {
            full_geneset_name <- colnames(gse39582_ssgsea_boxplot_data())[2]
            if (is.null(input$custom_gene_set_name) | input$custom_gene_set_name == ""){
                plot_title <- NULL
            }
            else{
                plot_title <- str_trunc(paste(stri_wrap(input$custom_gene_set_name, width = 30, whitespace_only = TRUE), collapse = "\n"), width = 85, side = "right", ellipsis = "...")
            }
            
        }
        
        # Set colours to use
        if (input$ssgsea_var_1 == "CMS"){
            cols <- cms_cols
        }
        else if (input$ssgsea_var_1 == "PDS"){
            cols <- pds_cols
        }
        else if (input$ssgsea_var_1 == "iCMS"){
            cols <- icms_cols
        }
        
        # Set comparisons for statistical testing
        if (input$ssgsea_var_1 == "CMS"){
            stat_comparisons <- list(c("CMS2", "CMS1"), 
                                     c("CMS2", "CMS3"),
                                     c("CMS2", "CMS4"))
        }
        else if (input$ssgsea_var_1 == "PDS"){
            stat_comparisons <- list(c("PDS2", "PDS3"), 
                                     c("PDS1", "PDS3"),
                                     c("PDS1", "PDS2"))
        }
        else if (input$ssgsea_var_1 == "iCMS"){
            stat_comparisons <- list(c("iCMS2", "iCMS3"))
        }
        
        
        ggplot(gse39582_ssgsea_boxplot_data(),
               aes_string(x = input$ssgsea_var_1, y = full_geneset_name, 
                          col = input$ssgsea_var_1, fill = input$ssgsea_var_1)) +
            geom_violin(width=0.8, alpha = 0.3, colour = NA) +
            geom_beeswarm(alpha = 0.5, size = 3) +
            geom_boxplot(width=0.1, color="white", alpha=0.2) +
            scale_fill_manual(name=input$ssgsea_var_1, values = cols, guide = 'none') +
            scale_colour_manual(name=input$ssgsea_var_1, values = cols, guide = 'none') +
            stat_compare_means(comparisons = stat_comparisons,
                               method = 'wilcox.test',
                               ## Add pairwise comparisons p-value
                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001,
                                                                0.01, 0.05, 1),
                                                  symbols = c("****", "***", "**",
                                                              "*", "ns"))) +
            xlab(NULL) +
            ylab('ssGSEA score') +
            ggtitle(plot_title) +
            box_theme.1
        
    })
    
    output$gse39582_ssgsea_boxplot <- renderPlot({
        gse39582_ssgsea_boxplot_obj()
    })
    
}
