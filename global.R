
## R packages
library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(GSVA)
library(stringi)
library(stringr)

## Load data ----
## FOCUS
focus_express <- fread(file.path('.', 'data', 'focus_gene_expression.txt'))
focus_clino <- fread(file.path('.','data', 'focus_clinical_data.txt'),
                                    stringsAsFactors = TRUE)


## GSE39582
gse39582_express <- fread(file.path('.', 'data', 'gse39582_gene_expression.txt'))
gse39582_clino <- fread(file.path('.', 'data', 'gse39582_clinical_data.txt'),
                                    stringsAsFactors = TRUE)


## Gene sets
gene_sets_list <- readRDS(file.path('.','data', 'hallmark_selected_c2_genesets_msigdbr_v7_4_1.rds'))


#### Relevel the PDS, CMS and iCMS calls ####

#### FOCUS #####

# Make PDS a factor variable with Mixed as last level
focus_clino$PDS_call <- factor(focus_clino$PDS_call,
                               levels = c(
                                   "PDS1",
                                   "PDS2",
                                   "PDS3",
                                   "Mixed"
                               ))

# Make CMS a factor variable with unknown (UNK) as last level
focus_clino$CMS_RF_0_4 <- factor(focus_clino$CMS_RF_0_4,
                                 levels = c(
                                     "CMS1",
                                     "CMS2",
                                     "CMS3",
                                     "CMS4",
                                     "UNK"
                                 ))

# Make CMS a factor variable with unknown (UNK) as last level
focus_clino$CMS_RF_0_5 <- factor(focus_clino$CMS_RF_0_5,
                                 levels = c(
                                     "CMS1",
                                     "CMS2",
                                     "CMS3",
                                     "CMS4",
                                     "UNK"
                                 ))

# Make iCMS a factor variable with unknown (UNK) as last level
focus_clino$iCMS_final_prediction <- 
    factor(focus_clino$iCMS_final_prediction,
           levels = c(
               "iCMS2",
               "iCMS3",
               "UNK"
           ))


#### GSE39582 #####

# Make PDS a factor variable with Mixed as last level
gse39582_clino$PDS_call <- factor(gse39582_clino$PDS_call,
                                  levels = c(
                                      "PDS1",
                                      "PDS2",
                                      "PDS3",
                                      "Mixed"
                                  ))

# Make CMS a factor variable with unknown (UNK) as last level
gse39582_clino$CMS <- factor(gse39582_clino$CMS,
                                    levels = c(
                                        "CMS1",
                                        "CMS2",
                                        "CMS3",
                                        "CMS4",
                                        "UNK"
                                    ))

# Make iCMS a factor variable with unknown (UNK) as last level
gse39582_clino$iCMS_final_prediction <- 
    factor(gse39582_clino$iCMS_final_prediction,
           levels = c(
               "iCMS2",
               "iCMS3",
               "UNK"
           ))


## Plot settings ----
# Define colours to use in plots for PDS subtypes
pds_cols <- c('PDS1' = '#920000',
              'PDS2' = '#332288',
              'PDS3' = '#D55E00',
              'Mixed' = 'grey50')

# Define colours to use in plots for CMS subtypes
cms_cols <- c("CMS1" = "#FFA54D",
              "CMS2" = "#0070B0",
              "CMS3" = "#CF75A8",
              "CMS4" = "#009C75",
              "UNK" = "grey50")

# Define colours to use in plots for iCMS subtypes
icms_cols <- c("iCMS2" = "#6A3D9A",
               "iCMS3" = "#FF8000",
               "UNK" = "grey50")

### basic theme - boxplot 
box_theme.1 <- theme(axis.title.x = element_blank(),
                     axis.text.x = element_text(size = 18, colour = "black"),
                     axis.title.y = element_text(size = 20),
                     axis.text.y = element_text(size = 18, colour = "black"),
                     axis.ticks.y = element_line(linewidth = 0.5, colour = "black"),
                     axis.ticks.length = unit(0.2, "cm"),
                     axis.line.y = element_line(linewidth = 0.5, colour = "black"),
                     axis.ticks.x = element_line(linewidth = 0.5, colour = "black"),
                     axis.line.x.bottom = element_line(linewidth = 0.5, colour = "black"),
                     panel.background = element_rect(fill = 'white'),
                     plot.margin = unit(c(4,2,0,2),'mm'),
                     plot.title = element_text(face = 'bold', size = 18,
                                               margin = margin(15,0,0,0),
                                               vjust = 4, hjust = 0.5, lineheight = 1),
                     legend.position = "none"
                     )
