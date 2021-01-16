library(shiny)
library(Seurat)
library(shinythemes)
library(shinycustomloader)
library(shinyjs)
library(igraph)
library(shinybusy)
library(XVector)
library(shinyWidgets)
library(shinyalert)
library(shinyFiles)
library(memisc)
library(AUCell)
library(dorothea)
library(networkD3)

library(bcellViper)
library(dplyr)
library(viper)
#library(V8)
jscode <- "shinyjs.toTop = function() {window.scrollTo(0, 0);}" ### for to top button


shinyUI(fluidPage(
  theme=shinytheme("cyborg"),
  tags$style(HTML("
        .tabs-above > .nav > li[class=active] > a {
           background-color: #000;
           color: #FFFF;
        }")),
  
  #setBackgroundImage(
  #  src = "download.jpg"
 # ),
 #themeSelector(),
 useShinyalert(),
  
  navbarPage(em("EcTracker"),
                   tabPanel("Home" , 
                            icon = icon("home"),
                          
                           HTML('<center><img src="aeae.gif", width="58%", height= "58%" ></center>'),
                         
                           
                          
                          column(8, align="center", offset = 2 ,print("Developed by: "), tags$u("Vishakha Gautam"), print("&" ), tags$u("Siddhant Kalra"),
                          
                          br(),print("For more details:"),
                          tags$head(tags$style(HTML("a {color: blue}"))),
                          tags$a(href="https://ahuja-lab.in/", "The Ahuja Lab", style = "color:yellow"))
                            
                          #tags$img(src='RStudio-Ball_2.png', heigth=400, width=400,)
                            
                            ),
                   
             
             tabPanel("Analysis" , icon = icon("angle-double-right"),
                     
                      #h5("Visualize the Feature Plot"),
                     navlistPanel( 
                     tabPanel(tags$b("1. Upload") , 
                               
                                  mainPanel(fileInput("file","Upload Expression Matrix"),
                                            
                                            h5("OR"),br(),
                                            fileInput("File10","Upload 10X Data Files (barcode.tsv, genes.tsv & matrix.mtx) ", multiple = TRUE),
                                            
                                            checkboxInput("testme", " Transdifferentiation (Sample Data)", 
                                                          value = FALSE),
                                            # checkboxInput("testme2","Sample Dataset 2", value = FALSE),
                                            
                                            br(),br(),
                                           # tags$b( "Select Parameters"),
                                            #checkboxInput(inputId = 'header', label = 'Header', value = FALSE),
                                           # checkboxInput(inputId = "stringAsFactors", "stringAsFactors", FALSE),
                                            #radioButtons(inputId = "rowname", label = "Row Names", choices=c(Yes=1, No=NULL), 1),
                                            #radioButtons(inputId = 'sep', label = 'Separator', choices = c(Comma=',',Semicolon=';',Tab='\t', Space=''), selected = ','),
                                            #actionButton("action000","sample")
                                           h6("Information about Data"),
                                           withLoader(tableOutput("dat"),type="image",loader = "ectracker_popup.gif"),
                                           h6("Dimensions"),
                                           withLoader(verbatimTextOutput("dimen"),type="image",loader = "ectracker_popup.gif"),
                                            
                                            
                                  )
                                  
                                ),
                       tabPanel(tags$b("2. scData Processing"),
                                
                                
                                actionButton("action1","One-Click Analysis*",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("paper-plane")),
                                br(),br(),
                               # print("Results are in different tabs")  ,
                                h6("Object"),
                               ## changes the color of action button
                             withLoader(verbatimTextOutput("seu_object"),type="image",loader = "ectracker_popup.gif"),br(),
                                actionButton("EX1","Explanation",icon = icon("book-open")),
                                withLoader(verbatimTextOutput("ex_object"),type="image",loader = "ectracker_popup.gif"),
                       
                       br(),br(),
                    
                                h6("Variable Feature Plot"),
                                
                               withLoader( plotOutput("p1"),type="image",loader = "ectracker_popup.gif"),
                                
                                br(), downloadButton(outputId="dndPlot","Download Plot"),    
                     actionButton("EX2","Explanation", icon = icon("book-open")),
                     
                     withLoader(verbatimTextOutput("ex_variable"),type="image",loader = "ectracker_popup.gif"),
                     
                     
                      #downloadButton("p1_down",label="Download as PDF"),
                      #submitButton("Next"),
                      #h5("Normalize Data"),
                    # actionButton("action2","Normalize Data"),
                     # withLoader(verbatimTextOutput("data_norm"),type="image",loader = "ectracker_popup.gif"),
                     # actionButton("EX3","Explanation",icon = icon("book-open")),
                      #withLoader(verbatimTextOutput("ex_normalization"),type="image",loader = "ectracker_popup.gif"),
                    
                       
                    
                      #verbatimTextOutput("variable_features"),
                    
                    br(),br(),
                      h6("Variable Features"),
                      withLoader(plotOutput("plot_vf"), type="image", loader = "ectracker_popup.gif"),
                     br(), downloadButton(outputId="vfPlot","Download Plot"),
                     actionButton("EX4","Explanation",icon = icon("book-open")),
                     
                     withLoader(verbatimTextOutput("ex_vfeature"),type="image",loader = "ectracker_popup.gif"),
                    
                     
                      
                    
                     #actionButton("action4","Scale Data"),
                    #  verbatimTextOutput("scale_data"),
                     # actionButton("EX5","Explanation", icon = icon("book-open")),
                      
                     # withLoader(verbatimTextOutput("ex_scale"),type="image",loader = "ectracker_popup.gif"),
                      
                      
                      br(),br(),
                      h6("PCA"),
                      withLoader(plotOutput("pca"),type="image",loader = "ectracker_popup.gif"),
                      br(),downloadButton(outputId="pcaPlot","Download Plot"),
                      actionButton("EX6","Explanation", icon = icon("book-open")),
                       withLoader(verbatimTextOutput("ex_pca"),type="image",loader = "ectracker_popup.gif"),
                    
                      
                     br(),br(),
                    
                      h6("Heat Map"),
                      withLoader(plotOutput("heatmap"),type="image",loader = "ectracker_popup.gif"),
                      br(),
                      
                      downloadButton(outputId="heatmapPlot","Download Plot"), actionButton("EX7","Explanation", icon = icon("book-open")),
                      
                      withLoader(verbatimTextOutput("ex_heatmap"),type="image",loader = "ectracker_popup.gif"),
                      
                      
                     br(), 
                     br(),
                     h6("JackStraw Plot"),
                      
                     withLoader(plotOutput("jackstraw"),type="image",loader = "ectracker_popup.gif"),
                      br(),
                    
                     downloadButton(outputId="jackstrawPlot","Download Plot"),actionButton("EX8","Explanation", icon = icon("book-open")),
                     
                     withLoader(verbatimTextOutput("ex_jackstraw"),type="image",loader = "ectracker_popup.gif"),
                      
                     br(),  
                     br(),
                     h6("Elbow Plot"),
                      withLoader(plotOutput("elbow"),type="image",loader = "ectracker_popup.gif"),
                     br(),
                    downloadButton(outputId="elbowPlot","Download Plot"),actionButton("EX9","Explanation", icon = icon("book-open")),
                    
                    withLoader(verbatimTextOutput("ex_elbow"),type="image",loader = "ectracker_popup.gif"),
                    br(),br(),

                    
                    
                      #actionButton("action9","Neighbours"),
                     # withLoader(verbatimTextOutput("neigh"),type="image",loader = "ectracker_popup.gif"),
                      
                      h6("UMAP"),
                      withLoader(plotOutput("umap"),type="image",loader = "ectracker_popup.gif"),
                      br(),
                      downloadButton(outputId="umapPlot","Download Plot"),
                      actionButton("EX10","Explanation", icon = icon("book-open")),
                      
                      withLoader(verbatimTextOutput("ex_umap"),type="image",loader = "ectracker_popup.gif"),
                    actionButton("toTop", "Top", icon = icon("arrow-alt-circle-up")),
                    
                    ),
                      
                    
                    tabPanel(tags$b("3. DEG Analysis"),
                             #textInput("log_fc","ENTER THRESHOLD VALUE FOR logFC (1 IS RECOMMENDED)"),
                             sliderInput("log_fc", "Select Fold Change Threshold (log)",
                                         min = 0, max = 20,
                                         value = 1),     
                      #actionButton("action11","Click here for Marker Identification",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),br(),
                    # withLoader(verbatimTextOutput("markers"),type="image",loader = "ectracker_popup.gif"),
                     downloadLink('downloadData', 'Download'),
                      br(),br(),
                     h6("Markers Table"),
                    withLoader(tableOutput("markers_display"),type="image",loader = "ectracker_popup.gif"),
                     useShinyjs(),
                     br(), 
                    fileInput("meta","Upload Meta File"),
                    checkboxInput("testme3", " Transdifferentiation (Meta file)", 
                                  value = FALSE),br(),
                    actionButton("m1","Table with Meta File "), 
                    withLoader(verbatimTextOutput("merge_meta"),type="image",loader = "ectracker_popup.gif"),
                    downloadLink('downloadData_meta', 'Download'),
                     br(),
                   br(),
                    h6("Generate Cluster Table"),
                      withLoader(verbatimTextOutput("cluster_ident"),type="image",loader = "ectracker_popup.gif"),
                     downloadLink('downloadData_cluster_ident', 'Download'),br(),br(),
                     extendShinyjs(text = jscode,functions = c("toTop")),
                    actionButton("toTop1", "Top", icon = icon("arrow-alt-circle-up")),
                    
                     
                     
                    ),
                   tabPanel(tags$b("4. CellEnrich "),
                            
                            tabsetPanel(type = "tabs",
                                        tabPanel(" AUCell (Method I)",
                                                 br(),br(),br(),
                                                # actionButton("au_main","One-Click Analysis",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                 
                                                 #actionButton("cell_geneset","Create Genesets"),
                                                 
                                                 # withLoader(verbatimTextOutput("geneset_making"),type="image",loader = "ectracker_popup.gif"),
                                                 
                                                 br(), # selectInput("aucell","Select Signature",choices = c( "Adrenal Gland","APC","Monocyte","Astrocyte","AT2_cell",	"B_cell",	"B_cell_plasmocyte",	"Basal_cell",	"CB_CD34+",	"Chondrocyte",	"Dendritic_cell",	"Endothelial_cell_APC",	"Endothelial_cell_to_mesenchymal",	"Endothelial_cell",	"Enterocyte", 	"Enterocyte_progenitor",	"Epithelial_cell_intermediated",	"Epithelial_cell",	"Erythroid_cell",	 "Erythroid_progenitor_cell",	"Fasciculata_cell",	"Fetal_acinar_cell",	"Fetal_chondrocyte",	"Fetal_endocrine_cell",	"Fetal_enterocyte",	"Fetal_epithelial_progenitor",	"Fetal_fibroblast",	"Fetal_skeletal_muscle","Fetal_mesenchymal_progenitor","Fetal_neuron", "Fetal_stromal_cell",	"Fibroblast",	"Gastric_chief_cell",	"Gastric_endocrine_cell",	"Goblet_cell",	"Hepatocyte_endodermal_cell",	"hESC87",	"Immature_sertoli_cell",	"Intercalated_cell",	"Intermediated_cell","Kidney_intercalated_cell",	"Loop_of_Henle",	"M2_Macrophage",	"Neutrophil",	"Macrophage",	"Mast_cell",	"Myeloid_cell","Mesothelial",	"Neutrophil (RPS high)",	"Oligodendrocyte",	"Pancreas exocrine cell",	"Primordial germ cell","Proximal_tubule_progenitor",	"Proliferating T cell",	"Sinusoidal_endothelial_cell",	"stromal_cell",	"smooth_muscle","Stratified_epithelial_cell",	"T_cell",	"Thyroid follicular cell","Ventricle cardiomyocyte","Ureteric bud cell")),
                                                 #selectInput("aucell","Select Signature",choices=c("Chemosensory - Signature","Astrocyte68","Antigen presenting cell (RPS high)41","Adrenal gland inflammatory cell102","AT2 cell30","b_cell","b_cell_plasmocyte","Basal_cell","CB CD34+23","Chondrocyte99","Dendritic_cell ","Endothelial_cell  ","Endothelial cell (endothelial to mesenchymal transition)66","Endothelial cell (APC)8","Enterocyte_progenitor","Enterocyte","Epithelial cell (intermediated)60","Epithelial_cell","Erythroid_cell ","Erythroid progenitor cell (RP high)12","Fasciculata_cell ","fetal B cells","fetal Basophil_Mast","fetal CD1C+ DCs","fetal CEPs","fetal CLEC9A+ DCs","fetal Collecting duct lineage","fetal Connecting tubule lineage","fetal Distal tubule lineage","fetal EBMPs","fetal EEPs","fetal Enteroendocrine cells_intestine_pancreas","fetal Enteroendocrine cells_stomach","fetal Erythroblast","fetal ETDs","fetal HSCs","fetal HSPCs","fetal IL1B+ Microglia","fetal ILC 3","fetal Islet beta cells","fetal Islet delta cells","fetal LEC","fetal Loop of Henle lineage","fetal Macrophages","fetal Meg progenitors","fetal Meg","fetal Megakaryoblasts","fetal Microglia","fetal Nephron progenitor cell","fetal NK cells","fetal pDCs","fetal Perivascular macrophages","fetal Phagocytic macrophages","fetal Plasma cells","fetal Podocyte lineage","fetal Proximal tubule lineage","fetal PTPRC+ Microglia","fetal Pulmonary neuroendocrine cells","fetal Renal vesicle","fetal S100A9+ DCs","fetal T cells","fetal TMEM119+ Microglia","fetal TRAF1+ APCs","fetal Ureteric tip","fetal Ureteric trunk","fetalAntigen-presenting macrophages","fetalEndocardium","VEC-adrenal","VEC-brain","VEC-kidney","VEC-liver","VEC-lung","VEC-placenta","VEC-spleen","Fetal acinar cell71","Fetal chondrocyte43","Fetal endocrine cell88","Fetal enterocyte 15","Fetal epithelial progenitor1","Fetal fibroblast17","Fetal skeletal muscle cell64","Fetal_mesenchymal_progenitor","Fetal_neuron ","Fetal_stromal_cell","Fibroblast","Gastric chief cell50","Gastric endocrine cell84","Goblet_cell ","Hepatocyte_Endodermal cell16","hESC87","Immature sertoli cell (Pre-Sertoli cell)92","Intercalated cell56","Intermediated cell97","Kidney intercalated cell101","Loop of Henle25","M2 Macrophage47","Neutrophil (RPS high)40","Neutrophil","Macrophage","Mast cell86","Myeloid cell93","Oligodendrocyte28","Pancreas exocrine cell44","Primordial germ cell62","Proximal tubule progenitor83","Proliferating T cell52","Sinusoidal endothelial cell61","Stromal_cell","Smooth_muscle_cell","Stratified_epithelial_cell","T_cell","Thyroid follicular cell32","Ventricle cardiomyocyte74","Ureteric bud cell77")),
                                                 # selectInput("aucell","Select Signature",choices=c("Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult b_cell","Adult b_cell_plasmocyte","Adult basal_cell","Adult CB CD34+23","Adult chondrocyte","Adult chemosensory ","Adult dendritic_cell ","Adult endothelial_cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte_progenitor","Adult enterocyte","Adult epithelial cell (intermediated)","Adult epithelial_cell","Adult erythroid_cell ","Adult erythroid progenitor cell (RP high)","Adult fasciculata_cell ","Fetal B cells","Fetal Basophil_Mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal Collecting duct lineage","Fetal Connecting tubule lineage","Fetal Distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal Enteroendocrine cells_intestine_pancreas","Fetal Enteroendocrine cells_stomach","Fetal Erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ Microglia","Fetal ILC 3","Fetal Islet beta cells","Fetal Islet delta cells","Fetal LEC","Fetal Loop of Henle lineage","Fetal Macrophages","Fetal Meg progenitors","Fetal Meg","Fetal Megakaryoblasts","Fetal Microglia","Fetal Nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal Perivascular macrophages","Fetal Phagocytic macrophages","Fetal Plasma cells","Fetal Podocyte lineage","Fetal Proximal tubule lineage","Fetal PTPRC+ Microglia","Fetal Pulmonary neuroendocrine cells","Fetal Renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ Microglia","Fetal TRAF1+ APCs","Fetal Ureteric tip","Fetal Ureteric trunk","Fetal Antigen-presenting macrophages","Fetal Endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal fibroblast","Fetal skeletal muscle cell","Fetal_mesenchymal_progenitor","Fetal_neuron ","Fetal_stromal_cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet_cell ","Adult hepatocyte_Endodermal cell","Adult hESC87","Adult immature sertoli cell (Pre-Sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 Macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal_cell","Adult smooth_muscle_cell","Adult stratified_epithelial_cell","Adult T_cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")),
                                                 # selectInput("aucell","Select Signature",choices=c("Devloping heart atrial cardiomyocyte","Devloping heart ventricular cardiomyocyte","Ectodermal cell differentiation","Ectodermal development","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Myofibroblast differentiation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult b_cell","Adult b_cell_plasmocyte","Adult basal_cell","Adult CB CD34+23","Adult chondrocyte","Adult chemosensory ","Adult dendritic_cell ","Adult endothelial_cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte_progenitor","Adult enterocyte","Adult epithelial cell (intermediated)","Adult epithelial_cell","Adult erythroid_cell ","Adult erythroid progenitor cell (RP high)","Adult fasciculata_cell ","Fetal B cells","Fetal Basophil_Mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal Collecting duct lineage","Fetal Connecting tubule lineage","Fetal Distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal Enteroendocrine cells_intestine_pancreas","Fetal Enteroendocrine cells_stomach","Fetal Erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ Microglia","Fetal ILC 3","Fetal Islet beta cells","Fetal Islet delta cells","Fetal LEC","Fetal Loop of Henle lineage","Fetal Macrophages","Fetal Meg progenitors","Fetal Meg","Fetal Megakaryoblasts","Fetal Microglia","Fetal Nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal Perivascular macrophages","Fetal Phagocytic macrophages","Fetal Plasma cells","Fetal Podocyte lineage","Fetal Proximal tubule lineage","Fetal PTPRC+ Microglia","Fetal Pulmonary neuroendocrine cells","Fetal Renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ Microglia","Fetal TRAF1+ APCs","Fetal Ureteric tip","Fetal Ureteric trunk","Fetal Antigen-presenting macrophages","Fetal Endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal fibroblast","Fetal skeletal muscle cell","Fetal_mesenchymal_progenitor","Fetal_neuron ","Fetal_stromal_cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet_cell ","Adult hepatocyte_Endodermal cell","Adult hESC87","Adult immature sertoli cell (Pre-Sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 Macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal_cell","Adult smooth_muscle_cell","Adult stratified_epithelial_cell","Adult T_cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")),
                                                 selectInput("aucell","Select Signature",choices=c("Ectodermal cell differentiation","Ectodermal development","Embryonic_stem_cell","Adult b_cell","Adult b_cell_plasmocyte","Fetal fibroblast","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult basal_cell","Adult CB CD34+23","Devloping heart atrial cardiomyocyte","Devloping heart ventricular cardiomyocyte","Adult chondrocyte","Adult chemosensory ","Adult dendritic_cell ","Adult endothelial_cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte_progenitor","Adult enterocyte","Adult epithelial cell (intermediated)","Adult epithelial_cell","Adult erythroid_cell ","Adult erythroid progenitor cell (RP high)","Adult fasciculata_cell ","Fetal B cells","Fetal Basophil_Mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal Collecting duct lineage","Fetal Connecting tubule lineage","Fetal Distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal Enteroendocrine cells_intestine_pancreas","Fetal Enteroendocrine cells_stomach","Fetal Erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ Microglia","Fetal ILC 3","Fetal Islet beta cells","Fetal Islet delta cells","Fetal LEC","Fetal Loop of Henle lineage","Fetal Macrophages","Fetal Meg progenitors","Fetal Meg","Fetal Megakaryoblasts","Fetal Microglia","Fetal Nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal Perivascular macrophages","Fetal Phagocytic macrophages","Fetal Plasma cells","Fetal Podocyte lineage","Fetal Proximal tubule lineage","Fetal PTPRC+ Microglia","Fetal Pulmonary neuroendocrine cells","Fetal Renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ Microglia","Fetal TRAF1+ APCs","Fetal Ureteric tip","Fetal Ureteric trunk","Fetal Antigen-presenting macrophages","Fetal Endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal skeletal muscle cell","Fetal_mesenchymal_progenitor","Fetal_neuron ","Fetal_stromal_cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet_cell ","Adult hepatocyte_Endodermal cell","Adult hESC87","Adult immature sertoli cell (Pre-Sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 Macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal_cell","Adult smooth_muscle_cell","Adult stratified_epithelial_cell","Adult T_cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")),
                                                 
                                                 br(), h6("Cell Ranking"), 
                                                 withLoader(plotOutput("cell_rank"),type="image",loader = "ectracker_popup.gif"),
                                                 br(), actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_cellrank"),type="image",loader = "ectracker_popup.gif"),
                                                 
                                                 # br(), br(),# actionButton("cell_Auc","Cell AUC"),
                                                 # withLoader(verbatimTextOutput("cell_auc"),type="image",loader = "ectracker_popup.gif"),
                                                 # br(), actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 # withLoader(verbatimTextOutput("ex_AUC"),type="image",loader = "ectracker_popup.gif"),
                                                 
                                                 br(), br(),   h6(" Cell Assignment"),
                                                 withLoader(plotOutput("au_cell"),type="image",loader = "ectracker_popup.gif"),br(),
                                                 actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_assign"),type="image",loader = "ectracker_popup.gif"),
                                                 
                                                 
                                                 # actionButton("au1","Cell Assignment Score"),
                                                 # withLoader(verbatimTextOutput("au_score"),type="image",loader = "ectracker_popup.gif"),
                                                 # actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 
                                                 # br(),br(), actionButton("au4","Heat Map"),
                                                 #  withLoader(plotOutput("hmap"),type="image",loader = "ectracker_popup.gif"),
                                                 # br(), actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 #  withLoader(verbatimTextOutput("ex_heat"),type="image",loader = "ectracker_popup.gif"),
                                                 
                                                 
                                                 br(), br(),  h6("UMAP"),
                                                 withLoader(plotOutput("tsne"),type="image",loader = "ectracker_popup.gif"),
                                                 br(),  actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_uma"),type="image",loader = "ectracker_popup.gif"),
                                                 
                                              br(),   br(),  h6("Plot for Signatures"),
                                                 withLoader(plotOutput("c_tsne"),type="image",loader = "ectracker_popup.gif"),
                                                 br(),# downloadButton(outputId="cell_umap","Download Plot"),
                                                  actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_umappp"),type="image",loader = "ectracker_popup.gif"),
                                                downloadButton(outputId="cell_assign","Download All Plots"),                                                 
                                                downloadLink(outputId="cell_assign_score","Download Cell Assignment Table"),
                                                
                                                
                                                 #  br(),br(),   actionButton("ngene","Detected Number of Genes"),
                                                 #      withLoader(plotOutput("n_gene"),type="image",loader = "ectracker_popup.gif"),
                                                 #      br(), downloadButton(outputId="cell_geneno","Download Plot"),
                                                br(), br(),  actionButton("toTop3", "Top", icon = icon("arrow-alt-circle-up")),
                                                 
                                                 
                                        ),
                                        
                                        tabPanel(" Stouffer Score (Method II)",
                                                # actionButton("stou_main","One-Click Analysis",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                 # selectInput("st_sc","Select Signature",choices=c("Ectodermal cell differentiation","Ectodermal development","Adult b_cell","Adult b_cell_plasmocyte","Fetal fibroblast","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult basal_cell","Adult CB CD34+23","Devloping heart atrial cardiomyocyte","Devloping heart ventricular cardiomyocyte","Adult chondrocyte","Adult chemosensory ","Adult dendritic_cell ","Adult endothelial_cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte_progenitor","Adult enterocyte","Adult epithelial cell (intermediated)","Adult epithelial_cell","Adult erythroid_cell ","Adult erythroid progenitor cell (RP high)","Adult fasciculata_cell ","Fetal B cells","Fetal Basophil_Mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal Collecting duct lineage","Fetal Connecting tubule lineage","Fetal Distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal Enteroendocrine cells_intestine_pancreas","Fetal Enteroendocrine cells_stomach","Fetal Erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ Microglia","Fetal ILC 3","Fetal Islet beta cells","Fetal Islet delta cells","Fetal LEC","Fetal Loop of Henle lineage","Fetal Macrophages","Fetal Meg progenitors","Fetal Meg","Fetal Megakaryoblasts","Fetal Microglia","Fetal Nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal Perivascular macrophages","Fetal Phagocytic macrophages","Fetal Plasma cells","Fetal Podocyte lineage","Fetal Proximal tubule lineage","Fetal PTPRC+ Microglia","Fetal Pulmonary neuroendocrine cells","Fetal Renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ Microglia","Fetal TRAF1+ APCs","Fetal Ureteric tip","Fetal Ureteric trunk","Fetal Antigen-presenting macrophages","Fetal Endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal skeletal muscle cell","Fetal_mesenchymal_progenitor","Fetal_neuron ","Fetal_stromal_cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet_cell ","Adult hepatocyte_Endodermal cell","Adult hESC87","Adult immature sertoli cell (Pre-Sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 Macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal_cell","Adult smooth_muscle_cell","Adult stratified_epithelial_cell","Adult T_cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")),
                                                 br(),br(),br(),br(),
                                                 selectInput("st_sc","Select Signature",choices=c("Devloping heart atrial cardiomyocyte","Devloping heart ventricular cardiomyocyte","Embryonic_stem_cell","Ectodermal cell differentiation","Ectodermal development","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Myofibroblast differentiation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult b_cell","Adult b_cell_plasmocyte","Adult basal_cell","Adult CB CD34+23","Adult chondrocyte","Adult chemosensory ","Adult dendritic_cell ","Adult endothelial_cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte_progenitor","Adult enterocyte","Adult epithelial cell (intermediated)","Adult epithelial_cell","Adult erythroid_cell ","Adult erythroid progenitor cell (RP high)","Adult fasciculata_cell ","Fetal B cells","Fetal Basophil_Mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal Collecting duct lineage","Fetal Connecting tubule lineage","Fetal Distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal Enteroendocrine cells_intestine_pancreas","Fetal Enteroendocrine cells_stomach","Fetal Erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ Microglia","Fetal ILC 3","Fetal Islet beta cells","Fetal Islet delta cells","Fetal LEC","Fetal Loop of Henle lineage","Fetal Macrophages","Fetal Meg progenitors","Fetal Meg","Fetal Megakaryoblasts","Fetal Microglia","Fetal Nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal Perivascular macrophages","Fetal Phagocytic macrophages","Fetal Plasma cells","Fetal Podocyte lineage","Fetal Proximal tubule lineage","Fetal PTPRC+ Microglia","Fetal Pulmonary neuroendocrine cells","Fetal Renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ Microglia","Fetal TRAF1+ APCs","Fetal Ureteric tip","Fetal Ureteric trunk","Fetal Antigen-presenting macrophages","Fetal Endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal fibroblast","Fetal skeletal muscle cell","Fetal_mesenchymal_progenitor","Fetal_neuron ","Fetal_stromal_cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet_cell ","Adult hepatocyte_Endodermal cell","Adult hESC87","Adult immature sertoli cell (Pre-Sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 Macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal_cell","Adult smooth_muscle_cell","Adult stratified_epithelial_cell","Adult T_cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")),
                                                 br(),
                                                 
                                                 h6("Stouffer Plot"),
                                                 withLoader(plotOutput("stouffer_sc_plot"),type="image",loader = "ectracker_popup.gif"),
                                                 br(), downloadButton(outputId="cell_stoudownload","Download Plot"),
                                                 
                                                 actionButton("ex12","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_stouffer"),type="image",loader = "ectracker_popup.gif"),
                                                 
                                                 br(),   br(),h6("One Sided Wilcoxon P-Value"),
                                                 withLoader(uiOutput("Stouffer_sc"),type="image",loader = "ectracker_popup.gif"), br(),
                                                downloadLink('stouffer_table_cell', 'Download Table'),br(),br(),
                                                 
                                                 
                                                 
                                                 
                                        ) 
                                        
                                        
                                        
                                        
                            )),
                   
                   tabPanel( tags$b("5. TissueEnrich"),
                             tabsetPanel(type = "tabs",
                                         
                                         
                                         
                                         tabPanel("Hypergeometric method (Method I)",
                                                  br(),
                                                  sliderInput("obs", "Select Cluster Number ",
                                                              min = 0, max = 20,
                                                              value = 0),
                                                  #numericInput("obs", "ENTER CLUSTER NUMBER", 0, min = 0, max = 1000),
                                                  
                                                  # actionButton("action24","Generate Specific Cluster Table"),br(),
                                                  
                                                  
                                                  # withLoader(verbatimTextOutput("cluster_name"),type="image",loader = "ectracker_popup.gif"),
                                                  downloadLink('downloadData_cluster_name', 'Download Cluster Table'),br(),br(),
                                                 # actionButton("tissue_main","One-Click Analysis",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                  
                                                  # br(),br(), actionButton("resetcluster", "Reset Cluster Number"),
                                                  
                                                  br(),br(),   numericInput("dataset", "Enter Dataset Number: 1: HPA ; 2: GTEx", 1, min = 1, max = 2), br(),
                                                  
                                                  h6("Generate Tissue Enrich Table"),
                                                  
                                                  withLoader(verbatimTextOutput("tissue_detail"),type="image",loader = "ectracker_popup.gif"),
                                                  
                                                  downloadLink('downloadtissuedata ', 'Download'),
                                                  
                                                  br(),  
                                                  br(),
                                                  
                                                 h6("Fold Change Plot"),
                                                  
                                                  withLoader(plotOutput("plot_fold"),type="image",loader = "ectracker_popup.gif"), 
                                                  br(),
                                                  downloadButton(outputId="tissue_fold1","Download Plot"),
                                                  
                                                  
                                                  br(),br(),
                                                  
                                                  
                                                  h6(" Log10(PValue)Plot"),
                                                  
                                                  withLoader(plotOutput("plot_log"),type="image",loader = "ectracker_popup.gif"),
                                                  br(),   downloadButton(outputId="tissue_log1","Download Plot"),
                                                  
                                                  
                                                  br(),br(),
                                                  
                                                  
                                                  h6("Tissue Specific Visualization"),
                                                  br(), uiOutput("t_or"),
                                                  uiOutput("tt_or"),br(),br(),
                                                  br(), plotOutput("tor_umap"),br(),
                                                  downloadButton(outputId="torPlot","Download Plot"),
                                                  actionButton("EX101","Explanation", icon = icon("book-open")),br(),
                                                  # withLoader(verbatimTextOutput("ex_umapp"),type="image",loader = "ectracker_popup.gif"), br(),
                                                  br(),br(),    
                                                  h6("Genomic Location"),br(),br(),
                                                  #  textInput("usergene","Enter Transcript Name",""),
                                                  withLoader(plotOutput("de_karyotype"),type="image",loader = "ectracker_popup.gif"),
                                                  br(), 
                                                 downloadButton(outputId="karyoplott","Download Plot"),
                                                 actionButton("EX2","Explanation", icon = icon("book-open")),
                                                 
                                                 withLoader(verbatimTextOutput("ex_karyo"),type="image",loader = "ectracker_popup.gif"),
                                                 
                                                 
                                                 br(),
                                                 
                                                 
                                                 
                                                 
                                                 actionButton("toTop4", "Top", icon = icon("arrow-alt-circle-up")),
                                       
                                                   ),
                                         #tabPanel( "Gene Regulaorty network",
                                         #  actionButton("action6000","Regulatory Network"),
                                         
                                         
                                         
                                         
                                         #withLoader(uiOutput("grn"),type="image",loader = "ectracker_popup.gif"),
                                                           # br(),downloadButton(outputId="grnPlot","Download Plot")
                                         # withLoader(verbatimTextOutput("grn_table"),type="image",loader = "ectracker_popup.gif"),
                                         
                                         #),
                                         
                                         tabPanel("Stouffer Score (Method II)",
                                                  
                                                  tabsetPanel( type= "tabs", 
                                                               
                                                               tabPanel("HPA Data",
                                                                        
                                                                        br(),br(),
                                                                      #  actionButton("hpa1_main","One-Click Analysis",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                                        
                                                                        selectInput("hpa11","Select Tissue",choices=c("Lymph Node",	"Tonsil",	"Appendix",	"Spleen",	"Bone Marrow"	,"Esophagus"	,"Skin",	"Colon",	"Rectum",	"Duodenum",	"Small Intestine",	"Stomach",	"Adipose Tissue",	"Lung",	"Placenta",	"Gallbladder",	"Urinary Bladder","Endometrium",	"Smooth Muscle","Fallopian Tube",	"Thyroid Gland",	"Ovary","Prostate",	"Kidney",	"Adrenal Gland",	"Brain",	"Salivary Gland","Pancreas",	"Skeletal Muscle",	"Liver",	"Heart Muscle")),
                                                                        
                                                                        
                                                                        h6("Stouffer Score "),
                                                                        withLoader(plotOutput("hpa_stouffer"),type="image",loader = "ectracker_popup.gif"),
                                                                        br(),  
                                                                      downloadButton(outputId="hpaplot","Download Plot"),
                                                                      actionButton("ex12","Explanation", icon = icon("book-open")),
                                                                      withLoader(verbatimTextOutput("ex_stouffer_hpa"),type="image",loader = "ectracker_popup.gif"),
                                                                      
                                                                      
                                                                      br(),br(),
                                                                        h6("One Sided Wilcoxon P-Value"),
                                                                        withLoader(uiOutput("hpa_set2"),type="image",loader = "ectracker_popup.gif"),
                                                                      downloadLink('stouffer_table_hpa', 'Download Table'),br(),br(),
                                                                        actionButton("toTop5", "Top", icon = icon("arrow-alt-circle-up")),
                                                                        
                                                                        
                                                               ),
                                                               
                                                               
                                                               br(), tabPanel("GTEx Data",br(),br(),
                                                                              selectInput("gtex11","Select Tissue",choices=c("Adrenal Gland","Adipose Tissue","Brain","Colon","Esophagus","Fallopian Tube","Bladder","Heart", "Muscle","Kidney","Liver","Lung","Ovary","Prostate","Salivary Gland","Skin","Small Intestine","Spleen","Stomach","Thyroid","Breast","Nerve","Blood Vessel","Uterus","Cervix/Uterine","Pituitary","Vagina","Pancreas")),
                                                                              
                                                                              br(),h6("Stouffer Score"),
                                                                              withLoader(plotOutput("gtex_stouffer"),type="image",loader = "ectracker_popup.gif"),
                                                                              br(), 
                                                                              downloadButton(outputId="gtexplot","Download Plot"),
                                                                              actionButton("ex12","Explanation", icon = icon("book-open")),
                                                                              withLoader(verbatimTextOutput("ex_stouffer_gtex"),type="image",loader = "ectracker_popup.gif"),
                                                                              
                                                                              br(),br(),
                                                                              
                                                                              h6("One Sided Wilcoxon P-Value"),
                                                                              withLoader(uiOutput("gtex_pstouffer"),type="image",loader = "ectracker_popup.gif"),
                                                                              downloadLink('stouffer_table_gtex', 'Download Table'),br(),br(),
                                                                              
                                                                              
                                                                              actionButton("toTop100", "Top", icon = icon("arrow-alt-circle-up")),
                                                               ) )             
                                         )
                                         
                                         
                             )),
                   
                   
                   tabPanel(tags$b("6. Gene Regulatory Network"),
                            tabsetPanel(
                              tabPanel("Cell Cluster",
                                       br(),br(),
                                       sliderInput("grn_obs", "Select Cluster Number ",
                                                   min = 0, max = 20,
                                                   value = 0),
                                       h6("Generate Specific Cluster Table"),br(),
                                       
                                       
                                       withLoader(verbatimTextOutput("grn_clust"),type="image",loader = "ectracker_popup.gif"),  
                                       
                              ),
                              
                              
                              tabPanel("CellEnrich",
                                       
                                      # actionButton("cell_grn_main","One-Click Analysis",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                       br(),br(),
                                      h6("Gene Regulatory Network for CellEnrich (Interaction between Transcription Factors and Target Genes)"),
                                      withLoader(simpleNetworkOutput("cell_grn_1"),type="image",loader = "ectracker_popup.gif"),
                                      h6("Gene Regulatory Network for CellEnrich (Interaction between Target Genes and Tissue Signatures)"),
                                      withLoader(simpleNetworkOutput("cell_grn"),type="image",loader = "ectracker_popup.gif"),
                                      

                                       actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                       withLoader(verbatimTextOutput("ex_cellgrn"),type="image",loader = "ectracker_popup.gif"),
                                       
                                       
                                       br(),br(),
                                       h6("Gene Regulatory Table"),
                                       withLoader(uiOutput("cell_grn_table"),type="image",loader = "ectracker_popup.gif"),br(),
                                       actionButton("ex_cellgrn","Explanation", icon = icon("book-open")), downloadLink('download_cell_gt', 'Download Table'),br(),br(),

                                       withLoader(verbatimTextOutput("ex_cellgrn_table"),type="image",loader = "ectracker_popup.gif"),
                                       br(), actionButton("toTop6", "Top", icon = icon("arrow-alt-circle-up")),
                                       
                              ),tabPanel("TissueEnrich",
                                         
                                         br(),br(),#actionButton("grn_main","One-Click Analysis",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                         
                                         h6("Gene Regulatory Network for HPA Datasets (Interaction between Transcription Factors and Target Genes)"),
                                         withLoader(simpleNetworkOutput("hpa_grn"),type="image",loader = "ectracker_popup.gif"),
                                         h6("Gene Regulatory Network for HPA Datasets (Interaction between Target Genes and Tissue Signatures)"),
                                         withLoader(simpleNetworkOutput("hpa_grn_1"),type="image",loader = "ectracker_popup.gif"),
                                         
                                         br(),
                                         
                                        # downloadButton(outputId="hpa_grndownload","Download Plot"),
                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                         withLoader(verbatimTextOutput("ex_hpagrn"),type="image",loader = "ectracker_popup.gif"),
                                         
                                         h6("Gene Regulatory Table (HPA)"),
                                         withLoader(uiOutput("hpa_grn_table"),type="image",loader = "ectracker_popup.gif"),
                                         br(),
                                         downloadLink('download_hpa_gt', 'Download Table'),br(),br(),
                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                         withLoader(verbatimTextOutput("ex_hpagrn_table"),type="image",loader = "ectracker_popup.gif"),br(),br(),
                                         
                                         h6("Gene Regulatory Network for GTEx data (Interaction between Transcription Factors and Target Genes)"),
                                         withLoader(simpleNetworkOutput("gtex_grn_1"),type="image",loader = "ectracker_popup.gif"),
                                         h6("Gene Regulatory Network for GTEx data (Interaction between Target Genes and Tissue Signatures)"),
                                         withLoader(simpleNetworkOutput("gtex_grn"),type="image",loader = "ectracker_popup.gif"),
                                         br(),
                                         
                                         
                                         #downloadButton(outputId="gtex_grndownload","Download Plot"),
                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                         withLoader(verbatimTextOutput("ex_gtexgrn"),type="image",loader = "ectracker_popup.gif"),
                                         
                                         br(), h6("Gene Regulatory Table (GTEx)"),
                                         withLoader(uiOutput("gtex_grn_table"),type="image",loader = "ectracker_popup.gif"),
                                         br(),br(),br(),

                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),downloadLink('download_gtex_gt', 'Download Table'),
                                         withLoader(verbatimTextOutput("ex_gtexgrn_table"),type="image",loader = "ectracker_popup.gif"),
                                         actionButton("toTop7", "Top", icon = icon("arrow-alt-circle-up")),
                              )
                              
                              
                            )
                            
                   )               
                   
                   
                   
                   
                   
                    )
                      
                      
             ),
             
            
             
             
             
         
             
            
             
            
             tabPanel("Tutorial", icon = icon("angle-double-right"),
                     h6( "Check Out the Video Tutorial here",align="center"),
                    tags$iframe(width="1500", height="650", src="atlast_final.mp4", frameborder="0", allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"),
                   h6("Check out the tutorial page link here:"), tags$a(href="final_documentation.html", "Tutorial.html", style = "color:yellow"),
                  # h5("Tutorial",align="center")
                   br(),br(),br(),br(),br(),
                  #  downloadButton(outputId="tut_download","Download Plot"),
                      
                      ),
             
            
             
             
             
             tabPanel("Contact Us", icon = icon("angle-double-right"),
                    #h4(  "Contact Us",align="left", style = "color:pink"),
   tabsetPanel( type = "tabs",
               
                         
                tabPanel("Developer",
                         HTML('<left><img src="vishakha.jpg", width="20%", height= "20%" ></left>'),
h6("Vishakha"), 
h6("Ph.D Scholar"),
h6("Email: vishakhag@iiitd.ac.in
"),
               ), 

tabPanel("Developer",
         HTML('<left><img src="siddhant.jpg", width="20%", height= "20%" ></left>'),
         h6("Siddhant"), 
         h6("M.tech student"),
         h6("Email: siddhant18241@iiitd.ac.in
")
          ),
tabPanel("Web Designer",
         HTML('<left><img src="aayushi.jpg", width="20%", height= "20%" ></left>'),
         h6("Aayushi Mittal"), 
         h6("Ph.D Scholar"),
         h6("Email: aayushim@iiitd.ac.in
")   
),
tabPanel("Testing",
         HTML('<left><img src="sanjay.jpg", width="20%", height= "20%" ></left>'),
         
         h6("Sanjay kumar Mohanty"), 
         h6("Ph.D Scholar"),
         h6("Email: sanjaym@iiitd.ac.in
")   ),
tabPanel("Collaborator",
         HTML('<left><img src="sengupta.PNG", width="20%", height= "20%" ></left>'),
         
        # h6("Sengupta Lab") ,
       br(),br(),  tags$a(href="https://www.debarka.com/", "Sengupta Lab", style = "color:yellow")
         
         
         
         
         )





))
#tabPanel("Final Summary Report",
       #  downloadButton("report", "Generate report")
#)
             
             
             
           ##  tabPanel("Other Receptors",  icon = icon("angle-double-right"),
                      #sidebarLayout(
                       # sidebarPanel(
                        #  uiOutput("other_or")
                          
                        #),
                        #mainPanel(plotOutput("other_umap"))
                      #)
                      #tableOutput("f_or"))
            # ),
             
            # tabPanel("User Input", icon = icon("angle-double-right"),
             #      sidebarLayout(
              #          sidebarPanel(("Enter gene name"),
               #           textInput("usergene","","")
                #        ),
                 #       mainPanel(plotOutput("user_umap"))
                  #    )
                      #tableOutput("f_or"))
             #),
             #tabPanel("GSVA*", icon = icon("angle-double-right"),
                    #  sidebarLayout(
                       # sidebarPanel(fileInput("csv_file","Upload the file"),
                              #       checkboxInput("testme1", "Sample Dataset"), 
                                                   
                                #     actionButton("action30"," GSVA"),br(),br(),
                              #       textInput("gene_name_gsva","Enter Transcript Name"),
                          #           textInput("gene_name_gsva_freq","Enter Frequency"),
                                     
                        #             actionButton("action300","PLOT GSVA")
                                     
                      #  ),
                      #  mainPanel(
                      #    print("It is an optional function. It will take around 6 hours or more"),
                       #   verbatimTextOutput("dimen_gsva"),
                       #   withLoader(tableOutput("gsva"),type="image",loader = "ectracker_popup.gif"),
                          
                        #  withLoader(plotOutput("gsva_plot"),type="image",loader = "ectracker_popup.gif")
                          
                       # ) 
                     # ))
)
)
)
