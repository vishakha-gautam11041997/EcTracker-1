

#library(V8)
jscode <- "shinyjs.toTop = function() {window.scrollTo(0, 0);}" ### for to top button


shinyUI(fluidPage(
  theme=shinytheme("cyborg"),
  tags$style(HTML("
        .tabs-above > .nav > li[class=active] > a {
           background-color: #000;
           color: #FFFF;
        }")),


 useShinyalert(),

  navbarPage(em("EcTracker"),
                   tabPanel("Home" ,
                            icon = icon("home"),

                           HTML('<center><img src="aeae.gif", width="58%", height= "58%" ></center>'),



                          column(8, align="center", offset = 2 ,print("Developed by: "), tags$u("Vishakha Gautam"), print("&" ), tags$u("Siddhant Kalra"),

                          br(),print("For more details:"),
                          tags$head(tags$style(HTML("a {color: blue}"))),
                          tags$a(href="https://ahuja-lab.in/", "The Ahuja Lab", style = "color:yellow"))


                            ),


             tabPanel("Analysis" , icon = icon("angle-double-right"),

                     navlistPanel(
                     tabPanel(tags$b("1. Upload") ,

                                  mainPanel(
                                    fileInput("file_batch","Upload matrix for Batch Correction"),
                                 em(h6("Note:" )),br(),
                                 em(h6("Takes atleast 15-20 minutes to run", style = "color:#FF0000;")),

                                    actionButton("batch1","Batch Correction"),
                                    withLoader(uiOutput("batch_corrected"),type="image",loader = "new_loader.gif"),

                                    downloadLink('downloadData_batch', 'Download File for further analysis', class = "butt"),
                                  tags$head(tags$style(".butt{color: #FFFF00;}")),
                                    br(),br(),


                                    fileInput("file","Upload Expression Matrix (Maximum Size Limit is 10 GB)"),
                                    #,br(),h6("Maximum Size Limit is 10GB"),

                                           br(), h5("OR"),br(),
                                            fileInput("File10","Upload 10X Data Files (barcode.tsv, genes.tsv & matrix.mtx) ", multiple = TRUE),

                                            checkboxInput("testme", "Sample Data (Transdifferentiation)",
                                                          value = FALSE),

                                            br(),br(),
                                             h6("Information about Data"),
                                           withLoader(tableOutput("dat"),type="image",loader = "new_loader.gif"),
                                           h6("Dimensions"),
                                           withLoader(verbatimTextOutput("dimen"),type="image",loader = "new_loader.gif"),
                                        # actionButton("actionsc1","Batch Correction method"),
                                       #  withloader()


                                  )

                                ),

                       tabPanel(tags$b("2. scData Processing"),

tabsetPanel(type = "tabs",

            tabPanel("dropClust (Faster)",
                    actionButton("dropclust1","One-Click Analysis",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("paper-plane")),
                    h6("Scatter Plot of PCA"),br(),

                     withLoader( plotOutput("drop_clust_plot"),type="image",loader = "Asset 1.jpg"),
                    downloadButton(outputId="drop_clust_scatter","Download Plot")

                        ),


            tabPanel("Seurat (Slower)",
                                actionButton("seurat1","One-Click Analysis*",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("paper-plane")),
                                br(),br(),
                                h6("Object"),
                             withLoader(verbatimTextOutput("seu_object"),type="image",loader = "Asset 1.jpg"),br(),
                                actionButton("EX1","Explanation",icon = icon("book-open")),
                                withLoader(verbatimTextOutput("ex_object"),type="image",loader = "Asset 1.jpg"),

                       br(),br(),

                                h6("Variable Feature Plot"),

                               withLoader( plotOutput("p1"),type="image",loader = "Asset 1.jpg"),

                                br(), downloadButton(outputId="dndPlot","Download Plot"),
                     actionButton("EX2","Explanation", icon = icon("book-open")),

                     withLoader(verbatimTextOutput("ex_variable"),type="image",loader = "Asset 1.jpg"),



                    br(),br(),
                      h6("Variable Features"),
                      withLoader(plotOutput("plot_vf"), type="image", loader = "Asset 1.jpg"),
                     br(), downloadButton(outputId="vfPlot","Download Plot"),
                     actionButton("EX4","Explanation",icon = icon("book-open")),

                     withLoader(verbatimTextOutput("ex_vfeature"),type="image",loader = "Asset 1.jpg"),




                      br(),br(),
                      h6("PCA"),
                      withLoader(plotOutput("pca"),type="image",loader = "Asset 1.jpg"),
                      br(),downloadButton(outputId="pcaPlot","Download Plot"),
                      actionButton("EX6","Explanation", icon = icon("book-open")),
                       withLoader(verbatimTextOutput("ex_pca"),type="image",loader = "Asset 1.jpg"),


                     br(),br(),

                      h6("Heat Map"),
                      withLoader(plotOutput("heatmap"),type="image",loader = "Asset 1.jpg"),
                      br(),

                      downloadButton(outputId="heatmapPlot","Download Plot"), actionButton("EX7","Explanation", icon = icon("book-open")),

                      withLoader(verbatimTextOutput("ex_heatmap"),type="image",loader = "Asset 1.jpg"),


                     br(),
                     br(),
                     h6("JackStraw Plot"),

                     withLoader(plotOutput("jackstraw"),type="image",loader = "Asset 1.jpg"),
                      br(),

                     downloadButton(outputId="jackstrawPlot","Download Plot"),actionButton("EX8","Explanation", icon = icon("book-open")),

                     withLoader(verbatimTextOutput("ex_jackstraw"),type="image",loader = "Asset 1.jpg"),

                     br(),
                     br(),
                     h6("Elbow Plot"),
                      withLoader(plotOutput("elbow"),type="image",loader = "Asset 1.jpg"),
                     br(),
                    downloadButton(outputId="elbowPlot","Download Plot"),actionButton("EX9","Explanation", icon = icon("book-open")),

                    withLoader(verbatimTextOutput("ex_elbow"),type="image",loader = "Asset 1.jpg"),
                    br(),br(),



                      h6("UMAP"),
                      withLoader(plotOutput("umap"),type="image",loader = "Asset 1.jpg"),
                      br(),
                      downloadButton(outputId="umapPlot","Download Plot"),
                      actionButton("EX10","Explanation", icon = icon("book-open")),

                      withLoader(verbatimTextOutput("ex_umap"),type="image",loader = "Asset 1.jpg"),
                    actionButton("toTop", "Top", icon = icon("arrow-alt-circle-up")),

                    ))


),


                    tabPanel(tags$b("3. DEG Analysis"),
                             tabsetPanel(type = "tabs",
                                         tabPanel("dropClust",

                                                  h6("Differentially Expressed Genes"),br(),
                                                  withLoader( uiOutput("drop_clust_markers"),type="image",loader = "Asset 1.jpg"),
                                                  downloadLink('drop_clust_de', 'Download', class = "butt"),

br(),br(),

                                                  #("Heat Map of the genes"),br(),
                                                 # withLoader( plotOutput("drop_clust_heatmap"),type="image",loader = "Asset 1.jpg"),
                                                  ),

                                         tabPanel("Seurat",  sliderInput("log_fc", "Select Fold Change Threshold (log)",
                                         min = 0, max = 20,step=0.25,
                                         value = 1),
                             selectInput("test_use","Choose the test for DE genes",choice=c("Wilcoxon","Bimod (Likelihood-ratio test for single cell gene expression)","ROC analysis","Student's t-test","Negative binomial generalized linear model","Poisson generalized linear model","Logistic Regression ","MAST","DESeq2")),
                             br(),
                               downloadLink('downloadData', 'Download', class = "butt"),
                      br(),br(),
                     h6("Markers Table"),
                    withLoader(tableOutput("markers_display"),type="image",loader = "Asset 1.jpg"),
                     useShinyjs(),
                     br(),
                    fileInput("meta","Upload Meta File"),
                    checkboxInput("testme3", " Transdifferentiation (Meta file)",
                                  value = FALSE),br(),
                    actionButton("m1","Table with Meta File "),
                    withLoader(verbatimTextOutput("merge_meta"),type="image",loader = "Asset 1.jpg"),
                    downloadLink('downloadData_meta', 'Download', class = "butt"),
                     br(),
                   br(),
                   h6("Feature Plot"),
                   textInput("feature_plot","Enter Gene Name"),
                   withLoader(plotOutput("feature_plot_result"),type="image",loader = "Asset 1.jpg"),
                   br(),
                   downloadButton(outputId="torPlot","Download Plot"),
                  #actionButton("EX101","Explanation", icon = icon("book-open")),br(),
                   br(), br(),
                   h6("Genomic Location"),br(),
                   withLoader(plotOutput("de_karyotype"),type="image",loader = "Asset 1.jpg"),
                   br(),
                   downloadButton(outputId="karyoplott","Download Plot"),
                   #actionButton("EX2","Explanation", icon = icon("book-open")),
                   
                  # withLoader(verbatimTextOutput("ex_karyo"),type="image",loader = "Asset 1.jpg"),
                   
                   br(),br(),
                    h6("Generate Cluster Table"),
                      withLoader(verbatimTextOutput("cluster_ident"),type="image",loader = "Asset 1.jpg"),
                     downloadLink('downloadData_cluster_ident', 'Download', class = "butt"),br(),br(),
                     extendShinyjs(text = jscode,functions = c("toTop")),
                    actionButton("toTop1", "Top", icon = icon("arrow-alt-circle-up")),



                    ))),





                   tabPanel(tags$b("4. CellEnrich "),

                            tabsetPanel(type = "tabs",
                                        tabPanel(" AUCell (Method I)",
                                                 br(),br(),br(),

                                                 br(),
                                                 actionButton("action1","One-Click Enrichment",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("paper-plane")),

                                                 selectInput("aucell","Select Signature",choices=c("Ectodermal cell differentiation","Ectodermal development","Embryonic_stem_cell","Adult b_cell","Adult b_cell_plasmocyte","Fetal fibroblast","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult basal_cell","Adult CB CD34+23","Devloping heart atrial cardiomyocyte","Devloping heart ventricular cardiomyocyte","Adult chondrocyte","Adult chemosensory ","Adult dendritic_cell ","Adult endothelial_cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte_progenitor","Adult enterocyte","Adult epithelial cell (intermediated)","Adult epithelial_cell","Adult erythroid_cell ","Adult erythroid progenitor cell (RP high)","Adult fasciculata_cell ","Fetal B cells","Fetal Basophil_Mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal Collecting duct lineage","Fetal Connecting tubule lineage","Fetal Distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal Enteroendocrine cells_intestine_pancreas","Fetal Enteroendocrine cells_stomach","Fetal Erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ Microglia","Fetal ILC 3","Fetal Islet beta cells","Fetal Islet delta cells","Fetal LEC","Fetal Loop of Henle lineage","Fetal Macrophages","Fetal Meg progenitors","Fetal Meg","Fetal Megakaryoblasts","Fetal Microglia","Fetal Nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal Perivascular macrophages","Fetal Phagocytic macrophages","Fetal Plasma cells","Fetal Podocyte lineage","Fetal Proximal tubule lineage","Fetal PTPRC+ Microglia","Fetal Pulmonary neuroendocrine cells","Fetal Renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ Microglia","Fetal TRAF1+ APCs","Fetal Ureteric tip","Fetal Ureteric trunk","Fetal Antigen-presenting macrophages","Fetal Endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal skeletal muscle cell","Fetal_mesenchymal_progenitor","Fetal_neuron ","Fetal_stromal_cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet_cell ","Adult hepatocyte_Endodermal cell","Adult hESC87","Adult immature sertoli cell (Pre-Sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 Macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal_cell","Adult smooth_muscle_cell","Adult stratified_epithelial_cell","Adult T_cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")),

                                                 br(), h6("Cell Ranking"),
                                                 withLoader(plotOutput("cell_rank"),type="image",loader = "Asset 1.jpg"),
                                                 br(), actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_cellrank"),type="image",loader = "Asset 1.jpg"),


                                                 br(), br(),   h6(" Cell Assignment"),
                                                 withLoader(plotOutput("au_cell_R"),type="image",loader = "Asset 1.jpg"),br(),
                                                 actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_assign"),type="image",loader = "Asset 1.jpg"),


                                                 # br(), br(),  h6("UMAP"),
                                                 # withLoader(plotOutput("tsne"),type="image",loader = "Asset 1.jpg"),
                                                 # br(),  actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 # withLoader(verbatimTextOutput("ex_uma"),type="image",loader = "Asset 1.jpg"),

                                              br(),   br(),  h6("UMAP for Signatures"),
                                                 withLoader(plotOutput("c_tsne"),type="image",loader = "Asset 1.jpg"),
                                                 br(),
                                                  actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_umappp"),type="image",loader = "Asset 1.jpg"),
                                                downloadButton(outputId="cell_assign","Download All Plots"),
                                                downloadLink(outputId="cell_assign_score","Download Cell Assignment Table", class = "butt"),


                                                 br(), br(),  actionButton("toTop3", "Top", icon = icon("arrow-alt-circle-up")),


                                        ),

                                        tabPanel(" Stouffer Based Enrichment Analysis (Method II)",
                                                 br(),br(),br(),br(),
                                                 selectInput("st_sc","Select Signature",choices=c("Devloping heart atrial cardiomyocyte","Devloping heart ventricular cardiomyocyte","Embryonic_stem_cell","Ectodermal cell differentiation","Ectodermal development","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Myofibroblast differentiation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult b_cell","Adult b_cell_plasmocyte","Adult basal_cell","Adult CB CD34+23","Adult chondrocyte","Adult chemosensory ","Adult dendritic_cell ","Adult endothelial_cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte_progenitor","Adult enterocyte","Adult epithelial cell (intermediated)","Adult epithelial_cell","Adult erythroid_cell ","Adult erythroid progenitor cell (RP high)","Adult fasciculata_cell ","Fetal B cells","Fetal Basophil_Mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal Collecting duct lineage","Fetal Connecting tubule lineage","Fetal Distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal Enteroendocrine cells_intestine_pancreas","Fetal Enteroendocrine cells_stomach","Fetal Erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ Microglia","Fetal ILC 3","Fetal Islet beta cells","Fetal Islet delta cells","Fetal LEC","Fetal Loop of Henle lineage","Fetal Macrophages","Fetal Meg progenitors","Fetal Meg","Fetal Megakaryoblasts","Fetal Microglia","Fetal Nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal Perivascular macrophages","Fetal Phagocytic macrophages","Fetal Plasma cells","Fetal Podocyte lineage","Fetal Proximal tubule lineage","Fetal PTPRC+ Microglia","Fetal Pulmonary neuroendocrine cells","Fetal Renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ Microglia","Fetal TRAF1+ APCs","Fetal Ureteric tip","Fetal Ureteric trunk","Fetal Antigen-presenting macrophages","Fetal Endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal fibroblast","Fetal skeletal muscle cell","Fetal_mesenchymal_progenitor","Fetal_neuron ","Fetal_stromal_cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet_cell ","Adult hepatocyte_Endodermal cell","Adult hESC87","Adult immature sertoli cell (Pre-Sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 Macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal_cell","Adult smooth_muscle_cell","Adult stratified_epithelial_cell","Adult T_cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")),
                                                 br(),

                                                 h6("Stouffer Plot"),
                                                 withLoader(plotOutput("stouffer_sc_plot"),type="image",loader = "Asset 1.jpg"),
                                                 br(), downloadButton(outputId="cell_stoudownload","Download Plot"),

                                                 actionButton("ex12","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_stouffer"),type="image",loader = "Asset 1.jpg"),

                                                 br(),   br(),h6("One Sided Wilcoxon P-Value"),
                                                 withLoader(uiOutput("Stouffer_sc"),type="image",loader = "Asset 1.jpg"), br(),
                                                downloadLink('stouffer_table_cell', 'Download Table', class = "butt"),br(),br(),




                                        )




                            )),

                   tabPanel( tags$b("5. TissueEnrich"),
                             tabsetPanel(type = "tabs",



                                         tabPanel("Hypergeometric method (Method I)",
                                                  br(),
                                                  sliderInput("obs", "Select Cluster Number ",
                                                              min = 0, max = 20,
                                                              value = 1),
                                                   downloadLink('downloadData_cluster_name', 'Download Cluster Table', class = "butt"),br(),br(),

                                                  br(),br(),   numericInput("dataset", "Enter Dataset Number: 1: HPA ; 2: GTEx", 1, min = 1, max = 2), br(),

                                                  h6("Generate TissueEnrich Table"),

                                                  withLoader(verbatimTextOutput("tissue_detail"),type="image",loader = "Asset 1.jpg"),

                                                  downloadLink('downloadtissuedata ', 'Download', class = "butt"),

                                                  br(),
                                                  br(),

                                                 h6("Fold Change Plot"),

                                                  withLoader(plotOutput("plot_fold"),type="image",loader = "Asset 1.jpg"),
                                                  br(),
                                                  downloadButton(outputId="tissue_fold1","Download Plot"),


                                                  br(),br(),


                                                  h6(" Log10 (P-Value) Plot"),

                                                  withLoader(plotOutput("plot_log"),type="image",loader = "Asset 1.jpg"),
                                                  br(),   downloadButton(outputId="tissue_log1","Download Plot"),


                                              


                                                
                                                 br(),
                                                 h6("Tissue-Specific Genes"),

                                                 uiOutput("t_or"),
                                                 uiOutput("tt_or"),br(),br(),
                                                 br(), #plotOutput("tor_umap"),br(),
                                                 #downloadButton(outputId="torPlot","Download Plot"),
                                                # actionButton("EX101","Explanation", icon = icon("book-open")),br(),
                                                # br(),br(),


                                                 actionButton("toTop4", "Top", icon = icon("arrow-alt-circle-up")),

                                                   ),


                                         tabPanel("Stouffer Based Enrichment Analysis (Method II)",

                                                  tabsetPanel( type= "tabs",

                                                               tabPanel("HPA Data",

                                                                        br(),br(),
                                                                        selectInput("hpa11","Select Tissue",choices=c("Lymph Node-HPA",	"Tonsil-HPA",	"Appendix-HPA",	"Spleen-HPA",	"Bone Marrow-HPA"	,"Esophagus-HPA"	,"Skin-HPA",	"Colon-HPA",	"Rectum-HPA",	"Duodenum-HPA",	"Small Intestine-HPA",	"Stomach-HPA",	"Adipose Tissue-HPA",	"Lung-HPA",	"Placenta-HPA",	"Gallbladder-HPA",	"Urinary Bladder-HPA","Endometrium-HPA",	"Smooth Muscle-HPA","Fallopian Tube-HPA",	"Thyroid Gland-HPA",	"Ovary-HPA","Prostate-HPA",	"Kidney-HPA",	"Adrenal Gland-HPA",	"Brain-HPA",	"Salivary Gland-HPA","Pancreas-HPA",	"Skeletal Muscle-HPA",	"Liver-HPA",	"Heart Muscle-HPA")),

                                                                        #selectInput("hpa11","Select Tissue",choices=c("Lymph Node",	"Tonsil",	"Appendix",	"Spleen",	"Bone Marrow"	,"Esophagus"	,"Skin",	"Colon",	"Rectum",	"Duodenum",	"Small Intestine",	"Stomach",	"Adipose Tissue",	"Lung",	"Placenta",	"Gallbladder",	"Urinary Bladder","Endometrium",	"Smooth Muscle","Fallopian Tube",	"Thyroid Gland",	"Ovary","Prostate",	"Kidney",	"Adrenal Gland",	"Brain",	"Salivary Gland","Pancreas",	"Skeletal Muscle",	"Liver",	"Heart Muscle")),


                                                                        h6("Stouffer Score "),
                                                                        withLoader(plotOutput("hpa_stouffer"),type="image",loader = "Asset 1.jpg"),
                                                                        br(),
                                                                      downloadButton(outputId="hpaplot","Download Plot"),
                                                                      actionButton("ex12","Explanation", icon = icon("book-open")),
                                                                      withLoader(verbatimTextOutput("ex_stouffer_hpa"),type="image",loader = "Asset 1.jpg"),


                                                                      br(),br(),
                                                                        h6("One Sided Wilcoxon P-Value"),
                                                                        withLoader(uiOutput("hpa_set2"),type="image",loader = "Asset 1.jpg"),
                                                                      downloadLink('stouffer_table_hpa', 'Download Table', class = "butt"),br(),br(),
                                                                        actionButton("toTop5", "Top", icon = icon("arrow-alt-circle-up")),


                                                               ),


                                                               br(), tabPanel("GTEx Data",br(),br(),
                                                                              selectInput("gtex11","Select Tissue",choices=c("Adrenal Gland","Adipose Tissue","Brain","Colon","Esophagus","Fallopian Tube","Bladder","Heart", "Muscle","Kidney","Liver","Lung","Ovary","Prostate","Salivary Gland","Skin","Small Intestine","Spleen","Stomach","Thyroid","Breast","Nerve","Uterus","Cervix/Uterine","Pituitary","Vagina","Pancreas")),

                                                                              br(),h6("Stouffer Score"),
                                                                              withLoader(plotOutput("gtex_stouffer"),type="image",loader = "Asset 1.jpg"),
                                                                              br(),
                                                                              downloadButton(outputId="gtexplot","Download Plot"),
                                                                              actionButton("ex12","Explanation", icon = icon("book-open")),
                                                                              withLoader(verbatimTextOutput("ex_stouffer_gtex"),type="image",loader = "Asset 1.jpg"),

                                                                              br(),br(),

                                                                              h6("One Sided Wilcoxon P-Value"),
                                                                              withLoader(uiOutput("gtex_pstouffer"),type="image",loader = "Asset 1.jpg"),
                                                                              downloadLink('stouffer_table_gtex', 'Download Table', class = "butt"),br(),br(),


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
                                                   value = 1),
                                       h6("Generate Specific Cluster Table"),br(),


                                       withLoader(verbatimTextOutput("grn_clust"),type="image",loader = "Asset 1.jpg"),

                              ),


                              tabPanel("CellEnrich",

                                       br(),br(),
                                      h6("Gene Regulatory Network for CellEnrich (Interaction between Transcription Factors and Target Genes)"),
                                      withLoader(simpleNetworkOutput("cell_grn_1"),type="image",loader = "Asset 1.jpg"),
                                      h6("Gene Regulatory Network for CellEnrich (Interaction between Target Genes and Tissue Signatures)"),
                                      withLoader(simpleNetworkOutput("cell_grn"),type="image",loader = "Asset 1.jpg"),


                                       actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                       withLoader(verbatimTextOutput("ex_cellgrn"),type="image",loader = "Asset 1.jpg"),


                                       br(),br(),
                                       h6("Gene Regulatory Table"),
                                       withLoader(uiOutput("cell_grn_table"),type="image",loader = "Asset 1.jpg"),br(),
                                       actionButton("ex_cellgrn","Explanation", icon = icon("book-open")),
                                      downloadLink('download_cell_gt', 'Download Table', class = "butt"),br(),br(),

                                       withLoader(verbatimTextOutput("ex_cellgrn_table"),type="image",loader = "Asset 1.jpg"),
                                       br(), actionButton("toTop6", "Top", icon = icon("arrow-alt-circle-up")),

                              ),tabPanel("TissueEnrich",

                                         br(),br(),

                                         h6("Gene Regulatory Network for HPA Datasets (Interaction between Transcription Factors and Target Genes)"),
                                         withLoader(simpleNetworkOutput("hpa_grn"),type="image",loader = "Asset 1.jpg"),
                                         h6("Gene Regulatory Network for HPA Datasets (Interaction between Target Genes and Tissue Signatures)"),
                                         withLoader(simpleNetworkOutput("hpa_grn_1"),type="image",loader = "Asset 1.jpg"),

                                         br(),

                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                         withLoader(verbatimTextOutput("ex_hpagrn"),type="image",loader = "Asset 1.jpg"),
                                         br(),br(),

                                         h6("Gene Regulatory Table (HPA)"),
                                         withLoader(uiOutput("hpa_grn_table"),type="image",loader = "Asset 1.jpg"),
                                         br(),
                                         downloadLink('download_hpa_gt', 'Download Table', class = "butt"),br(),br(),
                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                         withLoader(verbatimTextOutput("ex_hpagrn_table"),type="image",loader = "Asset 1.jpg"),br(),br(),

                                         h6("Gene Regulatory Network for GTEx data (Interaction between Transcription Factors and Target Genes)"),
                                         withLoader(simpleNetworkOutput("gtex_grn_1"),type="image",loader = "Asset 1.jpg"),
                                         h6("Gene Regulatory Network for GTEx data (Interaction between Target Genes and Tissue Signatures)"),
                                         withLoader(simpleNetworkOutput("gtex_grn"),type="image",loader = "Asset 1.jpg"),
                                         br(),


                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                         withLoader(verbatimTextOutput("ex_gtexgrn"),type="image",loader = "Asset 1.jpg"),

                                         br(),br(), h6("Gene Regulatory Table (GTEx)"),
                                         withLoader(uiOutput("gtex_grn_table"),type="image",loader = "Asset 1.jpg"),
                                         br(),

                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),downloadLink('download_gtex_gt', 'Download Table', class = "butt"),
                                         withLoader(verbatimTextOutput("ex_gtexgrn_table"),type="image",loader = "Asset 1.jpg"), br(),br(),
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
                   br(),br(),br(),br(),br(),

                      ),





             tabPanel("Contact Us", icon = icon("angle-double-right"),
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

       br(),br(),  tags$a(href="https://www.debarka.com/", "Sengupta Lab", style = "color:yellow")




         )





))

)
)
)
