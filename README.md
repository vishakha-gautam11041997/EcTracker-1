# EcTracker
EcTracker: Tracks silent genes activation in scRNA-seq data
## Workflow
<img src="main/www/workflow.PNG"> <br/>
###
**KEY POINTS**<br/> 
1. EcTracker possesses two distinct modules i.e. CellEnrich and TissueEnrich, to allow detection of cell-type or tissue-specific genes in the user-supplemented scRNA-sequencing dataset.
2. CellEnrich utilizes the widely used AUCell package and leverages on the bonafide fetal and adult tissue markers curated from the large-scale single-cell atlases.
3. EcTracker provides an enhanced visualization method for tracking the chromosomal locations of the selected genesets in a cluster-wise manner. 
4. EcTracker implements DoRothEA and its associated statistical method VIPER that allow users to identify central transcription factors driving the expression of ectopic transcripts present in the user datasets.

**Installation of the libraries** <br/>
All the libraries can be installed by two main commands: <br/>
1. install.packages(“Package_name”) <br/>
2. install.packages("BiocManager") <br/>
   BiocManager::install("Package_name")<br/>

**List of libraries**<br/>
karyoploteR (version: 1.14.1), org.Hs.eg.db (version:3.11.4), TxDb.Hsapiens.UCSC.hg19.knownGene (version:3.2.2), GenomicFeatures (version:1.40.1), GSEABase (version:1.50.1), graph (version:1.66.0), annotate (version:1.66.0), AnnotationDbi (version:1.50.3), GenomicRanges (version:1.40.0), GenomeInfoDb (version:1.24.2), viper (version: 1.22.0), dplyr (version:1.0.3), bcellViper (version:1.24.0), networkD3 (version:0.4), dorothea (version:1.0.1), AUCell (version:1.10.0), shiny (version:1.6.0), shinyFiles (version:0.9.0), shinyalert (version:2.0.0), shinyWidgets (version:0.5.6), shinybusy (version:0.2.2), igraph (version:1.2.6), shinyjs (version:2.0.0), shinycustomloader (version:0.9.0), shinythemes (version:1.2.0), Seurat (version:3.9.9.9008), ggplot2 (version:3.3.3)<br/>                                                     


**How to run files locally** <br/>
1. Download the github repository <br/>
2. Install all the required packages <br/>
3. Put all the files and folders in one directory.<br/>
4. Download the required files from drive link <br/> https://drive.google.com/drive/folders/121fOvDvO7MHIrUOtMBz8-o2EiXXSilXp?usp=sharing <br/>
5. Download the sample file from the link https://drive.google.com/drive/folders/1Hk__Muaww1aizAp1K81LYr90wtk7aBYT?usp=sharing <br/>
6. "www" folder should be there in order to get videos and loader gifs. <br/>
7. Set the working directory by command setwd("folder name") <br/>
8. Open the ui.R and server.R in RStudio <br/>
9. User can directly click on Run App or write command shiny::runApp()

