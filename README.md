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
All the libraries can be installed by two main commands <br/>
1. install.packages(“Package_name”) <br/>
2. install.packages("BiocManager") <br/>
   BiocManager::install("Package_name")<br/>

**How to run files locally** <br/>
1. Download the github repository <br/>
2. Insatll all the required packages <br/>
3. Put all the files and folders in one directory.<br/>
4. Run the ui.R and server.R <br/>
5. Once user have installed all the libraries Run App option will appear. User can press that and EcTracker will work in users server.


