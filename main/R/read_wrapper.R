
#' @importFrom utils read.csv untar unzip
read_data<-function(x=x, format="txt", header=T, sep=',', quote=NULL){
  
  fileEncoding = ""
  curDir = getwd()
  
  tmpdir <- file.path(dirname(x),"temp")
  
  # tar.gz
  if(any(grep(".tar.gz", basename(x)))){
    
      untar(x, compressed = "gzip", exdir = tmpdir)    
      files = list.files(path = tmpdir, full.names = TRUE, recursive = TRUE)
      
      
    } else if(any(grep(".zip", basename(x)))){
      
      tryCatch(
        {
          unzip(x,  exdir = tmpdir, overwrite = T)
          files = list.files(path = tmpdir, full.names = TRUE, recursive = TRUE)
        },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(e)
      }
      )
      
    }else if(format == "txt" && 
             any(grep("(.csv|.txt|.tsv)", basename(x)))){
      
      data =  input_csv(x, header, sep, quote)
      base::unlink(file.path(tmpdir,"/"), recursive = T)
      
      return(data)
      
    }else{stop("Invalid file format:")}
    

  if(format == "txt" && 
     any(grep("(.csv|.txt|.tsv)", files))){
    
    x = files[grep("(.csv|.txt|.tsv)", files)]
    
  }else if (format=="10x" &&
            any(grep("matrix.mtx", files)) && 
            any(grep("genes.tsv", files)) && 
            any(grep("barcodes.tsv", files))){
    

    x = dirname(files[grep("matrix.mtx", files)])
    
  }else {
    stop("Invalid file format.")
  }
  

  
  data = switch(format,
                "10x" = read10X(x),
                "txt" = input_csv(x, header, sep, quote))
  
  base::unlink(tmpdir, recursive = T)
  
  setwd(dir = curDir)
  
  return(data)
}



# u <- "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
# 
# 
# u.data<-read_data(u,format = "10x")
# 
# tz <-"E:/Projects/dropClust/data/pbmc3k_filtered_gene_bc_matrices.tar.gz"
# 
# tz.data<-read_data(tz,format = "10x")
# 
# 
# zip <-"E:/Projects/dropClust/data/formats/hg19.zip"
# 
# zip.data<-read_data(zip,format = "10x")
# 
# 
# csv <-"E:/Projects/batcheffect/Datasets/ColorectalTumor_CV1000_Processed_Data.csv"
# 
# csv.data<-read_data(csv,format = "txt")
# 
# csv.zp <-"E:/Projects/dropClust/data/formats/ColorectalTumor_CV1000_Processed_Data.zip"
# 
# csv.zp.data<-read_data(csv.zp,format = "txt")
# 
# 

# zp.10x = "C:/Projects/dropClust/data/pbmc3k.zip"
# 
# zz = read_data(zp.10x, format = "10x")

# x = "C:/Users/hp/AppData/Local/Temp/RtmpUVP6BR/7249d4e2363455da7da5fcb9/0.zip"

 # zz = read_data("C:/Users/hp/AppData/Local/Temp/RtmpUVP6BR/7249d4e2363455da7da5fcb9/0.zip", format = "txt")
