part1_RaMS <- function(FN1,ESImode,MS2_filter){
  # tt <- Sys.time()
  # setwd('D:/bio_inf/LipidIN/LipidIN 4-level hierarchical library v1.0.4/demo pos')
  # FN1 <- 'RQC-1-ddms.mzML';ESImode <- 'p';MS2_filter <- 0.1
  library(RaMS)
  library(dplyr)
  cat("\033[32mgrabMzmlData reading MS2!\n\033[0m")
  a <- grabMzmlData(FN1,grab_what=c("MS2"))
  a <- a$MS2
  # a$fragmz <- round(a$fragmz,4)
  # a <- a[which(a$int>0),]
  a$row_id <- 1:nrow(a)
  nrow.count <- a %>% group_by(rt,premz) %>% summarise(row_numbers=list(row_id),.groups='drop')
  nrow.count <- nrow.count %>% mutate(row_numbers=sapply(row_numbers,function(x) paste(x,collapse=", ")))
  rr <- lapply(1:nrow(nrow.count),function(i){
    rrr <- unlist(nrow.count[i,3]) %>% strsplit(", ") %>% unlist() %>% as.numeric()
    c <- a[rrr,]
    count.1 <- list(
      num=i,
      PrecursorMZ=as.numeric(c[1,2]),
      PrecursorIntensity=1,
      mode=ESImode,
      rt=as.numeric(c[1,1]),
      MS2mz=data.frame(da.temp.mz=c[,3],da.temp.intensity=c[,4])
    )
    aa <- count.1$MS2mz
    colnames(aa) <- c('da.temp.mz','da.temp.intensity')
    aa <- aa[which(aa$da.temp.intensity>MS2_filter*max(count.1$MS2mz[,2])),]
    count.1$MS2mz <- aa
    return(count.1)
  })
  neg.list <- do.call(list,rr)
  file3 <- paste(gsub('.mzML','',FN1),'_treated.rda',sep='')
  save(neg.list,file=file3)
  cat("\033[32mDONE!!!!!!\n \033[0m")
  # Sys.time()-tt
}


