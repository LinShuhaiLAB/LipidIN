preprocessing_RaMS <- function(FN1,ESImode,MS2_filter){
  cat("\033[32mgrabMzmlData reading MS2!\n\033[0m")
  a <- grabMzmlData(FN1,grab_what=c("MS2"))
  a <- a$MS2
  a$kind <- paste(a$rt,a$premz)
  b <- unique(a$kind)
  func <- function(i){
    c <- a[which(a$kind==b[i]),]
    count.1 <- list(
      num=i,
      PrecursorMZ=as.numeric(c[1,2]),
      PrecursorIntensity=1,
      mode=ESImode,
      rt=as.numeric(c[1,1]),
      MS2mz=data.frame(da.temp.mz=c[,3],da.temp.intensity=c[,4])
    )
    aa <- count.1$MS2mz
    aa <- aa[which(aa$da.temp.intensity>MS2_filter*max(count.1$MS2mz[,2])),]
    count.1$MS2mz <- aa
    return(count.1)
  }
  cat("\033[32mBegin MS filtration!\n\033[0m")
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  assign("a",a,envir=env)
  assign("b",b,envir=env)
  assign("ESImode",ESImode,envir=env)
  assign("MS2_filter",MS2_filter,envir=env)
  clusterExport(cl,"env")
  da2 <- parLapply(cl,1:length(b),func)
  stopCluster(cl)
  cat("\033[32mWriting .rda!\n\033[0m")
  neg.list <- do.call(list,da2)
  file3 <- paste(gsub('.mzML','',FN1),'_treated.rda',sep='')
  save(neg.list,file=file3)
  cat("\033[32mDONE!!!!!!\n \033[0m")
}


