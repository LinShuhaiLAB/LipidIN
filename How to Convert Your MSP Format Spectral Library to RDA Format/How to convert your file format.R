# How to convert your file format:
# Modify your parameters here ---------------------------------------------
# install.packages("pbapply")
FN1 <- 'D:/xxxxxxxx'
# The address where the working directory is located.
FN2 <- 'MoNA-export-All_LC-MS-MS_Orbitrap.msp'
# The name of the MSP file for the spectral library.
Meta_name <- 'Name: '
# Metabolite or lipid name prompt.
Formula_label <- 'Formula: '
# Formula prompt.
PrecursorMZ_label <- 'PrecursorMZ: '
# Precursor Ion prompt.
Precursor_type_label <- 'Precursor_type'
# Amalgamation prompt.
Ion_mode_label <- 'Ion_mode: '
# Ionization mode prompt.
Num_Peaks_label <- 'Num Peaks: '
# Prompt for the number of secondary peaks.

# run following -----------------------------------------------------------
setwd(FN1)
library(pbapply)
a <- FN2
ppm1 <- 1
ppm2 <- 1
data <- readLines(a)
data <- c('',data)
data <- c(data,'')
start_lines <- grep("^$",data)
t1 <- pblapply(1:(length(start_lines)-1),function(ii){
  a <- start_lines[ii]
  b <- start_lines[ii+1]
  d <- data[a:b]
  if(length(d)>2){
    Name <- gsub(Meta_name,'',d[grep(Meta_name,d)])
    Formula <- gsub(Formula_label,'',d[grep(Formula_label,d)])
    MZ <- as.numeric(gsub(PrecursorMZ_label,'',d[grep(PrecursorMZ_label,d)]))
    Adduct <- gsub(Precursor_type_label,'',d[grep(Precursor_type_label,d)])
    ESImode <- substr(gsub(Ion_mode_label,'',d[grep(Ion_mode_label,d)]),1,1)
    if(length(grep(Num_Peaks_label,d))==1&length(MZ)>0){
      m2 <- d[(grep(Num_Peaks_label,d)+1):(length(d)-1)]
      if(length(Name)>0&length(Formula)>0&is.na(MZ)==F&
         length(Adduct)>0&length(ESImode)>0&length(m2)>0){
        up <- MZ*(1+ppm1/(10^(6)))
        down <- MZ*(1-ppm1/(10^(6)))
        data.frame(Name=Name,Formula=Formula,MZ=MZ,CCS=0,
                   Adduct=Adduct,rt=0,ESImode=ESImode,up=up,down=down)
      }
    }
  }
})
t1 <- do.call(rbind,t1)
l1 <- pblapply(1:(length(start_lines)-1),function(ii){
  print(ii)
  a <- start_lines[ii]
  b <- start_lines[ii+1]
  d <- data[a:b]
  if(length(d)>2){
    Name <- gsub(Meta_name,'',d[grep(Meta_name,d)])
    Formula <- gsub(Formula_label,'',d[grep(Formula_label,d)])
    MZ <- as.numeric(gsub(PrecursorMZ_label,'',d[grep(PrecursorMZ_label,d)]))
    Adduct <- gsub(Precursor_type_label,'',d[grep(Precursor_type_label,d)])
    ESImode <- substr(gsub(Ion_mode_label,'',d[grep(Ion_mode_label,d)]),1,1)
    if(length(grep(Num_Peaks_label,d))==1&length(MZ)>0){
      m2 <- d[(grep(Num_Peaks_label,d)+1):(length(d)-1)]
      inf <- paste('DB#: ',Name,'-Fomula-',Formula,'-Add-',Adduct,sep='')
      if(length(Name)>0&length(Formula)>0&is.na(MZ)==F&
         length(Adduct)>0&length(ESImode)>0&length(m2)>0){
        ms2.spe <- data.frame(mz=as.numeric(gsub(' .*','',m2)),
                              intensity=as.numeric(gsub('.* ','',m2)),
                              up=0,down=0)
        ms2.spe <- ms2.spe[which(ms2.spe$intensity>=0.1*max(ms2.spe$intensity)),]
        ms2.spe <- ms2.spe[!duplicated(ms2.spe[,1]),]
        ms2.spe <- ms2.spe[order(ms2.spe[,1]),]
        ms2.spe$up <- (1+ppm2/(10^(6)))*ms2.spe$mz
        ms2.spe$down <- (1-ppm2/(10^(6)))*ms2.spe$mz
        list('inf'=inf,'ms2.spe'=ms2.spe)
      }
    }
  }
  
})
l1 <- Filter(Negate(is.null),l1)
t1$num <- c(1:nrow(t1))
t1 <- t1[order(t1$MZ),]
l1 <- l1[t1$num]
count.list=l1
count.data.frame=t1
save(count.list,count.data.frame,file=gsub('.msp','.rda',FN2))







