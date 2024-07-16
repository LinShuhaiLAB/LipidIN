MSDIAL.pretreat <- function(da) {
  # 去除unknown\RIKEN ---------------------------------------------------------------
  # da <- pos
  
  # 处理名字 --------------------------------------------------------------------
  d2 <- da$`Metabolite name`
  d2 <- gsub(' O-', '-O ', d2)
  d2 <- gsub(' P-', '-P ', d2)
  d2 <- gsub('-SN1', '', d2)
  d2 <- gsub('\\(methyl\\)', '', d2)
  d2 <- gsub('\\/N-', '_', d2)
  # 根据|划分
  da$total <- gsub('\\|.*', '', d2)
  d2 <- gsub('.*\\|', '', d2)
  d3 <- d2[str_detect(d2, '\\(FA')]# 处理有（FA 12:6）这种的
  # 将（FA和其他分开，统一放最后
  d4 <- gsub('.*\\(FA ', '', d3)
  d4 <- gsub('\\).*', '', d4)
  # 删除（FA xx）
  d3 <- gsub('\\(FA .*\\)', '', d3)
  # 在原本的名字后面加上-FA
  d3 <- gsub(' ', '-FA ', d3)
  d2[str_detect(d2, '\\(FA')] <- d3
  d2 <- gsub('AHexCer \\(O-', 'AHexCer-O ', d2)
  da$meta_name <- d2
  da$subclass <- gsub(' .*', '', d2)
  da$other <- gsub(' ', '', gsub('.* ', '', d2))
  da$other <- gsub('OH)', 'OH', da$other)
  da$other <- gsub('\\)', '_', da$other)
  da$other <- gsub('OH', 'OH)', da$other)
  da$other <- gsub('/', '_', da$other)
  # da$other <- paste("'",da$other,sep='')
  da <-
    da[, c(
      'Alignment ID',
      'Average Rt(min)',
      'Average Mz',
      'Metabolite name',
      'Adduct type',
      'meta_name',
      'total',
      'subclass',
      'other'
    )]
  colnames(da) <-
    c('Peak ID',
      'rt',
      'mz',
      'Title',
      'Adduct',
      'meta_name',
      'total',
      'subclass',
      'Chain')
  # da <- data.frame(da,'Total.C'=0,'Total Uns'=0)
  return(da)
}
MSDIAL.pretreat2 <- function(da) {
  # da <- d2
  da <- data.frame(da, Total.C = 0, 'Total Usa' = 0)
  da.total.chain <- da$Chain
  # 将每条链分开
  count <- matrix(0, 
                  nrow = length(da.total.chain), 
                  ncol = 16)
  colnames(count) <- c(
    'FA1',
    'C1',
    'U1',
    'other1',
    'FA2',
    'C2',
    'U2',
    'other2',
    'FA3',
    'C3',
    'U3',
    'other3',
    'FA4',
    'C4',
    'U4',
    'other4'
  )
  for (i in 1:length(da.total.chain)) {
    sp.chain <- unlist(strsplit(da.total.chain[i], '_'))
    for (j in 1:length(sp.chain)) {
      count[i, j * 4 - 3] <- sp.chain[j]
    }
  }
  count[, 5] <- gsub('\\(', ';\\(', count[, 5])# 针对第二列有部分（OH）前面没有；问题
  # 接下来按照每四列看一下X:Y;ZO
  for (ii in 1:4) {
    for (jj in 1:nrow(count)) {
      ddd <- count[jj, 4 * ii - 3]
      count[jj, 4 * ii - 2] <- unlist(strsplit(ddd, ':'))[1]
      count[jj, 4 * ii - 1] <-
        unlist(strsplit(unlist(strsplit(ddd, ':'))[2], ';'))[1]
      count[jj, 4 * ii] <-
        unlist(strsplit(unlist(strsplit(ddd, ':'))[2], ';'))[2]
    }
  }
  count <- as.data.frame(count)
  count[, c(2, 3, 6, 7, 10, 11, 14, 15)] <-
    apply(count[, c(2, 3, 6, 7, 10, 11, 14, 15)], 2, as.numeric)
  count[is.na(count)] <- 0
  final.count <-
    data.frame(
      'Total.C' = apply(count[, c(2, 6, 10, 14)], 1, sum),
      'Total Uns' = apply(count[, c(3, 7, 11, 15)], 1, sum),
      'other' = paste(count[, 4], count[, 8], count[, 12], count[, 16], sep =
                        '')
    )# 分别求和
  final.count$other <- gsub('0', '', 
                            final.count$other)
  m <- data.frame(
    'subclass' = da$subclass,
    'Total.C' = final.count$Total.C,
    'Total Uns' = paste(final.count$Total.Uns, 
                        final.count$other, sep =
                          ';'),
    'Title' = da$Title,
    'mz' = da$mz,
    'rt' = da$rt,
    'Adduct' = da$Adduct,
    'Chain' = da$Chain
  )
  return(m)
}
da.delect <- function(da) {
  # 依据得分作为剔除概率的剔除
  # da <- d
  # 计算rawmz和rt误差在指定范围内的点
  da$temp.cluster <- paste(rawmzcluster(da$X, 
                                        da$rawmz, 10 ^ (-5)),
                           rtcluster(da$X, 
                                     da$rt, 5 / 60),
                           sep = 'ab')
  da$scale.score <- 0
  k <- unique(da$temp.cluster)
  da <- do.call(rbind, lapply(1:length(k), function(i) {
    dd <- da[which(da$temp.cluster == k[i]),]
    if (length(unique(dd$score.match)) == 1) {
      dd$scale.score <- 1
    }
    if (length(unique(dd$score.match)) != 1) {
      dd$scale.score <-
        (dd$score.match - min(dd$score.match)) / (max(dd$score.match) - min(dd$score.match))
    }
    dn <- sample(1:nrow(dd), size = 1, prob = dd$scale.score)
    dd <- dd[as.numeric(dn),]
    return(dd)
  }))
  da <- da[,-(ncol(da) - 1)]
  return(da)
}
begin.point <- function(d, rep.time) {
  # rep.time <- 100
  # ratio <- 0.4
  # d <- d
  d.mz <- unique(rawmzcluster(d$X, d$rawmz, 10 ^ (-5)))
  if (length(d.mz) < 4) {
    # 点数少于4，报错
    warning("There are too few points to search and fit.")
  }
  if (length(d.mz) >= 4) {
    # 点数足够开始
    n <- max(4, ceiling(0.5 * length(d.mz)))# 抽取数量
    da.temp <- d
    da.temp$label <- rawmzcluster(d$X, d$rawmz, 10 ^ (-5))
    res <- lapply(1:rep.time, function(i) {
      # 用于重复rep.time次
      # 先对mz一样的随机保留一个
      da.temp1 <-
        do.call(rbind, lapply(1:length(d.mz), function(j) {
          dddd <- da.temp[which(da.temp$label == d.mz[j]),]# 提取出mz一样的
          ddddd <-
            sample(1:nrow(dddd),
                   size = 1,
                   prob = dddd$score.match)# 依据概率抽取一个
          return(dddd[ddddd,])
        }))
      da.temp1 <-
        da.temp1[sample(
          1:nrow(da.temp1),
          size = n,
          prob = da.temp1$score.match,
          replace = F
        ),]
      # 计算损失函数值
      a <- list()
      a[[1]] <- da.temp1$X
      a[[2]] <- target.function(da.temp1$rt, da.temp1$rawmz)
      names(a) <- c('num', 'score')
      return(a)
    })
    res.score <- do.call(c, lapply(1:length(res), function(i) {
      res[[i]][["score"]]
    }))
    if (max(res.score) == 0) {
      warning("More than 25% error points, unable to search and fit.")
    }
    else{
      res.score[which(res.score == 0)] <- 9999
      return(res[[which.min(res.score)]])
    }
  }
}
Monotonicity.judgment <- function(x, y) {
  da <- data.frame(x = x, y = y)
  da <- da[order(da$x),]
  error <- -5 * 10 ^ (-6) * mean(da$y)# 采用了5个ppm
  if (!is.unsorted(da$y) && all(diff(da$y) >= error)) {
    I1 <- 1
  }
  else{
    I1 <- 0
  }
  return(I1)
}
Smooth.judgment <- function(x, y) {
  test <-
    try(lm(y ~ poly(x, 2), data.frame(x = x, y = y)), silent = TRUE)
  if (('try-error' %in% class(test)) == FALSE) {
    I2 <- 1
  }
  else{
    I2 <- 0
  }
  return(I2)
}
SSE <- function(x, y) {
  fit <- lm(y ~ poly(x, 2), data.frame(x = x, y = y))
  a <- sqrt(sum(resid(fit) ^ 2))
  b <- mean(y) * length(x)
  return(1 + a / b)
}
target.function <- function(x, y) {
  I1 <- Monotonicity.judgment(x, y)
  I2 <- Smooth.judgment(x, y)
  if (I2 == 1) {
    I3 <- SSE(x, y)
  }
  else{
    I3 <- 0
  }
  return(I1 * I2 * I3)
}
calMS2 <- function(fm){
  # fm <- 'C6H101O2'
  # num.H <- 101
  pattern <- "H\\d+(?=[A-Z])"
  result <- gregexpr(pattern,fm,perl=TRUE)
  num.H <- as.numeric(gsub('H','',regmatches(fm,result)))
  if(num.H<=99){
    re <- calculateMolecularMass(fm)
  }
  if(num.H>99){
    fm2 <- fm
    f1 <- gsub('H.*','',fm2)
    f2 <- gsub('.*H','',fm2)
    f3 <- 'H99'
    f4 <- paste('H',num.H-99,sep='')
    re <- calculateMolecularMass(f1)+calculateMolecularMass(f2)+
      calculateMolecularMass(f3)+calculateMolecularMass(f4)
  }
  return(re)
}
all_search <- function(bp.point, other.point,target.function.tolerance) {
  # target.function.tolerance <- 2
  x.best <- bp.point$rt
  y.best <- bp.point$rawmz
  other.point <- rbind(cbind(other.point, label = 0),
                       cbind(bp.point, label = 'best.begin'))
  x <- other.point$rt
  y <- other.point$rawmz
  label <- rep('check', time = nrow(other.point))
  fit <- lm(y ~ poly(x, 2), data.frame(x = x.best, y = y.best))
  x.in.num <-
    which(other.point$rt <= max(x.best) &
            other.point$rt >= min(x.best))# 在起始点的范围内的点进行判断
  prep <- data.frame(
    x.num = x.in.num,
    x = other.point[x.in.num,]$rt,
    y = other.point[x.in.num,]$rawmz,
    predict(
      fit,
      newdata = data.frame(x = other.point[x.in.num,]$rt),
      interval = "prediction"
    )
  )
  if (length(which(prep$lwr == 'NaN')) > 0) {
    prep[which(prep$lwr == 'NaN'),]$lwr <- prep$fit - 0.1
  }
  if (length(which(prep$upr == 'NaN')) > 0) {
    prep[which(prep$upr == 'NaN'),]$upr <- prep$fit + 0.1
  }
  label[prep[which(prep$y <= prep$upr &
                     prep$y >= prep$lwr), 1]] <-
    'best.begin' # 置信区间内点命名
  mm <-
    data.frame(x = other.point$rt,
               y = other.point$rawmz,
               label = other.point$label)
  ####3.1.2全局搜索往下####
  # 从起始点的最小值开始
  for (xx in 1:nrow(other.point)) {
    x.best <- mm[which(label == 'best.begin'), 'x']
    y.best <- mm[which(label == 'best.begin'), 'y']
    x.min.data <- x.best[which(x.best == sort(x.best)[2])]
    y.min.data <- y.best[which(x.best == sort(x.best)[2])]
    # 计算最低点的切线的正弦值
    min.best.k <- sin(atan(
      2 * as.numeric(coef(fit))[2] * x.min.data +
        as.numeric(coef(fit))[3] * x.min.data
    ))
    # 寻找到符合的象限区间的点
    data.less <-
      data.frame(num = intersect(which(x < min(x.best)), which(y < min(y.best))),
                 x = x[intersect(which(x < min(x.best)), which(y <
                                                                 min(y.best)))],
                 y = y[intersect(which(x < min(x.best)), which(y <
                                                                 min(y.best)))])
    # 对于区间内的点求sin值
    if (nrow(data.less) > 0) {
      data.less <- do.call(rbind, lapply(1:nrow(data.less), function(i) {
        d = (x.min.data - data.less[i,]$x) / sqrt((x.min.data - data.less[i,]$x) ^
                                                    2 + (y.min.data - data.less[i,]$y) ^ 2)
        data.frame(data.less[i,], sin = d)
      }))
      # 进一缩小可行域在切线以下
      data.less <- data.less[which(data.less$sin < min.best.k),]
      # 贪婪计算可行域内每个点加入后对目标函数的改变
      if (nrow(data.less) > 0) {
        data.less <- do.call(rbind, lapply(1:nrow(data.less), function(i) {
          data.frame(data.less[i,], tf = target.function(c(x.best, data.less[i,]$x), c(y.best, data.less[i,]$y)))
        }))
        # 选择目标函数值最小并且满足设置的阈值
        if (data.less[which.min(data.less$tf),]$tf < target.function.tolerance) {
          label[data.less[which.min(data.less$tf),]$num] <- 'best.begin'
          # print(data.less[which.min(data.less$tf),]$num)
          mm <- data.frame(x = x,
                           y = y,
                           label = label)
        }
      }
    }
    if (nrow(data.less) == 0) {
      break
    }
  }
  ####3.1.3全局搜索往上###
  for (xx in 1:nrow(other.point)) {
    x.best <- mm[which(label == 'best.begin'), 'x']
    y.best <- mm[which(label == 'best.begin'), 'y']
    x.min.data <-
      x.best[which(x.best == sort(x.best, decreasing = TRUE)[2])]
    y.min.data <-
      y.best[which(x.best == sort(x.best, decreasing = TRUE)[2])]
    # 计算最低点的切线的正弦值
    min.best.k <- sin(atan(
      2 * as.numeric(coef(fit))[2] * x.min.data ^ 2 +
        as.numeric(coef(fit))[3] * x.min.data
    ))
    # 寻找到符合的象限区间的点
    data.less <-
      data.frame(num = intersect(which(x > max(x.best)), which(y > max(y.best))),
                 x = x[intersect(which(x > max(x.best)), which(y >
                                                                 max(y.best)))],
                 y = y[intersect(which(x > max(x.best)), which(y >
                                                                 max(y.best)))])
    if (nrow(data.less) > 0) {
      data.less <- do.call(rbind, lapply(1:nrow(data.less), function(i) {
        d = (data.less[i,]$x - x.min.data) / sqrt((x.min.data - data.less[i,]$x) ^
                                                    2 + (y.min.data - data.less[i,]$y) ^ 2)
        data.frame(data.less[i,], sin = d)
      }))
      # 进一缩小可行域在切线以下
      data.less <- data.less[which(data.less$sin < min.best.k),]
      # 贪婪计算可行域内每个点加入后对目标函数的改变
      if (nrow(data.less) > 0) {
        data.less <- do.call(rbind, lapply(1:nrow(data.less), function(i) {
          data.frame(data.less[i,], tf = target.function(c(x.best, data.less[i,]$x), c(y.best, data.less[i,]$y)))
        }))
        # 选择目标函数值最小并且满足设置的阈值
        if (data.less[which.min(data.less$tf),]$tf < target.function.tolerance) {
          label[data.less[which.min(data.less$tf),]$num] <- 'best.begin'
          # print(data.less[which.min(data.less$tf),]$num)
          mm <- data.frame(x = x,
                           y = y,
                           label = label)
        }
      }
    }
    if (nrow(data.less) == 0) {
      break
    }
  }
  # 将拟合曲线上的点包进来（只包括这些最优解起始点包含的范围）
  mm1 <- mm[which(mm$label %in% c('best.begin')),]
  fit <- lm(y ~ poly(x, 2), data.frame(x = mm1$x, y = mm1$y))
  prep <- data.frame(
    x = other.point$rt,
    y = other.point$rawmz,
    label = 0,
    predict(
      fit,
      newdata = data.frame(x = other.point$rt),
      interval = "prediction"
    )
  )
  if (length(which(prep$lwr == 'NaN')) > 0) {
    prep[which(prep$lwr == 'NaN'),]$lwr <- prep$fit - 0.1
  }
  if (length(which(prep$upr == 'NaN')) > 0) {
    prep[which(prep$upr == 'NaN'),]$upr <- prep$fit + 0.1
  }
  # label <- abs(prep$y-prep$fit)/prep$y*summary(fit)[["adj.r.squared"]]
  label <- abs(prep$y - prep$fit) / prep$y
  mm <- data.frame(rt = prep$x,
                   rawmz = prep$y,
                   label = label)
  mm <-
    merge(other.point[,-ncol(other.point)], mm, by = c('rt', 'rawmz'))
  return(mm)
}
LCI_nomultithread <- function(FNIII){
  cat("\033[32m",'loading R package and function.',"\n\033[0m")
  packages <- c('this.path','RaMS','Rcpp','tidyverse')
  installed_packages <- packages %in% rownames(installed.packages())
  if(all(installed_packages)){
    print("All required packages are installed.")
  } else{
    print("The following packages are missing:")
    print(packages[!installed_packages])
    
    # Install the missing packages
    install.packages(packages[!installed_packages])
  }
  library(this.path)
  library(RaMS)
  library(Rcpp)
  library(tidyverse)
  sourceCpp(paste(getwd(),'/EQ.cpp',sep=''))
  sourceCpp(paste(getwd(),'/calculateMolecularMass.cpp',sep=''))
  sourceCpp(paste(getwd(),'/rawmzcluster.cpp',sep=''))
  sourceCpp(paste(getwd(),'/removeRowsWithinError.cpp',sep=''))
  sourceCpp(paste(getwd(),'/sortMatrixByRow.cpp',sep=''))
  da <- read.csv(paste(gsub('.rda','',FNIII),'_NPG_processed.csv',sep=''))
  pos <- data.frame('Alignment ID' = da$Peak.num,'Average Rt(min)' = da$rt,
                    'Average Mz' = da$mz,'Metabolite name' = da$title,
                    'Adduct type' = da$Adduct,score.match = da$score.match,
                    score.ratio = da$score.ratio,check.names = F)
  d2 <- data.frame('Peak ID'=pos$`Alignment ID`,rt=pos$`Average Rt(min)`,
                   mz=pos$`Average Mz`,Title=pos$`Metabolite name`,
                   Adduct=pos$`Adduct type`,meta_name=pos$`Metabolite name`,
                   total=pos$`Metabolite name`,
                   subclass=gsub(' .*','',pos$`Metabolite name`),
                   Chain=gsub('.* ','',gsub('_CAH_.*','',pos$`Metabolite name`)))
  pos <- MSDIAL.pretreat2(d2)
  pos$Peak.num <- da$Peak.num
  pos$rawmz <- do.call(rbind,lapply(1:nrow(da),function(cm){
    calMS2(da[cm,]$Fomula)
  }))
  pos$Score <- da$Score
  pos$cluster <- paste(pos$subclass, pos$Total.Uns, pos$Adduct, sep = '')
  pos$score.match <- da$score.match
  pos$score.ratio <- da$score.ratio
  pos$X <- 1:nrow(pos)
  fc <- getwd()
  # 开始PHS 1 -----------------------------------------------------------------
  cat("\033[32m","LCI ECN","\n\033[0m")
  da <- pos
  da <- subset(da,select=-Peak.num)
  # 处理成list
  da.cluster <- unique(da$subclass)
  res <- do.call(list,lapply(1:length(da.cluster),function(i){
    da.temp <- da[which(da$subclass == da.cluster[i]), ]
    da.temp.cluster <- unique(da.temp$Total.Uns)
    res.temp <-
      do.call(list, lapply(1:length(da.temp.cluster), function(j) {
        da.temp[which(da.temp$Total.Uns == da.temp.cluster[j]), ]
      }))
    names(res.temp) <- da.temp.cluster
    return(res.temp)
  }))
  names(res) <- da.cluster
  m <- length(unique(paste(da$subclass,da$Total.Uns,sep='')))
  for_data <- do.call(rbind,lapply(1:length(res),function(ii){
    forfor <- lapply(1:length(res[[ii]]),function(fori){
      akb <- res[[ii]][[fori]]
      if(length(unique(paste(akb$Total.C,akb$Total.Uns)))>=3){
        data.frame(i=ii,j=fori)
      }
    })
    forfor <- do.call(rbind,forfor)
    return(forfor)
  }))
  da2 <- lapply(1:nrow(for_data),function(for_data_i){
    
    CI <- 0.05
    result.count <- NULL
    i <- for_data[for_data_i,1]
    j <- for_data[for_data_i,2]
    d <- res[[i]][[j]]
    d <- da.delect(d)# 去除总链长总不饱和度一样，rawmz一样，rt很接近的
    bp <- begin.point(d,100)
    result.count <- c()
    if (is.list(bp)==T){
      rrrreeesss <- lapply(1:10,function(ijk){
        d <- res[[i]][[j]]
        d <- da.delect(d)# 去除总链长总不饱和度一样，rawmz一样，rt很接近的
        bp <- begin.point(d,100)
        result <- NULL
        if(is.list(bp)==T){
          bp.point <- d[which(d$X %in% (bp[["num"]])),]
          # 删除和bp.point的rawmz一样的点
          other.point <- d[which(removeRowsWithinError(d,bp.point$rawmz,d$rawmz,10^(-6))==FALSE),]
          if (nrow(other.point)>0) {
            # 还有其他备选点
            result <- all_search(bp.point,other.point,target.function.tolerance=2)
          }
          if (nrow(other.point)==0) {
            # 没有其他备选点
            mm1 <- bp.point
            fit <- lm(y~poly(x,2),data.frame(x=mm1$rt,y=mm1$rawmz))
            prep <- data.frame(x=mm1$rt,y=mm1$rawmz,label=0, 
                               predict(fit,newdata=data.frame(x=mm1$rt), 
                                       interval="prediction",level=CI)
            )
            if (length(which(prep$lwr == 'NaN')) > 0) {
              prep[which(prep$lwr == 'NaN'), ]$lwr <- prep$fit - 0.1
            }
            if (length(which(prep$upr == 'NaN')) > 0) {
              prep[which(prep$upr == 'NaN'), ]$upr <- prep$fit + 0.1
            }
            label <- abs(prep$y - prep$fit) / prep$y
            result <- data.frame(rt = prep$x,rawmz=prep$y,label=label)
            result <- merge(bp.point,result,by=c('rt','rawmz'))
          }
          return(result)
        }
      })
      rppp <- do.call(rbind,rrrreeesss)
      if(is.null(rppp)==F){
        rppp <- rppp %>% group_by(X,Title,Adduct) %>% filter(label==min(label))
        result.count <- rbind(result.count,rppp)
      }
    }
    return(result.count)
  })
  result.count <- do.call(rbind,da2)
  if(is.null(result.count)==F){
    result <- result.count %>% group_by(Title) %>% filter(label==min(label))
    result.subclass <- unique(result[which(result$label<=0.1),]$subclass)# 最终鉴定结果中的
    op <- lapply(1:length(result.subclass),function(i) {
      # print(i)
      # i=9;j=1
      da.temp <- result[which(result$label<=0.1&result$subclass==result.subclass[i]),]# 能使用的点
      da.temp.cluster <- unique(da.temp$Total.Uns)
      r2 <- 0
      op <- c()
      # 选择拟合最好的作为中间模型
      for (j in 1:length(da.temp.cluster)) {
        # print(j)
        da.temp1 <- da.temp[which(da.temp$Total.Uns == da.temp.cluster[j]), ]
        if(length(unique(da.temp$Title))>=3){
          model <- lm(y ~ poly(x, 2), data.frame(x = da.temp$rt, y = da.temp$rawmz))# 得到模型
          rs <- summary(model)$adj.r.squared
          if (rs == 'NaN') {
            rs <- 1
          }
          if (r2 < rs) {
            r2 <- rs
            fit <- model
          }
        }
        
      }
      if (r2 > 0.6) {
        list.temp <- res[[result.subclass[i]]]# 得到该亚类的所有点
        list.name <- setdiff(names(list.temp), da.temp.cluster)# 提取所有未处理的不饱和度的名字
        # 对每个还没处理的分布饱和度进行处理
        if (length(list.name) > 0) {
          op <- do.call(rbind, lapply(1:length(list.name), function(k) {
            # print(c(i,j,k))
            dd <- list.temp[[list.name[k]]]
            # print(c(i,j,k))
            if (length(unique(dd$Total.C)) >= 2) {
              # 要求至少两个点
              # 依据得分的随机抽取起始点，得到var最小时的c值
              d.uc <- unique(dd$Total.C)
              c.var.rule <- 9999999999
              # 进行10次模拟
              for (mm in 1:100) {
                dd1 <- do.call(rbind, lapply(1:length(d.uc), function(jj) {
                  dddd <- dd[which(dd$Total.C == d.uc[jj]), ]# 提取出mz一样的
                  ddddd <- sample(1:nrow(dddd), size = 1, prob = dddd$score.match)# 依据概率抽取一个
                  return(dddd[ddddd, ])
                }))
                if(nrow(dd1)>1){
                  dd1.result <- predict(fit, data.frame(x = dd1$rt)) - dd1$rawmz
                }
                if(nrow(dd1)<=1){
                  dd1.result <- c(999.999,999)
                }
                c.var <- var(dd1.result)
                if (c.var < c.var.rule) {
                  c.var.rule <- c.var
                  c.result <- mean(dd1.result)
                }
                if (c.var >= c.var.rule) {
                  c.result <- 99999999999
                }
              }
              if(c.result!=99999999999){
                a <- predict(fit, data.frame(x = dd$rt)) - dd$rawmz - c.result
                label <- abs(a) / dd$rawmz
                op <- data.frame( rt = dd$rt, rawmz = dd$rawmz, X = dd$X, subclass = dd$subclass,
                                  Total.C = dd$Total.C, Total.Uns = dd$Total.Uns, Title = dd$Title,
                                  mz = dd$mz, Adduct = dd$Adduct, Chain = dd$Chain, cluster = dd$cluster,
                                  score.match = dd$score.match, score.ratio = dd$score.ratio, scale.score = 1,
                                  label = label)
                return(op)
              }
            }
          }))
        }
      }
      return(op)
    })
    op <- do.call(rbind,op)
    out.put <- rbind(result.count,op)
    write.csv(out.put,paste(gsub('.rda','',FNIII),'_PHSI.csv',sep=''))
  }
  # 开始PHS 2 -----------------------------------------------------------------
  cat("\033[32m","LCI ESCN+IUP","\n\033[0m")
  # da <- read.csv('N_B-30_treated_PHSI.csv')
  if(file.exists(paste(gsub('.rda','',FNIII),'_PHSI.csv',sep=''))==T){
    da <- read.csv(paste(gsub('.rda','',FNIII),'_PHSI.csv',sep=''))
    # 同一鉴定结果保留得分最小的
    da <- da %>% group_by(Title) %>% filter(label == min(label))
    da <- da[which(da$label < 0.1), ]
    da.SU <- separate(da, Chain, c('FA1', 'FA2', 'FA3', 'FA4'), sep = '_')
    # 构建不饱和度矩阵
    FA.ratio <- data.frame(
      ratio1 = gsub('.*:', '', da.SU$FA1),
      ratio2 = gsub('.*:', '', da.SU$FA2),
      ratio3 = gsub('.*:', '', da.SU$FA3),
      ratio4 = gsub('.*:', '', da.SU$FA4)
    )
    # 调整顺序
    FA.ratio[is.na(FA.ratio)] <- 'ZZZ'
    FA.ratio <- as.matrix(FA.ratio)
    da.SU[, c('FA1', 'FA2', 'FA3', 'FA4')] <- sortMatrixByRow(FA.ratio)
    da.SU$chain.cluster <- paste(da.SU$subclass, da.SU$Total.Uns, da.SU$FA1,
                                 da.SU$FA2, da.SU$FA3, da.SU$FA4, sep = '*')
    da.cluster <- unique(da.SU$chain.cluster)
    res <- do.call(list,lapply(1:length(da.cluster), function(i) {
      da.temp <- da.SU[which(da.SU$chain.cluster == da.cluster[i]), ]
      return(da.temp)
    }))
    names(res) <- da.cluster
    m <- length(res)
    da2 <- lapply(1:length(res),function(for_data_i){
      CI <- 0.05
      i <- for_data_i
      d <- res[[i]]
      d <- da.delect(d)# 去除总链长总不饱和度一样，rawmz一样，rt很接近的
      bp <- begin.point(d[, c('X', 'rt', 'rawmz')], 20)
      result.count <- c()
      if (is.list(bp)==T){
        rppp <- NULL
        ropp <- NULL
        ropp <- lapply(1:10,function(ijk){
          d <- res[[i]]
          d <- da.delect(d)# 去除总链长总不饱和度一样，rawmz一样，rt很接近的
          bp <- begin.point(d[, c('X', 'rt', 'rawmz')], 20)
          if (is.list(bp) == T){
            bp.point <- d[which(d$X %in% (bp[["num"]])), ]
            # 没有其他备选点
            mm1 <- bp.point
            fit <- lm(y ~ poly(x, 2), data.frame(x = mm1$rt, y = mm1$rawmz))
            prep <- data.frame(x = d$rt, y = d$rawmz, label = 0,
                               predict(fit, newdata = data.frame(x = d$rt),
                                       interval = "prediction", level = CI)
            )
            label <- abs(prep$y - prep$fit) / prep$y
            result <- data.frame(rt = prep$x, rawmz = prep$y, label = label)
            result <- merge(d, result, by = c('rt', 'rawmz'))
            return(result)
          }
        })
        rppp <- do.call(rbind,ropp)
        if(is.null(rppp)==F){
          rppp <- rppp %>% group_by(X,Title,Adduct) %>% filter(label.y == min(label.y))
          result.count <- rbind(result.count,rppp)
        }
      }
      return(result.count)
    })
    result.count <- do.call(rbind,da2)
    write.csv(result.count,paste(gsub('.rda','',FNIII),'_PHSII.csv',sep=''))
  }
  if(file.exists(paste(gsub('.rda','',FNIII),'_PHSI.csv',sep=''))==T&
     file.exists(paste(gsub('.rda','',FNIII),'_PHSII.csv',sep=''))==T){
    d1 <- pos
    d2 <- read.csv(paste(gsub('.rda','',FNIII),'_PHSI.csv',sep=''))
    d3 <- read.csv(paste(gsub('.rda','',FNIII),'_PHSII.csv',sep=''))
    d1 <- data.frame(X = d1$X,subclass = d1$subclass,Title = d1$Title,
                     mz = d1$mz,rt = d1$rt,rawmz = d1$rawmz,Adduct = d1$Adduct,
                     peak.num = d1$Peak.num,MS2.score1 = d1$score.match,
                     MS2.score2 = d1$score.ratio)
    d2 <- data.frame(X = d2$X,subclass = d2$subclass,Title = d2$Title,
                     mz = d2$mz,rt = d2$rt,rawmz = d2$rawmz,Adduct = d2$Adduct,
                     cluster1 = d2$cluster,MS2.score1 = d2$score.match,
                     MS2.score2 = d2$score.ratio,rule1.score = d2$label)
    d3 <- data.frame(X = d3$X,subclass = d3$subclass,Title = d3$Title,
                     mz = d3$mz,rt = d3$rt,awmz = d3$rawmz,Adduct = d3$Adduct,
                     cluster2 = paste(d3$cluster, d3$FA1, d3$FA2, d3$FA3, d3$FA4),
                     MS2.score1 = d3$score.match,MS2.score2 = d3$score.ratio,
                     rule1.score = d3$label.x,rule2.score = d3$label.y)
    merged <- merge(d1, d2, by = "X", all = TRUE)
    merged <- merge(merged, d3, by = 'X', all = T)
    merged <- data.frame(X = merged$X,peak.num = merged$peak.num,
                         subclass = merged$subclass.x,Title = merged$Title.x,
                         mz = merged$mz.x,rt = merged$rt.x,rawmz = merged$rawmz.x,
                         Adduct = merged$Adduct.x,MS2.score1 = merged$MS2.score1.x,
                         MS2.score2 = merged$MS2.score2.x,
                         rule1.score = merged$rule1.score.x,
                         rule2.score = merged$rule2.score,
                         cluster1 = merged$cluster1,
                         cluster2 = merged$cluster2)
    Standardization.score <- function(score, threshold) {
      # score <- da.test$rule1.score
      # threshold <- 0.1
      score <- data.frame(input = score, input.score = 0)
      score[is.na(score$input) == T, ]$input.score <-
        0.5#对于没有这个值的认为0.5
      score[which(score$input <= threshold), 2] <-
        (threshold - score[which(score$input <= threshold), ]$input) * 2.5 + 0.75
      score[which(score$input > threshold), 2] <-
        0.25 - (score[which(score$input > threshold), ]$input /
                  max(score[which(score$input >
                                    threshold), ]$input)) * 0.25
      return(score$input.score)
    }
    da.test <- merged
    da.test$rule1.score.Stand <- Standardization.score(da.test$rule1.score, 0.05)
    da.test$rule2.score.Stand <- Standardization.score(da.test$rule2.score, 0.05)
    ####合成第一部分结果####
    da.test$final.score <- da.test$MS2.score1+da.test$MS2.score2+
      da.test$rule1.score.Stand+da.test$rule2.score.Stand
    da.test$label <- as.numeric(1:nrow(da.test))
    da.test <- da.test %>% group_by(mz, rt) %>% mutate(label = label[1])
    da.test <- da.test[!duplicated(da.test),]
    write.csv(da.test,paste(gsub('.rda','',FNIII),'_part1_result.csv',sep=''))
  }
  if(file.exists(paste(gsub('.rda','',FNIII),'_PHSI.csv',sep=''))==F&
     file.exists(paste(gsub('.rda','',FNIII),'_PHSII.csv',sep=''))==F){
    d1 <- pos
    d1 <- data.frame(X = d1$X,subclass = d1$subclass,Title = d1$Title,
                     mz = d1$mz,rt = d1$rt,rawmz = d1$rawmz,Adduct = d1$Adduct,
                     peak.num = d1$Peak.num,MS2.score1 = d1$score.match,
                     MS2.score2 = d1$score.ratio)
    merged <- d1
    merged <- data.frame(X = merged$X,peak.num = merged$peak.num,
                         subclass = merged$subclass,
                         Title = merged$Title,
                         mz = merged$mz,rt = merged$rt,rawmz = merged$rawmz,
                         Adduct = merged$Adduct,MS2.score1 = merged$MS2.score1,
                         MS2.score2 = merged$MS2.score2,
                         rule1.score = 0.5,
                         rule2.score = 0.5,
                         cluster1 = 0.5,
                         cluster2 = 0.5)
    da.test <- merged
    da.test$rule1.score.Stand <- 0.5
    da.test$rule2.score.Stand <- 0.5
    ####合成第一部分结果####
    da.test$final.score <- da.test$MS2.score1+da.test$MS2.score2+
      da.test$rule1.score.Stand+da.test$rule2.score.Stand
    da.test$label <- as.numeric(1:nrow(da.test))
    da.test <- da.test %>% group_by(mz, rt) %>% mutate(label = label[1])
    da.test <- da.test[!duplicated(da.test),]
    write.csv(da.test,paste(gsub('.rda','',FNIII),'_part1_result.csv',sep=''))
  }
  
  
  # 整理文件
  setwd(gsub("\\/[^\\/]*$","",FNIII))
  Fk <- list.files()
  Fk <- Fk[grep(gsub('.rda','',gsub('.*\\/','',FNIII)),Fk)]
  d <- read.csv(Fk[grep('part1_result.csv',Fk)])
  d <- d %>% group_by(Title) %>% top_n(1,final.score)
  d <- d[,c("peak.num","subclass","Title","mz","rt","Adduct","MS2.score1",
            "MS2.score2","rule1.score.Stand","rule2.score.Stand","final.score")]
  file.remove(Fk[-grep('.rda',Fk)])
  d$top_n <- 'Top > 3'
  d1 <- d %>% group_by(peak.num) %>% top_n(3,final.score)
  d[which(d$Title%in%d1$Title),]$top_n <- 'Top 3'
  d1 <- d %>% group_by(peak.num) %>% top_n(2,final.score)
  d[which(d$Title%in%d1$Title),]$top_n <- 'Top 2'
  d1 <- d %>% group_by(peak.num) %>% top_n(1,final.score)
  d[which(d$Title%in%d1$Title),]$top_n <- 'Top 1'
  write.csv(d,gsub('.rda','_final_output.csv',Fk[grep('.rda',Fk)]),row.names=F)
}