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
part2 <- function(FNIII){
  cat("\033[32m",paste('Match files ',FNIII,sep=''),"\n\033[0m")
  load(FNIII)
  sample_mz <- purrr::map(.x=neg.list, .f=function(x){
    return(x$PrecursorMZ)
  }) %>% unlist()
  neg.list <- neg.list[order(sample_mz)]
  neg.list <- lapply(1:length(neg.list), function(i){
    neg.list[[i]][["num"]] <- as.character(neg.list[[i]][["num"]])
    return(neg.list[[i]])
  })
  # NPG-MS2 --------------------------------------------------------------------
  # 构造函数
  # compare <- new(Comparator)
  # # 加载数据
  # # par-1:sample_list , par-2:library_mz , par-3:ms1_ppm , par-4:ms2_ppm , par-5:library_list
  # compare$Load(neg.list,count.data.frame$MZ,ppm1,ppm2,count.list)
  compare$Load2(neg.list,ppm1,ppm2)
  # 比较谱图
  compare$Compare()
  # 输出csv文件
  compare$OutputCsv(paste(gsub('.rda','',FNIII),'_NPG_5ppm.csv',sep=''))
  # result processed --------------------------------------------------------
  cat("\033[32m",paste('Result processed ',FNIII,sep=''),"\n\033[0m")
  process_file(paste(gsub('.rda','',FNIII),'_NPG_5ppm.csv',sep=''))
  da <- read.csv(paste(gsub('.rda','',FNIII),'_NPG_5ppm_processed.csv',sep=''))
  colnames(da) <- c('Peak.num','mz','rt','intensity','title','score.match','score.ratio')
  if(length(grep('CH3COO',da$title))!=0){
    da <- da[-grep('CH3COO',da$title),]
  }
  d.temp <- da
  d.temp <- separate(d.temp, title, c('title', 'other'), sep = '-Fomula-')
  d.temp <- separate(d.temp, other, c('Fomula', 'CCS'), sep = '-CCS-')
  d.temp <- separate(d.temp, CCS, c('CCS', 'Adduct'), sep = '-Add-')
  da2 <- d.temp
  da2$Tmz <- 0
  da2$Tmz <- do.call(c,lapply(1:nrow(da2),function(cm){
    calMS2(da2[cm,]$Fomula)
  }))
  da <- da2[,-4]
  da <- da[!duplicated(da[,c(2,3,5)]),]
  write.csv(da,paste(gsub('.rda','',FNIII),'_NPG_processed.csv',sep=''))
}






