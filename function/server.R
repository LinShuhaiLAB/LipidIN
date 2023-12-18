server <- function(input, output) {
  options(warn = -1)# 取消warning
  options(shiny.maxRequestSize = 300000 * 1024 ^ 2)  # 设置最大上传文件大小为300000MB
  values <- reactiveValues(ESImode = NULL)
  message1 <-
    reactiveValues(message = 'Choose ESImode and Upload mzML file.')
  ####实时监视控制台####
  update_console <- function() {
    # 保证实时输出message变化的函数
    output$console <- renderPrint({
      cat(message1$message)
    })
  }
  observeEvent(values$a, {
    # 创建一个监视
    update_console()
  })
  autoInvalidate <- reactiveTimer(1000)# 实时更新
  observe({
    autoInvalidate()
    values$a <- values$a + 1
  })
  autoInvalidate <- reactiveTimer(1000)
  observe({
    autoInvalidate()
    values$a <- values$a + 1
  })
  ####选择ESImode####
  observeEvent(input$pn_input, {
    values$ESImode <- input$pn_input
    message1$message <-
      paste(values$ESImode, 'mode is Chosen.', sep = ' ')
  })
  ####选择峰处理####
  step2.file.name <- reactiveValues(file.name = NULL)
  observeEvent(input$file, {
    filename <- input$file$name
    step2.file.name$file.name <- filename
    message1$message <-
      paste(message1$message,
            'Reading the file: ',
            step2.file.name$file.name,
            sep = '\n')
    showNotification(
      "Now you can choose peaks processing or not.",
      type = c("error"),
      duration = 200,
      closeButton = T
    )
  })
  ####无参数读取峰####
  v2 <- reactiveValues(neg.list = NULL)
  observeEvent(input$start, {
    showNotification(
      "Start reading peaks without processing.",
      type = c("message"),
      duration = 200,
      closeButton = T
    )
    neg.list <- Peak.pick(step2.file.name$file.name)
    showNotification("Start pertreat.",
                     type = c("message"),
                     duration = 200)
    neg.list <- Peak.pick.pretreat(neg.list, values$ESImode)
    file3 <-
      paste(gsub('.mzML', '', step2.file.name$file.name),
            'treated.rda',
            sep = '')
    save(neg.list, file = file3)
    v2$neg.list <- neg.list
    message1$message <-
      paste(message1$message, 'Reading peaks Done!', sep = '\n')
    showNotification("Now you can start NPG-MS2.",
                     type = c("error"),
                     duration = 200)
  })
  ####有参数的情况下读取#####
  v1 <-
    reactiveValues(
      smooth_method = NULL,
      halfWindowSize = NULL,
      pickPeaks_method = NULL,
      pickfeature_method = NULL,
      mzerror = NULL
    )
  observeEvent(input$smooth_method, {
    v1$smooth_method <- input$smooth_method
  })
  observeEvent(input$halfWindowSize, {
    v1$halfWindowSize <- input$halfWindowSize
  })
  observeEvent(input$pickPeaks_method, {
    v1$pickPeaks_method <- input$pickPeaks_method
  })
  observeEvent(input$pickfeature_method, {
    v1$pickfeature_method <- input$pickfeature_method
  })
  observeEvent(input$mzerror, {
    v1$mzerror <- input$mzerror
  })
  observeEvent(input$start2, {
    showNotification(
      "Start reading peaks with processing.",
      type = c("message"),
      duration = 200
    )
    message1$message <- paste(
      message1$message,
      paste(
        'smooth method',
        v1$smooth_method,
        ', half Window Size',
        input$halfWindowSize,
        ', pick Peaks method',
        input$pickPeaks_method,
        ', pick feature method',
        input$pickfeature_method,
        ', ppm',
        v1$mzerror,
        sep = ': '
      ),
      sep = '\n'
    )
    neg.table <- Simple.pick(
      filename = step2.file.name$file.name,
      smooth_method = v1$smooth_method,
      halfWindowSize = input$halfWindowSize,
      pickPeaks_method = input$pickPeaks_method,
      method = input$pickfeature_method,
      ppm = v1$mzerror
    )
    showNotification(
      "Start reading peaks processed.",
      type = c("message"),
      duration = 200,
      closeButton = T
    )
    file2 <-
      paste(gsub('.mzML', '', step2.file.name$file.name),
            '_Center.mzML',
            sep = '')
    neg.data <- Peak.pick(file2)
    neg.list <- Peak.pick.pretreat(neg.data, values$ESImode)
    file3 <-
      paste(gsub('.mzML', '', step2.file.name$file.name),
            'treated.rda',
            sep = '')
    save(neg.list, file = file3)
    v2$neg.list <- neg.list
    message1$message <-
      paste(message1$message, 'Reading peaks Done!', sep = '\n')
    showNotification(
      "Now you can start PPG-MS2.",
      type = c("error"),
      duration = 200,
      closeButton = T
    )
  })
  
  ####根据参数读取峰####
  ####NPG-MS2####
  v3 <- reactiveValues(ppm = NULL)
  v51 <-
    reactiveValues(
      rda = NULL,
      s.ratio = NULL,
      NPG = NULL,
      part1 = NULL,
      topn = NULL,
      pos = NULL,
      res1 = NULL,
      res2 = NULL,
      cluster1 = NULL,
      cluster2 = NULL,
      da.result = NULL,
      count.data.frame = NULL,
      count.list = NULL
    )
  observeEvent(input$ppmerror, {
    v3$ppm <- input$ppmerror
  })
  observeEvent(input$library.file, {
    message1$message <- paste(message1$message,
                              paste('m/z tolerance', v3$ppm,
                                    sep = ': '),
                              sep = '\n')
    load(input$library.file$name)
    showNotification("Library loading ...",
                     type = c("message"),
                     duration = 200)
    v51$count.data.frame <- count.data.frame
    v51$count.list <- count.list
    SuperFastCompare(
      sample_List = v2$neg.list ,
      library_mz = count.data.frame$MZ ,
      ppm = v3$ppm ,
      library_list = count.list ,
      output_path =  "NPG_MS2.csv"
    )
    showNotification("NPG-MS2 match done!",
                     type = c("message"),
                     duration = 200)
    showNotification(
      "Now you can start MS2 result processing.",
      type = c("error"),
      duration = 200
    )
    message1$message <-
      paste(message1$message, 'NPG-MS2 match done!', sep = '\n')
  })
  ####结果处理####
  v4 <- reactiveValues(threshold = NULL)
  observeEvent(input$score, {
    v4$threshold <- input$score
  })
  observeEvent(input$NPG_MS2.file, {
    message1$message <- paste(message1$message,
                              paste('MS2 threshold', v4$threshold,
                                    sep = ': '),
                              sep = '\n')
    da <- read.csv(input$NPG_MS2.file$name)
    showNotification("Start data processing.",
                     type = c("message"),
                     duration = 100)
    colnames(da) <-
      c('Peak.num',
        'mz',
        'rt',
        'intensity',
        'title',
        'score.match',
        'score.ratio')
    # 处理正负离子问题
    if (values$ESImode == 'P') {
      da <- da[grep("\\+$", da$title), ]
    }
    if (values$ESImode == 'N') {
      da <- da[-grep("\\+$", da$title), ]
    }
    write.csv(da, 'NPG_MS2.csv')
    v51$NPG <- da
    da <-
      da[which(da$score.match >= v4$threshold |
                 da$score.ratio >= v4$threshold), ]
    output <- da
    d.temp <- output
    d.temp <- separate(d.temp, title, c('t', 'title'), sep = 'DB#: ')
    d.temp <- subset(d.temp, select = -t)
    d.temp <-
      separate(d.temp, title, c('title', 'other'), sep = '-Fomula-')
    d.temp <- separate(d.temp, other, c('Fomula', 'CCS'), sep = '-CCS-')
    d.temp <- separate(d.temp, CCS, c('CCS', 'Adduct'), sep = '-Add-')
    d.temp$Tmz <- 0
    showNotification("Start WM calculating.",
                     type = c("message"),
                     duration = 100)
    d.temp$Tmz <- calculateMolecularMass(d.temp$Fomula)
    da <- d.temp
    pos <- data.frame(
      'Alignment ID' = da$Peak.num,
      'Average Rt(min)' = da$rt,
      'Average Mz' = da$mz,
      'Metabolite name' = da$title,
      'Adduct type' = da$Adduct,
      score.match = da$score.match,
      score.ratio = da$score.ratio,
      check.names = F
    )
    # 处理不饱和度、分链信息 ------------------------------------------------------------
    showNotification(
      "Start Carbon chains processing.",
      type = c("message"),
      duration = 100
    )
    d2 <- MSDIAL.pretreat(pos)
    pos <- MSDIAL.pretreat2(d2)
    pos$Peak.num <- da$Peak.num
    pos$rawmz <- da$Tmz
    pos$Score <- da$Score
    pos$cluster <-
      paste(pos$subclass, pos$Total.Uns, pos$Adduct, sep = '')
    pos$score.match <- da$score.match
    pos$score.ratio <- da$score.ratio
    write.csv(pos, 'NPG_MS2_processed.csv')
    message1$message <-
      message1$message <-
      paste(message1$message, 'Standardization Done!', sep = '\n')
    showNotification("Now you can start RT based PHS.",
                     type = c("error"),
                     duration = 200)
  })
  ####RT based PHS####
  vpart1 <- reactiveValues(processed = NULL,
                           phs1 = NULL,
                           phs2 = NULL)
  observeEvent(input$processed.file, {
    da <- read.csv(input$processed.file$name)
    vpart1$processed <- da
    da <- subset(da, select = -Peak.num)
    # 处理成list
    showNotification("Start data processing.",
                     type = c("message"),
                     duration = 100)
    da.cluster <- unique(da$subclass)
    res <- do.call(list, lapply(1:length(da.cluster), function(i) {
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
    message1$message <-
      paste(message1$message, 'Converting format done!', sep = '\n')
    # 每个子子list中进行查找
    showNotification("Start PHS I.",
                     type = c("message"),
                     duration = 100)
    m <- length(unique(paste(da$subclass, da$Total.Uns, sep = '')))
    withProgress(message = 'Processing...', value = 0, {
      result.count <- c()
      for (i in 1:length(res)) {
        for (j in 1:length(res[[i]])) {
          # 在相同的亚类中
          incProgress(1 / m)
          d <- res[[i]][[j]]
          d <- da.delect(d)# 去除总链长总不饱和度一样，rawmz一样，rt很接近的
          bp <- begin.point(d, 100)
          if (is.list(bp) == T) {
            bp.point <- d[which(d$X %in% (bp[["num"]])), ]
            # 删除和bp.point的rawmz一样的点
            other.point <-
              d[which(removeRowsWithinError(d, bp.point$rawmz, d$rawmz, 10 ^ (-6)) ==
                        FALSE), ]
            if (nrow(other.point) > 0) {
              # 还有其他备选点
              result <-
                all_search(bp.point,
                           other.point,
                           target.function.tolerance = 2)
            }
            if (nrow(other.point) == 0) {
              # 没有其他备选点
              mm1 <- bp.point
              fit <-
                lm(y ~ poly(x, 2), data.frame(x = mm1$rt, y = mm1$rawmz))
              prep <- data.frame(
                x = mm1$rt,
                y = mm1$rawmz,
                label = 0,
                predict(
                  fit,
                  newdata = data.frame(x = mm1$rt),
                  interval = "prediction"
                )
              )
              if (length(which(prep$lwr == 'NaN')) > 0) {
                prep[which(prep$lwr == 'NaN'), ]$lwr <- prep$fit - 0.1
              }
              if (length(which(prep$upr == 'NaN')) > 0) {
                prep[which(prep$upr == 'NaN'), ]$upr <- prep$fit + 0.1
              }
              # label <- abs(prep$y-prep$fit)/prep$y*summary(fit)[["adj.r.squared"]]
              label <- abs(prep$y - prep$fit) / prep$y
              result <-
                data.frame(rt = prep$x,
                           rawmz = prep$y,
                           label = label)
              result <- merge(bp.point, result, by = c('rt', 'rawmz'))
            }
            result.count <- rbind(result.count, result)
          }
        }
      }
      m
    })# 进度条
    # 进度条
    message1$message <-
      paste(message1$message, 'PHS I done!', sep = '\n')
    # 对于数量不足的进行拟合迁移
    result <-
      result.count %>% group_by(Title) %>% filter(label == min(label))
    result.subclass <-
      unique(result[which(result$label <= 0.1), ]$subclass)# 最终鉴定结果中的
    op <- lapply(1:length(result.subclass), function(i) {
      da.temp <-
        result[which(result$label <= 0.1 &
                       result$subclass == result.subclass[i]), ]# 能使用的点
      da.temp.cluster <- unique(da.temp$Total.Uns)
      r2 <- 0
      op <- c()
      # 选择拟合最好的作为中间模型
      for (j in 1:length(da.temp.cluster)) {
        # print(c(i,j))
        # i <- 9
        # j <- 7
        da.temp1 <-
          da.temp[which(da.temp$Total.Uns == da.temp.cluster[j]), ]
        model <-
          lm(y ~ poly(x, 2),
             data.frame(x = da.temp$rt, y = da.temp$rawmz))# 得到模型
        rs <- summary(model)$adj.r.squared
        if (rs == 'NaN') {
          rs <- 1
        }
        if (r2 < rs) {
          r2 <- rs
          fit <- model
        }
      }
      if (r2 > 0) {
        list.temp <- res[[result.subclass[i]]]# 得到该亚类的所有点
        list.name <-
          setdiff(names(list.temp), da.temp.cluster)# 提取所有未处理的不饱和度的名字
        # 对每个还没处理的分布饱和度进行处理
        if (length(list.name) > 0) {
          op <- do.call(rbind, lapply(1:length(list.name), function(k) {
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
                  ddddd <-
                    sample(1:nrow(dddd),
                           size = 1,
                           prob = dddd$score.match)# 依据概率抽取一个
                  return(dddd[ddddd, ])
                }))
                dd1.result <-
                  predict(fit, data.frame(x = dd1$rt)) - dd1$rawmz
                c.var <- var(dd1.result)
                if (c.var < c.var.rule) {
                  c.var.rule <- c.var
                  c.result <- mean(dd1.result)
                }
              }
              a <-
                predict(fit, data.frame(x = dd$rt)) - dd$rawmz - c.result
              label <- abs(a) / dd$rawmz
              op <-
                data.frame(
                  rt = dd$rt,
                  rawmz = dd$rawmz,
                  X = dd$X,
                  subclass = dd$subclass,
                  Total.C = dd$Total.C,
                  Total.Uns = dd$Total.Uns,
                  Title = dd$Title,
                  mz = dd$mz,
                  Adduct = dd$Adduct,
                  Chain = dd$Chain,
                  cluster = dd$cluster,
                  score.match = dd$score.match,
                  score.ratio = dd$score.ratio,
                  scale.score = 1,
                  label = label
                )
              return(op)
            }
          }))
        }
      }
      return(op)
    })
    op <- do.call(rbind, op)
    out.put <- rbind(result.count, op)
    vpart1$phs1 <- out.put
    write.csv(out.put, 'PHS_I.csv')
    message1$message <-
      paste(message1$message, 'Migration lookup done!', sep = '\n')
    showNotification("Now you can start PSH II.",
                     type = c("error"),
                     duration = 200)
  })
  ####PHS II 分链####
  observeEvent(input$PHS_I.file, {
    da <- read.csv(input$PHS_I.file$name)
    # 同一鉴定结果保留得分最小的 -------------------------------------------------------------
    da <- da %>% group_by(Title) %>% filter(label == min(label))
    da <- da[which(da$label < 0.1), ]
    showNotification("Start data processing.",
                     type = c("message"),
                     duration = 100)
    # 获得分链矩阵 ------------------------------------------------------------------
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
    da.SU[, c('FA1', 'FA2', 'FA3', 'FA4')] <-
      sortMatrixByRow(FA.ratio)
    da.SU$chain.cluster <-
      paste(
        da.SU$subclass,
        da.SU$Total.Uns,
        da.SU$FA1,
        da.SU$FA2,
        da.SU$FA3,
        da.SU$FA4,
        sep = '*'
      )
    da.cluster <- unique(da.SU$chain.cluster)
    res <- do.call(list, lapply(1:length(da.cluster), function(i) {
      da.temp <- da.SU[which(da.SU$chain.cluster == da.cluster[i]), ]
      return(da.temp)
    }))
    names(res) <- da.cluster
    message1$message <-
      paste(message1$message, 'FA processing done!', sep = '\n')
    showNotification("Start PHS II.",
                     type = c("message"),
                     duration = 100)
    # 进行规则判断和得分 ---------------------------------------------------------------
    m <- length(res)
    withProgress(message = 'Processing...', value = 0, {
      result.count <- c()
      for (i in 1:length(res)) {
        incProgress(1 / m)
        d <- res[[i]]
        d <- da.delect(d)# 去除总链长总不饱和度一样，rawmz一样，rt很接近的
        bp <- begin.point(d[, c('X', 'rt', 'rawmz')], 20)
        if (is.list(bp) == T) {
          bp.point <- d[which(d$X %in% (bp[["num"]])), ]
          # 没有其他备选点
          mm1 <- bp.point
          fit <- lm(y ~ poly(x, 2), data.frame(x = mm1$rt, y = mm1$rawmz))
          prep <- data.frame(
            x = d$rt,
            y = d$rawmz,
            label = 0,
            predict(
              fit,
              newdata = data.frame(x = d$rt),
              interval = "prediction",
              level = 0.99
            )
          )
          label <- abs(prep$y - prep$fit) / prep$y
          result <- data.frame(rt = prep$x,
                               rawmz = prep$y,
                               label = label)
          result <- merge(d, result, by = c('rt', 'rawmz'))
          result.count <- rbind(result.count, result)
        }
      }
      m
    })
    vpart1$phs2 <- result.count
    write.csv(result.count, 'PHS_II.csv')
    message1$message <-
      paste(message1$message, 'PHS II done!', sep = '\n')
    showNotification(
      "Integrate the scores for the first part.",
      type = c("message"),
      duration = 200
    )
    d1 <- vpart1$processed
    d2 <- vpart1$phs1
    d3 <- vpart1$phs2
    d1 <-
      data.frame(
        X = d1$X,
        subclass = d1$subclass,
        Title = d1$Title,
        mz = d1$mz,
        rt = d1$rt,
        rawmz = d1$rawmz,
        Aadduct = d1$Adduct,
        peak.num = d1$Peak.num,
        MS2.score1 = d1$score.match,
        MS2.score2 = d1$score.ratio
      )
    d2 <-
      data.frame(
        X = d2$X,
        subclass = d2$subclass,
        Title = d2$Title,
        mz = d2$mz,
        rt = d2$rt,
        rawmz = d2$rawmz,
        Aadduct = d2$Adduct,
        cluster1 = d2$cluster,
        MS2.score1 = d2$score.match,
        MS2.score2 = d2$score.ratio,
        rule1.score = d2$label
      )
    d3 <-
      data.frame(
        X = d3$X,
        subclass = d3$subclass,
        Title = d3$Title,
        mz = d3$mz,
        rt = d3$rt,
        rawmz = d3$rawmz,
        Aadduct = d3$Adduct,
        cluster2 = paste(d3$cluster, d3$FA1, d3$FA2, d3$FA3, d3$FA4),
        MS2.score1 = d3$score.match,
        MS2.score2 = d3$score.ratio,
        rule1.score = d3$label.x,
        rule2.score = d3$label.y
      )
    merged <- merge(d1, d2, by = "X", all = TRUE)
    merged <- merge(merged, d3, by = 'X', all = T)
    merged <-
      data.frame(
        X = merged$X,
        peak.num = merged$peak.num,
        subclass = merged$subclass.x,
        Title = merged$Title.x,
        mz = merged$mz.x,
        rt = merged$rt.x,
        rawmz = merged$rawmz.x,
        Adduct = merged$Aadduct.x,
        MS2.score1 = merged$MS2.score1.x,
        MS2.score2 = merged$MS2.score2.x,
        rule1.score = merged$rule1.score.x,
        rule2.score = merged$rule2.score,
        cluster1 = merged$cluster1,
        cluster2 = merged$cluster2
      )
    showNotification(
      "Start score standardization.",
      type = c("message"),
      duration = 100
    )
    ####标准化得分####
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
    da.test$rule1.score.Stand <-
      Standardization.score(da.test$rule1.score, 0.1)
    da.test$rule2.score.Stand <-
      Standardization.score(da.test$rule2.score, 0.1)
    ####合成第一部分结果####
    
    da.test$final.score <-
      apply(da.test[, c("MS2.score1",
                        "MS2.score2",
                        "rule1.score.Stand",
                        "rule2.score.Stand")], 1, sum)
    da.test$label <- as.numeric(1:nrow(da.test))
    da.test <-
      da.test %>% group_by(mz, rt) %>% mutate(label = label[1])
    v51$part1 <- da.test
    write.csv(da.test, 'part1_result.csv')
    showNotification("Now you can start PG-MS2 match.",
                     type = c("error"),
                     duration = 100)
  })
  
  ####PG-MS2 match####
  v6 <- reactiveValues(strange = NULL)
  observeEvent(input$peak.rda.data, {
    showNotification("Start loading rda file.",
                     type = c("message"),
                     duration = 100)
    load(input$peak.rda.data$name)
    v51$rda <- neg.list
    showNotification("Loading rda file finished.",
                     type = c("message"),
                     duration = 100)
  })
  observeEvent(input$npgms2.data, {
    showNotification(
      "Start part1 result processing.",
      type = c("message"),
      duration = 100
    )
    da <- read.csv(input$npgms2.data$name, row.names = 1)
    colnames(da) <-
      c('Peak.num',
        'mz',
        'rt',
        'intensity',
        'title',
        'score.match',
        'score.ratio')
    da.tobe <- da[which(da$score.match > 0 & da$score.match < 0.6), ]
    output <- da.tobe
    d.temp <- output
    d.temp <- separate(d.temp, title, c('t', 'title'), sep = 'DB#: ')
    d.temp <- subset(d.temp, select = -t)
    d.temp <-
      separate(d.temp, title, c('title', 'other'), sep = '-Fomula-')
    d.temp <- separate(d.temp, other, c('Fomula', 'CCS'), sep = '-CCS-')
    d.temp <- separate(d.temp, CCS, c('CCS', 'Adduct'), sep = '-Add-')
    d.temp$Tmz <- 0
    d.temp$Tmz <- calculateMolecularMass(d.temp$Fomula)
    da <- d.temp
    pos <- data.frame(
      'Alignment ID' = da$Peak.num,
      'Average Rt(min)' = da$rt,
      'Average Mz' = da$mz,
      'Metabolite name' = da$title,
      'Adduct type' = da$Adduct,
      score.match = da$score.match,
      score.ratio = da$score.ratio,
      check.names = F
    )
    ####处理不饱和度、分链信息
    d2 <- MSDIAL.pretreat(pos)
    pos <- MSDIAL.pretreat2(d2)
    pos$peak.num <- da$Peak.num
    pos$rawmz <- da$Tmz
    pos$Score <- da$Score
    pos$cluster1 <- paste(pos$subclass, pos$Total.Uns, sep = '')
    pos$score.match <- da$score.match
    pos$score.ratio <- da$score.ratio
    da.SU <- separate(pos, Chain, c('FA1', 'FA2', 'FA3', 'FA4'), sep = '_')
    ####构建不饱和度矩阵
    FA.ratio <- data.frame(
      ratio1 = gsub('.*:', '', da.SU$FA1),
      ratio2 = gsub('.*:', '', da.SU$FA2),
      ratio3 = gsub('.*:', '', da.SU$FA3),
      ratio4 = gsub('.*:', '', da.SU$FA4)
    )
    FA.ratio[is.na(FA.ratio)] <- 'ZZZ'
    FA.ratio <- as.matrix(FA.ratio)
    da.SU[, c('FA1', 'FA2', 'FA3', 'FA4')] <-
      sortMatrixByRow(FA.ratio)
    pos$cluster2 <-
      paste(
        paste(da.SU$subclass, da.SU$Total.Uns, sep = ''),
        da.SU$FA1,
        da.SU$FA2,
        da.SU$FA3,
        da.SU$FA4,
        sep = ' '
      )
    pos$cluster2 <- gsub('ZZZ', 'NA', pos$cluster2)
    write.csv(pos, 'Low_score_processed.csv')
    rm(d.temp, d2, da, da.tobe, output, da.SU)
    v51$pos <- pos
  })
  observeEvent(input$part1.data, {
    showNotification(
      "Start part1 result processing.",
      type = c("message"),
      duration = 100
    )
    da.result <- read.csv(input$part1.data$name)
    da.result <-
      da.result[which(is.na(da.result$rule1.score) == F |
                        is.na(da.result$rule2.score) == F), ]
    da.result <-
      da.result[which(da.result$rule1.score.Stand > 0.5), ]# 筛选出优质点
    da.result$cluster1 <- gsub('\\[.*', '', da.result$cluster1)
    cluster1 <- unique(da.result$cluster1)
    da.result$cluster2 <- gsub("\\[.*?\\s", " ", da.result$cluster2)
    cluster12 <- unique(da.result$cluster2)
    showNotification(
      "Start part1-rules1 result processing.",
      type = c("message"),
      duration = 100
    )
    res1 <- lapply(1:length(cluster1), function(ii) {
      da.temp <- da.result[which(da.result$cluster1 == cluster1[ii]), ]
      if (length(unique(da.temp$mz)) >= 4) {
        # 规则一的拟合
        aa <-
          lm(y ~ poly(x, 2),
             data.frame(x = da.temp$rt, y = da.temp$rawmz))
        return(aa)
      }
    })
    names(res1) <- cluster1
    res1 <- res1[!sapply(res1, is.null)]
    cluster1 <- names(res1)
    showNotification(
      "Start part1-rules2 result processing finished.",
      type = c("message"),
      duration = 100
    )
    cluster2 <- do.call(c, lapply(1:length(cluster1), function(i) {
      da.temp <- da.result[which(da.result$cluster1 == cluster1[i]), ]
      if (length(unique(da.temp$mz)) >= 4) {
        # 规则一的拟合
        cluster2 <- na.omit(unique(da.temp$cluster2))
        if (length(cluster2) > 0) {
          # 规则二符合
          bb.num <- do.call(c, lapply(1:length(cluster2), function(j) {
            da.temp2 <- da.temp[which(da.temp$cluster2 == cluster2[j]), ]
            if (length(unique(da.temp2$rawmz)) >= 4) {
              j
            }
          }))
          cluster2[bb.num]
        }
      }
    }))# 找到所有可以的类
    res2 <- do.call(list, lapply(1:length(cluster2), function(j) {
      da.temp2 <- da.result[which(da.result$cluster2 == cluster2[j]), ]
      if (length(unique(da.temp2$rawmz)) >= 4) {
        bbb <- lm(y ~ poly(x, 2),
                  data.frame(x = da.temp2$rt, y = da.temp2$rawmz))
        bbb
      }
    }))
    names(res2) <- cluster2
    showNotification(
      "Part1 result processing finished.",
      type = c("message"),
      duration = 100
    )
    v51$res2 <- res2
    v51$res1 <- res1
    v51$cluster1 <- cluster1
    v51$cluster2 <- cluster2
    v51$da.result <- da.result
  })
  observeEvent(input$strength.ratio, {
    v51$s.ratio <- input$strength.ratio
  })
  observeEvent(input$topn, {
    v51$topn <- input$topn
  })
  observeEvent(input$start.part2, {
    showNotification("Start PG-MS2 match.",
                     type = c("message"),
                     duration = 100)
    pos <- v51$pos
    res1 <- v51$res1
    res2 <- v51$res2
    neg.list <- v51$rda
    cluster1 <- v51$cluster1
    cluster2 <- v51$cluster2
    da.result <- v51$da.result
    strength.ratio <- v51$s.ratio
    topn <- v51$topn
    count.list <- v51$count.list
    count.data.frame <- v51$count.data.frame
    # 进度条
    m <- nrow(pos)
    withProgress(message = 'MS2 spectral enhancing...', value = 0, {
      rrr <- lapply(1:nrow(pos), function(i) {
        incProgress(1 / m)
        if (pos[i, ]$cluster1 %in% cluster1) {
          # 如果他在类别里+如果有他的拟合方程
          k <-
            predict(res1[[pos[i, ]$cluster1]], data.frame(x = pos[i, ]$rt), interval =
                      "prediction")
          if (pos[i, ]$rawmz >= k[2] & pos[i, ]$rawmz <= k[3]) {
            # 在95%CI内
            neg.list.new <-
              neg.list[[paste('F1.S', sprintf("%04d", pos[i, ]$peak.num), sep = '')]]
            neg.list.new$score1 <- abs(k[1] - pos[i, ]$rawmz) / k[1]
            expand.score2 <- NA
            if (pos[i, ]$cluster2 %in% cluster2) {
              k <-
                predict(res2[[pos[i, ]$cluster2]], data.frame(x = pos[i, ]$rt), interval =
                          "prediction")
              if (pos[i, ]$rawmz >= k[2] & pos[i, ]$rawmz <= k[3]) {
                expand.score2 <- abs(k[1] - pos[i, ]$rawmz) / k[1]
              }
            }
            neg.list.new$score2 <- expand.score2
            neg.list.new$stength <- 0
            k <-
              da.result[which(da.result$cluster1 == pos[i, ]$cluster1), ]# 求所有同亚类+总不饱和度的点
            # 根据rawmz把差谱图的点塞进去，将其分为三个数据框
            # 数据框1，距离pos[i,]最近的一左一右的点
            nearest <-
              na.omit(rbind(
                data.frame(k[max(which(k$rawmz < pos[i, ]$rawmz)), ], kind = 'left'),
                data.frame(k[min(which(k$rawmz > pos[i, ]$rawmz)), ], kind =
                             'right')
              ))
            if (nrow(nearest) == 2 &
                length(unique(nearest$rawmz)) == 2) {
              # 既有左端点、又有右端点
              neg.list.new$stength <- 1
              # print(i)
              k21 <-
                neg.list[[paste('F1.S',
                                sprintf("%04d", nearest[1, ]$peak.num),
                                sep = '')]]$MS2mz# 左边谱图
              k22 <-
                neg.list[[paste('F1.S',
                                sprintf("%04d", nearest[2, ]$peak.num),
                                sep = '')]]$MS2mz# 右边谱图
              # 找到两个图共有的峰
              k21 <- data.frame(
                k21$da.temp.mz,
                up = (1 + 10 ^ (-3)) * k21$da.temp.mz,
                down = (1 - 10 ^ (-3)) * k21$da.temp.mz
              )# 计算一个谱图的上下限
              k20 <-
                neg.list[[paste('F1.S', sprintf("%04d", pos[i, ]$peak.num), sep = '')]]$MS2mz# 差谱图
              k2_12 <-
                do.call(rbind, lapply(1:nrow(k21), function(jj) {
                  # 找到21和22共有的
                  k22[which(k22$da.temp.mz <= k21[jj, ]$up &
                              k22$da.temp.mz >= k21[jj, ]$down), ]
                }))
              k2_12$up <- (1 + 10 ^ (-3)) * k2_12$da.temp.mz
              k2_12$down <- (1 - 10 ^ (-3)) * k2_12$da.temp.mz
              k2_12 <-
                do.call(rbind, lapply(1:nrow(k2_12), function(jj) {
                  # 排除和k20共有的
                  mm <-
                    which(k20$da.temp.mz <= k2_12[jj, ]$up &
                            k20$da.temp.mz >= k2_12[jj, ]$down)
                  if (length(mm) == 0) {
                    k2_12[jj, ]
                  }
                  else{
                    k2_12[jj, 1] <- k20[mm[which.max(k20[mm, ]$da.temp.intensity)], 1]
                    k2_12[jj, 2] <-
                      2 * (k20[which.max(k20[mm, ]$da.temp.intensity), 2] + strength.ratio *
                             k2_12[jj, 2]) / strength.ratio
                    k2_12[jj, ]
                  }
                }))
              # 对找到的峰进行加强
              if (nrow(k2_12) > 0) {
                # 说明有找到峰
                i1 <- 2 * 1.007825032 + 12
                index_left <-
                  1 / (abs(nearest[1, ]$rawmz - pos[i, ]$rawmz) / i1)
                index_right <-
                  1 / (abs(nearest[2, ]$rawmz - pos[i, ]$rawmz) / i1)
                strength.matrix <- data.frame(
                  da.temp.mz = k2_12[, 1],
                  da.temp.intensity = strength.ratio *
                    k2_12[, 2] * mean(index_left + index_right)
                )
                neg.list.new$MS2mz <-
                  rbind(neg.list.new$MS2mz, strength.matrix)
                neg.list.new$MS2mz <-
                  neg.list.new$MS2mz[order(neg.list.new$MS2mz$da.temp.mz), ]# 进行排序
                neg.list.new$MS2mz <- neg.list.new$MS2mz %>%
                  group_by(da.temp.mz) %>%
                  slice_max(da.temp.intensity)
              }
            }
            return(neg.list.new)
          }
        }
      })
      m
    })# 进度条
    rrr <- rrr[!sapply(rrr, is.null)]
    names(rrr) <- do.call(c, lapply(1:length(rrr), function(i) {
      rrr[[i]][["num"]]
    }))
    showNotification(
      "MS2 spectral enhancement finished.",
      type = c("message"),
      duration = 100
    )
    strange <- list()
    for (i in 1:length(rrr)) {
      if (rrr[[i]][["stength"]] == 1) {
        strange <- append(strange, list(rrr[[i]]))
      }
    }
    showNotification("Start PG-MS2 match.",
                     type = c("message"),
                     duration = 100)
    k <-
      process_numbers(do.call(c, lapply(1:length(strange), function(i) {
        strange[[i]][["num"]]
      })))
    names(strange) <- k
    strange <- lapply(1:length(strange), function(i) {
      strange[[i]][["num"]] <- k[i]
      strange[[i]]
    })
    v6$strange <- strange
    SuperFastCompare(
      sample_List = strange,
      library_mz = count.data.frame$MZ ,
      ppm = v3$ppm ,
      library_list = count.list ,
      output_path =  "PG_MS2.csv"
    )
    showNotification(
      "Now you can upload PG-MS2 result, to get final output.",
      type = c("error"),
      duration = 100
    )
  })
  ####整理结果####
  observeEvent(input$pg.data, {
    showNotification(
      "Start final output processing.",
      type = c("message"),
      duration = 100
    )
    da <- read.csv(input$pg.data$name)
    da <- da[!duplicated(da[, -1]), ]
    colnames(da) <-
      c('Peak.num',
        'mz',
        'rt',
        'intensity',
        'title',
        'score.match',
        'score.ratio')
    da.tobe <-
      da[which(da$score.match >= v4$threshold |
                 da$score.ratio >= v4$threshold), ]
    output <- da.tobe
    d.temp <- output
    d.temp <- separate(d.temp, title, c('t', 'title'), sep = 'DB#: ')
    d.temp <- subset(d.temp, select = -t)
    d.temp <-
      separate(d.temp, title, c('title', 'other'), sep = '-Fomula-')
    d.temp <- separate(d.temp, other, c('Fomula', 'CCS'), sep = '-CCS-')
    d.temp <- separate(d.temp, CCS, c('CCS', 'Adduct'), sep = '-Add-')
    d.temp$Tmz <- 0
    d.temp$Tmz <- calculateMolecularMass(d.temp$Fomula)
    da <- d.temp
    pos <- data.frame(
      'Alignment ID' = da$Peak.num,
      'Average Rt(min)' = da$rt,
      'Average Mz' = da$mz,
      'Metabolite name' = da$title,
      'Adduct type' = da$Adduct,
      score.match = da$score.match,
      score.ratio = da$score.ratio,
      check.names = F
    )
    d2 <- MSDIAL.pretreat(pos)
    pos <- MSDIAL.pretreat2(d2)
    pos$peak.num <- da$Peak.num
    pos$rawmz <- da$Tmz
    pos$Score <- da$Score
    pos$cluster1 <- paste(pos$subclass, pos$Total.Uns, sep = '')
    pos$score.match <- da$score.match
    pos$score.ratio <- da$score.ratio
    da.SU <- separate(pos, Chain, c('FA1', 'FA2', 'FA3', 'FA4'), sep = '_')
    ####构建不饱和度矩阵
    FA.ratio <- data.frame(
      ratio1 = gsub('.*:', '', da.SU$FA1),
      ratio2 = gsub('.*:', '', da.SU$FA2),
      ratio3 = gsub('.*:', '', da.SU$FA3),
      ratio4 = gsub('.*:', '', da.SU$FA4)
    )
    FA.ratio[is.na(FA.ratio)] <- 'ZZZ'
    FA.ratio <- as.matrix(FA.ratio)
    da.SU[, c('FA1', 'FA2', 'FA3', 'FA4')] <-
      sortMatrixByRow(FA.ratio)
    pos$cluster2 <-
      paste(
        paste(da.SU$subclass, da.SU$Total.Uns, sep = ''),
        da.SU$FA1,
        da.SU$FA2,
        da.SU$FA3,
        da.SU$FA4,
        sep = ' '
      )
    pos$cluster2 <- gsub('ZZZ', 'NA', pos$cluster2)
    sd <-
      data.frame(
        peak.num = gsub('-.*', '', da.tobe$Peak.num),
        pc = da.tobe$Peak.num,
        subclass = pos$subclass,
        Title = pos$Title,
        mz = pos$mz,
        rt = pos$rt,
        rawmz = pos$rawmz,
        Adduct = pos$Adduct,
        MS2.score1 = pos$score.match,
        MS2.score2 = pos$score.ratio,
        rule1.score.Stand = 0,
        rule2.score.Stand = 0,
        final.score = 0
      )
    # 挑选出2倍topn的数量
    sd <-
      sd %>% group_by(peak.num) %>% top_n(v51$topn * 2, (MS2.score1 + MS2.score2))
    Standardization.score <- function(score, threshold) {
      # score <- sd$rule1.score.Stand
      # threshold <- 0.1
      score <-
        data.frame(input = score, input.score = 0.5)# 对于没有这个值的认为0.5
      # score[is.na(score$input)==T,]$input.score <- 0.5#对于没有这个值的认为0.5
      score[which(score$input <= threshold), 2] <-
        (threshold - score[which(score$input <= threshold), ]$input) * 2.5 + 0.75
      score[which(score$input > threshold), 2] <-
        0.25 - (score[which(score$input > threshold), ]$input /
                  max(score[which(score$input >
                                    threshold), ]$input)) * 0.25
      return(score$input.score)
    }
    # 将规则一、规则二得分在list中找回来
    strange <- v6$strange
    ss <- do.call(rbind, lapply(1:length(strange), function(i) {
      data.frame(pc = strange[[i]][["num"]],
                 s1 = strange[[i]][["score1"]],
                 s2 = strange[[i]][["score2"]])
    }))
    merged <- merge(sd, ss, by = 'pc')
    sd$rule1.score.Stand <- Standardization.score(merged$s1, 0.1)
    sd$rule2.score.Stand <- Standardization.score(merged$s2, 0.1)
    sd$final.score <- apply(sd[, 9:12], 1, sum)
    sd <- sd[, -2]
    sd <- sd %>% group_by(peak.num) %>% top_n(3, final.score)
    write.csv(sd, 'part2_result.csv')
  })
  
}
