# Data preprocessing -----------------------------------------------------------
setwd('F:/...')
d1 <- read.csv('statTarget_judged_neg/statTarget/shiftCor/After_shiftCor/shift_all_cor.csv',check.names=F)
d2 <- read.csv('statTarget_judged_pos/statTarget/shiftCor/After_shiftCor/shift_all_cor.csv',check.names=F)
cc <- intersect(d1$sample,d2$sample)
d1 <- d1[which(d1$sample%in%cc),]
d2 <- d2[which(d2$sample%in%cc),]
fix_pos_neg_stattarget <- function(d1,d2){
  d1 <- d1[-which(d1$class=='QC'),]
  d2 <- d2[-which(d2$class=='QC'),]
  d1.meta.name <- colnames(d1)[-c(1:2)]
  d2.meta.name <- colnames(d2)[-c(1:2)]
  meta.d2 <- setdiff(d2.meta.name,d1.meta.name)
  d2 <- d2[,c(1,which(colnames(d2)%in%meta.d2))]
  data <- merge(d1,d2,by='sample')
  return(data)
}
da <- fix_pos_neg_stattarget(d1,d2)
rm(d1,d2)
da[,-c(1:2)] <- log(da[,-c(1:2)])
da$sample <- gsub('_Center.*','',da$sample)
da$sample <- gsub('\\.','-',da$sample)
da$class <- ifelse(da$class=="HD",0,1)
d1 <- da[which(da$class==0),]
d2 <- da[which(da$class==1),]
da <- rbind(d1,d2)
cn <- colnames(da)
rm(d1,d2,cn,fix_pos_neg_stattarget)
dd <- da
dd$cohort <- '1393 cohort'
dd[grep('TEST',dd$sample),]$cohort <- '333 cohort'
dd$sample <- gsub('	TEST-','',dd$sample)
dd <- dd[,c(1:2,ncol(dd),3:(ncol(dd)-1))]


# Random forest model for feature selection ------------------------------------
library(randomForest)
ln <- colnames(da)
colnames(da) <- c('sample','class',paste('Lipid',1:(length(colnames(da))-2),sep='_'))
set.seed(123)
d15 <- da[-grep('TEST',da$sample),]
d35 <- da[grep('TEST',da$sample),]
rf <- randomForest(factor(class)~.,data=d15[,-1],max_depth=5,ntree=100,importance=T)
predict <- predict(rf,d15[,-1])
table(predict,d15$class)
predict <- predict(rf,d35[,-1])
table(d35$class,predict)# 0.7147147
vip_max <- data.frame(rf$importance)
vip_max$Lipid <- ln[-c(1:2)]
m <- vip_max[order(vip_max$MeanDecreaseGini,decreasing=TRUE),][c(1:10),]
d15 <- cbind(d15[,1:2],d15[,row.names(m)])
d35 <- cbind(d35[,1:2],d35[,row.names(m)])
rm(rf,vip_max,da,predict)

# LightGBM model ---------------------------------------------------------------
library(keras)
library(reticulate)
use_python(".../python.exe",required=TRUE)
use_condaenv("PCPETG",required=TRUE)
library(lightgbm)
library(data.table)
library(caret)
dt <- data.table(d15[,-c(1:2)])
dt[,label:=d15$class]
x_matrix <- as.matrix(dt[, -"label"])
y_vector <- dt$label
lgb_data <- lgb.Dataset(data=x_matrix,label=y_vector)
params <- list(objective="binary",metric="binary_logloss",boosting_type="gbdt",
               num_leaves=31,
               learning_rate=0.01,
               num_iterations=210,
               max_depth=5)
model <- lgb.train(params = params, data = lgb_data)
x_test_dt <- as.matrix(d35[,-c(1:2)])
predictions <- predict(model, newdata = x_test_dt)
predicted_classes <- ifelse(predictions >0.5,1,0)
table(d35$class,predicted_classes)




