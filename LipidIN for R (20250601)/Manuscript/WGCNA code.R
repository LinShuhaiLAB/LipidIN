# Data preprocessing -------------------------------------------------------------------
setwd('F:/...')
d1 <-
  read.csv(
    'statTarget_judged_neg/statTarget/shiftCor/After_shiftCor/shift_all_cor.csv',
    check.names = F
  )
d2 <-
  read.csv(
    'statTarget_judged_pos/statTarget/shiftCor/After_shiftCor/shift_all_cor.csv',
    check.names = F
  )
fix_pos_neg_stattarget <- function(d1, d2) {
  d1 <- d1[-which(d1$class == 'QC'), ]
  d2 <- d2[-which(d2$class == 'QC'), ]
  d1.meta.name <- colnames(d1)[-c(1:2)]
  d2.meta.name <- colnames(d2)[-c(1:2)]
  meta.d2 <- setdiff(d2.meta.name, d1.meta.name)
  d2 <- d2[, c(1, which(colnames(d2) %in% meta.d2))]
  data <- merge(d1, d2, by = 'sample')
  return(data)
}
da <-
  fix_pos_neg_stattarget(d1, d2) # Combined positive and negative ion modes
rm(d1, d2)
da <-
  da[, -grep('FFA', colnames(da))] # Delete internal labeling FFA-d3 18:0
d1 <- da[which(da$class == 'BC'), ] # Deletion of disease groups
sample.name <- lapply(3:ncol(d1), function(i) {
  d11 <- d1[, i]
  lower_limit <- quantile(d11, 0.25)
  upper_limit <- quantile(d11, 0.75)
  which(d11 >= lower_limit & d11 <= upper_limit)
})
sample.name <- as.data.frame(table(do.call(c, sample.name)))
sample.name <- sample.name[which(sample.name$Freq > 0.8 * nrow(d1)), 1]
d1 <-
  d1[sample.name, ] # Filtering of lipids according to the 80% repetition principle
rm(da, sample.name)

# WGCNA_BC -------------------------------------------------------------------
library(WGCNA)
library(genefilter)
library(DT)
library(RColorBrewer)
library(gplots)
d1 <- t(d1)
colnames(d1) <- d1[1, ]
data <- d1[-c(1:2), ]
data <- apply(data, 2, as.numeric)
row.names(data) <- row.names(d1)[-c(1:2)]
info <-
  read.csv('.../BC_info.csv', check.names = F) # Loading information on clinical indicators
colnames(data) <-
  gsub('\\.', '-', gsub('_.*', '', gsub('.*_QC', '', colnames(data))))
data <- data[, which(colnames(data) %in% info$Label)]
cutoff <- 0.001
datExpr <- as.matrix(data)
rownames(datExpr) <- row.names(data)
datExpr <-
  varFilter(
    datExpr,
    var.func = IQR,
    var.cutoff = cutoff,
    filterByQuantile = TRUE
  )
dim(data)
dim(datExpr)
datExpr <- t(datExpr) # data transposition
gsg <-
  goodSamplesGenes(datExpr, verbose = 3) # Check for missing or outlying values
any(do.call(c, gsg) == F)
sampleTree <-
  hclust(dist(datExpr, method = 'manhattan'), method = "average")
sampleTree
par(mar = c(1, 1, 1, 1))
plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "",
  cex.lab = 1.2,
  cex.axis = 1.2,
  cex.main = 1.5
)
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
power <- 10
net <- blockwiseModules(
  datExpr,
  checkMissingData = TRUE,
  maxBlockSize = 40000,
  power = power,
  TOMType = "signed",
  saveTOMs = FALSE,
  corType = "pearson",
  saveTOMFileBase = "blockTOM",
  minModuleSize = 20,
  reassignThreshold = 0,
  mergeCutHeight = 0.15,
  networkType = "signed",
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  verbose = 3
)
mergedColors <- labels2colors(net$colors)
table(mergedColors)
plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
sampleMatrix <- as.data.frame(diag(x = 1, nrow = nSamples))
rownames(sampleMatrix) <- rownames(datExpr)
colnames(sampleMatrix) <- rownames(datExpr)
geneTraitSignificance <-
  as.data.frame(cor(datExpr, sampleMatrix, use = "p"))
geneTraitColor <-
  as.data.frame(numbers2colors(
    geneTraitSignificance,
    signed = TRUE,
    colors = colorRampPalette(c("blue", "white", "red"))(100)
  ))
names(geneTraitColor) <- colnames(sampleMatrix)
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
# Recalculate MEs with color labels
MEList <- moduleEigengenes(datExpr, colors = mergedColors)
MEs <- MEList$eigengenes
MET <- orderMEs(MEs)
nGenes <- ncol(datExpr)
if (nGenes > 400) {
  nSelect <- 400
} else {
  nSelect <- nGenes
}
dissTOM = 1 - TOMsimilarityFromExpr(datExpr, power = power)
select = sample(nGenes, size = nSelect)

selectTOM = dissTOM[select, select]

selectTree = hclust(as.dist(selectTOM), method = 'average')
selectColors = mergedColors[select]

plotDiss = selectTOM ^ power

diag(plotDiss) = NA

color <-
  colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
TOMplot(
  dissim = (1 - plotDiss),
  dendro = selectTree,
  Colors = selectColors,
  col = color,
  main = "Network heatmap plot, selected genes"
)
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
sampleMatrix <- as.data.frame(diag(x = 1, nrow = nSamples))
rownames(sampleMatrix) <- rownames(datExpr)
colnames(sampleMatrix) <- rownames(datExpr)
moduleSampleCor <- cor(MET, sampleMatrix, use = "p")
moduleSamplePvalue <- corPvalueStudent(moduleSampleCor, nSamples)
textMatrix <- paste(signif(moduleSampleCor, 2),
                    "\n(",
                    signif(moduleSamplePvalue, 2),
                    ")",
                    sep = "")
dim(textMatrix) <- dim(moduleSampleCor)
MEList <- moduleEigengenes(datExpr, colors = mergedColors)
MEs <- MEList$eigengenes
MET <- orderMEs(MEs)
dataTrait <- info[, -c(1:3)]
row.names(dataTrait) <- info$Label
dataTrait <-
  dataTrait[which(row.names(dataTrait) %in% colnames(sampleMatrix)), ]
moduleTraitCor = cor(MET, dataTrait, use = "p")
moduleTraitCor1 = cor(MET, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2),
                   "\n(",
                   signif(moduleTraitPvalue, 1),
                   ")",
                   sep = "")
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(moduleTraitCor),
  yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  xLabelsAngle = 45,
  cex.lab = 1,
  colorLabels = FALSE,
  colors = blueWhiteRed(100),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  yColorWidth = 0.03,
  xColorWidth = 0.05,
  main = paste("Module-Trait relationship")
)
write.csv(moduleTraitCor, 'BC correlation matrix.csv', row.names = T)
write.csv(moduleTraitPvalue, 'BC p value.csv', row.names = T)
moduleColors <- labels2colors(net$colors)
TOM <- TOMsimilarityFromExpr(datExpr, power = power)
modules <- moduleColors
allModules <- unique(moduleColors)
allModules
Genes <- colnames(datExpr)
inModule <- is.finite(match(moduleColors, modules))
modGenes <- Genes[inModule]
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modGenes, modGenes)

cyt <- exportNetworkToCytoscape(
  modTOM,
  edgeFile = "Edge.all.txt",
  nodeFile = "Node.all.txt",
  weighted = TRUE,
  threshold = 0.2,
  nodeNames = modGenes,
  altNodeNames = modGenes,
  nodeAttr = moduleColors[inModule]
)
ADmatrix <- abs(cor(datExpr, use = "p")) ^ power
ALLdegree <- intramodularConnectivity(ADmatrix, modules)
allIDs <- colnames(datExpr)
all_fpkm <- data[allIDs, ]
all_Connectivity <- ALLdegree[allIDs, ]
out_all <- data.frame(
  geneID = rownames(all_fpkm),
  moduleColors = moduleColors,
  all_Connectivity,
  all_fpkm,
  check.names = F
)
write.csv(out_all, 'BC cluster.csv', row.names = F)
