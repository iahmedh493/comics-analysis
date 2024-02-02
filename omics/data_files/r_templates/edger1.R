
require("edgeR")
require("limma")
require("RColorBrewer")
require("mixOmics")
require("RColorBrewer")
require("Glimma")
# ******* prepare datasets ********
rawCountTable <- read.table("$COUNT_MATRIX", header=TRUE, sep=",", row.names=1)
sampleInfo <- read.table("$DESIGN_MATRIX", header=TRUE, sep=",", row.names=1)

sampleInfo$batch <- as.factor(sampleInfo$batch)
sampleInfo$condition <- as.factor(sampleInfo$condition)
group <- as.factor(sampleInfo$condition)
all(colnames(rawCountTable) %in% rownames(sampleInfo))

all(colnames(rawCountTable) == rownames(sampleInfo))

x <- DGEList(rawCountTable, remove.zeros = TRUE)

x$samples$condition <- relevel(sampleInfo$condition, ref = "$REF")
x$samples$genotype <- sampleInfo$batch
x$samples$group <- group

lane <- sampleInfo$batch
x$samples$lane <- lane

geneid <- rownames(x)
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"),
                keytype="ENSEMBL")

genes <- genes[!duplicated(genes$ENSEMBL), ]
x$genes <- genes

cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
table(rowSums(x$counts==0)==9)

keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
samplenames <- colnames(x)

par(mar=c(1,1,1,1))
lcpm.cutoff <- log2(10/M + 2/L)

nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))

### plot1

plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
dev.copy(jpeg,filename="filteredVsRaw.jpg");
dev.off()
####

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data", ylab="Log-cpm")
dev.copy(jpeg,filename="normalizeEffect.jpg");
dev.off()

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data", ylab="Log-cpm")
dev.copy(jpeg,filename="Boxplots-of-log-CPM.jpg");
dev.off()
###
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")
dev.copy(jpeg,filename="MDS-plots-of-log-CPM.jpg");
dev.off()
#####
glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), groups=x$samples[,c(2,4)], launch=FALSE)
# saving path

design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  Untreatedvstreated = untreated-treated,
  levels = colnames(design))
contr.matrix



