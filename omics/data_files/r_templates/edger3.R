
require("edgeR")
require("limma")
require("RColorBrewer")
require("mixOmics")
# ******* prepare datasets ********
rawCountTable <- read.table("$COUNT_MATRIX", header=TRUE, sep=",", row.names=1)
sampleInfo <- read.table("$DESIGN_MATRIX", header=TRUE, sep=",", row.names=1)

sampleInfo$batch <- as.factor(sampleInfo$batch)
sampleInfo$condition <- as.factor(sampleInfo$condition)

all(colnames(rawCountTable) %in% rownames(sampleInfo))

all(colnames(rawCountTable) == rownames(sampleInfo))

dgeFull <- DGEList(rawCountTable, remove.zeros = TRUE)

dgeFull$samples$condition <- relevel(sampleInfo$condition, ref = "$REF")
dgeFull$samples$genotype <- sampleInfo$batch

# ********* end of data preparation *************
# Data exploration and quality assessment
pseudoCounts <- log2(dgeFull$counts + 1)

d <- colnames(pseudoCounts)
for (sample in d) {
  
  name1 = paste(sample,".png", sep = "")
  fname = paste("/Users/ibrahimahmed/projects/GUI/result_dir/edger_result/", name1, sep = "")
  png(filename = fname)
  hist(pseudoCounts[ , sample], main = "", xlab = "counts")
  dev.off()
}

par(mar = c(8,4,1,2))
#png(filename = "box_plot.png")
boxplot(pseudoCounts, col = "gray", las = 3, cex.names = 1)
dev.copy(jpeg,filename="/Users/ibrahimahmed/projects/GUI/result_dir/edger_result/plot.jpg");
dev.off()

limma::plotMA(pseudoCounts[ ,1:2], xlab = "M", ylab = "A", main = "")
abline(h = 0, col = "red")
dev.copy(jpeg,filename="/Users/ibrahimahmed/projects/GUI/result_dir/edger_result/plot2.jpg");
dev.off()

colConditions <- brewer.pal(3, "Set2")
colConditions <- colConditions[match(sampleInfo$condition,
                                     levels(sampleInfo$condition))]
pchGenotypes <- c(8, 15, 16)[match(sampleInfo$batch,
                                   levels(sampleInfo$batch))]
plotMDS(pseudoCounts, pch = pchGenotypes, col = colConditions)
legend("topright", lwd = 2, col = brewer.pal(3, "Set2")[1:2], 
       legend = levels(sampleInfo$condition))
legend("bottomright", pch = c(8, 15, 16), 
       legend = levels(sampleInfo$batch))
dev.copy(jpeg,filename="/Users/ibrahimahmed/projects/GUI/result_dir/edger_result/plot3.jpg");
dev.off()

sampleDists <- as.matrix(dist(t(pseudoCounts)))

cimColor <- colorRampPalette(rev(brewer.pal(9, "Reds")))(16)
cim(sampleDists, color = cimColor, symkey = FALSE, row.cex = 0.7, col.cex = 0.7)
dev.copy(jpeg,filename="/Users/ibrahimahmed/projects/GUI/result_dir/edger_result/plot4.jpg");
dev.off()
# ********* end of exploration ***************************************
dgeFull <- calcNormFactors(dgeFull, method="TMM")
normCounts <- cpm(dgeFull)
pseudoNormCounts <- cpm(dgeFull, log = TRUE, prior.count = 1)
par(mar = c(8,4,1,2))
boxplot(pseudoNormCounts, col = "gray", las = 3, cex.names = 1)
dev.copy(jpeg,filename="/Users/ibrahimahmed/projects/GUI/result_dir/edger_result/plot5.jpg");
dev.off()

plotMDS(pseudoNormCounts, pch = pchGenotypes, col = colConditions)
legend("topright", lwd = 2, col = brewer.pal(3, "Set2")[1:2], 
       legend = levels(sampleInfo$condition))
legend("bottomright", pch = c(8, 15, 16), 
       legend = levels(sampleInfo$batch))
dev.copy(jpeg,filename="/Users/ibrahimahmed/projects/GUI/result_dir/edger_result/plot5.jpg");
dev.off()

# Differential analysis
# create DGEList
dgeFull.group <- DGEList(rawCountTable, remove.zeros = TRUE, 
                         group = dgeFull$samples$condition)

dgeFull.group$samples$norm.factors <- dgeFull$samples$norm.factors
# Estimate dispersion
dgeFull.group <- estimateCommonDisp(dgeFull.group)
dgeFull.group <- estimateTagwiseDisp(dgeFull.group)
dgeFull.group

plotBCV(dgeFull.group)
dev.copy(jpeg,filename="/Users/ibrahimahmed/projects/GUI/result_dir/edger_result/plot8.jpg");
dev.off()
#Perform the test
dgeExactTest <- exactTest(dgeFull.group)

resExactTest <- topTags(dgeExactTest, n = nrow(dgeExactTest$table))

par(mfrow = c(1,2))
hist(resExactTest$table$PValue, xlab = "p-value", main = "raw p-values")
hist(resExactTest$table$FDR, xlab = "p-value", main = "adjusted p-values")
dev.copy(jpeg,filename="/Users/ibrahimahmed/projects/GUI/result_dir/edger_result/plot0.jpg");
dev.off()

selectedET <- resExactTest$table$FDR < 0.05 & abs(resExactTest$table$logFC) > 1
selectedET <- resExactTest$table[selectedET, ]

selectedET$updown <- factor(ifelse(selectedET$logFC > 0, "up", "down"))
write.table(selectedET, file = "/Users/ibrahimahmed/projects/GUI/result_dir/edger_result/tomatoDEG.csv", sep = ",")

#dgeFilt <- HTSFilter(dgeFull)$filteredData
#dgeFilt$samples$group <- dgeFilt$samples$condition
#plotSmear(dgeFilt, de.tags = rownames(selectedFilt))