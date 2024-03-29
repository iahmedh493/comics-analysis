require("DESeq2")
require("dplyr")
require("org.Hs.eg.db")
require("tibble")
require("clusterProfiler")
require("enrichplot")
require("ggplot2")
require("heatmaply")

count_data <- read.table("$COUNT_MATRIX", header=TRUE, sep= ",", row.names = 1)
design_matrix <- read.table("$DESIGN_MATRIX", header=TRUE, sep= ",")
dds <- DESeqDataSetFromMatrix(count_data, design_matrix, $FORMULA)
normdds <- DESeq2::estimateSizeFactors(dds)
lcpm <- log2(DESeq2::counts(normdds, normalized = FALSE) + 0.5)
dds_res <- DESeq(dds)

