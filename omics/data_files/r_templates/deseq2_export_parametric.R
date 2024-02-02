
res <- results(dds_res, contrast=c$CONTRAST)
res_ordered <- res[order(res$padj),]

df1 <- as.data.frame(res_ordered, stringsAsFactors=FALSE)
df1 <- tibble::rownames_to_column(df1, "X")
up <- dplyr::filter(df1, df1$log2FoldChange > 1 & df1$padj < 0.05)
up1 <- up[,1]
dn <- dplyr::filter(df1, df1$log2FoldChange > -1 & df1$padj < 0.05)
dn1 <- dn[,1]
total <- c(up1, dn1)
bkgd.genes2.entrez <- clusterProfiler::bitr(total,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db)
df = bkgd.genes2.entrez[!duplicated(bkgd.genes2.entrez$ENSEMBL),]
lcpm <- as.data.frame(lcpm, stringsAsFactors=FALSE)
lcpm <- tibble::rownames_to_column(lcpm, "X")
filtered <- lcpm %>% filter(lcpm$X %in% df$ENSEMBL)
genes <- df %>% filter(df$ENSEMBL %in% lcpm$X)
gene <- genes[,2]
row.names(filtered) <- gene
filtered <- dplyr::select(filtered, -X)
df <- log(filtered)
df <- na.omit(filtered)
df <- as.matrix(df)

#jpeg(file="/Users/ibrahimahmed/projects/GUI/filename.jpg")
#heatmap(df)
heatmaply(df[1:20,], file = "heatmaply_plot.html")

organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

original_gene_list <- df1$log2FoldChange
names(original_gene_list) <- df1$X
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

pp <- gseaplot2(gse, geneSetID = 1:3)
pp
ggsave("/Users/ibrahimahmed/projects/GUI/result_dir/DESeq2_result/clusterProfiler_GSEA_IVANOVA_small.pdf")
#ggsave(sprintf("/Users/ibrahimahmed/projects/GUI/GSEAps.pdf", p), width = 8, height = 6, onefile = T)
print(">>>><<><><><><><><><><")
print(total)
write.csv(as.data.frame(res_ordered),file="$OUTFILE_NAME")