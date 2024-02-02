
#require("DOSE")
#require("GO.db")
#require("GSEABase")
#require("org.Hs.eg.db")
###require("clusterProfiler")
#require("dplyr")
#require("tidyr")
#require("ggplot2")
#require("stringr")
#require("RColorBrewer")
#require("rWikiPathways")
#require("RCy3")
##require("pacman")
#require("stringApp")

if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways", update = FALSE)
}
library(rWikiPathways)
library(rWikiPathways)


load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
  print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}

library(clusterProfiler)
library(org.Hs.eg.db)

lung.expr <- read.table("$COUNT_MATRIX", header=TRUE, sep=",")
#lung.expr <- read.csv(system.file("$COUNT_MATRIX", package="rWikiPathways"),stringsAsFactors = FALSE)
up.genes <- lung.expr[lung.expr$log2FC > 1 & lung.expr$adj.P.Value < 0.05, 1]
dn.genes <- lung.expr[lung.expr$log2FC < -1 & lung.expr$adj.P.Value < 0.05, 1]
bkgd.genes <- lung.expr[,1]

up.genes.entrez <- clusterProfiler::bitr(up.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
#cat("\n\nWhich column contains my new Entrez IDs?\n")
#head(up.genes.entrez)
#keytypes(org.Hs.eg.db)

dn.genes.entrez <- clusterProfiler::bitr(dn.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez <- clusterProfiler::bitr(bkgd.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
print(":::::::::I am here::::::::::")
print(":::::::::I am here::::::::::")
print(":::::::::I am here::::::::::")

egobp <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez[[2]],
  universe = bkgd.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

barplot(egobp, showCategory = 10)
dotplot(egobp, showCategory = 10)
goplot(egobp)

print(":::::::::I am here::::::::::")
#options(ggrepel.max.overlaps = Inf)
ggplot(egobp[1:20], aes(x=reorder(Description, -pvalue), y=Count, fill=-p.adjust)) +
  geom_bar(stat = "identity") +
  #max.overlaps = getOption("ggrepel.max.overlaps", default = 10)
  coord_flip() +
  scale_fill_continuous(low="blue", high="red") +
  labs(x = "", y = "", fill = "p.adjust") +
  theme(axis.text=element_text(size=11))

## extract a dataframe with results from object of type enrichResult
egobp.results.df <- egobp@result

## create a new column for term size from BgRatio
egobp.results.df$term.size <- gsub("/(\\d+)", "", egobp.results.df$BgRatio)

## filter for term size to keep only term.size => 3, gene count >= 5 and subset
egobp.results.df <- egobp.results.df[which(egobp.results.df[,'term.size'] >= 3 & egobp.results.df[,'Count'] >= 5),]
egobp.results.df <- egobp.results.df[c("ID", "Description", "pvalue", "qvalue", "geneID")]

## format gene list column
egobp.results.df$geneID <- gsub("/", ",", egobp.results.df$geneID)

## add column for phenotype
egobp.results.df <- cbind(egobp.results.df, phenotype=1)
egobp.results.df <- egobp.results.df[, c(1, 2, 3, 4, 6, 5)]

## change column headers
colnames(egobp.results.df) <- c("Name","Description", "pvalue","qvalue","phenotype", "genes")

egobp.results.filename <-file.path("$OUTFILE_NAME",paste("clusterprofiler_cluster_enr_results.txt",sep="_"))
write.table(egobp.results.df,egobp.results.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

ewp.up <- clusterProfiler::enrichWP(
  up.genes.entrez[[2]],
  universe = bkgd.genes.entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.1, #p.adjust cutoff; relaxed for demo purposes
)

ewp.up <- DOSE::setReadable(ewp.up, org.Hs.eg.db, keyType = "ENTREZID")
barplot(ewp.up, showCategory = 8)
dotplot(ewp.up, showCategory = 8)

ewp.dn <- enrichWP(
  dn.genes.entrez[[2]],
  #universe = bkgd.genes[[2]],  #hint: comment out to get any results for demo
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.1, #p.adjust cutoff; relaxed for demo purposes
)

ewp.dn <- setReadable(ewp.dn, org.Hs.eg.db, keyType = "ENTREZID")
dotplot(ewp.dn, showCategory = 8)

lung.expr$fcsign <- sign(lung.expr$log2FC)
lung.expr$logfdr <- -log10(lung.expr$P.Value)
lung.expr$sig <- lung.expr$logfdr/lung.expr$fcsign
sig.lung.expr.entrez<-merge(lung.expr, bkgd.genes.entrez, by.x = "GeneID", by.y = "ENSEMBL")
gsea.sig.lung.expr <- sig.lung.expr.entrez[,8]
names(gsea.sig.lung.expr) <- as.character(sig.lung.expr.entrez[,9])
gsea.sig.lung.expr <- sort(gsea.sig.lung.expr,decreasing = TRUE)

gwp.sig.lung.expr <- clusterProfiler::gseWP(
  gsea.sig.lung.expr,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff
  organism = "Homo sapiens"
)

gwp.sig.lung.expr.df = as.data.frame(gwp.sig.lung.expr)
gwp.sig.lung.expr.df[which(gwp.sig.lung.expr.df$NES > 1),] #pathways enriched for upregulated lung cancer genes
gwp.sig.lung.expr.df[which(gwp.sig.lung.expr.df$NES < -1),] #pathways enriched for downregulated lung cancer genes

findPathwayNamesByText("lung cancer")

lc.pathways <- findPathwaysByText("lung cancer")  #quotes inside query to require both terms
human.lc.pathways <- lc.pathways %>%
  dplyr::filter(species == "Homo sapiens") # just the human lung cancer pathways
human.lc.pathways$name # display the pathway titles

lc.wpids <- human.lc.pathways$id
lc.wpids

ewp.up.wpids <- ewp.up$ID
ewp.up.wpids

#lung.expr <- read.csv("/Users/ibrahimahmed/October23/project1/mygenelist.txt", stringsAsFactors = FALSE)
#head(lung.expr)

