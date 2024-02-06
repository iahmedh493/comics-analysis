require("DESeq2")
require("limma")
require("Glimma")
require("edgeR")
require("EnsDb.Hsapiens.v79")

#files <- "/Users/ibrahimahmed/projects/GUI/counts_data.csv"
#meta_file <- "/Users/ibrahimahmed/projects/GUI/sample_info.csv"
#compare
#data <- read.csv("$COUNT_MATRIX", header=T, row.names="X")
#meta <- read.csv("$DESIGN_MATRIX", header=T, row.names = "cellLine")
data <- read.table("$COUNT_MATRIX", header=TRUE, sep= ",", row.names = 1)
meta <- read.table("$DESIGN_MATRIX", header=TRUE, sep= ",", row.names = 1)


geneid <- row.names(data)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneid, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
dup <- geneIDs1$GENEID[duplicated(geneIDs1$GENEID)]
geneIDs1[geneIDs1$GENEID %in% dup,][1:10,]
geneIDs1

mat <- match(geneid, geneIDs1$GENEID)
geneIDs1 <- geneIDs1[mat,]
geneIDs1[geneIDs1$GENEID %in% dup,][1:5,]

group<- as.factor(as.vector(meta['condition'])$condition)
count_matrix <- as.matrix(data)
d <- DGEList(counts = data, group=group)
d$samples$group <- group
d$samples$group <- group 
lane <- as.factor(rep(c("N61311","N052611","N080611", "N061011"), c(2,2,2,2))) 
d$samples$lane <- lane 
d$genes = geneIDs1
keep <- rowSums(cpm(d)>2) >= 4
d <- d[keep,]
cpm <- cpm(d)
lcpm <- cpm(d, log=TRUE)

glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), groups=d$samples$group, launch=FALSE, path = "$OUTFILE_NAME", folder="exploration", html="MDS-Plot.html")
#glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), groups=d$samples$group, launch=FALSE, path = "$OUTFILE_NAME", folder = "exploration", html="MDS-Plot.html")

design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))

#contr.matrix  <- makeContrasts($CONTRAST, levels = design)
contr.matrix <- makeContrasts(
  TreatedvUntreated = $CONTRAS, #treated-untreated, 
  levels = colnames(design))

v <- limma::voom(d, design)
vfit <- limma::lmFit(v, design)
vfit <- limma::contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)
output_dir <- "$OUTFILE_NAME"
save_path <- file.path(output_dir, "exploration", "volicano-plot.html")
#glimmaVolcano(efit, dge = d)
#htmlwidgets::saveWidget(glimmaVolcano(efit, dge = d), selfcontained = TRUE, "/Users/ibrahimahmed/projects/GUI/result_dir/exploration/volicano-plot.html")
htmlwidgets::saveWidget(glimmaVolcano(efit, dge = d), selfcontained = TRUE, save_path)






