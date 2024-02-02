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
glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), groups=d$samples$group, launch=FALSE, path = "/Users/ibrahimahmed/projects/GUI/result_dir", folder = "exploration", html="MDS-Plot.html")
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

#glimmaVolcano(efit, dge = d)
htmlwidgets::saveWidget(glimmaVolcano(efit, dge = d), selfcontained = TRUE, "/Users/ibrahimahmed/projects/GUI/result_dir/exploration/volicano-plot.html")


plotSA(efit)

#glXYPlot(x=vfit$coef, y=vfit$load, xlab="logFC", ylab="logodds",
#status=df, counts=vm, group=groups, anno=vfit$genes)

#glimmaVolcano(efit, dge = vfit)


summary(decideTests(efit))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)


head(tfit$genes$SYMBOL[de.common], n=20)

vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))
write.fit(tfit, dt, file="results.txt")

basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)
head(basal.vs.ml)
#glimmaVolcano
#glimmaVolcano(tfit, dge=v, html="vol-plot.html", status.colours=c("blue", "darkgrey", "red"))
#htmlwidgets::saveWidget(glMA, file="glimmaV2Example.html")

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13))

glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENTREZID", counts=lcpm, groups=group, launch=FALSE)



library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
   labRow=v$genes$SYMBOL[i], labCol=group,
   col=mycol, trace="none", density.info="none", 
   margin=c(8,6), lhei=c(2,10), dendrogram="column")
png("heatmap.2.png")
dev.off()
# png("heatmap.2.png")
load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p1.rdata"))
#load("/Users/ibrahimahmed/October23/SandwichApp/mouse_c2_v5p1.rdata")
idx <- ids2indices(Mm.c2,id=rownames(v)) 
cam.BasalvsLP <- camera(v,idx,design,contrast=contr.matrix[,1]) 
head(cam.BasalvsLP,5)

cam.LPvsML <- camera(v,idx,design,contrast=contr.matrix[,3]) 
head(cam.LPvsML,5)

barcodeplot(efit$t[,3], index=idx$LIM_MAMMARY_LUMINAL_MATURE_UP, index2=idx$LIM_MAMMARY_LUMINAL_MATURE_DN, main="LPvsML")

