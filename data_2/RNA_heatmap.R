library(gplots)
library(Heatplus)
library(vegan)
library(RColorBrewer)
library(circlize)


library(Heatplus)
library(vegan)
library(RColorBrewer)

MitDat <- read.csv("KK_Mitochondrial Energy Met.csv",header = T,stringsAsFactors = TRUE)
MitDat2 <- MitDat[-1,-1]
MitDat2  <- apply(as.matrix(MitDat2), 2, as.numeric)
row.names(MitDat2) <- MitDat$X[-1]

data.dist <- vegdist(MitDat2, method = "euclidean")
row.clus <- hclust(data.dist, "aver")





pdf("heatmap.pdf",height=12,width=15)
plot(annHeatmap2(data.matrix(MitDat2),
col = colorRampPalette(c("#d73027", "#fdae61","#ffffbf","#abd9e9","#2c7bb6"),
                       space = "rgb")(33), breaks = 33,
dendrogram = list(Row = list(dendro = as.dendrogram(row.clus))), legend = 2,
labels = list(Col = list(nrow = 10),cex=1.75)))
dev.off()



gg <- all.data
gg$class <- row.names(all.data)
chordDiagram(as.matrix(all.data),directional = 1,direction.type="arrows")
mat <- as.matrix(all.data)

df = data.frame(from = rep(rownames(mat), times = ncol(mat)),
                to = rep(colnames(mat), each = nrow(mat)),
                value = as.vector(mat),
                stringsAsFactors = FALSE)
grid.col = c(EV = "grey", shPDK1 = "grey", shPDK1_2 = "grey",rep("grey",25))
col_fun = colorRamp2(range(mat), c("#FFEEEE", "#FF0000"), transparency = 0.5)
chordDiagram(mat,  grid.col = "grey",column.col = c("red","blue","green"))


#tested for i in [1:5]
i=1.25

#--- Define par parameters ---
#par(mar = rep(0,4), cex = 2)
par(mar=rep(i*2,4), mgp=c(i*3,i*2, i*1.5), las=2,cex=0.45*i)

chordDiagram(mat,  grid.col = "grey",column.col = c("red","blue","green"),
direction.type = c("diffHeight", "arrows"),
link.arr.col = "grey", link.arr.length = 0.2,annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.75))
    circos.axis(h = "top", labels.cex = 1, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

circos.clear()


us.dend <- hclust(dist(scale(mat)))

EV <- mat[, 'EV']
shPDK1 <- mat[, 'shPDK1']
shPDK1_2 <- mat[, 'shPDK1_2']

## generate color maps
EV.cmap <- makecmap(EV, n = 25,colFn = jet,breaks = prettyLog)
shPDK1.cmap <- makecmap(shPDK1, n = 25, colFn = jet,breaks = prettyLog)
shPDK1_2.cmap <- makecmap(shPDK1_2, n = 25, colFn = jet,breaks = prettyLog )

us.mat <- data.frame(shPDK1_2 = cmap(shPDK1_2, EV.cmap ),
                     shPDK1 = cmap(shPDK1, EV.cmap ),
                     EV = cmap(EV, EV.cmap))

par(mar = c(4,5,4,4)+0.1)  # make space for color keys
dendromat(us.dend, us.mat)

#vkey(shPDK1_2.cmap, 'shPDK1_2',y=-0.1)
vkey(shPDK1.cmap, '',y=0.35)
#vkey(EV.cmap, 'EV',y=0.75)
