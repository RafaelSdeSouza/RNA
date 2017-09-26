library(gplots)  
library(Heatplus)
library(vegan)
library(RColorBrewer)

all.data <- read.csv("grant.csv")

all.data <- all.data[-c(10,12),]
row.names(all.data) <- all.data$Gene.name
all.data  <- all.data[,-1]
colnames(all.data) <- c("EV","shPDK1","shPDK1_2")


heatmap(as.matrix(all.data), Rowv = NA, Colv = NA)

data.dist <- vegdist(all.data, method = "bray")
row.clus <- hclust(data.dist, "aver")
heatmap(as.matrix(all.data), Rowv = as.dendrogram(row.clus), Colv = NA, margins = c(10, 5))





plot(annHeatmap2(as.matrix(all.data),

                 # the colours have to be fiddled with a little more than with the other functions. With 50 breaks in the heatmap densities, there needs to be 51 colours in the scale palette. Error messages from the plotting should help to determine how many colours are needed for a given number of breaks

                 col = colorRampPalette(c("blue2", "red2"), space = "rgb")(24), breaks = 24,

                 # dendrogram controls how dendrograms are made and can be calculatated witin the function. The internal calculation can be overridden and externally calculated dendrograms can be displayed.

                 dendrogram = list(Row = list(dendro = as.dendrogram(row.clus))), legend = 3, # this puts the colour-scale legend on the plot. The number indicates the side on which to plot it (1 = bottom, 2 = left, 3 = top, 4 = right)

                 labels = list(Col = list(nrow = 12)) # gives more space for the Genus names
                 ))
