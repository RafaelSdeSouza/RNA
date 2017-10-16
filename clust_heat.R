library(gplots)
library(Heatplus)
library(vegan)
library(RColorBrewer)
library(circlize)
require(genefilter)
require(DESeq2)

all.data <- read.csv("grant.csv")

all.data <- all.data[-c(10,12),]
row.names(all.data) <- all.data$Gene.name
all.data  <- all.data[,-1]
colnames(all.data) <- c("EV","shPDK1","shPDK1_2")


colData(as.matrix(all.data))
DESeqDataSetFromMatrix(as.matrix(all.data))


sampleDists <- dist( t( assay(all.data) ) )
as.matrix( sampleDists )[ 1:3, 1:3 ]




countData <- as.matrix(round(all.data))
treatment <- factor(c("EV","shPDK1","shPDK1_2"))
rna <- row.names(all.data) 
dds <- DESeqDataSetFromMatrix(countData, colData=DataFrame(treatment),design =  ~ rna+treatment  )
rld <- rlog(dds)
sampleDists <- dist( t( assay(rld) ) )

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$treatment, 
                                     rld$patient, sep="-" )
colnames(sampleDistMatrix) <- NULL   
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

