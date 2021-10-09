

library("DESeq")

# Human Samples
countData <- as.matrix(read.csv("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountHuman_v1/RCountH/Homo_sapiens.GRCh38.104.geneonly.csv", row.names="gene_id"))
#We need to compare :
#  A1 with B1.
countData1<-countData[,-c(2,4)]
head(countData1)

#Putting the condition.
dataDesign=data.frame(row.names=colnames(countData1),condition=c("C","D"))
dataDesign


cds = newCountDataSet(countData1,dataDesign)
conds <- factor( c( "C", "D") )
cds <- newCountDataSet( countData1, conds )
head( counts(cds) )

cds <- estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )


cds <- estimateDispersions( cds , method='blind')

 str( fitInfo(cds) )
 plotDispEsts <- function( cds )
    {
      plot(
        rowMeans( counts( cds, normalized=TRUE ) ),
        fitInfo(cds)$perGeneDispEsts,
        pch = '.', log="xy" )
      xg <- 10^seq( -.5, 5, length.out=300 )
       lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
     }
 plotDispEsts( cds )
 #Figure legend :
 #Figure 1: Empirical (black dots) and fitted (red lines) dispersion values plotted against mean
 #expression strength.
 
 head( fData(cds) )     
       
       
#Calling differential expression
 res <- nbinomTest( cds, "C", "D" )
 head(res)
 
