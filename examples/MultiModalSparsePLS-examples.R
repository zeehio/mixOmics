\dontrun{
library(MultiAssayExperiment)
data("miniACC")
miniACC <- MatchedAssayExperiment(miniACC)
ncomp <- 2
keep_features <- list(RNASeq2GeneNorm = rep(20, ncomp), 
                      gistict= rep(10, ncomp), 
                      miRNASeqGene= rep(30, ncomp))
keep_features <- NULL
res <- MultiModalSparsePLS(data = miniACC, 
                           formula = RNASeq2GeneNorm ~ gistict + miRNASeqGene,
                           keep_features = keep_features, ncomp = ncomp,
                           design = 'null', scale = TRUE)

## sample plot
plotIndiv(res, pch = 16, group = miniACC$pathologic_stage, legend = TRUE)
}
