\dontrun{
library(MultiAssayExperiment)
data("miniACC")
miniACC <- MatchedAssayExperiment(miniACC)
ncomp <- 2
## mixOmics needs unique feature names across modalities
miniACC <- make_unique_feature_names(miniACC)
keep_features <- list(RNASeq2GeneNorm = rep(20, ncomp), 
                      gistict= rep(10, ncomp), 
                      miRNASeqGene= rep(30, ncomp))

res <- MultiModalSparsePLSDA(data = miniACC, 
                             formula = vital_status ~
                                 RNASeq2GeneNorm + gistict + miRNASeqGene,
                             keep_features = keep_features, ncomp = ncomp,
                             design = 'full', scale = TRUE)

## sample plot
plotIndiv(res, pch = 16, group = miniACC$pathologic_stage, legend = TRUE)
circosPlot(res, cutoff = 0.8, comp = 1, block.labels.adj = -0.1, showIntraLinks = TRUE)
}
