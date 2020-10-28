\dontrun{
library(MultiAssayExperiment)
data("miniACC")
miniACC <- MatchedAssayExperiment(miniACC)
ncomp <- 2
## mixOmics needs unique feature names across modalities
miniACC <- make_unique_feature_names(miniACC)
test_keep_features <- list(RNASeq2GeneNorm = c(20, 30), 
                           gistict= c(20, 30), 
                           miRNASeqGene= c(20, 30, 40))

tune.res <- TuneMultiModalSparsePLSDA(data = miniACC, 
                                      formula = vital_status ~ 
                                          RNASeq2GeneNorm + gistict + miRNASeqGene,
                                      test_keep_features = test_keep_features, 
                                      ncomp = ncomp, design = 'full', scale = TRUE, 
                                      BPPARAM = BiocParallel::MulticoreParam(parallel::detectCores()-1))

## sample plot
plot(tune.res)
tune.res$choice.keepX
}
