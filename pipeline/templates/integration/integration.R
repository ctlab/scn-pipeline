
## INTEGRATION

whole.features <- SelectIntegrationFeatures(object.list = whole, nfeatures = 2000)
whole <- PrepSCTIntegration(object.list = whole, anchor.features = whole.features,
                            verbose = FALSE)
whole.anchors <- FindIntegrationAnchors(object.list = whole, normalization.method = "SCT",
                                        anchor.features = whole.features, verbose = FALSE)
whole.integrated <- IntegrateData(anchorset = whole.anchors, normalization.method = "SCT",
                                  verbose = FALSE)