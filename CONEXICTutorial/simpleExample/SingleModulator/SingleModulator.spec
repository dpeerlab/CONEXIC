score NormalGamma alpha=2 lambda=1 
data GeneExpressionFile ../../data/TCGA_GBM_Agilent_Exp.txt regulators=candidate.List CNVData=../../data/biolearn.gene.matrix
moduleInitiation RegCopyNumberClustering  amplified_list=candidate.List.AMP deleted_list=candidate.List.DEL AllowSelfRegulation=FALSE noUpdown 
ClusterOnly
