score NormalGamma alpha=2 lambda=1 
data GeneExpressionFile ../../data/TCGA_GBM_Agilent_Exp.txt regulators=../../data/GBMPeakGenesSep10.ALL CNVRegulators=../../data/GBMPeakGenesSep10.AMP,../../data/GBMPeakGenesSep10.DEL CNVData=../../data/biolearn.gene.matrix
moduleInitiation RegCopyNumberClustering  amplified_list=../../data/GBMPeakGenesSep10.AMP deleted_list=../../data/GBMPeakGenesSep10.DEL amplified_list_cnv=../../data/GBMPeakGenesSep10.AMP deleted_list_cnv=../../data/GBMPeakGenesSep10.DEL first_output=SingleModulatorBoot.1.step1  second_output=SingleModulatorBoot.1.step2  rejected_list=SingleModulatorBoot.1.rejected_genes  permutations=10000 pvaluethreshold=0.001 WelchTTestthreshold=0.05 AllowSelfRegulation=FALSE MinClusterSize=20 noUpdown 
ClusterOnly
Sample 136
NumRuns 1
