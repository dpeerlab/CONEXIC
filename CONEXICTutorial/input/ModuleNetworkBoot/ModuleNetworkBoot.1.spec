score NormalGamma alpha=2 lambda=1
Prior LeafPenalty 20
Prior RegulatorPenalty 15
Constraint LeafMaximum 6
algorithm GreedyHillClimbing local
data GeneExpressionFile ../../data/TCGA_GBM_Agilent_Exp.txt regulators=../../data/GBMPeakGenesSep10.ALL otherRegulators=../../data/GBMPeakGenesSep10.continuous.CN.matrix
assignmentAlgorithm GeronemoIteration likelihoodcutoff=2 
choiceTest GeronemoTest
Sample 136
NumRuns 1
