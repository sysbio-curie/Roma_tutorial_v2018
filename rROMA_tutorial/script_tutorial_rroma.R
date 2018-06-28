library(rRoma)
library(readr)

###################
### Load data

samples <- read.table(file = "sample_annotation.txt",header=TRUE) ### read sample annotation
expr <- as.matrix(read.table("matrix_gene_expression.txt",header=T, row.names=1)) ## Read expression matrix
expr <- expr[, intersect(samples$Sample, colnames(expr))]  ### Select samples with available annotation

###################
### Set sample annotation as factor variable

Group <- as.character(samples$Group)
names(Group) <- samples$Sample
Group <- factor(Group, levels = c("G1", "G2", "G3", "G4"))

###################
#### Testing different signatures for a given pathway (ex: wnt)

wntGMT <- ReadGMTFile("Unsigned_wnt_path.gmt", SearchString = NULL, Mode = "ANY")

Data.wnt <- rRoma.R(ExpressionMatrix = expr, ModuleList = wntGMT, FixedCenter = TRUE, 
                    MaxGenes = 1000, PCSignMode="CorrelateAllWeightsByGene", 
                    PCAType = "DimensionsAreSamples")

###################
### Explore table of results 

Data.wnt$ModuleMatrix  ### Two WNT signatures (WNT_pthw_Metastasis,wnt_IPA) 
                       ### perform better than the others

###################
#### Testing many signatures from a given DB (ex: MsigDB Hallmarsks )

AllHall <- SelectFromMSIGdb("HALLMARK")
#selGMT <- SelectFromMSIGdb(SearchString = c("xenobiotic", "keg"), Mode = "ALL")
#Hallv6 <- ReadGMTFile("h.all.v6.1.symbols.gmt", SearchString = NULL, Mode = "ANY")

Data.hall <- rRoma.R(ExpressionMatrix = expr, ModuleList = AllHall, FixedCenter = TRUE, 
                     MaxGenes = 1000, PCSignMode="CorrelateAllWeightsByGene", 
                     PCAType = "DimensionsAreSamples")

#Data.hallv6 <- rRoma.R(ExpressionMatrix = expr, ModuleList = Hallv6, FixedCenter = TRUE, MaxGenes = 1000, PCSignMode="CorrelateAllWeightsByGene", PCAType = "DimensionsAreSamples")


###################
### Create plots to visualize Module actvity 

### Boxplot modules (real activity score vs sampling)
### Heatmap activity score (selected according to VarThr pv threshold, VarMode = "Wil" or "PPV" for Wilcoxon or permutation test)

pdf(file = "Plot_heatmap_activity.pdf")
AggData.FC <- Plot.Genesets(RomaData = Data.hall,
              Selected = SelectGeneSets(RomaData = Data.hall, VarThr = 1e-05,
                                        VarMode = "Wil", VarType = "Over"),
              GenesetMargin = 20, SampleMargin = 14, cluster_cols = TRUE,
              GroupInfo = Group, AggByGroupsFL = c("mean", "sd"),
              HMTite = "Overdispersed genesets (Fixed center)")

dev.off()

###################
### Statistical comparison across samples (differential analysis based on pathway activity)

pdf(file = "Plot_diff_activity.pdf")
CompareAcrossSamples(RomaData = Data.hall,
                     Selected = SelectGeneSets(RomaData = Data.hall, VarThr = 1e-05,
                                               VarMode = "Wil", VarType = "Over"),
                     Groups = Group)
dev.off()

###################
### Top contributing genes

png(file = "Heatmap_Top_genes.png")
GeneMat <- GetTopContrib(Data.hall,
                         Selected = SelectGeneSets(RomaData = Data.hall, VarThr = 1e-5,
                                                   VarMode = "Wil", VarType = "Over"),
                         nGenes = .1, OrderType = "Abs", Mode = "Wei", Plot = TRUE)


dev.off()


GeneMat <- GetTopContrib(Data.hall,
                         Selected = SelectGeneSets(RomaData = Data.hall, VarThr = 1e-5,
                                                   VarMode = "Wil", VarType = "Over"),
                         ExpressionMatrix = expr,
                         nGenes = .1, OrderType = "Abs", Mode = "Cor", Plot = TRUE)

pdf(file = "Plot_Top_genes.pdf")

PlotGeneWeight(RomaData = Data.hall, PlotGenes = 30,
               ExpressionMatrix = expr, LogExpression = FALSE,
               Selected = SelectGeneSets(RomaData = Data.hall, VarThr = 1e-5,
                                         VarMode = "Wil", VarType = "Over"),
               PlotWeigthSign = TRUE)

dev.off()


###################
### Visualization of ROMA score on ACSN maps

Data.NFC.CS <- rRoma.R(ExpressionMatrix = expr,
                       ModuleList = ReadGMTFile("http://acsn.curie.fr/files/survival_v1.1.gmt"),
                       FixedCenter = FALSE, MaxGenes = 1000,
                       UseParallel = TRUE, nCores = 8, ClusType = "PSOCK",
                       PCSignMode="CorrelateAllWeightsByGene")

PlotOnACSN(RomaData = Data.NFC.CS, SampleName = names(Group[Group == "G2"]),
           AggScoreFun = "median", FilterByWei = 30, 
           DispMode = c("Module", "Gene"), DataName = "Normal", 
           QTop = .99, QBottom = .1,
           Selected = NULL,
           MapURL = 'https://acsn.curie.fr/navicell/maps/survival/master/index.php',
           AggGeneFun = "median")





###################
### rROMA Dashboard

library(rRomaDash)
rRomaDash()
