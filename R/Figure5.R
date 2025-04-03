##### Figure 5 and Figure S5
### working pipeline
# 1. index ·
# 2. estimate tumor purity ·
# 3. TSNet tumor purity ·
# 4.1 add CAF, endo and epi and run BayesDebulk ·
# 4.2 modualize immune pathway in RNA and Pro and numberize the modules from 1-k separately
# 5. clustering using immune cell and Protein modual
# 6. check function, prognosis and ZDHHC11
# 7. making heatmap
rm(list = ls()); try(dev.off(), silent = T)
# packages loading
library(tidyverse); library(data.table)
library(ggpubr); library(cowplot); library(ggthemes); library(ggplot2); library(ggrepel); library(ggside)
library(maftools)
library(ComplexHeatmap); library(circlize)
library(openxlsx)

##### prep1. data loading #####
##### 1 omics
load('../DataMain.Rda')
load('ImmuneCluster.RData') # Output
load('Index of this part.RData')
# data for Table S5

#### 2 color
load('../Colors (ggsci).RData')
library(wesanderson)
scales::show_col(ColJournal$Science)
ColImmune <- c(ColJournal$Nature[5], ColJournal$Nature[1],
               ColJournal$Nature[4], ColJournal$Nature[9])
names(ColImmune) <- c('A','B','C','D')

### 3 Genesets
load('../Hallmark_KEGG_Reactome.RData')
ref.go <- clusterProfiler::read.gmt('../C5 GO CC symbol.txt') # GO-CC

# 4 COSMIC Cancer Gene Census
ref.cosmic <- fread('../COSMIC genes.csv')

# 5. others
# Kinase
library(KSEAapp)
ref.kinase = KSData
# Chr info
load('../ChrBandLength.RData')
# Gene coordination
GeneCoord <- rtracklayer::import('../Homo_sapiens.GRCh38.108.gtf') %>%
  as.data.frame() %>% select(1,2,3,7,12) %>% distinct(seqnames, gene_name, .keep_all = T) %>% na.omit()
GeneCoord <- mutate(GeneCoord, seqnames = paste0('chr', seqnames) %>% factor(levels = levels(Chr$Chr), ordered = T)) %>%
  na.omit() %>% arrange(seqnames, start)
GeneCoord = filter(GeneCoord, )
# add cytoband info
GeneCoord$ChrBand = apply(GeneCoord, 1, function(x){
  CurrentBand = filter(Chr, Chr == as.character(x[1]), Start <= as.numeric(x[2]), 
                       End >= as.numeric(x[3]))$ChrBand 
  CurrentBand = ifelse(length(CurrentBand) == 0, NA, as.character(CurrentBand))
}) %>%
  unlist() %>% as.character()
GeneCoord = filter(GeneCoord, ChrBand != '')
GeneCoord = as_tibble(GeneCoord)


##### Crucial: Construct immune subtypes using Wang Pei method #####
##### 1.0 Index #####
IndexDecov <- rbind(select(GBC_Main_Index, 1,2,4,6) %>% 
                      bioinfoamateur::dfm_change_colnames(c('Sample','WES','RNA','Pro')), 
                    mutate(GBC_Main_Index, Sample = paste0(Sample, '_P')) %>% select(1,3,5,7) %>% 
                      bioinfoamateur::dfm_change_colnames(c('Sample','WES','RNA','Pro')))
IndexDecov$Tumor <- ifelse(!grepl(fixed('_P'),IndexDecov$Sample), 'Yes', NA)
# retain samples having proteomic data
IndexDecov <- filter(IndexDecov, !is.na(Pro))

##### 1.1 Omic matrix prep #####
DecovRNA <- GBC_Main_RNA_logTPM[,na.omit(IndexDecov$RNA)]
DecovPro <- GBC_Main_Pro[,na.omit(IndexDecov$Pro)]
colnames(DecovRNA) = sapply(colnames(DecovRNA), function(x){filter(IndexDecov, RNA == x)$Sample})
colnames(DecovPro) = sapply(colnames(DecovPro), function(x){filter(IndexDecov, Pro == x)$Sample})

##### 1.2 estimating tumor purity using TSNet and calc immune&stroma score #####
##### 1.2.0 add origin tumor purity by HE
IndexDecov$TumoPurityImage <- sapply(IndexDecov$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Tumor_purity %>% as.numeric()}) %>% as.numeric()

##### 1.2.1 step 1 prior knowledge using ESTIMATE or Image
library(IOBR)
### RNA
TempEstimate_RNA <- deconvo_tme(eset = DecovRNA, method = "estimate")
IndexDecov$TumoPurityESTIMATE_RNA <- sapply(IndexDecov$Sample, function(x){filter(TempEstimate_RNA, ID == x)$TumorPurity_estimate %>% as.numeric()}) %>% as.numeric()
IndexDecov$TumoPurityESTIMATE_RNA[is.na(IndexDecov$Tumor)] = NA
### Pro
TempEstimate_Pro <- deconvo_tme(eset = DecovPro, method = "estimate")
IndexDecov$TumoPurityESTIMATE_Pro <- sapply(IndexDecov$Sample, function(x){filter(TempEstimate_Pro, ID == x)$TumorPurity_estimate %>% as.numeric()}) %>% as.numeric()
IndexDecov$TumoPurityESTIMATE_Pro[is.na(IndexDecov$Tumor)] = NA
### Mean
IndexDecov$TumorPurityInteg <- apply(IndexDecov[,6:8], 1, function(x){mean(as.numeric(as.vector(x)),
                                                                           na.rm = T, trim = 0)})
# cor.test
cor.test(IndexDecov$TumoPurityImage, IndexDecov$TumorPurityInteg)
cor.test(IndexDecov$TumoPurityESTIMATE_RNA, IndexDecov$TumorPurityInteg)

# conclusion: cor is high between ESTIMATE and image assessment, using image TP as the input for TSNet
#
rm(TempEstimate, TempEstimate_Pro, TempEstimate_RNA)

##### 1.2.2 run TSNet 
##### source
# source('./TSNet/TSNet_function.R') # see: https://github.com/petraf01/TSNet
# OutputPurity <- deNet_purity(exprM = t(DecovPro[,1:193]) %>% as.matrix(),
#                              purity = as.numeric(IndexDecov$TumorPurityInteg[1:193]))
# save(OutputPurity, file = 'PurityByTSNet.RData')
load('PurityByTSNet.RData')
# add to index
IndexDecov$TumorPurityTSNet <- c(OutputPurity$Epurity, rep(0,137))
cor.test(IndexDecov$TumoPurityImage, IndexDecov$TumorPurityTSNet)

##### 1.2.3 add stroma and immune score (based on proteomics)
TempxCell_Pro <- deconvo_tme(eset = DecovPro, method = "xcell", arrays = FALSE)
# add
IndexDecov$ImmuneScore <- sapply(IndexDecov$Sample, function(x){filter(TempxCell_Pro, ID == x)$ImmuneScore_xCell %>% as.numeric()}) %>% as.numeric()
IndexDecov$StromaScore <- sapply(IndexDecov$Sample, function(x){filter(TempxCell_Pro, ID == x)$StromaScore_xCell %>% as.numeric()}) %>% as.numeric()
# add score as sum cell frac
#rm 
rm(TempxCell_Pro, i)

##### 1.3 Making Input for BayesDebulk and run on server #####
##### 1.3.0 sourcing and func prep
# library(invgamma); library(truncnorm)
# source('./BayesDeBulk/BayesDeBulk.R')
# source('./BayesDeBulk/LM22_markers.R')
# load('./BayesDeBulk/LM22.rda')

##### 1.3.1 running LM22_marker list
# OutputMarkers <- LM22_markers(list(DecovRNA, DecovPro)) %>% as_tibble()
# add cell types: fibro and endo
# TempAddCellType <- c('Fibroblast','Endothelial','SmoothMuscle','Epithelial')
# TempAllCellType <- c(names(table(OutputMarkers$V1)), names(table(OutputMarkers$V2))) %>% unique()
# add old cell type marker higher than new cell type
# for (i in 1:length(TempAllCellType)) {
#   TempMarker <- filter(OutputMarkers, V1 == TempAllCellType[i])$marker.unique %>% unique()
#   TempOut <- expand.grid(V1 = TempAllCellType[i], 
#                          V2 = TempAddCellType,
#                          marker.unique = TempMarker) %>% as_tibble()
#   OutputMarkers = rbind(OutputMarkers, TempOut)
# }

# Fibroblast: CD36, PDGFRB, C5AR2, S100A4, CD70, PDPN, VIM, ITGA5, MME, PDGFRA, FAP, ACTA2,
# Endothelial: PECAM1, VEGFA, KDR, CD34, ITGB1, CD74
# Smooth muscle: MYH11, TPM1, MYL9, LMOD1, MYL6, GUCY1B1, GUCY1A1, MYL12B, MYL6B  
# Epithelial: CA2, AQP1, AKR1C4, UGT2B15, UGT2B10, ACOX2, FXYD2                   
# OutputMarkers <- rbind(OutputMarkers,
#                        expand.grid(V1 = 'Fibroblast', V2 = TempAllCellType,
#                                    marker.unique = c('CD36', 'PDGFRB', 'C5AR2', 'S100A4', 'CD70', 
#                                                      'PDPN', 'VIM', 'ITGA5', 'MME', 'PDGFRA', 'FAP', 'ACTA2')) %>% as_tibble(),
#                        expand.grid(V1 = 'Endothelial', V2 = TempAllCellType,
#                                    marker.unique = c('PECAM1', 'VEGFA', 'KDR', 'CD34', 'ITGB1', 'CD74')) %>% as_tibble(),
#                        expand.grid(V1 = 'SmoothMuscle', V2 = TempAllCellType,
#                                    marker.unique = c('MYH11', 'TPM1', 'MYL9', 'LMOD1', 'MYL6', 
#                                                      'GUCY1B1', 'GUCY1A1', 'MYL12B', 'MYL6B')) %>% as_tibble(),
#                        expand.grid(V1 = 'Epithelial', V2 = TempAllCellType,
#                                    marker.unique = c('CA2', 'AQP1', 'AKR1C4', 'UGT2B15', 'UGT2B10', 'ACOX2', 'FXYD2')) %>% as_tibble()) %>%
#   as.matrix()
# load('/Users/fuzile/Desktop/たいようけい/我参与的项目/主要参与/2023.08 GBC MO/6 写作/Code Ocean prep/Figure 5/BayesDB_All_On_Server.RData')
# write.xlsx(as.data.frame(OutputMarkers), file = 'Input markers for BayesDebulk.xlsx', rowNames = T, colNames = T, overwrite = T)

##### 1.3.2 running BayesDeBulk function [On server] and normalize
##### run on server
# OutputBDB <- BayesDeBulk(n.iter = 10000, burn.in = 1000,
#                          Y = list(t(scale(t(DecovRNA))) %>% as.data.frame(),
#                                   t(scale(t(DecovPro))) %>% as.data.frame()),
#                          markers = OutputMarkers)
# save(OutputBDB, file = 'BayesDeBulk_Output0311.RData')
# save(list = ls(), file = 'BayesDB_All_On_Server.RData')
load('BayesDeBulk_Output.RData') # downloaded from [server]
OutputBDB <- OutputBDB$cell.fraction %>% as.data.frame()
##### norm
# add purity
OutputBDB$Purity <- sapply(rownames(OutputBDB), function(x){filter(IndexDecov, Sample == x)$TumorPurityTSNet})
OutputBDB = OutputBDB[IndexDecov$Sample, ]
### norm to sum(1-TP)
OutputBDB_Normed_TP <- apply(OutputBDB[,1:26], 2, function(x){x*(1-OutputBDB$Purity)}) %>% t() %>% as.data.frame()
OutputBDB_Normed_TP = OutputBDB_Normed_TP[,IndexDecov$Sample]
# t
OutputBDB = as.data.frame(t(OutputBDB))

##### 1.4 Immune Module calc and clustering #####
##### 1.4.1 Signature processing 
SigImmune <- openxlsx::read.xlsx('./Immune Signatures.xlsx', sheet = 2, startRow = 5, rowNames = F, colNames = T) %>% as.data.frame()
SigImmune <- SigImmune[,c(1,3)]
SigImmune <- filter(SigImmune, !grepl('XCELL', Signature, ignore.case = T))
### separate
SigImmune <- separate_rows(SigImmune, Gene, sep = "/")
SigImmuneGSVA <- bioinfoamateur::dfm_change_colnames(SigImmune, c('term','gene')) %>%
  bioinfoamateur::Enrich_create_GSVA_ref_list() 
##### 1.4.2 calc score in RNA and Pro level
ImmModualRNA <- GSVA::gsva(as.matrix(DecovRNA), SigImmuneGSVA, method = 'ssgsea') %>% t() %>% as.data.frame()
ImmModualPro <- GSVA::gsva(as.matrix(DecovPro), SigImmuneGSVA, method = 'ssgsea') %>% t() %>% as.data.frame()
##### 1.4.3 calc distance
ImmModualRNA_Dis <- cor(ImmModualRNA, method = 'spearman')
ImmModualRNA_Dis = 1-ImmModualRNA_Dis
ImmModualPro_Dis <- cor(ImmModualPro, method = 'spearman')
ImmModualPro_Dis = 1-ImmModualPro_Dis
##### 1.4.4 clustering
# library(ConsensusClusterPlus)
# set.seed(0312)
# ImmModualClusterRNA <- ConsensusClusterPlus(d = ImmModualRNA_Dis, maxK = 10, reps = 10, pFeature = 0.8, distance="spearman",
#                                             clusterAlg = "pam", seed = 0312, title = 'ImmModualRNA', plot = 'png', verbose = T)
# ImmModualClusterPro <- ConsensusClusterPlus(d = ImmModualPro_Dis, maxK = 10, reps = 10, pFeature = 0.8, distance="spearman",
#                                             clusterAlg = "pam", seed = 0312, title = 'ImmModualPro', plot = 'png', verbose = T)
# save(ImmModualClusterRNA, ImmModualClusterPro,
#      file = 'ImmClusterSubtypingAndNaming.RData')

##### 1.4.5 extract cluster
load('ImmClusterSubtypingAndNaming.RData')
OutputImmuneModual <- tibble(Pathway = colnames(ImmModualRNA_Dis),
                             RNAModual = sapply(colnames(ImmModualRNA_Dis), function(x){as.numeric(ImmModualClusterRNA[[4]]$consensusClass[x])}),
                             RNAModualName = rep(1, 334),
                             ProModual = sapply(colnames(ImmModualPro_Dis), function(x){as.numeric(ImmModualClusterPro[[5]]$consensusClass[x])}),
                             ProModualName = rep(1, 334))
OutputImmuneModual$RNAModualName <- case_when(OutputImmuneModual$RNAModual == 1 ~ 'Lymphocytes',
                                              OutputImmuneModual$RNAModual == 2 ~ 'Myeloid',
                                              OutputImmuneModual$RNAModual == 3 ~ 'Stromal',
                                              OutputImmuneModual$RNAModual == 4 ~ 'Interferon')
OutputImmuneModual$ProModualName <- case_when(OutputImmuneModual$ProModual == 1 ~ 'Myeloid',
                                              OutputImmuneModual$ProModual == 2 ~ 'Stromal',
                                              OutputImmuneModual$ProModual == 3 ~ 'Lymphocytes1',
                                              OutputImmuneModual$ProModual == 4 ~ 'Lymphocytes2',
                                              OutputImmuneModual$ProModual == 5 ~ 'Interferon')
#
rm(ImmModualClusterPro, ImmModualClusterRNA, ImmModualPro_Dis, ImmModualRNA_Dis)

##### 1.5 Immune Module score calc #####
##### calc mean AND scale all (T + N)
# pro
ImmModualAvg_Pro <- as.data.frame(t(ImmModualPro)) %>% t() %>% scale() %>% t() %>% as.data.frame()
ImmModualAvg_Pro$Modual <- OutputImmuneModual$ProModualName
ImmModualAvg_Pro <- aggregate(ImmModualAvg_Pro[,-ncol(ImmModualAvg_Pro)], by = list(Modual = ImmModualAvg_Pro$Modual), mean) %>%
  column_to_rownames(var = 'Modual')
# rna
ImmModualAvg_RNA <- as.data.frame(t(ImmModualRNA)) %>% t() %>% scale() %>% t() %>% as.data.frame()
ImmModualAvg_RNA$Modual <- OutputImmuneModual$RNAModualName
ImmModualAvg_RNA <- aggregate(ImmModualAvg_RNA[,-ncol(ImmModualAvg_RNA)], by = list(Modual = ImmModualAvg_RNA$Modual), mean) %>%
  column_to_rownames(var = 'Modual')

##### separate T and N
### pro
ImmModualAvg_Pro_T <- ImmModualAvg_Pro[, filter(IndexDecov, Tumor == 'Yes')$Sample]
ImmModualAvg_Pro_N <- ImmModualAvg_Pro[, filter(IndexDecov, is.na(Tumor))$Sample]
Heatmap(ImmModualAvg_Pro_T, column_km = 4,
        col = colorRamp2(c(-2, 0, 2), c(ColColor$`Low-GreenBlue`[8], "white", ColColor$`High-OrangeRed`[8])))
Heatmap(ImmModualAvg_Pro_N,
        col = colorRamp2(c(-2, 0, 2), c(ColColor$`Low-GreenBlue`[8], "white", ColColor$`High-OrangeRed`[8])))
# mean
apply(ImmModualAvg_Pro_T, 1, sd); cat('\n');apply(ImmModualAvg_Pro_N, 1, sd)
### rna
ImmModualAvg_RNA_T <- ImmModualAvg_RNA[, filter(IndexDecov, Tumor == 'Yes', !is.na(RNA))$Sample]
ImmModualAvg_RNA_N <- ImmModualAvg_RNA[, filter(IndexDecov, is.na(Tumor), !is.na(RNA))$Sample]
Heatmap(ImmModualAvg_RNA_T, column_km = 4,
        col = colorRamp2(c(-2, 0, 2), c(ColColor$`Low-GreenBlue`[8], "white", ColColor$`High-OrangeRed`[8])))
Heatmap(ImmModualAvg_RNA_N,
        col = colorRamp2(c(-2, 0, 2), c(ColColor$`Low-GreenBlue`[8], "white", ColColor$`High-OrangeRed`[8])))
# mean
apply(ImmModualAvg_RNA_T, 1, sd, na.rm = T); cat('\n');apply(ImmModualAvg_RNA_N, 1, sd, na.rm = T)

##### supply the lacked sample in mRNA level
### func prep
ExtendMatrix <- function(OldDf, RefDf) {
  # find diff sample name
  TempSampleName = setdiff(colnames(RefDf), colnames(OldDf))
  # make new extended matrix
  ExtendedDf <- matrix(NA, nrow = nrow(OldDf), ncol = ncol(RefDf))
  ExtendedDf[, 1:ncol(OldDf)] <- as.matrix(OldDf)
  ExtendedDf <- as.data.frame(ExtendedDf)
  rownames(ExtendedDf) = rownames(OldDf)
  colnames(ExtendedDf) = c(colnames(OldDf), TempSampleName)
  # arrange sample name
  ExtendedDf = ExtendedDf[,colnames(RefDf)]
  #
  return(ExtendedDf)
}
#
ImmModualAvg_RNA_T = ExtendMatrix(ImmModualAvg_RNA_T, ImmModualAvg_Pro_T)
ImmModualAvg_RNA_N = ExtendMatrix(ImmModualAvg_RNA_N, ImmModualAvg_Pro_N)
# vis test with NAs (just try)
Heatmap(ImmModualAvg_RNA_N, cluster_columns = F, na_col = 'grey90',
        col = colorRamp2(c(-2, 0, 2), 
                         c(ColColor$`Low-GreenBlue`[8], "white", ColColor$`High-OrangeRed`[8])),
        height = unit(6,'cm'), width = unit(137/30,'cm'))

##### 1.6 Integrated clustering using BayesDeBulk output and Proteom-based immune modual score #####
library(CancerSubtypes)
# ##### Input prep
# ClusterInput_Cell <- OutputBDB_Normed_TP[,filter(IndexDecov, Tumor == 'Yes')$Sample]
# ClusterInput_Modual <- ImmModualPro[filter(IndexDecov, Tumor == 'Yes')$Sample, ] %>% t() %>% as.data.frame()
# # check dis
# data.checkDistribution(ClusterInput_Cell)
# dev.off()
# data.checkDistribution(ClusterInput_Modual)
# dev.off()
# # z-score norm
# # ClusterInput_Cell <- data.normalization(ClusterInput_Cell, type="feature_zscore", log2=FALSE)
# # ClusterInput_Modual <- data.normalization(ClusterInput_Modual, type="feature_zscore", log2=FALSE)
# ##### Clustering using km
# set.seed(0320)
# # 0313 ver
# ClusterOutput <- ExecuteCC(clusterNum = 4, d = list(CellFrac = ClusterInput_Cell,
#                                                     ImmModual = ClusterInput_Modual),
#                            maxK=10, reps = 50, clusterAlg = "km", title="ImmCluster")
# save(ClusterOutput, file = 'IntegClusterRawResult_0320.RData')
# load('IntegClusterRawResult.RData')
### check prognosis
# add
IndexDecov$OStime <- sapply(IndexDecov$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$OS_days/30 %>% as.numeric()})
IndexDecov$OS <- sapply(IndexDecov$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$OS_event %>% as.numeric()})
IndexDecov$OStime <- as.numeric(IndexDecov$OStime)
IndexDecov$OS <- as.numeric(IndexDecov$OS)
# check
survAnalysis(mainTitle="ImmCluster",
             IndexDecov$OStime[1:193],
             IndexDecov$OS[1:193],
             ClusterOutput$group,
             ClusterOutput$distanceMatrix,similarity=TRUE)
# add to index
IndexDecov$ImmCluster <- c(ClusterOutput$group, rep(NA, 137))

##### 1.6.1 prognosis #####
library(ggsurvfit)
survfit2(Surv(OStime, OS) ~ ImmCluster, data = IndexDecov[1:193,]) %>%
  ggsurvfit() +
  scale_ggsurvfit() + add_pvalue(caption = "Log-rank {p.value}") +
  add_risktable(risktable_group = "risktable_stats") +
  scale_ggsurvfit()


##### 1.6.2 paras with subtypes #####
plot_grid(ggboxplot(filter(IndexDecov, Tumor == 'Yes'),'ImmCluster','StromaScore') + 
            stat_compare_means(ref.group = 2, method = 't.test'),
          ggboxplot(filter(IndexDecov, Tumor == 'Yes'),'ImmCluster','ImmuneScore') + 
            stat_compare_means(ref.group = 2, method = 't.test'))

##### 1.7 retain important info #####
# save(IndexDecov, file = 'KeyVar_Imm.RData')
# load('KeyVar_Imm.RData')

##### 1.8 save all result (In the future, this part of analysis could start from this part) #####
# save(DecovPro, DecovRNA, ClusterInput_Cell, ClusterInput_Modual,
#      OutputBDB, OutputBDB_Normed_TP, ImmModualAvg_Pro, ImmModualPro,
#      IndexDecov,
#      file = 'AfterClusteringAll.RData')
# load('AfterClusteringAll.RData')

##### Figure 5A #####
VisIndexT <- IndexDecov %>% filter(Tumor == 'Yes') %>% arrange(ImmCluster)
##### Matrix
# cell and modual
ClusterInput_Cell <- OutputBDB_Normed_TP[,filter(IndexDecov, Tumor == 'Yes')$Sample]
VisInput_Cell <- ClusterInput_Cell[,VisIndexT$Sample] %>% t() %>% scale() %>% t() %>% as.data.frame()
VisInput_Modual <- ImmModualAvg_Pro_T[,VisIndexT$Sample] %>% t() %>% scale() %>% t() %>% as.data.frame()
VisInput_ModualRNA <- ImmModualAvg_RNA_T[,VisIndexT$Sample] %>% t() %>% scale() %>% t() %>% as.data.frame()
# markers
TempMarkers <- c('CD3E','CD8A','CD14','HLA-A','HLA-B','HLA-DRA','HLA-DQA1','GZMA','GZMB','PRF1','CCL5',
                 'CD274','TIGIT','TNFSF18','SIGLEC9','VSIR','VSIG4','IDO1','CEACAM1','CD276',   # tumor cell
                 'CTLA4','LAG3','HAVCR2','TNFRSF18','PDCD1','CD40','CD40LG')                    # immune cell
VisInput_MarkerRNA <- ExtendMatrix(DecovRNA[,1:134], DecovPro[,1:193])[TempMarkers, VisIndexT$Sample] %>% t() %>% scale %>% t() 
VisInput_MarkerRNA <- VisInput_MarkerRNA[rownames(VisInput_MarkerRNA) != 'NA',]
VisInput_MarkerPro <- DecovPro[TempMarkers, VisIndexT$Sample] %>% t() %>% scale %>% t() %>% na.omit()
##### Anno
### add info
# naming cluster
VisIndexT$ImmCluster_Name = case_when(VisIndexT$ImmCluster == 1 ~ 'A',      # immune - stroma - 
                                      VisIndexT$ImmCluster == 2 ~ 'B',      # immune + stroma - hot
                                      VisIndexT$ImmCluster == 3 ~ 'C',      # immune + stroma + 
                                      VisIndexT$ImmCluster == 4 ~ 'D') %>%  # immune - stroma + 
  factor(levels = c('A','B','C','D'), ordered = T)
## pathologic type
VisIndexT$pathology = sapply(VisIndexT$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Pathological_type}) %>% as.character()
VisIndexT$pathology = case_when(VisIndexT$pathology == 'adenocarcinoma' ~ 'AC',
                                VisIndexT$pathology == 'adenosquanmous carcinoma' ~ 'AS',
                                VisIndexT$pathology == 'neuroendocrine carcinoma' ~ 'NE',
                                !is.na(VisIndexT$pathology) ~ 'Others') %>%
  factor(levels = c('AC','AS','NE','Others'), ordered = T)
#
library(wesanderson)
ColPatho2 = c(wes_palette('Darjeeling1', 3, type = c("discrete")), 'grey70')
names(ColPatho2) = c('AC','AS','NE','Others')
# cholelithiasis
VisIndexT$Cholelithiasis = sapply(VisIndexT$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$cholelithiasis}) %>% as.numeric()
VisIndexT$Cholelithiasis = ifelse(VisIndexT$Cholelithiasis == 1, 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)
# LI
VisIndexT$LiverInvasion = sapply(VisIndexT$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Regional_invasion}) %>% as.numeric()
VisIndexT$LiverInvasion = ifelse(VisIndexT$LiverInvasion == 1, 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)
# Distal metastasis
VisIndexT$DistalMetastasis = sapply(VisIndexT$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Distal_metastasis}) %>% as.numeric()
VisIndexT$DistalMetastasis = ifelse(VisIndexT$DistalMetastasis == 1, 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)
# TMB
GBC_Main_Maf <- maftools::read.maf(maf = GBC_Main_Mutation)
TempTMB <- tmb(GBC_Main_Maf)
VisIndexT$TMB <- sapply(VisIndexT$Sample, function(x){filter(TempTMB, Tumor_Sample_Barcode == x)$total_perMB %>% as.numeric()}) %>% as.numeric()
rm(TempTMB)
# TP53 mut
VisIndexT$TP53Mut <- case_when(VisIndexT$Sample %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'TP53')$Tumor_Sample_Barcode ~ 'Yes',
                               VisIndexT$Sample %in% GBC_Main_Mutation$Tumor_Sample_Barcode ~ 'No',
                               !is.na(VisIndexT$Sample) ~ NA) %>%
  factor(levels = c('Yes','No'), ordered = T)
# ERBB2 Amp
# load('/Users/fuzile/Desktop/たいようけい/我参与的项目/主要参与/2023.08 GBC MO/3 新结果/0923 -/队列总体描述/20231209分型结果命名.RData')
VisIndexT$ERBB2Amp <- sapply(VisIndexT$Sample, function(x){filter(Input, Patient == x)$CNV_ERBB2}) %>% as.character()
VisIndexT$ERBB2Amp = ifelse(VisIndexT$ERBB2Amp == 'Amplification', 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)
rm(Input)
### make anno
AnnoSample <- columnAnnotation(`Immune clusters` = VisIndexT$ImmCluster_Name,
                               `Tumor purity` = VisIndexT$TumorPurityTSNet,
                               `Immune score` = VisIndexT$ImmuneScore,
                               `Stroma score` = VisIndexT$StromaScore,
                               Pathology = VisIndexT$pathology,
                               Cholelithiasis = VisIndexT$Cholelithiasis,
                               `Liver invasion` = VisIndexT$LiverInvasion,
                               `Distal metastasis` = VisIndexT$DistalMetastasis,
                               TMB = VisIndexT$TMB,
                               `TP53 Mutation` = VisIndexT$TP53Mut,
                               `ERBB2 Amplificatoin` = VisIndexT$ERBB2Amp,
                               col = list(
                                 `Immune clusters` = ColImmune,
                                 `Tumor purity` = colorRamp2(c(0, 0.4, 0.8), c("white", ColColor$`High-Orange`[4],  ColColor$`High-Orange`[9])),
                                 `Stroma score` = colorRamp2(c(0, 0.2, 0.4), c("white", ColColor$`Single-Brown`[4],  ColColor$`Single-Brown`[9])),
                                 `Immune score` = colorRamp2(c(0, 0.15, 0.3), c("white", ColColor$`High-Orange`[4],  ColColor$`High-Orange`[9])),
                                 Pathology = ColPatho2,
                                 Cholelithiasis = c('Yes' = 'Black', 'No' = 'white'),
                                 `Liver invasion` = c('Yes' = 'Black', 'No' = 'white'),
                                 `Distal metastasis` = c('Yes' = 'Black', 'No' = 'white'),
                                 TMB = colorRamp2(c(0, 3, 6), c("white", ColColor$`Low-Indigo`[4],  ColColor$`Low-Indigo`[9])),
                                 `TP53 Mutation` = c('Yes' = 'Black', 'No' = 'white'),
                                 `ERBB2 Amplificatoin` = c('Yes' = 'Black', 'No' = 'white')
                               ),
                               simple_anno_size = unit(0.3333, "cm"),
                               gp = gpar(col = "white", lwd = 0.2))

##### Vis
# cell and modual
HmCell <- Heatmap(VisInput_Cell, cluster_columns = F, name = 'Z-scored\nNormalized cell fractions',
                  col = colorRamp2(c(-2.5,0,2.5), c(ColColor$`Low-Blue`[9], 'white',
                                                    ColColor$`High-YellowOrange`[9])),
                  column_split = VisIndexT$ImmCluster, show_column_names = F, top_annotation = AnnoSample,
                  column_title = NULL, column_gap = unit(2,'mm'),
                  height = unit(26/2.5,'cm'), width = unit(193/15,'cm'))
HmModual <- Heatmap(VisInput_Modual, cluster_columns = F, name = 'Z-scored\nFeature levels',
                    column_split = VisIndexT$ImmCluster, show_column_names = F, column_gap = unit(2,'mm'),
                    col = colorRamp2(c(-3,0,3), c(ColColor$`Low-LightBlue`[9],'white',ColColor$`High-Red`[9])),
                    height = unit(5/2.5,'cm'), width = unit(193/15,'cm'))
HmModual_RNA <- Heatmap(VisInput_ModualRNA, cluster_columns = F, show_heatmap_legend = F,
                        column_split = VisIndexT$ImmCluster, show_column_names = F, column_gap = unit(2,'mm'),
                        col = colorRamp2(c(-3,0,3), c(ColColor$`Low-LightBlue`[9],'white',ColColor$`High-Red`[9])),
                        height = unit(6/2.5,'cm'), width = unit(193/15,'cm'))
# marker
HmMarker_RNA <- Heatmap(VisInput_MarkerRNA, cluster_columns = F, show_heatmap_legend = F,
                        column_split = VisIndexT$ImmCluster, show_column_names = F, column_gap = unit(2,'mm'),
                        col = colorRamp2(c(-3,0,3), c(ColColor$`Low-LightBlue`[9],'white',ColColor$`High-Red`[9])),
                        height = unit(26/2.5,'cm'), width = unit(193/15,'cm'))
HmMarker_Pro <- Heatmap(VisInput_MarkerPro, cluster_columns = F, show_heatmap_legend = F,
                        column_split = VisIndexT$ImmCluster, show_column_names = F, column_gap = unit(2,'mm'),
                        col = colorRamp2(c(-3,0,3), c(ColColor$`Low-LightBlue`[9],'white',ColColor$`High-Red`[9])),
                        height = unit(14/2.5,'cm'), width = unit(193/15,'cm'))
draw(HmCell %v% HmModual %v% HmModual_RNA %v% HmMarker_Pro %v% HmMarker_RNA)
draw(HmCell %v% HmModual %v% HmMarker_Pro,
     heatmap_legend_side = 'left', annotation_legend_side = 'left') 
# data for Table S6
write.xlsx(filter(IndexDecov, Tumor == 'Yes'), file = 'Metadata immune.xlsx', rowNames = T, colNames = T, overwrite = T)
write.xlsx(VisInput_Cell, file = 'immune cell.xlsx', rowNames = T, colNames = T, overwrite = T)
write.xlsx(VisInput_Modual, file = 'immune module.xlsx', rowNames = T, colNames = T, overwrite = T)
#
rm(HmCell, HmMarker_Pro, HmMarker_RNA, HmModual, HmModual_RNA, AnnoSample,
   VisInput_ModualRNA, VisIndexT, VisInput_Cell, VisInput_MarkerPro, VisInput_MarkerRNA, VisInput_Modual,
   ImmModualAvg_Pro_N, ImmModualAvg_Pro_T, ImmModualAvg_RNA_N, ImmModualAvg_RNA_T, ImmModualAvg_RNA,
   ImmModualPro_Dis, ImmModualRNA, ImmModualRNA_Dis,
   OutputBDB_Normed, OutputBDB_Normed_P, OutputBDB_Normed_T, OutputBDB_Normed_TP_P, OutputBDB_Normed_TP_T,
   ClusterInput_Cell, ClusterInput_Modual,
   IndexOld, Output, OutputPurity, Temp,
   Output, SigImmune, SigImmuneGSVA, TempMarkers, ClusterOutput)

########## Variables used for the following analysis ##########
##### Index
# old
IndexDecov
# new
Index <- IndexDecov
# Tumor only
IndexT <- filter(Index, Tumor == 'Yes')

##### Cell frac
OutputBDB = as.data.frame(t(OutputBDB))# without norm by tumor purity
OutputBDB_Normed_TP

##### Immune modual
ImmModualPro
ImmModualAvg_Pro

##### Figure S5A #####
SelectedCells = c('Dendritic.cells.activated','Macrophages.M1',
                  'Plasma.cells','T.cells.CD8','T.cells.regulatory..Tregs.',
                  'NK.cells.activated','Fibroblast','Endothelial','Neutrophils')
Input = tibble(Sample = Index$Sample,
               Group = c(as.character(Index$ImmCluster_Name[1:193]), rep('NAT',137))) %>%
  cbind(as_tibble(t(OutputBDB_Normed_TP[SelectedCells,Index$Sample]))) %>%
  as_tibble()
colnames(Input)[3:11] = c('Activated DCs','M1 Macrophages','Plasma cells','CD8+T cells','Tregs',
                          'Activated NKs','Fibroblasts','Endothelial','Neutrophils')
##### transform
Input <- reshape2::melt(Input, id.vars = c('Sample', 'Group'),
                        variable.name = 'Cell', value.name = 'Fraction') %>%
  as_tibble() %>%
  mutate(Group = factor(Group, levels = c('A','B','C','D','NAT'), ordered = T))
##### Vis
ggboxplot(Input, "Cell", "Fraction", color = "Group") +
  scale_color_manual(values = c(ColImmune, 'NAT' = ColJournal$Science[3])) +
  stat_compare_means(method = 'anova') + 
  xlab('') + ylab('Normalized cell fraction') 

##### Figuer S5B #####
Input = tibble(Sample = IndexT$Sample,
               OStime = IndexT$OStime, OS = IndexT$OS,
               ImmuneScore = IndexT$ImmuneScore, StromaScore = IndexT$StromaScore,
               TumorPurity = IndexT$TumorPurityTSNet)
Input = cbind(Input, as_tibble(t(OutputBDB_Normed_TP[,IndexT$Sample]))) %>% as_tibble()
##### output prep of Uni-Cox model HR and Pval
Output = tibble(expand.grid(Variables = colnames(Input)[4:32],
                            Cohort = c('All','Immunotherapy','Chemotherapy'),
                            HR = 1,
                            HR_Pval = 1,
                            LogRank_Pval = 1))
##### calc
library(survival); library(coin)
for (i in 1:nrow(Output)) {
  TempV = as.character(Output$Variables[i])
  TempDf <- cbind(Input[, c(1:3)], as.data.frame(Input[,TempV]))
  colnames(TempDf)[4] = 'V'
  if (as.character(Output$Cohort[i]) == 'All') {
    TempDf$HiLow = ifelse(TempDf$V >= median(TempDf$V), 'Hi', 'Low') %>% factor()
    # COX
    TempCOXModel <- coxph(Surv(OStime, OS) ~ V, data = TempDf)
    TempCOXModel <- summary(TempCOXModel)
    # Log-Rank
    TempLR <- logrank_test(Surv(OStime, OS) ~ HiLow, data = TempDf, type = 'logrank')
  }
  if (as.character(Output$Cohort[i]) == 'Immunotherapy') {
    TempDf <- filter(TempDf, Sample %in% filter(GBC_Main_Clinical, `Adjuvant.therapies_1` == 'immunotherapy')$Patient_ID)
    TempDf$HiLow = ifelse(TempDf$V >= median(TempDf$V), 'Hi', 'Low') %>% factor()
    #
    TempCOXModel <- coxph(Surv(OStime, OS) ~ V, data = TempDf)
    TempCOXModel <- summary(TempCOXModel)
    #
    TempLR <- logrank_test(Surv(OStime, OS) ~ HiLow, data = TempDf, type = 'logrank')
  }
  if (as.character(Output$Cohort[i]) == 'Chemotherapy') {
    TempDf <- filter(TempDf, Sample %in% filter(GBC_Main_Clinical, `Adjuvant.therapies_1` == 'chemotherapy')$Patient_ID)
    TempDf$HiLow = ifelse(TempDf$V >= median(TempDf$V), 'Hi', 'Low') %>% factor()
    #
    TempCOXModel <- coxph(Surv(OStime, OS) ~ V, data = TempDf)
    TempCOXModel <- summary(TempCOXModel)
    #
    TempLR <- logrank_test(Surv(OStime, OS) ~ HiLow, data = TempDf, type = 'logrank')
  }
  Output$HR[i] = as.numeric(TempCOXModel$coefficients[1,2])
  Output$HR_Pval[i] = as.numeric(TempCOXModel$coefficients[1,5])
  Output$LogRank_Pval[i] = pvalue(TempLR)
}
rm(i, TempDf, TempV, TempCOXModel, TempLR)
##### Vis-HR bubble
library(ggbreak)
Output$LogHR = log10(Output$HR)
Output$LogHR_FDR = -log10(Output$HR_Pval)
Output$Col_HR = case_when(Output$HR > 1 &  Output$HR_Pval <= 0.05 ~ 'Harm',
                          Output$HR < 1 &  Output$HR_Pval <= 0.05 ~ 'Protection',
                          !is.na(Output$HR) ~ 'Others')
Output = filter(Output, Variables != 'Epithelial')
plot_grid(ggscatter(filter(Output, Cohort == "Immunotherapy"), 'LogHR', 'LogHR_FDR', size = 'LogHR_FDR', color = 'Col_HR') +
            scale_color_manual(values = c('Harm' = 'red', 'Protection' = 'lightblue', 'Others' = 'grey90')) +
            xscale("log2", .format = TRUE) +
            ggrepel::geom_text_repel(data = filter(Output, Cohort == "Immunotherapy", Col_HR == 'Harm'),aes(label = Variables),
                                     box.padding = 0.5, max.overlaps = Inf) +
            theme(legend.position = 'right') + ggtitle('Adjuvant immunotherapy (N=34)'),
          ggscatter(filter(Output, Cohort == "Chemotherapy"), 'LogHR', 'LogHR_FDR', size = 'LogHR_FDR', color = 'Col_HR') +
            scale_color_manual(values = c('Harm' = 'red', 'Protection' = 'lightblue', 'Others' = 'grey90')) +
            xscale("-log2", .format = TRUE) +
            ggrepel::geom_text_repel(data = filter(Output, Cohort == "Chemotherapy", Col_HR == 'Protection'),aes(label = Variables),
                                     box.padding = 0.5, max.overlaps = Inf) +
            theme(legend.position = 'right') + ggtitle('Adjuvant Chemotherapy (N=115)'),
          nrow = 2) 
# data for Table S6
write.xlsx(Output, file = 'HR.xlsx', rowNames = T, colNames = T, overwrite = T)
#
rm(Output, Input)


##### Figure 5B #####
CM <- data.frame(A = ifelse(IndexT$ImmCluster_Name == 'A', 1, 0),
                 B = ifelse(IndexT$ImmCluster_Name == 'B', 1, 0),
                 C = ifelse(IndexT$ImmCluster_Name == 'C', 1, 0),
                 D = ifelse(IndexT$ImmCluster_Name == 'D', 1, 0))
PathwayPro <- GSVA::gsva(as.matrix(GBC_Main_Pro[,IndexT$Pro]), TotalPathwayGSVA, method = 'gsva')
## pro
DEA_PathPro <- bioinfoamateur::core_Differential_analysis_continuous(CM, PathwayPro,
                                                     log = T, method = 't.test', p.adj = T, show.belong = T)
DEA_PathPro$Pathway = rownames(DEA_PathPro)

##### Vis prep
SelectPathway <- tibble(
  Pathway = c(
    # cc
    'KEGG Cell cycle','KEGG DNA replication','HALLMARK_G2M_CHECKPOINT', # C1
    # metab
    'HALLMARK_FATTY_ACID_METABOLISM','KEGG Oxidative phosphorylation', # C1
    'HALLMARK_GLYCOLYSIS','KEGG N-Glycan biosynthesis','KEGG Lysine degradation', # C4
    # biologic process
    'KEGG Apoptosis','KEGG Autophagy',  # C2
    'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', # C3
    'HALLMARK_HYPOXIA','HALLMARK_ANGIOGENESIS', # C4
    # immune
    'REACTOME Antigen Presentation','REACTOME Signaling by Interleukins','REACTOME Cytokine Signaling in Immune system',
    'REACTOME Innate Immune System','REACTOME Adaptive Immune System','REACTOME Interferon gamma signaling',   # C2
    'HALLMARK_COMPLEMENT', # C4
    # signalings
    'REACTOME Signaling by NOTCH','REACTOME Signaling by Hedgehog','REACTOME Signaling by FGFR',   # C1
    'KEGG ErbB signaling pathway','KEGG NF-kappa B signaling pathway','KEGG VEGF signaling pathway','HALLMARK_PI3K_AKT_MTOR_SIGNALING', # C2
    'HALLMARK_TGF_BETA_SIGNALING','KEGG Hippo signaling pathway','REACTOME EPH-Ephrin signaling',   # C3
    'REACTOME Signaling by PDGF','HALLMARK_WNT_BETA_CATENIN_SIGNALING','KEGG AMPK signaling pathway','KEGG HIF-1 signaling pathway', # C4
    # cell junction
    'KEGG Tight junction','REACTOME Cell junction organization','KEGG Focal adhesion', # C3
    # ECM degradation
    'KEGG ECM-receptor interaction','REACTOME Activation of Matrix Metalloproteinases',
    'REACTOME Degradation of the extracellular matrix' # C4
  ),
  Category = c(rep('Cell cycle', 3), rep('Metabolism', 5), rep('Biologic processes', 5),
               rep('Immune response', 7), rep('Signalings', 14), rep('Cell junction', 3), rep('ECM', 3)))
##### Vis input prep
### value matrix
Input <- PathwayPro[SelectPathway$Pathway, IndexT$Pro] %>% t() %>% as.data.frame()
Input = aggregate(Input, by = list(Subtype = IndexT$ImmCluster_Name), mean) %>%
  column_to_rownames(var = 'Subtype') %>% t() %>% as.data.frame() %>%
  t() %>% scale() %>% t() %>% as.data.frame()
### sig matrix
InputSig <- Input
InputSig[!is.na(InputSig)] = 0
for(i in 1:nrow(InputSig)) {
  TempVec = DEA_PathPro[rownames(InputSig)[i], ] %>% as.vector()
  InputSig[i, 1] = ifelse(TempVec[1] <= 0.05, '*', NA)
  InputSig[i, 2] = ifelse(TempVec[3] <= 0.05, '*', NA)
  InputSig[i, 3] = ifelse(TempVec[5] <= 0.05, '*', NA)
  InputSig[i, 4] = ifelse(TempVec[7] <= 0.05, '*', NA)
}
##### Vis prep
# Anno
AnnoCluster <- columnAnnotation(Subtype = c('A','B','C','D'),
                                col = list(
                                  Subtype = ColImmune
                                ),
                                gp = gpar(col = "white"),
                                annotation_name_side = "left",
                                simple_anno_size = unit(0.5, "cm"))
AnnoPathway <- rowAnnotation(Category = SelectPathway$Category,
                             col = list(
                               Category = c('Cell cycle' = ColJournal$COSMICsignature[1],
                                            'Biologic processes' = ColJournal$COSMICsignature[2],
                                            'ECM' = ColJournal$COSMICsignature[3],
                                            'Cell junction' = ColJournal$COSMICsignature[4],
                                            'Immune response' = ColJournal$COSMICsignature[5],
                                            'Metabolism' = ColJournal$COSMICsignature[6],
                                            'Signalings' = ColJournal$Science[8])
                             ),
                             gp = gpar(col = "white"),
                             # annotation_name_side = "left",
                             simple_anno_size = unit(0.5, "cm"))
##### Vis
set.seed(0317)
Heatmap(Input, cluster_columns = F, name = 'NES', 
        col = colorRamp2(c(-1.5,-1,0,1,1.5), c(ColColor$`Low-BlueGreen`[9], ColColor$`Low-BlueGreen`[4], 
                                               "white", 
                                               ColColor$`High-Orange`[2], ColColor$`High-Orange`[9])),
        top_annotation = AnnoCluster, left_annotation = AnnoPathway,
        column_split = c('A','B','C','D'), column_title = NULL, 
        row_split = SelectPathway$Category, row_title = NULL, cluster_rows = F, 
        rect_gp = gpar(col = "white", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(!is.na(InputSig[i, j]))
            grid.text(InputSig[i, j], x, y, gp = gpar(fontsize = 15), just = c("centre","center"))
        },
        height = unit(40/2,'cm'), width = unit(4,'cm')) %>%
  draw(heatmap_legend_side = 'left', annotation_legend_side = 'left')
# data for Table S6
write.xlsx(Input, file = 'Immune cluster function.xlsx', rowNames = T, colNames = T, overwrite = T)
#
rm(SelectPathway, Input, InputSig, AnnoCluster, AnnoPathway, i)

##### Figure 5C #####
library(ggsurvfit)
survfit2(Surv(OStime, OS) ~ ImmCluster_Name, data = filter(Index, Tumor == 'Yes')) %>%
  ggsurvfit(lwd = 1) +
  add_quantile(y_value = 0.5, linetype = "dotted", color = "grey30", 
               linewidth = 0.5) + 
  scale_color_manual(values = ColImmune) +
  scale_ggsurvfit() + 
  annotate('text', x = 40, y = 0.75, label = 'Log-Rank P = 0.042',
           color = 'black', family = 'arial') +
  xlab('Time (momth)') + ylab('Overall survival') +
  theme(axis.text = element_text(colour = 'black', family = 'arial'),
        legend.position = 'right') 

##### Figure 5D #####
CM <- data.frame(C = ifelse(filter(IndexT, ImmCluster_Name %in% c('C','D'))$ImmCluster_Name == 'C', 1, 0),
                 D = ifelse(filter(IndexT, ImmCluster_Name %in% c('C','D'))$ImmCluster_Name == 'D', 1, 0))
DEA_Pro <- bioinfoamateur::core_Differential_analysis_continuous(CM, DecovPro[,filter(IndexT, ImmCluster_Name %in% c('C','D'))$Sample], 
                                                 log = T, p.adj = T)
DEA_Pro$Gene = rownames(DEA_Pro)
DEA_Phos <- bioinfoamateur::core_Differential_analysis_continuous(CM, GBC_Main_ProLevelPhos_knn[,filter(IndexT, ImmCluster_Name %in% c('C','D'))$Pro], 
                                                  log = T, p.adj = T)
DEA_Phos$Gene = rownames(DEA_Phos)
##### GSEA
SelectGene = tibble(Gene = c('ITGAM','ITGB2','ITGA2B','ITGA3',
                             'CMAS','GNE','ST3GAL1','C1GALT1C1',
                             'IDO1','VSIR','NT5E','CD276'),
                    Category = c(rep('Integrin',4), rep('N-Glycan',4), rep('Checkpoints',4)))
##### Input prep
InputCD <- tibble(Sample = filter(IndexT, ImmCluster_Name %in% c('C','D'))$Sample,
                  Pro = filter(IndexT, ImmCluster_Name %in% c('C','D'))$Pro,
                  Subtype = filter(IndexT, ImmCluster_Name %in% c('C','D'))$ImmCluster_Name)
InputCD1 <- DecovPro[filter(SelectGene, Category == 'Integrin')$Gene,InputCD$Sample] %>% 
  t() %>% as.data.frame() %>% mutate(Subtype = filter(IndexT, ImmCluster_Name %in% c('C','D'))$ImmCluster_Name) %>%
  reshape2::melt(id.var = 'Subtype', variable.name = 'Gene', value.name = 'Abundance') %>% as_tibble()
InputCD2 <- DecovPro[filter(SelectGene, Category == 'N-Glycan')$Gene,InputCD$Sample] %>% 
  t() %>% as.data.frame() %>% mutate(Subtype = filter(IndexT, ImmCluster_Name %in% c('C','D'))$ImmCluster_Name) %>%
  reshape2::melt(id.var = 'Subtype', variable.name = 'Gene', value.name = 'Abundance') %>% as_tibble()
InputCD3 <- DecovPro[filter(SelectGene, Category == 'Checkpoints')$Gene,InputCD$Sample] %>% 
  t() %>% as.data.frame() %>% mutate(Subtype = filter(IndexT, ImmCluster_Name %in% c('C','D'))$ImmCluster_Name) %>%
  reshape2::melt(id.var = 'Subtype', variable.name = 'Gene', value.name = 'Abundance') %>% as_tibble()
##### Vis
plot_grid(ggboxplot(InputCD1, 'Gene', 'Abundance', color = 'Subtype') + 
            scale_color_manual(values = ColImmune) +
            stat_pvalue_manual(compare_means(
              Abundance ~ Subtype, data = InputCD1, group.by = "Gene",
              method = "t.test"), x = "Gene", y.position = 3.25,
              label = "p.signif",
              position = position_dodge(0.8)) +
            xlab('') + ggtitle('Integrins') +
            theme(legend.position = 'none'),
          ggboxplot(InputCD2, 'Gene', 'Abundance', color = 'Subtype') + 
            scale_color_manual(values = ColImmune) +
            stat_pvalue_manual(compare_means(
              Abundance ~ Subtype, data = InputCD2, group.by = "Gene",
              method = "t.test"), x = "Gene", y.position = 3.25,
              label = "p.signif",
              position = position_dodge(0.8)) +
            xlab('') + ggtitle('N-Glycan related') +
            theme(legend.position = 'none'),
          ggboxplot(InputCD3, 'Gene', 'Abundance', color = 'Subtype') + 
            scale_color_manual(values = ColImmune) +
            stat_pvalue_manual(compare_means(
              Abundance ~ Subtype, data = InputCD3, group.by = "Gene",
              method = "t.test"), x = "Gene", y.position = 3.25,
              label = "p.signif",
              position = position_dodge(0.8)) +
            xlab('') + ggtitle('Immune checkpoints') +
            theme(legend.position = 'none'),
          ncol = 1)  
#
rm(InputCD, InputCD1, InputCD2, InputCD3, SelectGene, DEA_Pro, DEA_Phos)

##### Figure S5C #####
library(ggnewscale)
library(ggside)
##### Input
Input <- tibble(Sample = filter(IndexDecov, Tumor == 'Yes')$Sample,
                `ImmuneScore` = scale(filter(IndexDecov, Tumor == 'Yes')$ImmuneScore %>% as.numeric()) %>% as.numeric(),
                # `StromaScore` = scale(filter(IndexDecov, Tumor == 'Yes')$StromaScore %>% as.numeric()) %>% as.numeric(),
                StromaScore = filter(IndexDecov, Tumor == 'Yes')$StromaScore %>% as.numeric(),
                TumorPurity = filter(IndexDecov, Tumor == 'Yes')$TumorPurityTSNet,
                Subtype = filter(IndexDecov, Tumor == 'Yes')$ImmCluster_Name
)
### use rank
Input$ImmRank <- rank(Input$ImmuneScore)
Input$StrRank <- rank(Input$StromaScore)
##
mean(Input$ImmuneScore)
mean(Input$StromaScore)
##### Vis
# Col
scales::show_col(ColJournal$Science)
ColImmune <- c(ColJournal$Nature[5], ColJournal$Nature[1],
               ColJournal$Nature[4], ColJournal$Nature[9])
names(ColImmune) <- c('A','B','C','D')
# Vis
ggplot(Input, aes(x = ImmRank, y = StrRank, color = TumorPurity)) +
  geom_point() +
  scale_x_continuous(limits = c(-200,400)) + 
  # scale_y_continuous(limits = c(-50,250)) +
  scale_color_gradient(low = ColColor$`High-Orange`[1], high = ColColor$`High-Orange`[9]) +
  theme_classic2() + 
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'left') + 
  geom_vline(xintercept = 108, linetype = 2, lwd = 1, color = 'grey') +
  geom_hline(yintercept = 111, linetype = 2, lwd = 1, color = 'grey') +
  new_scale_color() + new_scale_fill() +
  stat_ellipse(aes(fill = Subtype), geom = "polygon", alpha = 0.3, level = 0.75) +
  geom_xsideboxplot(aes(y = Subtype, color = Subtype), orientation = "y") +
  geom_ysideboxplot(aes(x = Subtype, color = Subtype), orientation = "x") +
  scale_fill_manual(values = ColImmune) +
  scale_color_manual(values = ColImmune) +
  xlab('Immune score (Ranking)') + ylab('Stroma score (Ranking)')
# Vis 0510
ggplot(Input, aes(x = ImmRank, y = StrRank, color = Subtype)) +
  geom_point() +
  scale_x_continuous(limits = c(-200,400)) + 
  # scale_y_continuous(limits = c(-50,250)) +
  scale_color_manual(values = ColImmune) + 
  theme_classic2() + 
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'left') + 
  geom_vline(xintercept = 108, linetype = 2, lwd = 1, color = 'grey') +
  geom_hline(yintercept = 111, linetype = 2, lwd = 1, color = 'grey') +
  new_scale_color() + new_scale_fill() +
  # stat_ellipse(aes(fill = Subtype), geom = "polygon", alpha = 0.3, level = 0.75) +
  geom_xsideboxplot(aes(y = Subtype, color = Subtype), orientation = "y") +
  geom_ysideboxplot(aes(x = Subtype, color = Subtype), orientation = "x") +
  scale_fill_manual(values = ColImmune) +
  scale_color_manual(values = ColImmune) +
  xlab('Immune score (Ranking)') + ylab('Stroma score (Ranking)') 
#
rm(Input)
# A: Imm- Str- Tumor+ 
# B: Imm+ Str- Tumor-
# C: Imm± Str+ Tumor-
# D: Imm+ Str+ Tumor+

##### Figure S5D #####
Input <- tibble(Sample = IndexT$Sample,
                WES = IndexT$WES,
                Pro = IndexT$Pro,
                Subtype = IndexT$ImmCluster_Name)
##### add TMB
GBC_Main_Maf <- maftools::read.maf(maf = GBC_Main_Mutation)
TempTMB <- tmb(GBC_Main_Maf)
Input$TMB <- sapply(Input$Sample, function(x){filter(TempTMB, Tumor_Sample_Barcode == x)$total_perMB %>% as.numeric()}) %>% as.numeric()
rm(TempTMB)

##### add CIN
InputCIN = GBC_Main_CNV
InputCIN$ChrBand = sapply(rownames(InputCIN), function(x){
  Band = filter(GeneCoord, gene_name == x)$ChrBand
  Band = ifelse(length(Band) == 0, NA, as.character(Band))
})
InputCIN = na.omit(InputCIN)
# calc band mean and band weight
InputCIN = aggregate(InputCIN[,-ncol(InputCIN)], by = list(ChrBand = InputCIN$ChrBand), mean) %>%
  column_to_rownames(var = 'ChrBand')
InputCalcIndex = tibble(ChrBand = rownames(InputCIN),
                        Chr = paste0('chr',str_split_fixed(rownames(InputCIN), '[pq]', 2)[,1]))
InputCalcIndex$BandLength = sapply(InputCalcIndex$ChrBand, function(x){filter(Chr, ChrBand == x)$Length})
InputCalcIndex$ChrLength = sapply(InputCalcIndex$Chr, function(x){filter(Chr, Chr == x)$Length %>% sum()})
InputCalcIndex$BandW = InputCalcIndex$BandLength / InputCalcIndex$ChrLength
# Calc CIN
Output = tibble(Sample = names(apply(InputCIN, 2, function(x){sum(abs(x) * InputCalcIndex$BandW)})),
                CIN = apply(InputCIN, 2, function(x){sum(abs(x) * InputCalcIndex$BandW)}))
# add
Input$CIN = sapply(Input$WES, function(x){filter(Output, Sample == x)$CIN}) %>% as.numeric()
# rm
rm(InputCIN, InputCalcIndex, Output)
##### Vis
plot_grid(ggboxplot(Input, 'Subtype', 'TMB', color = 'Subtype', add = 'jitter') +
            scale_color_manual(values = ColImmune) +
            stat_compare_means(method = 'anova') + xlab('') +
            theme(legend.position = 'none'),
          ggboxplot(Input, 'Subtype', 'CIN', color = 'Subtype', add = 'jitter') +
            scale_color_manual(values = ColImmune) +
            stat_compare_means(method = 'anova') + xlab('') +
            theme(legend.position = 'none') + 
            stat_pvalue_manual(compare_means(CIN ~ Subtype, data = Input, method = "t.test") %>%
                                 mutate(y.position = c(19:24)), label = "p.signif"),
          nrow = 1) 
#
rm(Input, Output)

##### Figure 5E #####
Input18q <- InputBandAvg %>%
  mutate(Arm = paste0(str_split_fixed(rownames(InputBandAvg),'[pq]', 2)[,1] %>% as.numeric(),
                      ifelse(grepl('p', rownames(InputBandAvg)), 'p', 'q')))
# agg
Input18q = aggregate(Input18q[,-ncol(Input18q)], by = list(Arm = Input18q$Arm), mean)
# 
Input18q = reshape2::melt(Input18q, id.var = c('Arm'), variable.name = 'Sample', value.name = 'ArmCNV') %>%
  as_tibble() %>% filter(Arm == '18q')
Input18q$Subtype = sapply(Input18q$Sample, function(x){filter(IndexGenome, WES == x)$ImmCluster_Name})
##### Vis
ggboxplot(Input18q,'Subtype','ArmCNV', color = 'Subtype') + 
  scale_color_manual(values = ColImmune) +
  stat_pvalue_manual(compare_means(ArmCNV ~ Subtype, data = Input18q, method = "t.test") %>%
                       mutate(y.position = c(0.6,0.7,0.8,0.9,1,1.1)), label = "p.signif") +
  xlab('') + ylab('Chr18q CNV') +
  theme(legend.position = 'none')
# data for Table S6
write.xlsx(Input18q, file = 'Chr18q-cluster.xlsx', rowNames = T, colNames = T, overwrite = T)
#
Input18q$Del18q <- ifelse(Input18q$ArmCNV < -0.25, 'Del', 'Others') %>%
  factor(levels = c('Del','Others'), ordered = T)
table(Input18q$Del18q)

##### Figure S5E #####
# data for Table S6
ORA <- clusterProfiler::enricher(filter(GeneCoord, grepl('18q',ChrBand))$gene_name %>% unique(), TERM2GENE = TotalPathway)@result %>%
  filter(pvalue <= 0.05) %>% arrange(pvalue) %>%
  mutate(LogPval = -log10(pvalue))
write.xlsx(ORA, file = 'Chr18q gene function.xlsx', rowNames = T, colNames = T, overwrite = T)
# ORA
ORA <- clusterProfiler::enricher(filter(GeneCoord, grepl('18q',ChrBand))$gene_name %>% unique(), TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.25) %>% select(ID, pvalue) %>% arrange(pvalue) %>%
  mutate(LogPval = -log10(pvalue))
##### Vis
ggbarplot(ORA, 'ID', 'LogPval', 
          color = ColColor$`Single-Brown`[9], fill = ColColor$`Single-Brown`[9]) +
  scale_y_continuous(expand = c(0,0)) +
  xlab('') +
  coord_flip()
#
rm(ORA)

##### Figure 5F #####
Input18q$ImmuneScore = sapply(Input18q$Sample, function(x){filter(IndexGenome, WES == x)$ImmuneScore}) %>% as.numeric()
Input18q$StromaScore = sapply(Input18q$Sample, function(x){filter(IndexGenome, WES == x)$StromaScore}) %>% as.numeric()
##### Vis
ggboxplot(Input18q, 'Del18q', 'ImmuneScore', color = 'Del18q') +
  stat_compare_means(method = 't.test') +
  ylab('xCell Immune Score') + xlab('Chr 18q') +
  scale_color_manual(values = c('Del' = 'blue', 'Others' = 'grey')) +
  theme(legend.position = 'none') 

##### Figure S5F #####
InputCellFrac <- tibble(Sample = IndexGenome$Sample,
                        WES = IndexGenome$WES)
InputCellFrac$Del18q <- sapply(InputCellFrac$WES, function(x){filter(Input18q, Sample == x)$Del18q})
##### add cell fracs
InputCellFrac <- cbind(InputCellFrac, t(OutputBDB_Normed_TP[, InputCellFrac$Sample])) %>% as_tibble()
# melt
InputCellFrac <- reshape2::melt(InputCellFrac, id.var = c('Sample','WES','Del18q'),
                                value.name = 'CellFrac', variable.name = c('CellType')) %>% as_tibble()
##### Vis 1st
ggboxplot(InputCellFrac, 'CellType', 'CellFrac', color = 'Del18q') +
  stat_pvalue_manual(compare_means(
    CellFrac ~ Del18q, data = InputCellFrac, group.by = "CellType",
    method = "t.test"), x = "CellType", y.position = 0.4,
    label = "p.signif",
    position = position_dodge(0.8)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
##### filter
InputCellFrac <- filter(InputCellFrac, CellType %in% c('B.cells.naive','T.cells.CD8','T.cells.CD4.memory.activated',
                                                       'T.cells.follicular.helper','NK.cells.activated','Macrophages.M1',
                                                       'Dendritic.cells.resting'))
##### Vis final
ggboxplot(InputCellFrac, 'CellType', 'CellFrac', color = 'Del18q') +
  scale_color_manual(values = c('Del' = 'blue', 'Others' = 'grey')) +
  stat_pvalue_manual(compare_means(
    CellFrac ~ Del18q, data = InputCellFrac, group.by = "CellType",
    method = "t.test"), x = "CellType", y.position = 0.4,
    label = "p.signif",
    position = position_dodge(0.8)) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = 'none') +
  xlab('') + ylab('Cell fraction') 
# 
rm(InputCellFrac)

##### Figure 5G #####
Input18q$Pro <- sapply(Input18q$Sample, function(x){filter(IndexGenome, WES == x)$Pro}) %>% as.character()
##### DEA
CM <- data.frame(Del = ifelse(as.character(Input18q$Del18q) == 'Del', 1, 0),
                 Others = ifelse(as.character(Input18q$Del18q) != 'Del', 1, 0))
DEA_Pro <- bioinfoamateur::core_Differential_analysis_continuous(CM, GBC_Main_Pro[,Input18q$Pro],
                                                 log = T,p.adj = T)
DEA_Pro$Protein = rownames(DEA_Pro)
DEA_Phos <- bioinfoamateur::core_Differential_analysis_continuous(CM, GBC_Main_Phos_knn[,Input18q$Pro],
                                                  log = T,p.adj = T)
DEA_Phos$Protein = str_split_fixed(rownames(DEA_Phos),':',2)[,1]

##### data for Table S6
GSEAPro <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = DEA_Pro$Protein, FC = DEA_Pro$`Del-logFC`))
GSEAPhos <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = DEA_Phos$Protein, FC = DEA_Phos$`Del-logFC`))
# 
GSEAPro <- clusterProfiler::GSEA(GSEAPro, TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.05) %>% arrange(desc(NES)) %>%
  mutate(LogPval = -log10(pvalue), Omic = 'Protein')
GSEAPhos <- clusterProfiler::GSEA(GSEAPhos, TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.05) %>% arrange(desc(NES)) %>%
  mutate(LogPval = -log10(pvalue), Omic = 'Phosphorylation')
write.xlsx(rbind(GSEAPro, GSEAPhos), file = 'GSEA 18q.xlsx', rowNames = T, colNames = T, overwrite = T)

##### GSEA input
GSEAPro <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = DEA_Pro$Protein, FC = DEA_Pro$`Del-logFC`))
GSEAPhos <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = DEA_Phos$Protein, FC = DEA_Phos$`Del-logFC`))
##### GSEA
GSEAPro <- clusterProfiler::GSEA(GSEAPro, TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.05) %>% select(ID, NES, pvalue) %>% arrange(desc(NES)) %>%
  mutate(LogPval = -log10(pvalue), Omic = 'Protein')
GSEAPhos <- clusterProfiler::GSEA(GSEAPhos, TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.05) %>% select(ID, NES, pvalue) %>% arrange(desc(NES)) %>%
  mutate(LogPval = -log10(pvalue), Omic = 'Phosphorylation')


##### Vis Prep (Del vs WT)
SelectPathway <- c('HALLMARK_MYC_TARGETS_V1','HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                   'KEGG DNA replication','REACTOME Mismatch Repair',
                   'HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_TNFA_SIGNALING_VIA_NFKB',
                   'REACTOME Toll-like Receptor Cascades','HALLMARK_IL6_JAK_STAT3_SIGNALING')
# re GSEA
GSEAPro <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = DEA_Pro$Protein, FC = DEA_Pro$`Del-logFC`))
GSEAPhos <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = DEA_Phos$Protein, FC = DEA_Phos$`Del-logFC`))
GSEAPro <- clusterProfiler::GSEA(GSEAPro, TERM2GENE = TotalPathway, pvalueCutoff = 1)@result %>%
  filter(ID %in% SelectPathway) %>% select(ID, NES, pvalue, p.adjust) %>% arrange(desc(NES)) %>%
  mutate(LogPval = -log10(pvalue), Omic = 'Protein')
GSEAPhos <- clusterProfiler::GSEA(GSEAPhos, TERM2GENE = TotalPathway, pvalueCutoff = 1)@result %>%
  filter(ID %in% SelectPathway) %>% select(ID, NES, pvalue, p.adjust) %>% arrange(desc(NES)) %>%
  mutate(LogPval = -log10(pvalue), Omic = 'Phosphorylation')
# Integ
GSEA <- rbind(GSEAPro, GSEAPhos) %>% as_tibble()
GSEA$ID = factor(GSEA$ID, levels = SelectPathway, ordered = T)
GSEA = filter(GSEA, p.adjust <= 0.05)
##### Vis
ggbarplot(GSEA, 'ID', 'NES', color = 'Omic', fill = 'Omic',
          position = position_dodge()) +
  coord_flip() + xlab('') + ylab('NES') +
  scale_color_manual(values = c('Protein' = ColJournal$Science[5], 'Phosphorylation' = ColJournal$Science[4])) +
  scale_fill_manual(values = c('Protein' = ColJournal$Science[5], 'Phosphorylation' = ColJournal$Science[4])) +
  scale_y_continuous(expand = c(0,0)) 

##### Figure S5G #####
IndexGenome = IndexDecov %>% filter(Tumor == 'Yes', !is.na(WES))
# prep
InputBand <- mutate(GBC_Main_CNV[, IndexGenome$WES], Gene = rownames(GBC_Main_CNV[, IndexGenome$WES]))
InputBand$Band <- sapply(InputBand$Gene, function(x){as.character(filter(GeneCoord, gene_name == x)$ChrBand)[1]}) %>% as.character()
InputBand <- filter(InputBand, Band != 'character(0)')
InputBand <- InputBand[, colnames(InputBand) != 'Gene']
# calc avg
InputBandAvg <- aggregate(InputBand[,-ncol(InputBand)], by = list(Band = InputBand$Band), mean) %>% 
  as.data.frame() %>% column_to_rownames(var = 'Band')
##### DEA of cytobands
CM <- data.frame(Hot = ifelse(IndexGenome$ImmCluster_Name == 'B', 1, 0),
                 Others = ifelse(IndexGenome$ImmCluster_Name != 'B', 1, 0))
DEABand <- bioinfoamateur::core_Differential_analysis_continuous(CM, InputBandAvg, log = T, p.adj = T, show.belong = T, method = 't.test')
##### 
Input18q <- as.data.frame(t(InputBandAvg[grepl('18q',rownames(InputBandAvg)), ]))
Input18q$Subtype = sapply(rownames(Input18q), function(x){filter(IndexGenome, WES == x)$ImmCluster_Name})
##### Vis
Input18q = reshape2::melt(Input18q, id.var = 'Subtype', variable.name = 'Band', value.name = 'BandCNV') %>% as_tibble()
ggboxplot(Input18q, 'Band', 'BandCNV', color = 'Subtype') +
  scale_color_manual(values = ColImmune) +
  xlab('') + ylab('Band-level CNV') + 
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust = 1)) 
#
rm(Input18q)





