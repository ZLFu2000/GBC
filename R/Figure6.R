##### Figure 6 and Figure S6
rm(list = ls()); try(dev.off(), silent = T)
# packages loading
library(tidyverse); library(data.table)
library(ggpubr); library(cowplot); library(ggthemes); library(ggplot2); library(ggrepel)
library(maftools)
library(ComplexHeatmap); library(circlize)
library(openxlsx)

##### prep1. data loading #####
##### 1 omics
load('../DataMain.Rda')
#
load('../Figure 5/KeyVar_Imm.RData')
IndexDecov$ImmCluster_Name = case_when(IndexDecov$ImmCluster == 1 ~ 'A',      # immune - stroma - 
                                       IndexDecov$ImmCluster == 2 ~ 'B',      # immune + stroma - 
                                       IndexDecov$ImmCluster == 3 ~ 'C',      # immune + stroma + 
                                       IndexDecov$ImmCluster == 4 ~ 'D') %>%  # immune - stroma + 
  factor(levels = c('A','B','C','D'), ordered = T)
table(IndexDecov$ImmCluster_Name)
# as new
index.raw = IndexDecov; rm(IndexDecov)

#### 2 color
load('../Colors (ggsci).RData')
library(wesanderson)
scales::show_col(ColJournal$Science)
# Col Immune
col.imm = c(ColJournal$Nature[5], ColJournal$Nature[1],
            ColJournal$Nature[4], ColJournal$Nature[9])
names(col.imm) <- c('A','B','C','D')
# Col multiomics cluster
col.mo = c(ColJournal$COSMICsignature[1], ColJournal$COSMICsignature[3],
           ColJournal$COSMICsignature[5], ColJournal$COSMICsignature[6])
names(col.mo) <- c('MO1','MO2','MO3','MO4')
# col pathology
library(wesanderson)
col.patho = c(wes_palette('Darjeeling1', 3, type = c("discrete")), 'grey70')
names(col.patho) = c('AC','AS','NE','Others')

### 3 Genesets
load('../Hallmark_KEGG_Reactome.RData')

# 4 COSMIC Cancer Gene Census
ref.cosmic <- fread('../COSMIC genes.csv')

# 5. others
# drugs
ref.drugs <- fread('./DGIdb.txt')

##### prep2. index #####
# Gender 
index.raw$Male = sapply(index.raw$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Sex}) %>% as.character()
index.raw$Male = ifelse(index.raw$Tumor == 'Yes', index.raw$Male, NA)
index.raw$Male = ifelse(index.raw$Male == 'Male', 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)
# Age
index.raw$Age = sapply(index.raw$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Age}) %>% as.numeric()
## pathologic type
index.raw$pathology = sapply(index.raw$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Pathological_type}) %>% as.character()
index.raw$pathology = case_when(index.raw$pathology == 'adenocarcinoma' ~ 'AC',
                                index.raw$pathology == 'adenosquanmous carcinoma' ~ 'AS',
                                index.raw$pathology == 'neuroendocrine carcinoma' ~ 'NE',
                                !is.na(index.raw$pathology) ~ 'Others') %>%
  factor(levels = c('AC','AS','NE','Others'), ordered = T)
#
library(wesanderson)
ColPatho2 = c(wes_palette('Darjeeling1', 3, type = c("discrete")), 'grey70')
names(ColPatho2) = c('AC','AS','NE','Others')
# cholelithiasis
index.raw$Cholelithiasis = sapply(index.raw$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$cholelithiasis}) %>% as.numeric()
index.raw$Cholelithiasis = ifelse(index.raw$Cholelithiasis == 1, 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)
# LI
index.raw$LiverInvasion = sapply(index.raw$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Regional_invasion}) %>% as.numeric()
index.raw$LiverInvasion = ifelse(index.raw$LiverInvasion == 1, 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)
# Distal metastasis
index.raw$DistalMetastasis = sapply(index.raw$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Distal_metastasis}) %>% as.numeric()
index.raw$DistalMetastasis = ifelse(index.raw$DistalMetastasis == 1, 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)
# TNM
index.raw$TNM = sapply(index.raw$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Tumor_stage2}) %>% as.numeric()
index.raw$TNM = case_when(index.raw$TNM == 1 ~ 'I', index.raw$TNM == 2 ~ 'II',
                          index.raw$TNM == 3 ~ 'III', index.raw$TNM == 4 ~ 'IV') %>%
  factor(levels = c('I','II','III','IV'), ordered = T)
# TMB
GBC_Main_Maf <- maftools::read.maf(maf = GBC_Main_Mutation)
TempTMB <- tmb(GBC_Main_Maf)
index.raw$TMB <- sapply(index.raw$Sample, function(x){filter(TempTMB, Tumor_Sample_Barcode == x)$total_perMB %>% as.numeric()}) %>% as.numeric()
rm(TempTMB)
# TP53 mut
index.raw$TP53Mut <- case_when(index.raw$Sample %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'TP53')$Tumor_Sample_Barcode ~ 'Yes',
                               index.raw$Sample %in% GBC_Main_Mutation$Tumor_Sample_Barcode ~ 'No',
                               !is.na(index.raw$Sample) ~ NA) %>%
  factor(levels = c('Yes','No'), ordered = T)
# ERBB3 mut
index.raw$ERBB3Mut <- case_when(index.raw$Sample %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB3')$Tumor_Sample_Barcode ~ 'Yes',
                                index.raw$Sample %in% GBC_Main_Mutation$Tumor_Sample_Barcode ~ 'No',
                                !is.na(index.raw$Sample) ~ NA) %>%
  factor(levels = c('Yes','No'), ordered = T)
# ERBB2 Amp
index.raw$ERBB2Amp <- sapply(index.raw$Sample, function(x){filter(Input, Patient == x)$CNV_ERBB2 %>% as.character()}) %>% as.character()
index.raw$ERBB2Amp = ifelse(index.raw$ERBB2Amp == 'Amplification', 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)
rm(Input)
# MYC Amp
index.raw$MYCAmp <- sapply(index.raw$Sample, function(x){filter(Input, Patient == x)$CNV_MYC %>% as.character()}) %>% as.character()
index.raw$MYCAmp = ifelse(index.raw$MYCAmp == 'Amplification', 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)
rm(Input)

##### retain tumor only
index.raw.T = filter(index.raw, Tumor == 'Yes')

########## Mission 1. add MO cluster (N = 129) ##########
##### 1.1 data call back [wasted] #####

##### 1.2 redo MO cluster: prep #####
index.mutual <- filter(index.raw.T, !is.na(WES), !is.na(RNA), !is.na(Pro))
# order are the same
##### Mut
mo.mut <- pmin(table(GBC_Main_Mutation$Hugo_Symbol, GBC_Main_Mutation$Tumor_Sample_Barcode), 1) %>%
  as.data.frame()
mo.mut <- dcast(mo.mut, formula = Var1 ~ Var2)
rownames(mo.mut) <- mo.mut$Var1; mo.mut <- dplyr::select(mo.mut, -Var1)
mo.mut = mo.mut[, index.mutual$Sample]
mo.mut[1:5,1:5]
##### CNV
mo.cnv <- GBC_Main_CNV[, index.mutual$WES]
colnames(mo.cnv) = index.mutual$Sample
##### RNA
mo.rna <- GBC_Main_RNA_logTPM[, index.mutual$RNA]
colnames(mo.rna) = index.mutual$Sample
##### Protein
mo.pro <- GBC_Main_Pro[, index.mutual$Pro]
colnames(mo.pro) = index.mutual$Sample
##### Phos
mo.phos <- GBC_Main_Phos_knn[, index.mutual$Pro]
colnames(mo.phos) = index.mutual$Sample

##### 1.3 feature selection and integ as list containing all 5 omic data #####
library(MOVICS)
##### mut
temp.mut <- getElites(dat = mo.mut,
                      method = 'freq', elite.pct = 0.1)
temp.mut = temp.mut$elite.dat
##### CNV
temp.cnv <- getElites(dat = mo.cnv,
                      method = 'mad', elite.pct = 0.05)
temp.cnv = temp.cnv$elite.dat
##### RNA
temp.rna <- getElites(dat = mo.rna,
                      method = 'mad', elite.pct = 0.05)
temp.rna = temp.rna$elite.dat
##### Pro
temp.pro <- getElites(dat = mo.pro,
                      method = 'mad', elite.pct = 0.05)
temp.pro = temp.pro$elite.dat
##### Phos
temp.phos <- getElites(dat = mo.phos,
                       method = 'mad', elite.pct = 0.05)
temp.phos = temp.phos$elite.dat

##### integ
mo.data <- list(Mutation = temp.mut, 
                CNV = temp.cnv, 
                RNA = temp.rna, 
                Protein = temp.pro, 
                Phosphorylation = temp.phos)
rm(temp.cnv, temp.mut, temp.phos, temp.pro, temp.rna)

##### 1.4 select best cluster N #####
movics.optmN <- getClustNum(data        = mo.data,
                            is.binary   = c(T,F,F,F,F), # necessary
                            try.N.clust = 2:6,        # try cluster number from 2 to 6
                            fig.name    = "1_Optm_N")

##### 1.5 Clustering and saving #####
# using all methods
set.seed(1120)
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", 'iClusterBayes',
                                            "LRAcluster", "ConsensusClustering", 
                                            "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 4,
                         type        = c("binomial", "gaussian", "gaussian", "gaussian", "gaussian"))

##### 1.6 check MO cluster results #####
##### consistency?
moic.consist <- getConsensusMOIC(moic.res.list = moic.res.list,
                                 fig.name      = "2_Consensus_heatmap",
                                 distance      = "euclidean",
                                 linkage       = "average")
##### sil
getSilhouette(sil      = moic.consist$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "2_SILHOUETTE",
              height   = 5.5,
              width    = 5)

##### 1.7 extract cluster #####
movics.cluster = moic.consist$clust.res
table(movics.cluster$clust)
##### add to index
index.mutual$MOcluster = case_when(movics.cluster$clust == 1 ~ 'MO1',
                                   movics.cluster$clust == 2 ~ 'MO2',
                                   movics.cluster$clust == 3 ~ 'MO3',
                                   movics.cluster$clust == 4 ~ 'MO4') %>%
  factor(levels = c('MO1','MO2','MO3','MO4'), ordered = T)

##### 1.8 survival obv in primary clustering #####
#
temp.surv.info = index.mutual[,c(13,14)] %>% as.data.frame()
rownames(temp.surv.info) = index.mutual$Sample
colnames(temp.surv.info) = c('futime','fustat')
temp.surv.info$futime = temp.surv.info$futime * 30
#
surv.brca <- compSurv(moic.res         = moic.consist,
                      surv.info        = temp.surv.info,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      fig.name        = "2_OS_movics")
#
rm(temp.surv.info, surv.brca)

##### 1.9 exam clinipatho paras #####
##### make metadata input
temp.clini.info = index.mutual %>% as.data.frame() %>% 
  column_to_rownames(var = 'Sample')
temp.clini.info <- select(temp.clini.info, futime = OStime, fustat = OS,
                          ImmCluster_Name, pathology, Cholelithiasis, 
                          LiverInvasion, DistalMetastasis, TNM, 
                          TP53Mut, ERBB3Mut, ERBB2Amp, MYCAmp,
                          # numeric
                          Age, TumorPurityTSNet, ImmuneScore, StromaScore, TMB)
temp.clini.info$ImmCluster_Name = as.character(temp.clini.info$ImmCluster_Name)
temp.clini.info$pathology = as.character(temp.clini.info$pathology)
temp.clini.info$Cholelithiasis = as.character(temp.clini.info$Cholelithiasis)
temp.clini.info$LiverInvasion = as.character(temp.clini.info$LiverInvasion)
temp.clini.info$DistalMetastasis = as.character(temp.clini.info$DistalMetastasis)
temp.clini.info$TNM = as.character(temp.clini.info$TNM)
temp.clini.info$TP53Mut = as.character(temp.clini.info$TP53Mut)
temp.clini.info$ERBB3Mut = as.character(temp.clini.info$ERBB3Mut)
temp.clini.info$ERBB2Amp = as.character(temp.clini.info$ERBB2Amp)
temp.clini.info$MYCAmp = as.character(temp.clini.info$MYCAmp)

##### calc signif
moic.clinipatho <- compClinvar(moic.res      = moic.consist,
                               var2comp      = temp.clini.info, # data.frame needs to summarize (must has row names of samples)
                               # exactVars = c('ImmuneScore'),
                               strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                               doWord        = TRUE, # generate .docx file in local path
                               tab.name      = "2_clinipatho_features")
moic.clinipatho = as.data.frame(moic.clinipatho)
##### check word
View(fread('2_clinipatho_features.txt'))
#
rm(temp.clini.info)
# Record
# Sig: immCluster(Tumor purity, immune score, stroma score), TNM stage, ERBB2 Amp

##### 1.10 diff mutations #####
moic.mutdif <- compMut(moic.res     = moic.consist,
                       mut.matrix   = mo.data$Mutation, # binary somatic mutation matrix
                       doWord       = F,    # generate table in .docx format
                       doPlot       = TRUE, # draw OncoPrint
                       freq.cutoff  = 0.05,    # keep those genes that mutated in at least 5% of samples
                       p.adj.cutoff = 0.05,  # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                       innerclust   = TRUE, # perform clustering within each subtype
                       width        = 6, 
                       height       = 6,
                       fig.name     = "2_mutation_diff")

##### 1.11 TMB between clusters #####
moic.tmbdif <- compTMB(moic.res     = moic.consist,
                       maf          = mutate(GBC_Main_Mutation, Tumor_Seq_Allele1 = 'C'),
                       rmDup        = TRUE, # remove duplicated variants per sample
                       rmFLAGS      = FALSE, # keep FLAGS mutations
                       exome.size   = 38, # estimated exome size
                       test.method  = "parametric", # statistical testing method
                       fig.name     = "2_TMB")

##### rm of Mission 1 #####
rm(moic.tmbdif, moic.mutdif)
detach("package:MOVICS", unload = TRUE)
save(index.mutual, file = 'MO129_1123.rda')

########## Mission 2. training (N = 64) by ML models ##########
library(xgboost)
library(pROC)
##### 2.1 test the performance of xgBoost on our data #####
##### input mtx prep
ml.cv.data <- cbind(t(mo.mut), t(mo.cnv), t(mo.rna), t(mo.pro), t(mo.phos))
### mtx to xgb.DMatrix
ml.cv.data <- xgb.DMatrix(data = as.matrix(ml.cv.data),
                          label = as.numeric(index.mutual$MOcluster) - 1)
### 10X CV
set.seed(1123)
cv_results <- xgb.cv(
  params = list(
    objective = "multi:softprob", num_class = 4, 
    eval_metric = "mlogloss",     
    eta = 0.1,                    
    max_depth = 6,                 
    subsample = 0.8,              
    colsample_bytree = 0.8        
  ),
  data = ml.cv.data,
  nrounds = 100,                
  nfold = 10,                   
  verbose = TRUE, nthread = 8,
  prediction = TRUE             
)

### draw ROC
# temp.pred.labels <- cv_results$pred; colnames(temp.pred.labels) = c(1,2,3,4)
# #  ROC 
# temp.roc.data <- multiclass.roc(factor(as.numeric(index.mutual$MOcluster), 
#                                        levels = c(1,2,3,4), ordered = T),
#                                 temp.pred.labels)
# temp.auc <- temp.roc.data$auc # Macro-AUC = 0.95
# 
# ### draw roc
# temp.roc.data <- temp.roc.data$rocs
# temp.vis.roc <- data.frame()
# # 
# for (combo in names(temp.roc.data)) {
#   # 
#   current_roc <- temp.roc.data[[combo]]
#   
#   # roc1
#   roc1 <- data.frame(
#     Sensitivity = current_roc[[1]]$sensitivities,
#     Specificity = 1 - current_roc[[1]]$specificities,
#     Combination = combo,
#     Curve = "1st"
#   )
#   
#   # roc2
#   roc2 <- data.frame(
#     Sensitivity = current_roc[[2]]$sensitivities,
#     Specificity = 1 - current_roc[[2]]$specificities,
#     Combination = combo,
#     Curve = "2nd"
#   )
#   
#   # combine
#   temp.vis.roc <- bind_rows(temp.vis.roc, roc1, roc2)
#   #
#   rm(current_roc, roc1, roc2)
# }
# # vis
# ggplot(temp.vis.roc, aes(x = Specificity, y = Sensitivity, color = Combination, linetype = Curve)) +
#   geom_line(size = 1) +
#   theme_minimal() +
#   labs(
#     title = paste("Multiclass ROC Curve (XGBoost)\nMacro-AUC:",
#                   round(as.numeric(temp.auc), 3)),     
#     x = "1 - Specificity",
#     y = "Sensitivity",
#     color = "Class Combination",
#     linetype = "Curve"
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     legend.position = "right",
#     text = element_text(family = "Arial", color = "black"), 
#     panel.grid = element_blank(), 
#     panel.background = element_rect(fill = "white", color = "white"), 
#     axis.line = element_line(color = "black"), 
#     axis.ticks = element_line(color = "black"), 
#     axis.text = element_text(color = "black")  
#   )

##### 2.2 Feature selection #####
library(randomForest)
# 
perform_feature_selection <- function(data, labels, method = "threshold", 
                                      threshold = NULL, cumulative_rate = 0.9, 
                                      is_binary = FALSE) {
  # 1
  if (is_binary) {
    # binary
    p_values <- apply(data, 1, function(feature) {
      chisq.test(table(feature, labels))$p.value
    })
  } else {
    # continuous
    p_values <- apply(data, 1, function(feature) {
      summary(aov(feature ~ labels))[[1]][["Pr(>F)"]][1]
    })
  }
  
  # sig?
  significant_features <- names(p_values)[p_values < 0.05]
  significant_data <- data[significant_features, ]
  
  # 2 feature selection
  rf_model <- randomForest(t(significant_data), as.factor(labels), importance = TRUE)
  importance_scores <- importance(rf_model)  
  # select
  selected_features <- rownames(importance_scores)[
    importance_scores[, "MeanDecreaseAccuracy"] > 0 & 
      importance_scores[, "MeanDecreaseGini"] > 0
  ]
  
  
  # 
  selected_data <- significant_data[selected_features, ]
  # return(selected_data)
  # return(significant_data)
}

# 
ml.fs.mut <- perform_feature_selection(mo.mut, index.mutual$MOcluster, method = "threshold", is_binary = TRUE)
ml.fs.cnv <- perform_feature_selection(mo.cnv, index.mutual$MOcluster, method = "threshold", is_binary = F)
ml.fs.rna <- perform_feature_selection(mo.rna, index.mutual$MOcluster, method = "threshold", is_binary = F)
ml.fs.pro <- perform_feature_selection(mo.pro, index.mutual$MOcluster, method = "threshold", is_binary = F)
ml.fs.phos <- perform_feature_selection(mo.phos, index.mutual$MOcluster, method = "threshold", is_binary = F)
# change rownames
rownames(ml.fs.mut) <- paste("MUT", rownames(ml.fs.mut), sep = ":")
rownames(ml.fs.cnv) <- paste("CNV", rownames(ml.fs.cnv), sep = ":")
rownames(ml.fs.rna) <- paste("RNA", rownames(ml.fs.rna), sep = ":")
rownames(ml.fs.pro) <- paste("PRO", rownames(ml.fs.pro), sep = ":")
rownames(ml.fs.phos) <- paste("PHO", rownames(ml.fs.phos), sep = ":")

##### 2.3 normalize #####
# check
if (!all(colnames(ml.fs.mut) == colnames(ml.fs.cnv) &&
         colnames(ml.fs.cnv) == colnames(ml.fs.rna) &&
         colnames(ml.fs.rna) == colnames(ml.fs.pro) &&
         colnames(ml.fs.pro) == colnames(ml.fs.phos))) {
  stop("Column names (sample names) are not consistent across matrices!")
}
# scale
zscore_standardize <- function(matrix) {
  t(scale(t(matrix), center = TRUE, scale = TRUE))  
}
ml.fs.nor.cnv <- zscore_standardize(ml.fs.cnv)
ml.fs.nor.rna <- zscore_standardize(ml.fs.rna)
ml.fs.nor.pro <- zscore_standardize(ml.fs.pro)
ml.fs.nor.phos <- zscore_standardize(ml.fs.phos)
# rbind
ml.fs.nor <- rbind(ml.fs.mut, ml.fs.nor.cnv, ml.fs.nor.rna, ml.fs.nor.pro, ml.fs.nor.phos)
#
cat("Integrated data dimensions: ", dim(ml.fs.nor), "\n")
cat("Number of missing values in integrated data: ", sum(is.na(ml.fs.nor)), "\n")

##### 2.4 select best params using 10xCV #####
# 1
ml.xgb.dmtx <- xgb.DMatrix(data = t(as.matrix(ml.fs.nor)), 
                           label = as.numeric(index.mutual$MOcluster) - 1)

# net para
param_grid <- expand.grid(
  eta = c(0.01, 0.05, 0.1),           
  max_depth = c(4, 6),          
  subsample = c(0.8),          
  colsample_bytree = c(0.8),    
  nrounds = c(75, 125)           
)

best_params <- NULL
best_logloss <- Inf

# 
for (i in 1:nrow(param_grid)) {
  # 
  params <- list(
    objective = "multi:softprob",
    num_class = 4,
    eval_metric = "mlogloss",
    eta = param_grid$eta[i],
    max_depth = param_grid$max_depth[i],
    subsample = param_grid$subsample[i],
    colsample_bytree = param_grid$colsample_bytree[i]
  )
  cat('Current params:\n')
  params
  
  # 
  cv_results <- xgb.cv(
    params = params,
    data = ml.xgb.dmtx,               
    nrounds = param_grid$nrounds[i],
    nfold = 10,                      
    verbose = T,
    prediction = FALSE,
    nthread = 9                     
  )
  
  # logloss
  min_logloss <- min(cv_results$evaluation_log$test_mlogloss_mean)
  
  # 
  if (min_logloss < best_logloss) {
    best_logloss <- min_logloss
    best_params <- params
    best_nrounds <- param_grid$nrounds[i]
  }
}

# 
# record: eta = 0.05, depth = 6, nround = 125, loss = 0.3961896
cat("Best parameters:\n")
print(best_params)
cat("Best number of rounds: ", best_nrounds, "\n")
cat("Best logloss: ", best_logloss, "\n")
rm(best_nrounds, best_logloss)


##### 2.5 model training #####
# training
set.seed(1127)
ml.model <- xgb.train(
  params = best_params,
  data = ml.xgb.dmtx,
  nrounds = 125,
  verbose = TRUE
)
# saving
save(ml.model, file = 'Model_xgBoost_1124.Rda')

##### 2.6 apply model on 129 samples to get ROC curve #####
##### 2.6.1 predict
temp.predict <- predict(ml.model, ml.xgb.dmtx)  # xgb.DMatrix 
# to mtx
temp.predict <- matrix(temp.predict, nrow = 129, byrow = TRUE)
colnames(temp.predict) = c(1,2,3,4)
# macro AUC = 1, no meaning 

##### 2.7 applying model on the remaining 64 samples [prediction] #####
##### 2.7.1 input prep #####
##### index
index.remain <- filter(index.raw.T, !Sample %in% index.mutual$Sample)
### mut
mo.remain.mut <- pmin(table(GBC_Main_Mutation$Hugo_Symbol, GBC_Main_Mutation$Tumor_Sample_Barcode), 1) %>%
  as.data.frame() %>% filter(Var2 %in% index.remain$Sample)
mo.remain.mut <- dcast(mo.remain.mut, Var1 ~ Var2, value.var = "Freq", fill = 0) %>% as.data.frame()
mo.remain.mut <- column_to_rownames(mo.remain.mut, var = 'Var1')
head(mo.remain.mut, c(6,6))
# add missing sample
for (sample in setdiff(index.remain$Sample, colnames(mo.remain.mut))) {mo.remain.mut[[sample]] <- NA}
mo.remain.mut <- mo.remain.mut[, index.remain$Sample]
# 0-1 mtx
mo.remain.mut <- as.data.frame(lapply(mo.remain.mut, as.numeric), row.names = rownames(mo.remain.mut))
# rowname order
mo.remain.mut <- mo.remain.mut[str_split_fixed(rownames(ml.fs.mut), 'MUT:', 2)[,2], ]

### cnv
mo.remain.cnv <- GBC_Main_CNV[str_split_fixed(rownames(ml.fs.nor.cnv), 'CNV:', 2)[,2], na.omit(index.remain$WES)]
colnames(mo.remain.cnv) = sapply(colnames(mo.remain.cnv),function(x){filter(index.remain, WES == x)$Sample})
# add missing sample
for (sample in setdiff(index.remain$Sample, colnames(mo.remain.cnv))) {mo.remain.cnv[[sample]] <- NA}
mo.remain.cnv <- mo.remain.cnv[, index.remain$Sample]

### rna
mo.remain.rna <- GBC_Main_RNA_logTPM[, na.omit(index.remain$RNA)]
colnames(mo.remain.rna) = sapply(colnames(mo.remain.rna),function(x){filter(index.remain, RNA == x)$Sample})
# add missing sample
for (sample in setdiff(index.remain$Sample, colnames(mo.remain.rna))) {mo.remain.rna[[sample]] <- NA}
mo.remain.rna <- mo.remain.rna[str_split_fixed(rownames(ml.fs.nor.rna), 'RNA:', 2)[,2], index.remain$Sample]

### pro
mo.remain.pro <- GBC_Main_Pro[, index.remain$Pro]
colnames(mo.remain.pro) = sapply(colnames(mo.remain.pro),function(x){filter(index.remain, Pro == x)$Sample})
# add missing sample
for (sample in setdiff(index.remain$Sample, colnames(mo.remain.pro))) {mo.remain.pro[[sample]] <- NA}
mo.remain.pro <- mo.remain.pro[str_split_fixed(rownames(ml.fs.nor.pro), 'PRO:', 2)[,2], index.remain$Sample]

### phos
mo.remain.phos <- GBC_Main_Phos_knn[, index.remain$Pro]
colnames(mo.remain.phos) = sapply(colnames(mo.remain.phos),function(x){filter(index.remain, Pro == x)$Sample})
# add missing sample
for (sample in setdiff(index.remain$Sample, colnames(mo.remain.phos))) {mo.remain.pro[[sample]] <- NA}
mo.remain.phos <- mo.remain.phos[str_split_fixed(rownames(ml.fs.nor.phos), 'PHO:', 2)[,2], index.remain$Sample]

##### na check
cat("Number of NA values after processing: ", sum(is.na(mo.remain.cnv)), "\n")

##### 2.7.2 prediction #####
#scale
mo.remain.cnv <- zscore_standardize(mo.remain.cnv)
mo.remain.rna <- zscore_standardize(mo.remain.rna)
mo.remain.pro <- zscore_standardize(mo.remain.pro)
mo.remain.phos <- zscore_standardize(mo.remain.phos)
# change row name
rownames(mo.remain.mut) <- paste("MUT", rownames(mo.remain.mut), sep = ":")
rownames(mo.remain.cnv) <- paste("CNV", rownames(mo.remain.cnv), sep = ":")
rownames(mo.remain.rna) <- paste("RNA", rownames(mo.remain.rna), sep = ":")
rownames(mo.remain.pro) <- paste("PRO", rownames(mo.remain.pro), sep = ":")
rownames(mo.remain.phos) <- paste("PHO", rownames(mo.remain.phos), sep = ":")
# integ
mo.remain.mut[is.na(mo.remain.mut)] = 0
mo.remain.cnv[is.na(mo.remain.cnv)] = 0
mo.remain.rna[is.na(mo.remain.rna)] = 0
mo.remain.pro[is.na(mo.remain.pro)] = 0
mo.remain.phos[is.na(mo.remain.phos)] = 0
ml.xgb.pred.dmtx <- rbind(mo.remain.mut, mo.remain.cnv, mo.remain.rna, mo.remain.pro, mo.remain.phos) %>% t()
ml.xgb.pred.dmtx <- xgboost::xgb.DMatrix(data = as.matrix(ml.xgb.pred.dmtx))
ml.xgb.pred.dmtx
# save model and input
# save(mo.remain.mut, mo.remain.cnv, mo.remain.rna, mo.remain.pro, mo.remain.phos, ml.xgb.pred.dmtx, ml.model,
#      file = 'Prediction_Input_and_Model.Rda')
# predict
ml.predict <- predict(ml.model, ml.xgb.pred.dmtx)         # prediction
ml.predict <- matrix(ml.predict, nrow = 64, byrow = TRUE)
ml.predict
ml.predict.final <- max.col(ml.predict)  # 
table(ml.predict.final)
table(index.mutual$MOcluster)
# add tag
index.remain$MOcluster <- case_when(ml.predict.final == 1 ~ 'MO1',
                                    ml.predict.final == 2 ~ 'MO2',
                                    ml.predict.final == 3 ~ 'MO3', 
                                    ml.predict.final == 4 ~ 'MO4') %>%
  factor(levels = levels(index.mutual$MOcluster), ordered = T)
# integ index
index.final <- rbind(index.mutual, index.remain) %>% arrange(Sample)
index.final$MOcluster %>% table()
# rm
rm(mo.remain.mut, mo.remain.cnv, mo.remain.phos, mo.remain.pro, mo.remain.rna,
   best_params, cv_results, param_grid, params, i, min_logloss, ml.predict.final, ml.xgb.dmtx, ml.xgb.pred.dmtx, sample, mutation_rate,
   ml.fs.cnv, ml.fs.mut, ml.fs.nor, ml.fs.phos, ml.fs.rna, ml.fs.pro,
   ml.fs.nor.cnv, ml.fs.nor.phos, ml.fs.nor.pro, ml.fs.nor.rna,
   ml.predict)

##### 2.7.3 feature contribution (reload model) #####
##### load model
library(xgboost)
load('Prediction_Input_and_Model.Rda')
##### view importance first
ml.xgb.pred.dmtx <- rbind(mo.remain.mut, mo.remain.cnv, mo.remain.rna, mo.remain.pro, mo.remain.phos) %>% t()
ml.importance <- xgb.model.dt.tree(feature_names = ml.model$feature_names, model = ml.model)
# add cluster of trees
ml.importance <- ml.importance %>%
  mutate(cluster = (Tree %% 4) + 1)  
# importance
ml.importance <- ml.importance %>%
  filter(Feature != "Leaf") %>%           
  group_by(cluster, Feature) %>%         
  summarise(
    Freq = n(),                           
    AvgGain = mean(Quality, na.rm = TRUE) 
  ) %>%
  mutate(Importance = Freq * AvgGain) %>% 
  arrange(cluster, desc(Importance))      
# extract feature
ml.importance$Omic <- str_split_fixed(ml.importance$Feature, fixed(':'), 2)[,1]
ml.importance$Gene <- str_split_fixed(ml.importance$Feature, fixed(':'), 2)[,2]
# function
View(clusterProfiler::enricher(unique(filter(ml.importance, cluster == 4)$Gene),
                               TERM2GENE = TotalPathway)@result)
##### vis
ml.top.features <- ml.importance %>%
  group_by(cluster) %>%
  top_n(10, Importance) %>%
  arrange(cluster, desc(Importance))
ml.top.features$cluster = case_when(ml.top.features$cluster == 1 ~ 'MO1',
                                    ml.top.features$cluster == 2 ~ 'MO2',
                                    ml.top.features$cluster == 3 ~ 'MO3',
                                    ml.top.features$cluster == 4 ~ 'MO4')

ggplot(ml.top.features, aes(x = reorder(Feature, Importance), y = Importance, fill = as.factor(cluster))) +
  geom_bar(stat = "identity", show.legend = FALSE) +  
  scale_fill_manual(values = col.mo) +
  coord_flip() +
  facet_wrap(~cluster, scales = "free_y") +
  theme(
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA),  
    panel.grid = element_blank(),                                
    axis.line = element_line(color = "black"),                   
    axis.ticks = element_line(color = "black"),                  
    axis.text = element_text(family = "Arial", color = "black", size = 10),  
    axis.title = element_text(family = "Arial", color = "black", size = 12), 
    strip.text = element_text(family = "Arial", color = "black", size = 12), 
    plot.title = element_text(family = "Arial", color = "black", size = 14, hjust = 0.5) 
  )  +
  labs(title = "Feature importance",
       x = "Feature", y = "Importance",
       fill = "Cluster")
# data for Table S7
write.xlsx(ml.importance, file = 'ML feature importance.xlsx', rowNames = T, colNames = T, overwrite = T)

##### 2.7.4 crucial feature function of each cluster #####
ml.feature.ora <- rbind(MO1 = clusterProfiler::enricher(unique(filter(ml.importance, cluster == 1)$Gene),
                                                        TERM2GENE = TotalPathway)@result %>% filter(p.adjust <= 0.05) %>% arrange(desc(p.adjust)) %>% mutate(Cluster = 'MO1'),
                        MO2 = clusterProfiler::enricher(unique(filter(ml.importance, cluster == 2)$Gene),
                                                        TERM2GENE = TotalPathway)@result %>% filter(p.adjust <= 0.05) %>% arrange(desc(p.adjust)) %>% mutate(Cluster = 'MO2'),
                        MO3 = clusterProfiler::enricher(unique(filter(ml.importance, cluster == 3)$Gene),
                                                        TERM2GENE = TotalPathway)@result %>% filter(p.adjust <= 0.05) %>% arrange(desc(p.adjust)) %>% mutate(Cluster = 'MO3'),
                        MO4 = clusterProfiler::enricher(unique(filter(ml.importance, cluster == 4)$Gene),
                                                        TERM2GENE = TotalPathway)@result %>% filter(p.adjust <= 0.05) %>% arrange(desc(p.adjust)) %>% mutate(Cluster = 'MO4'))
# data for Table S7
write.xlsx(ml.feature.ora, file = 'ORA ML feature.xlsx', rowNames = T, colNames = T, overwrite = T)
#
ml.feature.ora <- select(ml.feature.ora, ID, `p.adjust`, Cluster)
##### Vis prep
ml.feature.ora <- ml.feature.ora %>%
  arrange(Cluster, p.adjust) %>%
  mutate(UniqueID = paste0(ID, " (", Cluster, ")"))  

ml.feature.ora <- ml.feature.ora %>% mutate(log_p = -log10(p.adjust))

ml.feature.ora <- ml.feature.ora %>%
  arrange(Cluster, desc(log_p)) %>%
  mutate(UniqueID = factor(UniqueID, levels = rev(UniqueID))) 
##### Vis
ggplot(ml.feature.ora, aes(x = UniqueID, y = log_p, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = col.mo) +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip() + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),                                
    axis.line = element_line(color = "black"),                   
    axis.ticks = element_line(color = "black"),                  
    axis.text = element_text(family = "Arial", color = "black", size = 10),
    axis.title = element_text(family = "Arial", color = "black", size = 12),
    legend.title = element_text(family = "Arial", color = "black", size = 10),
    legend.text = element_text(family = "Arial", color = "black", size = 9)
  ) +
  labs(title = NULL, x = NULL, y = "-log10(FDR)", fill = "Cluster") 

#
rm(mo.remain.mut, mo.remain.cnv, mo.remain.phos, mo.remain.pro, mo.remain.rna,
   best_params, cv_results, param_grid, params, i, min_logloss, ml.predict.final, ml.xgb.dmtx, ml.xgb.pred.dmtx, sample, mutation_rate,
   ml.fs.cnv, ml.fs.mut, ml.fs.nor, ml.fs.phos, ml.fs.rna, ml.fs.pro,
   ml.fs.nor.cnv, ml.fs.nor.phos, ml.fs.nor.pro, ml.fs.nor.rna,
   ml.predict, ml.top.features, ml.importance)


# ### Plan B
# xgboost::xgb.importance(feature_names = ml.model$feature_names, model = ml.model) %>% View()
# ### calc importance
# ml.importance <- lapply(1:4, function(class_idx) {
#   # 提取第 class_idx 类别的树集合
#   booster <- xgboost::xgb.model.dt.tree(feature_names = ml.model$feature_names, model = ml.model)
#   # 计算特征重要性
#   importance <- booster[Feature != "Leaf", .(Gain = sum(Quality), Frequency = .N), by = Feature]
#   # 3. 排序特征重要性
#   importance <- importance[order(-Gain)]
#   
#   importance[, Class := paste("Class", class_idx)]
#   return(importance)
# })

########## Mission 3. Analysis using final cluster info ##########
##### 3.1 Integ all samples and check clinipatho features #####
table(index.final$MOcluster, index.final$ImmCluster_Name)
chisq.test(table(index.final$MOcluster, index.final$ImmCluster_Name))

##### 3.2 add other features (TNB, CIN, etc) #####
### add CIN
# Chr info
load('/Users/fuzile/Desktop/たいようけい/基本文件/常用特殊数据集/ChrBandLength.RData')
# Gene Location 
GeneCoord <- rtracklayer::import('/Users/fuzile/Desktop/たいようけい/基本文件/常用特殊数据集/Homo_sapiens.GRCh38.108.gtf') %>%
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
# calc
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
##### add to Input
index.final$CIN = sapply(index.final$WES, function(x){as.numeric(filter(Output, Sample == x)$CIN)}) %>% as.numeric()
# Input$CINGroup = ifelse(Input$CIN >= median(Input$CIN, na.rm = T), 'CIN High','CIN Low') %>%
#   factor(levels = c('CIN High','CIN Low'), ordered = T)
#
rm(InputCIN, InputCalcIndex, Output)

### crucial mut and CNV manually into heatmap index
## mut
index.final$KRASMut <- case_when(index.final$Sample %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'KRAS')$Tumor_Sample_Barcode ~ 'Yes',
                                 index.final$Sample %in% GBC_Main_Mutation$Tumor_Sample_Barcode ~ 'No',
                                 !is.na(index.final$Sample) ~ NA) %>%
  factor(levels = c('Yes','No'), ordered = T)
index.final$ELF3Mut <- case_when(index.final$Sample %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ELF3')$Tumor_Sample_Barcode ~ 'Yes',
                                 index.final$Sample %in% GBC_Main_Mutation$Tumor_Sample_Barcode ~ 'No',
                                 !is.na(index.final$Sample) ~ NA) %>%
  factor(levels = c('Yes','No'), ordered = T)
## CNV
# Func
SeparatedCNV <- function(vec) {
  vec = case_when(vec == -2 ~ 'Deletion',
                  vec == -1 ~ 'Shallow Deletion',
                  vec == 0 ~ 'Diploid', 
                  vec == 1 ~ 'Gain',
                  vec == 2 ~ 'Amplification') %>%
    factor(levels = c('Amplification','Gain', 'Diploid', 'Shallow Deletion', 'Deletion'), ordered = T)
  return(vec)
}
# categroy CNV
GBC_Main_CNV_cate <- fread('/Users/fuzile/Desktop/たいようけい/我参与的项目/主要参与/2023.08 GBC MO/1 原始数据/原始数据整理/0905/整理结果 0908/Main cohort/CNV/all_thresholded.by_genes.txt') %>%
  select(-2,-3) %>% column_to_rownames(var = 'Gene Symbol')
colnames(GBC_Main_CNV_cate) = GBC_Main_Index$WES_T[match(colnames(GBC_Main_CNV_cate), GBC_Main_Index$Sample)]
GBC_Main_CNV_cate = GBC_Main_CNV_cate[, colnames(GBC_Main_CNV)]
# add
index.final$`CNV_ERBB2` <- sapply(index.final$WES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV_cate['ERBB2',x]), NA)})
index.final$`CNV_ERBB2` <- SeparatedCNV(index.final$`CNV_ERBB2`)
index.final$`CNV_MYC` <- sapply(index.final$WES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV_cate['MYC',x]), NA)})
index.final$`CNV_MYC` <- SeparatedCNV(index.final$`CNV_MYC`)
table(index.final$`CNV_MYC`);table(index.final$`CNV_ERBB2`)

### TNB
GBC.tnb <- readRDS('Res_Neos.rds')
GBC.tnb.df <- data.frame(Sample = character(), Total_Aff_nM = numeric(), stringsAsFactors = FALSE)
# 遍历 GBC.tnb 列表，计算每个样本的新抗原负荷总和
for (sample in names(GBC.tnb)) {
  # 检查样本中是否存在 `Aff(nM) BindLevel` 列
  if ("Aff(nM) BindLevel" %in% colnames(GBC.tnb[[sample]])) {
    # 计算当前样本 `Aff(nM)` 列的总和
    total_affinity <- sum(GBC.tnb[[sample]]$`Aff(nM) BindLevel`, na.rm = TRUE)
    # 将结果存储到数据框中
    GBC.tnb.df <- rbind(GBC.tnb.df, data.frame(Sample = sample, Total_Aff_nM = total_affinity))
  }
  rm(total_affinity, sample)
}
GBC.tnb.df <- as_tibble(GBC.tnb.df)
GBC.tnb.df$LogTNB <- log10(GBC.tnb.df$Total_Aff_nM + 1)
# add to index
index.final$TNB <- sapply(index.final$Sample, function(x){filter(GBC.tnb.df, Sample == x)$LogTNB %>% as.numeric()}) %>% as.numeric()
#
rm(GBC.tnb.df, GBC.tnb)

# add cohort
index.final$Cohort <- ifelse(index.final$Sample %in% index.mutual$Sample, 'MOcluster', 'xgBoost')

##### 3.2.5 save final index [containing all features, including TNB, TMB, CIN and MOcluster] #####
save(index.final, file = 'FinalFinalIndex_1201.Rda')
# data for Table S7
write.xlsx(index.final[,c(1,29,36)], file = 'Tag.xlsx', rowNames = T, colNames = T, overwrite = T)

##### 3.3 KM plot #####
library(ggsurvfit)
library(survival)
index.final[,c('OStime','OS','MOcluster')]
### vis prep
temp.fit <- survfit(Surv(time = index.final$OStime, event = index.final$OS) ~ MOcluster, data = index.final)
# calc p
temp.lrp <- survdiff(Surv(time = index.final$OStime, event = index.final$OS) ~ MOcluster, data = index.final)
temp.lrp <- 1 - pchisq(temp.lrp$chisq, length(temp.lrp$n) - 1)
# vis
library(survminer)
ggsurvplot(
  temp.fit,
  data = index.final,
  palette = col.mo,              
  conf.int = F,             
  risk.table = T,            
  pval = TRUE,                   
  pval.coord = c(40, 0.2),       
  xlab = "Time (months)", 
  ylab = "Survival Probability",
  legend.title = "Cluster",
  legend.labs = c("MO1", "MO2", "MO3", "MO4"),
  ggtheme = theme_classic(base_family = "Arial") + 
    theme(
      panel.grid = element_blank(),             
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(family = "Arial", color = "black", size = 10),
      axis.title = element_text(family = "Arial", color = "black", size = 12),
      legend.text = element_text(family = "Arial", color = "black", size = 9),
      legend.title = element_text(family = "Arial", color = "black", size = 10),
      plot.title = element_text(family = "Arial", color = "black", size = 14, hjust = 0.5)
    )
) 
#
rm(temp.fit, temp.lrp)

##### 3.4 Heatmap input prep and DEA #####
### mut
load('Prediction_Input_and_Model.Rda')
hm.mtx.mut <- pmin(table(GBC_Main_Mutation$Hugo_Symbol, GBC_Main_Mutation$Tumor_Sample_Barcode), 1) %>%
  as.data.frame()
hm.mtx.mut <- dcast(hm.mtx.mut, formula = Var1 ~ Var2)
rownames(hm.mtx.mut) <- hm.mtx.mut$Var1; hm.mtx.mut <- dplyr::select(hm.mtx.mut, -Var1)
hm.mtx.mut <- hm.mtx.mut[str_split_fixed(rownames(mo.remain.mut), 'MUT:',2)[,2], ]
hm.mtx.mut <- hm.mtx.mut[, filter(index.final, !is.na(WES))$Sample]

# mutation selection by freq
mutation_rate <- rowSums(hm.mtx.mut) / nrow(hm.mtx.mut)
temp.vis.mut <- data.frame(Gene = rownames(hm.mtx.mut), p.value = NA, 
                           MO1 = NA, MO2 = NA, MO3 = NA, MO4 = NA,
                           Highest_Frequency_Cluster = NA,
                           Highest_Frequency = NA) %>% as.tibble()
temp.tag <- sapply(colnames(hm.mtx.mut), function(x){filter(index.final, Sample == x)$MOcluster})
# 
for (i in 1:nrow(hm.mtx.mut)) {
  gene = rownames(hm.mtx.mut)[i]
  gene_mut <- hm.mtx.mut[i, ] %>% as.numeric()
  # 
  freq_by_cluster <- tapply(gene_mut, temp.tag, function(x) sum(x) / length(x))
  # add each freq
  temp.vis.mut$MO1[i] = as.numeric(freq_by_cluster)[1]
  temp.vis.mut$MO2[i] = as.numeric(freq_by_cluster)[2]
  temp.vis.mut$MO3[i] = as.numeric(freq_by_cluster)[3]
  temp.vis.mut$MO4[i] = as.numeric(freq_by_cluster)[4]
  
  # 
  highest_cluster <- names(which.max(freq_by_cluster))
  highest_frequency <- max(freq_by_cluster, na.rm = TRUE)
  # for chisq.test
  group_table <- table(temp.tag, gene_mut)
  if (all(dim(group_table) > 1)) {
    temp.vis.mut$p.value[i] <- fisher.test(group_table)$p.value
  } else {
    temp.vis.mut$p.value[i] <- NA
  }
  # add freq
  temp.vis.mut$Highest_Frequency_Cluster[temp.vis.mut$Gene == gene] <- highest_cluster
  temp.vis.mut$Highest_Frequency[temp.vis.mut$Gene == gene] <- highest_frequency
  rm(i, gene, gene_mut, group_table, highest_cluster, highest_frequency)
}
# adjust pval
temp.vis.mut$padj <- p.adjust(temp.vis.mut$p.value, method = 'fdr')
# add total freq
temp.vis.mut$Freq <- rowSums(hm.mtx.mut) / nrow(hm.mtx.mut)
# select mutation
hm.mtx.mut.select <- hm.mtx.mut[arrange(filter(temp.vis.mut, Freq >= 0.05),Highest_Frequency_Cluster)$Gene,]


### cnv
hm.mtx.cnv <- GBC_Main_CNV[, na.omit(index.final$WES)]
#
rm(mo.remain.mut, mo.remain.cnv, mo.remain.phos, mo.remain.pro, mo.remain.rna,
   ml.model)

### others
# RNA
hm.mtx.rna <- GBC_Main_RNA_logTPM[, na.omit(index.final$RNA)]
# Protein
hm.mtx.pro <- GBC_Main_Pro[, na.omit(index.final$Pro)]
# Phos
hm.mtx.phos <- GBC_Main_Phos_knn[, na.omit(index.final$Pro)]

##### change to uniform colnames
colnames(hm.mtx.cnv) = sapply(colnames(hm.mtx.cnv), function(x){filter(index.final, WES == x)$Sample})
colnames(hm.mtx.rna) = sapply(colnames(hm.mtx.rna), function(x){filter(index.final, RNA == x)$Sample})
colnames(hm.mtx.pro) = sapply(colnames(hm.mtx.pro), function(x){filter(index.final, Pro == x)$Sample})
colnames(hm.mtx.phos) = sapply(colnames(hm.mtx.phos), function(x){filter(index.final, Pro == x)$Sample})

##### DEA
# CM
CM.hm <- data.frame(MO1 = ifelse(index.final$MOcluster == 'MO1', 1, 0),
                    MO2 = ifelse(index.final$MOcluster == 'MO2', 1, 0),
                    MO3 = ifelse(index.final$MOcluster == 'MO3', 1, 0),
                    MO4 = ifelse(index.final$MOcluster == 'MO4', 1, 0))
rownames(CM.hm) = index.final$Sample
# DEA
hm.DEA.cnv <- bioinfoamateur::core_Differential_analysis_continuous(CM.hm[colnames(hm.mtx.cnv),], hm.mtx.cnv, log = T, p.adj = F, show.belong = T)
hm.DEA.rna <- bioinfoamateur::core_Differential_analysis_continuous(CM.hm[colnames(hm.mtx.rna),], hm.mtx.rna, log = T, p.adj = T, show.belong = T)
hm.DEA.pro <- bioinfoamateur::core_Differential_analysis_continuous(CM.hm[colnames(hm.mtx.pro),], hm.mtx.pro, log = T, p.adj = T, show.belong = T)
hm.DEA.phos <- bioinfoamateur::core_Differential_analysis_continuous(CM.hm[colnames(hm.mtx.phos),], hm.mtx.phos, log = T, p.adj = T, show.belong = T)
# add feature name
hm.DEA.cnv$Feature <- rownames(hm.DEA.cnv)
hm.DEA.rna$Feature <- rownames(hm.DEA.rna)
hm.DEA.pro$Feature <- rownames(hm.DEA.pro)
hm.DEA.phos$Feature <- rownames(hm.DEA.phos)
# check features
list(CNV = table(hm.DEA.cnv$belong),
     RNA = table(hm.DEA.rna$belong),
     Pro = table(hm.DEA.pro$belong),
     Phos = table(hm.DEA.phos$belong))
# additional check for cnv
filter(hm.DEA.cnv, Feature %in% ref.COSMIC$`Gene Symbol`) %>% View()

##### get final input matrix
# func prep
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
} # extend samples to the full cohort sample
# input prep
hm.input <- list(Mut = ExtendMatrix(hm.mtx.mut.select, hm.mtx.pro),
                 CNV = ExtendMatrix(hm.mtx.cnv[filter(hm.DEA.cnv, Feature %in% ref.cosmic$`Gene Symbol`)$Feature, ], hm.mtx.pro),
                 RNA = ExtendMatrix(hm.mtx.rna, hm.mtx.pro)[rownames(hm.DEA.rna), ],
                 Pro = ExtendMatrix(hm.mtx.pro, hm.mtx.pro)[rownames(hm.DEA.pro), ],
                 Phos = ExtendMatrix(hm.mtx.phos, hm.mtx.pro)[rownames(hm.DEA.phos), ])

##### 3.5 Heatmap Vis #####
library(ComplexHeatmap)
library(circlize)
### Order
hm.vis = arrange(index.final, MOcluster)

##### Annotation: Up to Down
# MOcluster-ImmuneCluster-Age-Gender-Pathology-TNM-LI-DisMet-TMB
hm.anno.sample <- columnAnnotation(`MO clusters` = hm.vis$MOcluster,
                                   `Immune clusters` = hm.vis$ImmCluster_Name,
                                   Method = hm.vis$Cohort,
                                   Age = hm.vis$Age,
                                   Gender = ifelse(hm.vis$Male == 'Yes', 'Male', 'Female'),
                                   Pathology = hm.vis$pathology,
                                   `TNM stage` = hm.vis$TNM,
                                   `Liver invasion` = hm.vis$LiverInvasion,
                                   `Distal metastasis` = hm.vis$DistalMetastasis,
                                   TMB = hm.vis$TMB,
                                   TNB = hm.vis$TNB,
                                   `TP53 Mutation` = hm.vis$TP53Mut,
                                   `ELF3 Mutation` = hm.vis$ELF3Mut,
                                   CIN = hm.vis$CIN,
                                   `ERBB2 CNV` = hm.vis$CNV_ERBB2,
                                   `MYC CNV` = hm.vis$CNV_MYC,
                                   col = list(
                                     `MO clusters` = col.mo,
                                     `Immune clusters` = col.imm,
                                     `Method` = c('MOcluster' = ColJournal$Nature[1],
                                                  'xgBoost' = ColJournal$Nature[2]),
                                     Age = colorRamp2(c(25,62,86), c(ColColor$`High-Orange`[1], ColColor$`High-Orange`[4],ColColor$`High-Orange`[8])),
                                     Gender = c('Male' = 'black', 'Female' = 'white'),
                                     Pathology = col.patho,
                                     `TNM stage` = c('I' = 'grey90','II' = 'grey70','III' = 'grey30','IV' = 'black'),
                                     `Liver invasion` = c('Yes' = 'Black', 'No' = 'white'),
                                     `Distal metastasis` = c('Yes' = 'Black', 'No' = 'white'),
                                     TMB = colorRamp2(c(0, 5, 11.54), c("white", ColColor$`Low-Indigo`[4],  ColColor$`Low-Indigo`[9])),
                                     TNB = colorRamp2(c(0, 3.5, 5.528), c("white", ColColor$`Low-Indigo`[2],  ColColor$`Low-Indigo`[9])),
                                     CIN = colorRamp2(c(0.41, 3.5, 6), c("white", ColColor$`Low-Indigo`[2],  ColColor$`Low-Indigo`[9])),
                                     `TP53 Mutation` = c('Yes' = 'black', 'No' = 'white'),
                                     `ELF3 Mutation` = c('Yes' = 'black', 'No' = 'white'),
                                     `ERBB2 CNV` = c('Amplification' = '#a9006a', 'Gain' = 'pink','Diploid' = 'white', 'Shallow Deletion' = '#8a90a4', 'Deletion' = '#162049'), # 淡紫色六位代码：#d480b4
                                     `MYC CNV` = c('Amplification' = '#a9006a', 'Gain' = 'pink','Diploid' = 'white', 'Shallow Deletion' = '#8a90a4', 'Deletion' = '#162049')
                                   ),
                                   simple_anno_size = unit(0.3333, "cm"),
                                   gp = gpar(col = "white", lwd = 0.2),
                                   na_col = 'grey90')
##### Draw heatmap
# mutation vis
hm.plot.mut <- Heatmap(hm.input$Mut[, hm.vis$Sample], cluster_columns =  F, cluster_rows = F, name = 'Mutation',
                       top_annotation = hm.anno.sample, row_title = 'Mutation', na_col = 'grey90',
                       show_column_names = F, show_row_names = F, column_split = hm.vis$MOcluster, column_gap = unit(2,'mm'),
                       col = colorRamp2(c(0,1), c('white','black')),
                       rect_gp = gpar(col = "white", lwd = 0.2),
                       height = unit(0.33*3, 'cm'), width = unit(193/15, 'cm'))
# other vis
hm.plot.cnv <- Heatmap(t(scale(t(hm.input$CNV[, hm.vis$Sample]))), cluster_columns = F, cluster_rows = F, name = 'Feature Z-score',
                       show_column_names = F, show_row_names = F, column_split = hm.vis$MOcluster, column_gap = unit(2,'mm'),
                       row_title = 'CNV',na_col = 'grey90',
                       # rect_gp = gpar(col = "white", lwd = 0.2),
                       col = colorRamp2(c(-1.5,0,1.5), c(ColColor$`Low-Blue`[9],'white',ColColor$`High-Red`[9])),
                       height = unit(0.33*6, 'cm'), width = unit(193/15, 'cm'))
hm.plot.rna <- Heatmap(t(scale(t(hm.input$RNA[, hm.vis$Sample]))), cluster_columns = F, cluster_rows = F, name = 'Feature Z-score',
                       show_column_names = F, show_row_names = F,  row_title = 'mRNA',na_col = 'grey90',
                       column_split = hm.vis$MOcluster, column_gap = unit(2,'mm'),
                       rect_gp = gpar(col = "white", lwd = 0.2),
                       col = colorRamp2(c(-1.5,0,1.5), c(ColColor$`Low-Blue`[9],'white',ColColor$`High-Red`[9])),
                       height = unit(4, 'cm'), width = unit(193/15, 'cm'))
hm.plot.pro <- Heatmap(t(scale(t(hm.input$Pro[, hm.vis$Sample]))), cluster_columns = F, cluster_rows = F, show_heatmap_legend = F,
                       show_column_names = F, show_row_names = F,  row_title = 'Protein',na_col = 'grey90',
                       column_split = hm.vis$MOcluster, column_gap = unit(2,'mm'),
                       rect_gp = gpar(col = "white", lwd = 0.2),
                       col = colorRamp2(c(-1.5,0,1.5), c(ColColor$`Low-Blue`[9],'white',ColColor$`High-Red`[9])),
                       height = unit(4, 'cm'), width = unit(193/15, 'cm'))
hm.plot.phos <- Heatmap(t(scale(t(hm.input$Phos[, hm.vis$Sample]))), cluster_columns = F, cluster_rows = F, show_heatmap_legend = F,
                        show_column_names = F, show_row_names = F,  row_title = 'Phosphosite',na_col = 'grey90',
                        column_split = hm.vis$MOcluster, column_gap = unit(2,'mm'),
                        rect_gp = gpar(col = "white", lwd = 0.2),
                        col = colorRamp2(c(-1.5,0,1.5), c(ColColor$`Low-Blue`[9],'white',ColColor$`High-Red`[9])),
                        height = unit(4, 'cm'), width = unit(193/15, 'cm'))
draw(hm.plot.mut %v% hm.plot.cnv %v% hm.plot.rna %v% hm.plot.pro %v% hm.plot.phos,
     heatmap_legend_side = 'left', annotation_legend_side = 'left') 
draw(hm.plot.mut %v% hm.plot.cnv)

##### 3.5.5 find features of each omic each cluster [Annotations beside the heatmap] #####
# search method
arrange(filter(input.dea$Phos, 
               str_split_fixed(Feature, fixed(':'), 2)[,1] %in% filter(TotalPathway, 
                                                                       Pathway %in% c('KEGG N-Glycan biosynthesis','REACTOME Sialic acid metabolism','REACTOME O-linked glycosylation',
                                                                                      'HALLMARK_GLYCOLYSIS','REACTOME ERBB2 Activates PTK6 Signaling',
                                                                                      'REACTOME Serotonin Neurotransmitter Release Cycle'))$Gene,
               belong == 4), desc(`MO4-logFC`))$Feature
##### mutation
# C3: CTNNB1
# C4: CHD4
##### CNV
# C1: MYC, BRAF, BCL2
# C4: KRAS, PTK6, FGFR4
##### mRNA 
# C1: KRT6A, S100A8
# C2: SFRP2(Wnt), FBLN1(cell adhesion and migration along protein fibers)
# C3: A(regulate expression of CYPs), ALDH3A1
# C4: MUC4(immune barrier, siglac ligand), ERBB2
##### Protein
# C1: MCM2(cell cycle), PSMB5(immunoproteasome)
# C2: JAK2, AKT3
# C3: AMACR(regulating FFA biosynthesis), PRDX5(relieve ROS)
# C4: MUC5AC, GALNT4(PD-L1 Glycan)
##### Phosphorylation
# C1: (pro-cell cycle) CDK1:T14, CDC23:T596
# C2: VIM:S325(pro EMT and migration), GAB1:S645(activate PI3K/AKT and MAPK through binding GRB2)
# C3: PDHA1:S293(pro-OXPHOS), BCKDHA:S347(BCAA metabolism)
# C4: no clear functions

##### 3.6 Function #####
library(GSVA)
# Prolevel-phos
hm.mtx.prophos <- GBC_Main_ProLevelPhos_knn[, na.omit(index.final$Pro)]
colnames(hm.mtx.prophos) = sapply(colnames(hm.mtx.prophos), function(x){filter(index.final, Pro == x)$Sample})
# GSVA
gsva.rna <- gsva(as.matrix(hm.mtx.rna), TotalPathwayGSVA, method = 'gsva') %>% t() %>% scale() %>% t() %>% as.data.frame()
gsva.pro <- gsva(as.matrix(hm.mtx.pro), TotalPathwayGSVA, method = 'gsva') %>% t() %>% scale() %>% t() %>% as.data.frame()
gsva.phos <- gsva(as.matrix(hm.mtx.prophos), TotalPathwayGSVA, method = 'gsva') %>% t() %>% scale() %>% t() %>% as.data.frame()
# DEA
# rownames(CM.hm) = index.final$Sample
gsva.DEA.rna <- bioinfoamateur::core_Differential_analysis_continuous(CM.hm[colnames(gsva.rna),], gsva.rna, log = T, p.adj = F, show.belong = T)
gsva.DEA.pro <- bioinfoamateur::core_Differential_analysis_continuous(CM.hm[colnames(gsva.pro),], gsva.pro, log = T, p.adj = F, show.belong = T)
gsva.DEA.phos <- bioinfoamateur::core_Differential_analysis_continuous(CM.hm[colnames(gsva.phos),], gsva.phos, log = T, p.adj = F, show.belong = T)
# add feature name
gsva.DEA.rna$Pathway <- rownames(gsva.DEA.rna)
gsva.DEA.pro$Pathway <- rownames(gsva.DEA.pro)
gsva.DEA.phos$Pathway <- rownames(gsva.DEA.phos)
# check belong
list(RNA = table(gsva.DEA.rna$belong),
     Protein = table(gsva.DEA.pro$belong),
     ProPhos = table(gsva.DEA.phos$belong))
##### pathway selection
pathway.selection <- c('KEGG Cell cycle','KEGG Mismatch repair','REACTOME DNA Replication',
                       'KEGG Autophagy','REACTOME Antigen processing','KEGG NF-kappa B signaling pathway', # MO1: proliferation - immune hot
                       'KEGG JAK-STAT signaling pathway','KEGG cGMP-PKG signaling pathway',
                       'KEGG Ras signaling pathway','KEGG ErbB signaling pathway',
                       'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','REACTOME Extracellular matrix organization', # MO2: signaling transduction - CAF dominant
                       'HALLMARK_FATTY_ACID_METABOLISM','HALLMARK_BILE_ACID_METABOLISM',
                       'HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_MYC_TARGETS_V1',
                       'KEGG Valine, leucine and isoleucine degradation','KEGG Peroxisome', # MO3: Metabolism - immune desert
                       'KEGG N-Glycan biosynthesis','REACTOME Sialic acid metabolism','REACTOME O-linked glycosylation',
                       'HALLMARK_GLYCOLYSIS','REACTOME ERBB2 Activates PTK6 Signaling',
                       'REACTOME Serotonin Neurotransmitter Release Cycle') # MO4: ERBB2-Driven - immune inhibition
##### Making input
input.path <- expand.grid(Pathway = pathway.selection, Omics = c('mRNA','Proteomics','Phosphoproteomics'),
                          Cluster = c('MO1','MO2','MO3','MO4')) %>% as.tibble()
input.path$Pathway = as.character(input.path$Pathway)
input.path$Omics = as.character(input.path$Omics)
input.path$Cluster = as.character(input.path$Cluster)
input.path$Score = 0
for (i in 1:nrow(input.path)) {
  # RNA
  if (as.character(input.path$Omics[i]) == 'mRNA') {
    input.path$Score[i] = gsva.rna[input.path$Pathway[i],filter(index.final, MOcluster == input.path$Cluster[i], !is.na(RNA))$Sample] %>% 
      as.numeric() %>% mean()
  }
  # Pro
  if (as.character(input.path$Omics[i]) == 'Proteomics') {
    input.path$Score[i] = gsva.pro[input.path$Pathway[i],filter(index.final, MOcluster == input.path$Cluster[i])$Sample] %>% 
      as.numeric() %>% mean()
  }
  # Phos
  if (as.character(input.path$Omics[i]) == 'Phosphoproteomics') {
    input.path$Score[i] = gsva.phos[input.path$Pathway[i],filter(index.final, MOcluster == input.path$Cluster[i])$Sample] %>% 
      as.numeric() %>% mean()
  }
}
# Transform
input.path$Pathway = factor(input.path$Pathway, levels = pathway.selection, ordered = T)
input.path$Omics = factor(input.path$Omics, levels = c('mRNA','Proteomics','Phosphoproteomics'), ordered = T)
input.path$RowName = paste0(input.path$Cluster,'_',input.path$Omics)
input.path$RowName = factor(input.path$RowName, levels = c('MO1_mRNA','MO1_Proteomics','MO1_Phosphoproteomics',
                                                           'MO2_mRNA','MO2_Proteomics','MO2_Phosphoproteomics',
                                                           'MO3_mRNA','MO3_Proteomics','MO3_Phosphoproteomics',
                                                           'MO4_mRNA','MO4_Proteomics','MO4_Phosphoproteomics'),
                            ordered = T)
input.path2 <- dcast(input.path, RowName ~ Pathway, value.var = 'Score') %>% 
  column_to_rownames(var = 'RowName') %>% t() %>% scale()
input.path2 = input.path2[pathway.selection,]
##### Vis
anno.cluster = rowAnnotation(Cluster = c(rep('MO1',6),rep('MO2',6),rep('MO3',6),rep('MO4',6)),
                             col = list(
                               Cluster = col.mo
                             ),gp = gpar(col = "white"),
                             show_annotation_name = F,
                             simple_anno_size = unit(1/1.5, "cm"))
anno.omics = columnAnnotation(Omics = rep(c('mRNA','Protein','Phosphoprotein'), 4),
                              col = list(
                                Omics = c('mRNA' = '#F1A42A', 'Protein' = '#4E185F', 'Phosphoprotein' = '#20706C')
                              ),gp = gpar(col = "white"),
                              show_annotation_name = F,
                              simple_anno_size = unit(1/1.5, "cm"))

Heatmap(t(scale(t(input.path2))), cluster_rows = F, cluster_columns = F, name = 'Pathway Activity', row_names_side = 'left',
        left_annotation = anno.cluster, top_annotation = anno.omics, show_column_names = F,
        col = circlize::colorRamp2(c(-2, 0, 2), c(ColColor$`Low-Indigo`[5], "white", ColColor$`High-Red`[5])),
        rect_gp = gpar(col = "white", lwd = 2), na_col = 'grey', column_names_rot = 30,
        column_split = c(rep('MO1',3),rep('MO2',3),rep('MO3',3),rep('MO4',3)), row_title = NULL,
        row_split = c(rep('MO1',6),rep('MO2',6),rep('MO3',6),rep('MO4',6)), 
        height = unit(16,'cm'), width = unit(8,'cm')) 
# data for Table S7
write.xlsx(as.data.frame(input.path2), file = 'Func heatmap.xlsx', rowNames = T, colNames = T, overwrite = T)

##### 3.7 rm current variables #####
input.raw <- list(Mut = ExtendMatrix(hm.mtx.mut, hm.mtx.pro),
                  CNV = ExtendMatrix(hm.mtx.cnv, hm.mtx.pro),
                  RNA = ExtendMatrix(hm.mtx.rna, hm.mtx.pro),
                  Pro = ExtendMatrix(hm.mtx.pro, hm.mtx.pro),
                  Phos = ExtendMatrix(hm.mtx.phos, hm.mtx.pro))
input.dea <- list(CNV = hm.DEA.cnv,
                  RNA = hm.DEA.rna,
                  Pro = hm.DEA.pro,
                  Phos = hm.DEA.phos)
rm(anno.cluster, anno.omics, hm.anno.sample, hm.DEA.cnv, hm.DEA.phos, hm.DEA.pro, hm.DEA.rna,hm.input,
   hm.mtx.cnv, hm.mtx.mut, hm.mtx.phos, hm.mtx.pro, hm.mtx.mut.select, hm.mtx.prophos, hm.mtx.rna,
   hm.vis, input.path, input.path2, mo.cnv, mo.rna, mo.pro, mo.phos, mo.mut, mo.cnv, mo.rna, mo.pro, mo.phos,
   i, ml.xgb.pred.dmtx, mutation_rate, temp.tag, pathway.selection,
   hm.plot.mut, hm.plot.cnv, hm.plot.rna, hm.plot.pro, hm.plot.phos, freq_by_cluster)
gc()

########## Mission 4. Genomic features ##########
# target: TP53 mut, ERBB3 mut, ERBB2 Amp, MYC Amp

##### 4.1 differential mutations of each clusters #####
##### prep: create a df just like temp.vis.mut
temp.bar.mut <- pmin(table(GBC_Main_Mutation$Hugo_Symbol, GBC_Main_Mutation$Tumor_Sample_Barcode), 1) %>% as.data.frame()
temp.bar.mut <- dcast(temp.bar.mut, formula = Var1 ~ Var2)
rownames(temp.bar.mut) <- temp.bar.mut$Var1; temp.bar.mut <- dplyr::select(temp.bar.mut, -Var1)
temp.bar.mut <- temp.bar.mut[, filter(index.final, !is.na(WES))$Sample]

# mutation selection by freq
temp.bar.mut.stat <- data.frame(Gene = rownames(temp.bar.mut), p.value = NA, 
                                MO1 = NA, MO2 = NA, MO3 = NA, MO4 = NA,
                                Highest_Frequency_Cluster = NA,
                                Highest_Frequency = NA) %>% as.tibble()
temp.tag <- sapply(colnames(temp.bar.mut), function(x){filter(index.final, Sample == x)$MOcluster})
# 
for (i in 1:nrow(temp.bar.mut.stat)) {
  gene = temp.bar.mut.stat$Gene[i]
  gene_mut <- temp.bar.mut[i, ] %>% as.numeric()
  # 
  freq_by_cluster <- tapply(gene_mut, temp.tag, function(x) sum(x) / length(x))
  # add each freq
  temp.bar.mut.stat$MO1[i] = as.numeric(freq_by_cluster)[1]
  temp.bar.mut.stat$MO2[i] = as.numeric(freq_by_cluster)[2]
  temp.bar.mut.stat$MO3[i] = as.numeric(freq_by_cluster)[3]
  temp.bar.mut.stat$MO4[i] = as.numeric(freq_by_cluster)[4]
  
  # 
  highest_cluster <- names(which.max(freq_by_cluster))
  highest_frequency <- max(freq_by_cluster, na.rm = TRUE)
  # for chisq.test
  group_table <- table(temp.tag, gene_mut)
  if (all(dim(group_table) > 1)) {
    temp.bar.mut.stat$p.value[i] <- fisher.test(group_table)$p.value
  } else {
    temp.bar.mut.stat$p.value[i] <- NA
  }
  # add freq
  temp.bar.mut.stat$Highest_Frequency_Cluster[temp.bar.mut.stat$Gene == gene] <- highest_cluster
  temp.bar.mut.stat$Highest_Frequency[temp.bar.mut.stat$Gene == gene] <- highest_frequency
  rm(i, gene, gene_mut, group_table, highest_cluster, highest_frequency)
}
# adjust pval
temp.bar.mut.stat$padj <- p.adjust(temp.bar.mut.stat$p.value, method = 'fdr')
# add total freq
temp.bar.mut.stat$Freq <- rowSums(temp.bar.mut) / ncol(temp.bar.mut)

##### 4.2 differential CNVs of each clusters #####
##### prep: create a df just like temp.vis.mut
temp.bar.cnv <- GBC_Main_CNV_cate
colnames(temp.bar.cnv) <- sapply(colnames(temp.bar.cnv), function(x){filter(index.final, WES == x)$Sample})
temp.bar.cnv <- temp.bar.cnv[, filter(index.final, !is.na(WES))$Sample]
# mutation selection by freq
temp.bar.cnv.stat <- data.frame(Gene = rownames(temp.bar.cnv), p.value = NA, 
                                MO1 = NA, MO2 = NA, MO3 = NA, MO4 = NA,
                                Type = NA, # CNV amp or del
                                Highest_Frequency_Cluster = NA,
                                Highest_Frequency = NA) %>% as.tibble()
temp.tag <- sapply(colnames(temp.bar.cnv), function(x){filter(index.final, Sample == x)$MOcluster})

# 
temp.tag.cluster <- unique(levels(temp.tag))

# calc params
for (gene in rownames(temp.bar.cnv)) {
  cnv_values <- temp.bar.cnv[gene, ]
  #
  deletion_count <- sum(cnv_values %in% c(-2, -1))
  amplification_count <- sum(cnv_values %in% c(1, 2))
  if (deletion_count > amplification_count) {
    gene_type <- "Deletion"
  } else if (amplification_count > deletion_count) {
    gene_type <- "Amplification"
  } else {
    gene_type <- "Balanced"
  }
  temp.bar.cnv.stat$Type[temp.bar.cnv.stat$Gene == gene] <- gene_type
  # 
  freq_values <- c()  # 
  for (i in 1:length(temp.tag.cluster)) {
    cluster = temp.tag.cluster[i]
    cluster_samples <- which(temp.tag == cluster)
    if (gene_type == "Deletion") {
      freq <- sum(cnv_values[cluster_samples] == -2) / length(cluster_samples)
    } else if (gene_type == "Amplification") {
      freq <- sum(cnv_values[cluster_samples] == 2) / length(cluster_samples)
    } else {
      freq <- 0
    }
    freq_values <- c(freq_values, freq)
    temp.bar.cnv.stat[temp.bar.cnv.stat$Gene == gene, (i+2)] <- as.numeric(freq)
  }
  # 
  group_table <- table(temp.tag, cnv_values %in% c(-2, 2))
  if (all(dim(group_table) > 1)) {
    temp.bar.cnv.stat$p.value[temp.bar.cnv.stat$Gene == gene] <- chisq.test(group_table)$p.value
  } else {
    temp.bar.cnv.stat$p.value[temp.bar.cnv.stat$Gene == gene] <- NA
  }
  # 
  temp.bar.cnv.stat$Highest_Frequency[temp.bar.cnv.stat$Gene == gene] <- max(freq_values, na.rm = TRUE)
  temp.bar.cnv.stat$Highest_Frequency_Cluster[temp.bar.cnv.stat$Gene == gene] <- as.character(temp.tag.cluster)[which.max(freq_values)]
  # rm
  rm(cnv_values, deletion_count, amplification_count, gene_type, freq_values, freq, group_table, i)
}
# 
head(temp.bar.cnv.stat)
# adjust pval
temp.bar.cnv.stat$padj <- p.adjust(temp.bar.cnv.stat$p.value, method = 'fdr')

#
rm(temp.bar.mut, temp.bar.cnv, temp.tag, temp.tag.cluster)

##### 4.3 Filtering genomic events #####
##### Mutation
### record
# Necessary: ELF3(in MO4, promoting glycolysis), SMAD4(in MO4, promoting glycolysis), CTNNB1(in MO3, promoting FFA metabolism)
# pick manually
filter(temp.bar.mut.stat, Gene %in% c('ELF3','SMAD4','CTNNB1'))
# pick by cluster
filter(temp.bar.mut.stat, Highest_Frequency_Cluster == 'MO3', p.value <= 0.05) %>%
  arrange(desc(Freq))
# pick by freq and cluster
filter(temp.bar.mut.stat, Freq >= 0.05, p.value <= 0.05, Highest_Frequency >= 0.1)

# final pick
temp.bar.mut.stat.filter <- filter(temp.bar.mut.stat, Gene %in% c('ELF3','SMAD4','CTNNB1')) %>%
  arrange(Highest_Frequency_Cluster)
temp.bar.mut.stat.filter
# recording
# ELF3: in MO3 and MO4, promoting glycolysis
# SMAD4: in MO4, promoting glycolysis
# CTNNB1: in MO3, promoting FFA metabolism

##### CNV
### record
# Necessary: ERBB2、MO3 metabolic genes
# pick manually
filter(temp.bar.cnv.stat, Gene %in% c('MYC','ERBB2','AMACR'))
# pick by freq and cluster
arrange(filter(temp.bar.cnv.stat, p.value <= 0.05, 
               Highest_Frequency >= 0.1, Highest_Frequency_Cluster == 'MO3'),
        desc(Highest_Frequency))$Gene
# final pick
temp.bar.cnv.stat.filter <- filter(temp.bar.cnv.stat, Gene %in% c('ERBB2','AMACR')) %>%
  arrange(Highest_Frequency_Cluster)
temp.bar.cnv.stat.filter


##### 4.4 barplot visualization #####
VisGenomeEvent <- function(gene, Freqdf){
  plot_data <- Freqdf %>%
    filter(Gene == gene) %>%
    select(Gene, MO1, MO2, MO3, MO4) %>%
    pivot_longer(cols = starts_with("M"), names_to = "Cluster", values_to = "Frequency") %>%
    mutate(Cluster = factor(Cluster, levels = c("MO4", "MO3", "MO2", "MO1"), ordered = T))
  # 
  ggplot(plot_data, aes(x = Frequency, y = Cluster)) +
    geom_bar(stat = "identity", fill = col.mo, width = 0.7) +
    geom_text(aes(label = scales::percent(Frequency, accuracy = 1)), 
              hjust = -0.2, size = 4.5, color = "black", family = "Arial") +
    labs(
      title = paste(gene, "mutation"),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      text = element_text(family = "Arial", color = "black"),
      plot.title = element_text(hjust = 0.5, face = "italic"),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),          
      axis.ticks = element_line(color = "black")          
    ) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1))
}
plot_grid(VisGenomeEvent(gene = 'CTNNB1', Freqdf = temp.bar.mut.stat.filter),
          VisGenomeEvent(gene = 'SMAD4', Freqdf = temp.bar.mut.stat.filter),
          VisGenomeEvent(gene = 'ELF3', Freqdf = temp.bar.mut.stat.filter),
          VisGenomeEvent(gene = 'ERBB2', Freqdf = temp.bar.cnv.stat.filter),
          nrow = 2)
#
rm(temp.bar.cnv.stat, temp.bar.mut.stat,
   temp.bar.cnv.stat.filter, temp.bar.mut.stat.filter, gene, freq_by_cluster, temp.vis.mut)

########## Mission 5. scores: NE, EMT, TMB, TNB, CIN, Immune score ##########
##### 5.1.1 calc NE and add to index.final #####
##### 5.1.1 geneset prep
ref.ne <- list(ref.7 = c('TFF3','ENO2','BEX1','MAP10','MYCN','KIF1A','RGS7'),
               ref.10 = c("SCG3", "CHGA", "CHGB", "CHRNB2", "PCSK1",
                          "ELAVL4", "ENO2", "SCN3A", "SYP", "NKX2-1"),
               ref.50 = c("BEX1", "ASCL1", "INSM1", "CHGA", "TAGLN3", "KIF5C", "CRMP1", "SCG3", "SYT4", "RTN1", "MYT1",
                          "SYP", "KIF1A", "TMSB15A", "SYN1", "SYT11", "RUNDC3A", "TFF3", "CHGB", "FAM57B", "SH3GL2",
                          "BSN", "SEZ6", "TMSB15B", "CELF3"),
               ref.70 = c("ASXL3", "CAND2", "ETV5", "GPX2", "JAKMIP2", "KIAA0408", "SOGA3", "TRIM9", "BRINP1", "C7orf76",
                          "GNAO1", "KCNB2", "KCND2", "LRRC16B", "MAP10", "NRSN1", "PCSK1", "PROX1", "RGS7", "SCG3",
                          "SEC11C", "SEZ6", "ST8SIA3", "SVOP", "SYT11", "AURKA", "DNMT1", "EZH2", "MYCN"))
##### 5.1.2 gsva
ref.ne.gsva <- data.frame(NE7 = c(ref.ne$ref.7, rep(NA, length(ref.ne$ref.70)-length(ref.ne$ref.7))),
                          NE10 = c(ref.ne$ref.10, rep(NA, length(ref.ne$ref.70)-length(ref.ne$ref.10))),
                          NE50 = c(ref.ne$ref.50, rep(NA, length(ref.ne$ref.70)-length(ref.ne$ref.50))),
                          NE70 = c(ref.ne$ref.70))
ssgsea.ne <- GSVA::gsva(as.matrix(GBC_Main_RNA_TPM), ref.ne.gsva, method = 'ssgsea') %>% 
  t() %>% scale() %>% t() %>% as.data.frame()
# add to index
index.final$NE7 <- sapply(index.final$RNA, function(x){ifelse(!is.na(x),unname(as.numeric(ssgsea.ne['NE7',x])), NA)})
index.final$NE10 <- sapply(index.final$RNA, function(x){ifelse(!is.na(x),unname(as.numeric(ssgsea.ne['NE10',x])), NA)})
index.final$NE50 <- sapply(index.final$RNA, function(x){ifelse(!is.na(x),unname(as.numeric(ssgsea.ne['NE50',x])), NA)})
index.final$NE70 <- sapply(index.final$RNA, function(x){ifelse(!is.na(x),unname(as.numeric(ssgsea.ne['NE70',x])), NA)})

##### 5.1.2 calc invasion score #####
ref.invasion <- c("ADAMTS1", "ALDH3A1", "ANGPT1", "ANGPTL4", "CASP8", "CCNE2", "CCR7", "CD44", "CD82", 
                  "CDH1", "CDH11", "CDH2", "CDH6", "CLDN7", "COL1A1", "COL4A2", "COL6A1", "CSF1", 
                  "CSF3", "CST7", "CTGF", "CTSB", "CTSD", "CTSK", "CTSL1", "CXCL1", "CXCL12", "CXCR4", 
                  "CXCR6", "DRG1", "EREG", "FGF8", "FLT1", "FLT4", "GPI", "GSN", "HGF", "HIF1A", 
                  "HMGB1", "ID1", "IGFBP7", "IL13RA2", "ISG20", "JAG1", "KISS1", "KLRC2", "KYNU", 
                  "LTBP1", "MAP2K4", "MAP2K5", "MAP2K7", "MCAM", "MET", "METAP2", "MMP1", "MMP10", 
                  "MMP11", "MMP13", "MMP14", "MMP2", "MMP7", "MYC", "NEDD9", "NF2", "NME1", "NME2", 
                  "NME4", "PAX5", "PDGFA", "PLAUR", "PTGS2", "RUNX1", "SERPINE1", "SERPINB5", "SOX4", 
                  "SPARC", "SPP1", "SRC", "SYK", "TFF1", "TGFB1", "TIMP1", "TIMP2", "TIMP3", "TIMP4", 
                  "TNC", "TP53", "VEGFA")
ref.invasion.gsva <- data.frame(Invasion = ref.invasion)
# data for table S7
write.xlsx(ref.invasion.gsva, file = 'Invasion geneset.xlsx', rowNames = T, colNames = T, overwrite = T)
# ssgsea
ssgsea.invasion <- GSVA::gsva(as.matrix(GBC_Main_Pro), ref.invasion.gsva, method = 'ssgsea') %>% 
  t() %>% scale() %>% t() %>% as.data.frame()
# add to index
index.final$InvasionScore <- sapply(index.final$Pro, function(x){ifelse(!is.na(x),unname(as.numeric(ssgsea.invasion['Invasion',x])), NA)})
# vis
temp.comp <- list(c('MO2','MO1'),c('MO3','MO1'),c('MO4','MO1'),
                  c('MO3','MO2'),c('MO2','MO4'),c('MO4','MO3'))
ggboxplot(index.final, 'MOcluster', 'InvasionScore') + 
  stat_compare_means(comparisons = temp.comp, label.y = c(6,5.5,5),
                     method = 'wilcox') +
  stat_compare_means(method = 'kruskal')

##### 5.1.3 extract EMT score #####
index.final$EMTScore <- as.numeric(gsva.pro['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',])
# vis
ggboxplot(index.final, 'MOcluster', 'EMTScore') + 
  stat_compare_means(comparisons = temp.comp, label.y = c(6,5.5,5),
                     method = 'wilcox') +
  stat_compare_means(method = 'kruskal')

##### 5.1.x obv primary (MO4 top?) #####
temp.comp <- list(c('MO4','MO1'),
                  c('MO4','MO2'),
                  c('MO4','MO3'))
plot_grid(ggboxplot(index.final, 'MOcluster', 'NE7') + 
            stat_compare_means(comparisons = temp.comp, label.y = c(6,5.5,5),
                               method = 'wilcox') +
            stat_compare_means(method = 'kruskal'),
          ggboxplot(index.final, 'MOcluster', 'NE10') + 
            stat_compare_means(comparisons = temp.comp, label.y = c(6,5.5,5),
                               method = 'wilcox') +
            stat_compare_means(method = 'kruskal'),
          ggboxplot(index.final, 'MOcluster', 'NE50') + 
            stat_compare_means(comparisons = temp.comp, label.y = c(6,5.5,5),
                               method = 'wilcox') +
            stat_compare_means(method = 'kruskal'),
          ggboxplot(index.final, 'MOcluster', 'NE70', color = 'MOcluster') + 
            stat_compare_means(comparisons = temp.comp, label.y = c(6,5.5,5),
                               method = 'wilcox') +
            stat_compare_means(method = 'kruskal') +
            scale_color_manual(values = col.mo),
          nrow = 1)
# conclusion: use NE7 as the final NE score

##### 5.2 vis all scores #####
temp.comp <- list(c('MO2','MO1'),c('MO3','MO1'),c('MO4','MO1'),
                  c('MO3','MO2'),c('MO2','MO4'),c('MO4','MO3'))
##### TMB and TNB
# TMB
ggboxplot(index.final, 'MOcluster', 'TMB', color = 'MOcluster') +
  stat_compare_means(method = 'anova') +
  scale_color_manual(values = col.mo) +
  labs(x = NULL, y = 'TMB') +
  theme(panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA),   
        panel.grid = element_blank(),                                
        axis.line = element_line(color = "black"),                   
        axis.ticks = element_line(color = "black"),                  
        axis.text = element_text(family = "Arial", color = "black", size = 10),  
        axis.title = element_text(family = "Arial", color = "black", size = 12), 
        strip.text = element_text(family = "Arial", color = "black", size = 12), 
        plot.title = element_text(family = "Arial", color = "black", size = 14, hjust = 0.5), 
        legend.position = 'none')  
# TNB
ggboxplot(index.final, 'MOcluster', 'TNB', color = 'MOcluster') +
  stat_compare_means(method = 'anova') +
  scale_color_manual(values = col.mo) +
  labs(x = NULL, y = 'Tumor neoantigen burden') +
  theme(panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA),  
        panel.grid = element_blank(),                               
        axis.line = element_line(color = "black"),                   
        axis.ticks = element_line(color = "black"),                  
        axis.text = element_text(family = "Arial", color = "black", size = 10),  
        axis.title = element_text(family = "Arial", color = "black", size = 12), 
        strip.text = element_text(family = "Arial", color = "black", size = 12), 
        plot.title = element_text(family = "Arial", color = "black", size = 14, hjust = 0.5), 
        legend.position = 'none')  

##### CIN
ggboxplot(index.final, 'MOcluster', 'CIN', color = 'MOcluster') +
  stat_compare_means(comparisons = temp.comp, label.y = c(seq(max(index.final$CIN, na.rm = T)+0.1,
                                                              max(index.final$CIN, na.rm = T)+5,
                                                              length.out = 6)),
                     method = 't.test') +
  stat_compare_means(method = 'anova') +
  scale_color_manual(values = col.mo) +
  labs(x = NULL, y = 'Chromosome instability index') +
  theme(panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA),  
        panel.grid = element_blank(),                                
        axis.line = element_line(color = "black"),                   
        axis.ticks = element_line(color = "black"),                  
        axis.text = element_text(family = "Arial", color = "black", size = 10),  
        axis.title = element_text(family = "Arial", color = "black", size = 12), 
        strip.text = element_text(family = "Arial", color = "black", size = 12), 
        plot.title = element_text(family = "Arial", color = "black", size = 14, hjust = 0.5), 
        legend.position = 'none')  

##### NE
ggboxplot(index.final, 'MOcluster', 'NE7', color = 'MOcluster') +
  stat_compare_means(comparisons = temp.comp, label.y = c(seq(max(index.final$NE7, na.rm = T)+0.1,
                                                              max(index.final$NE7, na.rm = T)+1.5,
                                                              length.out = 6)),
                     method = 't.test') +
  stat_compare_means(method = 'anova') +
  scale_color_manual(values = col.mo) +
  labs(x = NULL, y = 'NE score') +
  theme(panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA),   
        panel.grid = element_blank(),                                
        axis.line = element_line(color = "black"),                   
        axis.ticks = element_line(color = "black"),                  
        axis.text = element_text(family = "Arial", color = "black", size = 10),  
        axis.title = element_text(family = "Arial", color = "black", size = 12), 
        strip.text = element_text(family = "Arial", color = "black", size = 12), 
        plot.title = element_text(family = "Arial", color = "black", size = 14, hjust = 0.5),
        legend.position = 'none') 

##### EMT and Invasion
ggboxplot(index.final, 'MOcluster', 'EMTScore', color = 'MOcluster') +
  stat_compare_means(comparisons = temp.comp, label.y = c(seq(max(index.final$EMTScore, na.rm = T)+0.1,
                                                              max(index.final$EMTScore, na.rm = T)+1.5,
                                                              length.out = 6)),
                     method = 't.test') +
  stat_compare_means(method = 'anova') +
  scale_color_manual(values = col.mo) +
  labs(x = NULL, y = 'EMT score') +
  theme(panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA),  
        panel.grid = element_blank(),                                
        axis.line = element_line(color = "black"),                   
        axis.ticks = element_line(color = "black"),                 
        axis.text = element_text(family = "Arial", color = "black", size = 10),  
        axis.title = element_text(family = "Arial", color = "black", size = 12), 
        strip.text = element_text(family = "Arial", color = "black", size = 12), 
        plot.title = element_text(family = "Arial", color = "black", size = 14, hjust = 0.5), 
        legend.position = 'none') 
#
ggboxplot(index.final, 'MOcluster', 'InvasionScore', color = 'MOcluster') +
  stat_compare_means(comparisons = temp.comp, label.y = c(seq(max(index.final$InvasionScore, na.rm = T)+0.1,
                                                              max(index.final$InvasionScore, na.rm = T)+1.5,
                                                              length.out = 6)),
                     method = 't.test') +
  stat_compare_means(method = 'anova') +
  scale_color_manual(values = col.mo) +
  labs(x = NULL, y = 'Invasion score') +
  theme(panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA),   
        panel.grid = element_blank(),                                
        axis.line = element_line(color = "black"),                   
        axis.ticks = element_line(color = "black"),                  
        axis.text = element_text(family = "Arial", color = "black", size = 10),  
        axis.title = element_text(family = "Arial", color = "black", size = 12), 
        strip.text = element_text(family = "Arial", color = "black", size = 12), 
        plot.title = element_text(family = "Arial", color = "black", size = 14, hjust = 0.5), 
        legend.position = 'none') 

########## Mission 6. EMT explaining ##########
##### 6.1 calc partial and complete EMT and vis #####
ref.pEMT <- read.xlsx('../pEMT_genes.xlsx', colNames = T)
ssgsea.pEMT <- GSVA::gsva(as.matrix(input.raw$Pro), ref.pEMT, method = 'ssgsea') %>% 
  t() %>% scale() %>% t() %>% as.data.frame()
# add
index.obv <- select(index.final,1,2,3,4,15,29)
index.obv$pEMT <- sapply(index.obv$Sample, function(x){ifelse(!is.na(x),unname(as.numeric(ssgsea.pEMT['common_pMET_genes',x])), NA)})
index.final$pEMT = index.obv$pEMT
# vis
ggboxplot(index.obv, 'MOcluster', 'pEMT', color = 'MOcluster') +
  stat_compare_means(method = 'anova') +
  scale_color_manual(values = col.mo) +
  labs(x = NULL, y = 'pEMT') +
  theme(panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA),  
        panel.grid = element_blank(),                               
        axis.line = element_line(color = "black"),                  
        axis.ticks = element_line(color = "black"),                  
        axis.text = element_text(family = "Arial", color = "black", size = 10),  
        axis.title = element_text(family = "Arial", color = "black", size = 12), 
        strip.text = element_text(family = "Arial", color = "black", size = 12), 
        plot.title = element_text(family = "Arial", color = "black", size = 14, hjust = 0.5), 
        legend.position = 'none')  # score pEMT 4*4p

##### 6.2 find potential pathway associated with EMT: TGF-b #####
##### pathway selection
temp.pathway <- grep('signal', filter(gsva.DEA.pro, belong == 2)$Pathway, ignore.case = T, value = T) %>% unique()
temp.pathway <- grep('KEGG', temp.pathway, ignore.case = T, value = T) %>% unique()
##### calc cor
temp.cor <- gsva.pro[c('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', temp.pathway),] %>% t()
temp.cor.cor <- cor(temp.cor, method = 'spearman')
##### vis prep
index.obv <- select(index.final,1,2,3,4,15,29)
index.obv$EMTScore <- sapply(index.obv$Sample, function(x){as.numeric(gsva.pro['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',x])}) %>% as.numeric()
index.obv$TGFbScore <- sapply(index.obv$Sample, function(x){as.numeric(gsva.pro['KEGG TGF-beta signaling pathway',x])}) %>% as.numeric()
index.final$EMTScore = index.obv$EMTScore
index.final$TGFbScore = index.obv$TGFbScore
##### vis
ggplot(index.obv, aes(x = EMTScore, y = TGFbScore, color = MOcluster)) +
  geom_point(size = 3, alpha = 0.7) +  
  stat_ellipse(aes(fill = MOcluster), 
               geom = "polygon", alpha = 0.2, show.legend = FALSE) +  
  scale_color_manual(values = col.mo) + 
  scale_fill_manual(values = col.mo) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  
  stat_cor(method = "spearman", size = 5, label.x.npc = "left") + 
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"), 
    axis.text = element_text(family = "Arial", color = "black", size = 10),
    axis.title = element_text(family = "Arial", color = "black", size = 12),
    legend.title = element_text(family = "Arial", color = "black", size = 10),
    legend.text = element_text(family = "Arial", color = "black", size = 9)
  ) +
  labs(
    title = "Correlation between EMT and TGFb Scores with MOcluster",
    x = "EMT Score",
    y = "TGFb Score",
    color = "MOcluster"
  ) 

##### 6.3 pick crucial genes of partial and complete EMT #####
'MMP3' %in% rownames(input.raw$RNA); 'MMP3' %in% rownames(input.raw$Pro)
## epi
# CDH1 mRNA + Protein (E-cadherin)
# OCLN mRNA + Protein
# DSP mRNA + Protein
# TJP1 mRNA + Protein
## Mes
# CDH2 mRNA + Protein (N-cadherin)
# VIM mRNA + Protein
# TWIST1 mRNA
# MMP2 mRNA + Protein
# MMP3 mRNA + Protein
# MMP9 mRNA
# MMP10 mRNA
##### input prep
identical(colnames(input.raw$RNA), index.final$Sample)
input.emt <- data_frame(`CDH1 mRNA` = as.numeric(input.raw$RNA['CDH1',]),
                        `CDH1 Protein` = as.numeric(input.raw$Pro['CDH1',]),
                        `OCLN mRNA` = as.numeric(input.raw$RNA['OCLN',]),
                        `OCLN Protein` = as.numeric(input.raw$Pro['OCLN',]),
                        `DSP mRNA` = as.numeric(input.raw$RNA['DSP',]),
                        `DSP Protein` = as.numeric(input.raw$Pro['DSP',]),
                        `TJP1 mRNA` = as.numeric(input.raw$RNA['TJP1',]),
                        `TJP1 Protein` = as.numeric(input.raw$Pro['TJP1',]),
                        # M
                        `CDH2 mRNA` = as.numeric(input.raw$RNA['CDH2',]),
                        `CDH2 Protein` = as.numeric(input.raw$Pro['CDH2',]),
                        `VIM mRNA` = as.numeric(input.raw$RNA['VIM',]),
                        `VIM Protein` = as.numeric(input.raw$Pro['VIM',]),
                        `TWIST1 mRNA` = as.numeric(input.raw$RNA['TWIST1',]),
                        `MMP2 mRNA` = as.numeric(input.raw$RNA['MMP2',]),
                        `MMP2 Protein` = as.numeric(input.raw$Pro['MMP2',]),
                        `MMP3 mRNA` = as.numeric(input.raw$RNA['MMP3',]),
                        `MMP3 Protein` = as.numeric(input.raw$Pro['MMP3',]),
                        `MMP9 mRNA` = as.numeric(input.raw$RNA['MMP9',]),
                        `MMP10 mRNA` = as.numeric(input.raw$RNA['MMP10',]))
rownames(input.emt) <- index.final$Sample
input.emt <- as.data.frame(t(input.emt)) %>% t() %>% scale() %>% t() %>% as.data.frame()

##### Vis
temp.hm.anno <- columnAnnotation(`MO cluster` = index.final$MOcluster,
                                 col = list(
                                   `MO cluster` = col.mo
                                 ),
                                 simple_anno_size = unit(0.5, "cm"),
                                 gp = gpar(col = "white", lwd = 0.2),
                                 na_col = 'grey90')
set.seed(1219)
temp.hm.anno.feature <- rowAnnotation(`Mechanism` = c(rep('Epithelial', 8), rep('Mesenchymal', 11)),
                                      col = list(
                                        Mechanism = c('Epithelial'=ColJournal$Science[1],
                                                      'Mesenchymal'=ColJournal$Science[2])
                                      ),
                                      simple_anno_size = unit(0.5, "cm"),
                                      gp = gpar(col = "white", lwd = 0.2),
                                      na_col = 'grey90')
# vis
Heatmap(input.emt, cluster_rows = F, cluster_columns = F, name = 'Scaled abundance',
        top_annotation = temp.hm.anno, show_column_names = F,
        left_annotation = temp.hm.anno.feature,
        column_split = index.final$MOcluster, column_gap = unit(2,'mm'),
        row_split = c(1,1,1,1,1,1,1,1,
                      2,2,2,2,2,2,2,2,2,2,2), row_gap = unit(2,'mm'),
        col = colorRamp2(c(-1.5,0,1.5), c(ColColor$`Low-Blue`[9],'white',ColColor$`High-Red`[9])),
        height = unit(0.5*13, 'cm'), width = unit(129/10, 'cm')) 
#
rm(input.emt, temp.pathway, temp.cor, temp.cor.cor)


########## Mission 7. MYC function ##########
##### 7.1 add MYC phos sites into index #####
index.obv <- select(index.final,1,2,3,4,15,29)
# MYC S62# MYC S62Sample
index.obv$MYC_S62 <- sapply(index.obv$Pro, function(x){as.numeric(GBC_Main_Phos_NA['MYC:S62',x])}) %>% as.numeric()
table(!is.na(index.obv$MYC_S62), index.obv$MOcluster)
# MYC S161
index.obv$MYC_S161 <- sapply(index.obv$Pro, function(x){as.numeric(GBC_Main_Phos_NA['MYC:S161',x])}) %>% as.numeric()
table(!is.na(index.obv$MYC_S161), index.obv$MOcluster)
# MYC S347
index.obv$MYC_S347 <- sapply(index.obv$Pro, function(x){as.numeric(GBC_Main_Phos_NA['MYC:S347',x])}) %>% as.numeric()
table(!is.na(index.obv$MYC_S347), index.obv$MOcluster)
# MYC S348
index.obv$MYC_S348 <- sapply(index.obv$Pro, function(x){as.numeric(GBC_Main_Phos_NA['MYC:S348',x])}) %>% as.numeric()
table(!is.na(index.obv$MYC_S348), index.obv$MOcluster)

##### 7.2 Vis #####
temp.comp <- list(c('MO2','MO1'),c('MO3','MO1'),c('MO4','MO1'),
                  c('MO3','MO2'),c('MO2','MO4'),c('MO4','MO3'))
# S62
ggboxplot(index.obv, 'MOcluster', 'MYC_S62', color = 'MOcluster') +
  stat_compare_means(method = 'kruskal') +
  scale_color_manual(values = col.mo) +
  labs(x = NULL, y = 'MYC S62') +
  theme(panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA),  
        panel.grid = element_blank(),                               
        axis.line = element_line(color = "black"),                   
        axis.ticks = element_line(color = "black"),                  
        axis.text = element_text(family = "Arial", color = "black", size = 10),  
        axis.title = element_text(family = "Arial", color = "black", size = 12), 
        strip.text = element_text(family = "Arial", color = "black", size = 12), 
        plot.title = element_text(family = "Arial", color = "black", size = 14, hjust = 0.5), 
        legend.position = 'none') 
# S161
ggboxplot(index.obv, 'MOcluster', 'MYC_S161', color = 'MOcluster') +
  stat_compare_means(method = 'kruskal') +
  scale_color_manual(values = col.mo) +
  labs(x = NULL, y = 'MYC S161') +
  theme(panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA),  
        panel.grid = element_blank(),                                
        axis.line = element_line(color = "black"),                   
        axis.ticks = element_line(color = "black"),                  
        axis.text = element_text(family = "Arial", color = "black", size = 10),  
        axis.title = element_text(family = "Arial", color = "black", size = 12), 
        strip.text = element_text(family = "Arial", color = "black", size = 12), 
        plot.title = element_text(family = "Arial", color = "black", size = 14, hjust = 0.5), 
        legend.position = 'none') # MYC S161 4*4p
# S347
ggboxplot(index.obv, 'MOcluster', 'MYC_S347', color = 'MOcluster') +
  stat_compare_means(method = 'kruskal') +
  scale_color_manual(values = col.mo) +
  labs(x = NULL, y = 'MYC S347') +
  theme(panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA),   
        panel.grid = element_blank(),                                
        axis.line = element_line(color = "black"),                   
        axis.ticks = element_line(color = "black"),                  
        axis.text = element_text(family = "Arial", color = "black", size = 10),  
        axis.title = element_text(family = "Arial", color = "black", size = 12), 
        strip.text = element_text(family = "Arial", color = "black", size = 12), 
        plot.title = element_text(family = "Arial", color = "black", size = 14, hjust = 0.5), 
        legend.position = 'none') 
# S348
ggboxplot(index.obv, 'MOcluster', 'MYC_S348', color = 'MOcluster') +
  stat_compare_means(method = 'kruskal') +
  scale_color_manual(values = col.mo) +
  labs(x = NULL, y = 'MYC S348') +
  theme(panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA),   
        panel.grid = element_blank(),                                
        axis.line = element_line(color = "black"),                   
        axis.ticks = element_line(color = "black"),                  
        axis.text = element_text(family = "Arial", color = "black", size = 10),  
        axis.title = element_text(family = "Arial", color = "black", size = 12), 
        strip.text = element_text(family = "Arial", color = "black", size = 12), 
        plot.title = element_text(family = "Arial", color = "black", size = 14, hjust = 0.5), 
        legend.position = 'none') # MYC S348 4*4p

#  data for Table S7
write.xlsx(index.final[,c(1,29,30,35,37,41:44)], file = 'Pathway activity mo.xlsx', rowNames = T, colNames = T, overwrite = T)

########## Mission 8. Drug response prediction ##########
library(oncoPredict)

##### 8.1 input prep #####
##### training set
# drug response
op.train.response <- openxlsx::read.xlsx('../GDSC2_fitted_dose_response_24Jul22.xlsx') # GDSC database
# exp and model
load('../Depmap数据集.RData') # Depmap database
op.train.exp <- DepmapRNAseq
op.train.model <- Model
rm(DepmapCNV,DepmapDependency,DepmapEffect,DepmapMutation, DepmapRNAseq, Model, AchillesDependency)
# confirm cell line
intersect(op.train.response$COSMIC_ID, op.train.model$COSMICID)
# exp + metadata
op.train.model = filter(op.train.model, ModelID %in% intersect(op.train.model$ModelID, colnames(op.train.exp)))
op.train.exp = op.train.exp[, intersect(op.train.model$ModelID, colnames(op.train.exp))]

# drug response to matrix
op.train.response.mtx = dcast(op.train.response, DRUG_NAME ~ SANGER_MODEL_ID, value.var = 'AUC', mean, na.rm = T) %>% 
  column_to_rownames(var = 'DRUG_NAME')

# no blood
op.train.model.filter <- filter(op.train.model, !OncotreeLineage %in% c('Myeloid','Lymphoid'),
                                ModelID %in% colnames(op.train.exp),
                                SangerModelID %in% op.train.response$SANGER_MODEL_ID | COSMICID %in% op.train.response$COSMIC_ID | CellLineName %in% op.train.response$CELL_LINE_NAME)
# 
op.train.response.mtx <- op.train.response.mtx[,op.train.model.filter$SangerModelID]
# 
op.train.exp = op.train.exp[,op.train.model.filter$ModelID]
colnames(op.train.exp) = op.train.model.filter$SangerModelID

##### 8.2 prediction #####
calcPhenotype(trainingExprData = op.train.exp[intersect(rownames(GBC_Main_RNA_logTPM), rownames(op.train.exp)), ] %>% as.matrix(),
              trainingPtype = op.train.response.mtx %>% t() %>% as.matrix(),
              testExprData = GBC_Main_RNA_logTPM[intersect(rownames(GBC_Main_RNA_logTPM), rownames(op.train.exp)), ] %>% as.matrix(),
              batchCorrect = 'standardize',
              powerTransformPhenotype = T,
              removeLowVaryingGenes = 0.2,
              removeLowVaringGenesFrom = 'rawData',
              minNumSamples = 10,
              selection = 1,
              printOutput = T,
              pcr = F,
              report_pc = F,
              cc = T,
              percent = 80,
              rsq = T)

##### 8.3 read output #####
op.result <- fread('calcPhenotype_Output/DrugPredictions.csv')
# format arrangement
op.result <- column_to_rownames(op.result, var = 'V1') %>% t() %>% as.data.frame()
colnames(op.result) <- sapply(colnames(op.result),
                              function(x){as.character(filter(index.final, RNA == x)$Sample)}) %>% unlist()

##### 8.4 dea #####
index.drug <- filter(index.final, Sample %in% colnames(op.result))
op.result <- op.result[, index.drug$Sample]
# AUC to sensitivity
op.result <- max(as.matrix(op.result)) - op.result
# data for Table S7
write.xlsx(op.result, file = 'Drug sensitivity.xlsx', rowNames = T, colNames = T, overwrite = T)
#
cm.drug <- data.frame(MO1 = ifelse(index.drug$MOcluster == 'MO1', 1, 0),
                      MO2 = ifelse(index.drug$MOcluster == 'MO2', 1, 0),
                      MO3 = ifelse(index.drug$MOcluster == 'MO3', 1, 0),
                      MO4 = ifelse(index.drug$MOcluster == 'MO4', 1, 0))
# dea
op.dea <- bioinfoamateur::core_Differential_analysis_continuous(cm.drug, op.result, method = 't.test',
                                                                p.adj = T, log = T, show.belong = T)

##### 8.5 pick drugs and draw sensitivity heatmap #####
##### drug
info.drug <- data.frame(Drug = c("Cisplatin", "Topotecan", "Gemcitabine",
                                 "CDK9_5038", "Wee1 Inhibitor", "Pevonedistat",           # MO1 sensitive
                                 "Dasatinib", "Dabrafenib", "Foretinib", "Temsirolimus",
                                 "WIKI4", "XAV939",                                       # MO2 sensitive Wnt inhibitor
                                 'POMHEX',           # ENO2 inhibitor
                                 'Dihydrorotenone',  # mito inhibitor
                                 'Leflunomide'      # indirectly inhibit Gln metabolism
),
Mechanism = c(rep('Chemotherapy', 3), rep('Cell cycle inhibition', 2), rep('Immune activation', 1),
              rep('Signal transduction', 4), rep('EMT', 2),
              rep('Metabolism', 3) ))

##### 8.6 Vis by heatmap #####
op.hm.input <- op.result[info.drug$Drug, index.drug$Sample] %>% t() %>% scale() %>% t() %>% as.data.frame()
op.hm.anno <- columnAnnotation(`MO cluster` = index.drug$MOcluster,
                               col = list(
                                 `MO cluster` = col.mo
                               ),
                               simple_anno_size = unit(0.5, "cm"),
                               gp = gpar(col = "white", lwd = 0.2),
                               na_col = 'grey90')
set.seed(1219)
op.hm.anno.drug <- rowAnnotation(`Mechanism` = info.drug$Mechanism,
                                 col = list(
                                   Mechanism = c('Chemotherapy'=ColJournal$Science[1],
                                                 'Cell cycle inhibition'=ColJournal$Science[2],
                                                 'Immune activation'=ColJournal$Science[3],
                                                 'Signal transduction'=ColJournal$Science[4],
                                                 'EMT'=ColJournal$Science[5],
                                                 'Metabolism'=ColJournal$Science[6])
                                 ),
                                 simple_anno_size = unit(0.5, "cm"),
                                 gp = gpar(col = "white", lwd = 0.2),
                                 na_col = 'grey90')

Heatmap(op.hm.input, cluster_rows = F, cluster_columns = F, name = 'Relative sensitivity',
        top_annotation = op.hm.anno, show_column_names = F,
        left_annotation = op.hm.anno.drug,
        column_split = index.drug$MOcluster, column_gap = unit(2,'mm'),
        row_split = c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3), row_gap = unit(2,'mm'),
        col = colorRamp2(c(-2,-1.5,0,1.5,2), c(ColColor$`Low-Indigo`[9], ColColor$`Low-Indigo`[7],
                                               'white',ColColor$`High-Orange`[7], ColColor$`High-Orange`[9])),
        height = unit(0.5*15, 'cm'), width = unit(129/12, 'cm')) 
#
rm(op.train.exp, op.train.model, op.train.model.filter, op.train.response, op.train.response.mtx,
   op.dea, op.hm.input, op.hm.anno, op.hm.anno.drug, op.result)

########## Mission 9. Drug target expression ##########
##### 9.1 input prep and dea ####
dt.pro <- intersect(filter(ref.drugs, !is.na(drug_claim_name), !is.na(drug_claim_primary_name), !is.na(drug_name))$gene_name, rownames(GBC_Main_Pro))
dt.mtx <- input.raw$Pro[dt.pro, index.final$Sample]
dt.dea <- bioinfoamateur::core_Differential_analysis_continuous(CM.hm, dt.mtx, log = T, p.adj = F, method = 't.test', show.belong = T)
dt.dea$Target <- rownames(dt.dea)

##### 9.2 select feature note #####
'HK2' %in% rownames(input.raw$RNA); 'HK2' %in% rownames(input.raw$Pro)
'NAE1' %in% dt.pro
# C1(Prolif-AP-immune hot):
# C2(Signal-EMT-CAF dominant):
# C3(Metabolism-MYC-immune desert):
# C4(Glycosylation-ERBB amp-immune inhibition)
dt.input <- data_frame(# C1(Prolif-AP-immune hot): cell cycle inhibitors + ICIs
  # ICI
  `PD-1 mRNA` = as.numeric(input.raw$RNA['PDCD1',]),
  `PD-L1 mRNA` = as.numeric(input.raw$RNA['CD274',]),
  `CTLA-4 mRNA` = as.numeric(input.raw$RNA['CTLA4',]),
  `TIGIT mRNA` = as.numeric(input.raw$RNA['TIGIT',]),
  `TIM-3 mRNA` = as.numeric(input.raw$RNA['HAVCR2',]),
  `LAG3 mRNA` = as.numeric(input.raw$RNA['LAG3',]),
  # `CD47 mRNA` = as.numeric(input.raw$RNA['CD47',]),
  `CD47 Protein` = as.numeric(input.raw$Pro['CD47',]),
  # Cell cycle
  # `CDK6 mRNA` = as.numeric(input.raw$RNA['CDK6',]),
  `CDK6 Protein` = as.numeric(input.raw$Pro['CDK6',]),
  # `ATM mRNA` = as.numeric(input.raw$RNA['ATM',]),
  `ATM Protein` = as.numeric(input.raw$Pro['ATM',]),
  # `WEE1 mRNA` = as.numeric(input.raw$RNA['WEE1',]),
  `AURKA mRNA` = as.numeric(input.raw$RNA['AURKA',]),
  
  # C2(Signal-EMT-CAF dominant):
  # signal
  # `PDGFRA mRNA` = as.numeric(input.raw$RNA['PDGFRA',]),
  `PDGFRA Protein` = as.numeric(input.raw$Pro['PDGFRA',]),
  # `PDGFRB mRNA` = as.numeric(input.raw$RNA['PDGFRB',]),
  `PDGFRB Protein` = as.numeric(input.raw$Pro['PDGFRB',]),
  # `FGFR1 mRNA` = as.numeric(input.raw$RNA['FGFR1',]),
  `FGFR1 Protein` = as.numeric(input.raw$Pro['FGFR1',]),
  # `JAK2 mRNA` = as.numeric(input.raw$RNA['JAK2',]),
  `JAK2 Protein` = as.numeric(input.raw$Pro['JAK2',]),
  # `ABL1 mRNA` = as.numeric(input.raw$RNA['ABL1',]),
  `ABL1 Protein` = as.numeric(input.raw$Pro['ABL1',]),
  # `AKT3 mRNA` = as.numeric(input.raw$RNA['AKT3',]),
  `AKT3 Protein` = as.numeric(input.raw$Pro['AKT3',]),
  # Wnt EMT
  # `RAC1 mRNA` = as.numeric(input.raw$RNA['RAC1',]),
  `RAC1 Protein` = as.numeric(input.raw$Pro['RAC1',]),
  # `AXL mRNA` = as.numeric(input.raw$RNA['AXL',]),
  `AXL Protein` = as.numeric(input.raw$Pro['AXL',]),
  
  # C3(Metabolism-MYC-immune desert):
  # MYC
  # `BRD4 mRNA` = as.numeric(input.raw$RNA['BRD4',]),
  `BRD4 Protein` = as.numeric(input.raw$Pro['BRD4',]),
  # `HSP90AA1 mRNA` = as.numeric(input.raw$RNA['HSP90AA1',]),
  `HSP90AA1 Protein` = as.numeric(input.raw$Pro['HSP90AA1',]),
  # FFA and AA Metabolism
  # `CPT2 mRNA` = as.numeric(input.raw$RNA['CPT2',]),
  `CPT2 Protein` = as.numeric(input.raw$Pro['CPT2',]),
  # `ACAT1 mRNA` = as.numeric(input.raw$RNA['ACAT1',]),
  `ACAT1 Protein` = as.numeric(input.raw$Pro['ACAT1',]),
  # `NDUFS2 mRNA` = as.numeric(input.raw$RNA['NDUFS2',]),
  `NDUFS2 Protein` = as.numeric(input.raw$Pro['NDUFS2',]),
  # `HSD17B10 mRNA` = as.numeric(input.raw$RNA['HSD17B10',]),
  `HSD17B10 Protein` = as.numeric(input.raw$Pro['HSD17B10',]),
  # `HSPA9 mRNA` = as.numeric(input.raw$RNA['HSPA9',]),
  `HSPA9 Protein` = as.numeric(input.raw$Pro['HSPA9',]),
  # `BCAT2 mRNA` = as.numeric(input.raw$RNA['BCAT2',]),
  `BCAT2 Protein` = as.numeric(input.raw$Pro['BCAT2',]),
  # `BCKDHA mRNA` = as.numeric(input.raw$RNA['BCKDHA',]),
  `BCKDHA Protein` = as.numeric(input.raw$Pro['BCKDHA',]),
  # `ACADSB mRNA` = as.numeric(input.raw$RNA['ACADSB',]),
  `ACADSB Protein` = as.numeric(input.raw$Pro['ACADSB',]),
  # `AMACR mRNA` = as.numeric(input.raw$RNA['AMACR',]),
  `AMACR Protein` = as.numeric(input.raw$Pro['AMACR',]),
  
  # C4(Glycosylation-ERBB amp-immune inhibition):
  # glycolysis
  # `HK2 mRNA` = as.numeric(input.raw$RNA['HK2',]),
  `HK2 Protein` = as.numeric(input.raw$Pro['HK2',]),
  # `PFKFB3 mRNA` = as.numeric(input.raw$RNA['PFKFB3',]),
  `PFKFB3 Protein` = as.numeric(input.raw$Pro['PFKFB3',]),
  # glycosylation
  # `MGAT2 mRNA` = as.numeric(input.raw$RNA['MGAT2',]),
  `MGAT2 Protein` = as.numeric(input.raw$Pro['MGAT2',]),
  # `GNE mRNA` = as.numeric(input.raw$RNA['GNE',]),
  `GNE Protein` = as.numeric(input.raw$Pro['GNE',]),
  # `FUT3 mRNA` = as.numeric(input.raw$RNA['FUT3',]),
  `FUT3 Protein` = as.numeric(input.raw$Pro['FUT3',]),
  # Other AA Metabolism and Redox
  # `SOD2 mRNA` = as.numeric(input.raw$RNA['SOD2',]),
  `SOD2 Protein` = as.numeric(input.raw$Pro['SOD2',]),
  # `ASL mRNA` = as.numeric(input.raw$RNA['ASL',]),
  `ASL Protein` = as.numeric(input.raw$Pro['ASL',]),
  # ERBB2
  # `ERBB2 mRNA` = as.numeric(input.raw$RNA['ERBB2',]),
  `ERBB2 Protein` = as.numeric(input.raw$Pro['ERBB2',]))
rownames(dt.input) <- index.final$Sample
dt.input <- as.data.frame(t(dt.input)) %>% t() %>% scale() %>% t() %>% as.data.frame()
##### Vis
temp.hm.anno <- columnAnnotation(`MO cluster` = index.final$MOcluster,
                                 col = list(
                                   `MO cluster` = col.mo
                                 ),
                                 simple_anno_size = unit(0.5, "cm"),
                                 gp = gpar(col = "white", lwd = 0.2),
                                 na_col = 'grey90')
# vis
Heatmap(dt.input, cluster_rows = F, cluster_columns = F, name = 'Scaled abundance',
        top_annotation = temp.hm.anno, show_column_names = F,
        column_split = index.final$MOcluster, column_gap = unit(2,'mm'),
        row_split = c(rep('MO1',10),rep('MO2',8),
                      rep('MO3',11),rep('MO4',8)), row_gap = unit(2,'mm'),
        col = colorRamp2(c(-1.5,0,1.5), c(ColColor$`Low-Blue`[9],'white',ColColor$`High-Red`[9])),
        height = unit(0.5*30, 'cm'), width = unit(129/10, 'cm')) 
#
rm(input.emt, temp.pathway, temp.cor, temp.cor.cor)



########## 10. Sankey plot ##########
##### Input prep
sk.input <- tibble(MOcluster = index.final$MOcluster,
                   Immcluster = index.final$ImmCluster_Name)
#
sk.input <- tibble(Circ = names(table(paste(sk.input$MOcluster, sk.input$Immcluster))),
                   Freq = table(paste(sk.input$MOcluster, sk.input$Immcluster)))
sk.input$MOcluster = str_split_fixed(sk.input$Circ, ' ', 2)[,1] %>%
  factor(levels = c('MO1','MO2','MO3','MO4'), ordered = T)
sk.input$Immcluster = str_split_fixed(sk.input$Circ, ' ', 2)[,2] %>%
  factor(levels = c('A','B','C','D'), ordered = T)
sk.input$Treatment <- case_when(sk.input$MOcluster == 'MO1' ~ 'AAAAAAAAA',
                                sk.input$MOcluster == 'MO2' ~ 'BBBBBBBBB',
                                sk.input$MOcluster == 'MO3' ~ 'CXXXXXXXXX',
                                sk.input$MOcluster == 'MO4' ~ 'DDDDDDDDDD')
##### sankey plot
library(ggalluvial)
library(ggthemes)
ggplot(sk.input,
       aes(y = Freq, axis1 = Treatment, axis2 = MOcluster, axis3 = Immcluster)) +
  geom_alluvium(aes(fill = MOcluster), width = 1/12) +
  geom_stratum(width = 1/12, fill = "white", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("MOcluster", "Immcluster"), expand = c(.05, .05)) +
  # scale_fill_brewer(type = "qual", palette = "Set1") +
  scale_fill_manual(values = col.mo) + 
  ggtitle("Sample distribution") +
  theme_classic2() + theme(axis.text.x = element_text(color = 'black'),
                           axis.text.y = element_text(color = 'black')) +
  coord_flip() # 8*12 Sankey


##### ex. any feature obv #####
index.obv <- select(index.final,1,2,3,4,15,29)
##### ICIs 
# CD40
index.obv$CD40_Pro <- sapply(index.obv$Sample, function(x){as.numeric(input.raw$Pro['CD40',x])}) %>% as.numeric()
ggboxplot(index.obv, 'MOcluster', 'CD40_Pro') + stat_compare_means()
# CD276
index.obv$CD276_Pro <- sapply(index.obv$Sample, function(x){as.numeric(input.raw$Pro['CD276',x])}) %>% as.numeric()
ggboxplot(index.obv, 'MOcluster', 'CD276_Pro') + stat_compare_means()
# CEACAM1
# CD274 mRNA
index.obv$CD274_mRNA <- sapply(index.obv$Sample, function(x){as.numeric(input.raw$RNA['CD274',x])}) %>% as.numeric()
ggboxplot(index.obv, 'MOcluster', 'CD274_mRNA') + stat_compare_means()

##### MYC [phosphosite all passed, mRNA not passed]
# MYC S62
index.obv$MYC_S62 <- sapply(index.obv$Pro, function(x){as.numeric(GBC_Main_Phos_NA['MYC:S62',x])}) %>% as.numeric()
ggboxplot(index.obv, 'MOcluster', 'MYC_S62') + stat_compare_means()
# MYC S161
index.obv$MYC_S161 <- sapply(index.obv$Pro, function(x){as.numeric(GBC_Main_Phos_NA['MYC:S161',x])}) %>% as.numeric()
ggboxplot(index.obv, 'MOcluster', 'MYC_S161') + stat_compare_means()
table(!is.na(index.obv$MYC_S161), index.obv$MOcluster)
# MYC S347
index.obv$MYC_S347 <- sapply(index.obv$Pro, function(x){as.numeric(GBC_Main_Phos_NA['MYC:S347',x])}) %>% as.numeric()
ggboxplot(index.obv, 'MOcluster', 'MYC_S347') + stat_compare_means()
table(!is.na(index.obv$MYC_S347), index.obv$MOcluster)
# MYC S348
index.obv$MYC_S348 <- sapply(index.obv$Pro, function(x){as.numeric(GBC_Main_Phos_NA['MYC:S348',x])}) %>% as.numeric()
ggboxplot(index.obv, 'MOcluster', 'MYC_S348') + stat_compare_means()
table(!is.na(index.obv$MYC_S348), index.obv$MOcluster)
# MYC mRNA
index.obv$MYC_mRNA <- sapply(index.obv$Sample, function(x){as.numeric(input.raw$RNA['MYC',x])}) %>% as.numeric()
ggboxplot(index.obv, 'MOcluster', 'MYC_mRNA') + stat_compare_means()

##### Clinical metadata #####
##### 1.1 clinical 
Clinical.metadata <- tibble(Patient = GBC_Main_Clinical$Patient_ID, 
                            Gender = GBC_Main_Clinical$Sex,
                            Age = GBC_Main_Clinical$Age,
                            `Tumor location` = case_when(GBC_Main_Clinical$Tumor_location == '1' ~ 'fundus',
                                                         GBC_Main_Clinical$Tumor_location == '2' ~ 'body',
                                                         GBC_Main_Clinical$Tumor_location == '3' ~ 'neck',
                                                         GBC_Main_Clinical$Tumor_location == '4' ~ 'cystic duct',
                                                         GBC_Main_Clinical$Tumor_location == '5' ~ 'more than 4 regions',
                                                         GBC_Main_Clinical$Tumor_location == '1,2,3' ~ 'fundus + body + neck',
                                                         GBC_Main_Clinical$Tumor_location == '1,2' ~ 'fundus + body',
                                                         GBC_Main_Clinical$Tumor_location == '2,3' ~ 'body + neck'),
                            `Histologic type` = GBC_Main_Clinical$Pathological_type,
                            `Number of tumor` = case_when(grepl(fixed(';'), GBC_Main_Clinical$Tumor_size_diameter) ~ 2,
                                                          !is.na(GBC_Main_Clinical$Tumor_size_diameter) ~ 1,
                                                          is.na(GBC_Main_Clinical$Tumor_size_diameter) ~ NA),
                            `Tumor size (cm)` = str_split_fixed(GBC_Main_Clinical$Tumor_size_diameter, 'cm', 2)[,1] %>% as.numeric(),
                            `Hepatic invasion` = case_when(GBC_Main_Clinical$Regional_invasion == 1 ~ 'yes',
                                                           GBC_Main_Clinical$Regional_invasion == 0 ~ 'no'),
                            `Lymph node metastasis` = case_when(GBC_Main_Clinical$Regional_lymph_node_metastasis == 1 ~ 'yes',
                                                                GBC_Main_Clinical$Regional_lymph_node_metastasis == 0 ~ 'no'),
                            `Distal metastasis` = case_when(GBC_Main_Clinical$Distal_metastasis == 1 ~ 'yes',
                                                            GBC_Main_Clinical$Distal_metastasis == 0 ~ 'no'),
                            `TNM` = GBC_Main_Clinical$TNM,
                            `TNM stage` = GBC_Main_Clinical$Tumor_stage,
                            `Vascular invasion` = case_when(GBC_Main_Clinical$Vascular_invasion == 1 ~ 'yes',
                                                            GBC_Main_Clinical$Vascular_invasion == 0 ~ 'no'),
                            `Perineural invasion` = case_when(GBC_Main_Clinical$Perineural_invasion == 1 ~ 'yes',
                                                             GBC_Main_Clinical$Perineural_invasion == 0 ~ 'no'),
                            `Cholelithiasis` = case_when(GBC_Main_Clinical$cholelithiasis == 1 ~ 'yes',
                                                         GBC_Main_Clinical$cholelithiasis == 0 ~ 'no'),
                            `CA19-9 (U/mL)` = as.numeric(GBC_Main_Clinical$CA199),
                            `Carcinoembryonic antigen (ng/mL)` = as.numeric(GBC_Main_Clinical$CEA),
                            `Adjuvant therapy` = GBC_Main_Clinical$Adjuvant.therapies,
                            `Survival (month)` = as.numeric(GBC_Main_Clinical$OS_days)/30,
                            `Status` = case_when(GBC_Main_Clinical$OS_event == 1 ~ 'dead',
                                                 GBC_Main_Clinical$OS_event == 0 ~ 'alive'),
                            `Tumor WES ID` = GBC_Main_Index$WES_T, `NAT WES ID` = GBC_Main_Index$WES_P,
                            `Tumor RNA-seq ID` = GBC_Main_Index$RNA_T, `NAT RNA-seq ID` = GBC_Main_Index$RNA_P,
                            `Tumor Proteomic&Phosphoproteomic ID` = GBC_Main_Index$ProPhos_T,
                            `NAT Proteomic&Phosphoproteomic ID` = GBC_Main_Index$ProPhos_P)

# Data for Table S1
write.xlsx(Clinical.metadata, file = '../Clinical.xlsx', rowNames = T, colNames = T, overwrite = T)

##### 1.2 mut
Omic.mut <- as.tibble(select(GBC_Main_Mutation, 
                             Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position, Variant_Classification, Variant_Type,
                             Reference_Allele, Tumor_Seq_Allele2, Hugo_Symbol, txChange, aaChange))
Omic.mut$`Methodology` = case_when(GBC_Main_Mutation$Type == 'Tonly' ~ 'Tumor only mode',
                                   GBC_Main_Mutation$Type == 'TN' ~ 'Tumor-NAT mode')
# Data for Table S1
write.xlsx(Omic.mut, file = '../Omics_mut.xlsx', rowNames = T, colNames = T, overwrite = T)

##### 1.3 cnv
# Data for Table S1
write.xlsx(GBC_Main_CNV, file = '../Omics_cnv.xlsx', rowNames = T, colNames = T, overwrite = T)

##### 1.3 rna
# Data for Table S1
write.xlsx(GBC_Main_RNA_logTPM, file = '../Omics_rna.xlsx', rowNames = T, colNames = T, overwrite = T)

##### 1.4 pro
# Data for Table S1
write.xlsx(GBC_Main_Pro, file = '../Omics_pro.xlsx', rowNames = T, colNames = T, overwrite = T)

##### 1.4 phos
# Data for Table S1
write.xlsx(GBC_Main_Phos_knn, file = '../Omics_phos.xlsx', rowNames = T, colNames = T, overwrite = T)

##### 1.5 Data table for CPTAC
# data for CPTAC table


