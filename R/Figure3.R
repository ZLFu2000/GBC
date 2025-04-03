##### Figure 3 and Figure S3
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
# TCGA Pan-cancer (XENA)
load('../TCGA_Pan-cancer.RData')
rm(Pan_CNV, Pan_CNV_cut, Pan_Mutation, Pan_survival)
Pan_Clinical <- fread('../TCGA pan metadata')
# select NE: PCPG, THCA, SKCA, GBC, ACC
table(Pan_Clinical$`cancer type abbreviation`)
Pan_Clinical <- filter(Pan_Clinical, `cancer type abbreviation` %in% c("ACC",'PCPG','SKCM','THCA'))
Pan_Clinical$`cancer type abbreviation` <- factor(Pan_Clinical$`cancer type abbreviation`, c("ACC",'PCPG','SKCM','THCA'), ordered = T)
# arrange
temp.mutual.sample <- intersect(Pan_Clinical$sample, colnames(Pan_RNA))
Pan_Clinical <- filter(Pan_Clinical, sample %in% temp.mutual.sample) %>% arrange(`cancer type abbreviation`)
# sample selection
Pan_RNA <- Pan_RNA[, Pan_Clinical$sample]

#### 2 color
load('../Colors (ggsci).RData')
library(wesanderson)
# without TAC
ColPatho = c(wes_palette('Darjeeling1', 3, type = c("discrete")))
names(ColPatho) = c('AC','AS','NE')
# with TAC
ColPatho2 = c(wes_palette('Darjeeling1', 4, type = c("discrete")), ColJournal$Nature[3])
names(ColPatho2) = c('AC','AS','NE','Adjacent')
# only NE and AC/AS
ColPatho3 = c('brown',wes_palette('Darjeeling1', 3, type = c("discrete"))[3])
names(ColPatho3) = c('AC_AS','NE')

### 3 Genesets
load('../Hallmark_KEGG_Reactome.RData')
# 4 COSMIC Cancer Gene Census
ref.cosmic <- fread('../COSMIC genes.csv')
# 5 Gene location
load('/Users/fuzile/Desktop/たいようけい/基本文件/常用特殊数据集/ChrBandLength.RData')
#
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

##### prep2. index #####
Index <- tibble(Patient = rep(GBC_Main_Index$Sample, 2),
                Sample = c(paste0(GBC_Main_Index$Sample,'_T'), paste0(GBC_Main_Index$Sample, '_P')),
                SampleWES = c(GBC_Main_Index$WES_T, GBC_Main_Index$WES_P),
                SampleRNA = c(GBC_Main_Index$RNA_T, GBC_Main_Index$RNA_P),
                SamplePro = c(GBC_Main_Index$ProPhos_T, GBC_Main_Index$ProPhos_P),
                SamplePhos = c(GBC_Main_Index$ProPhos_T, GBC_Main_Index$ProPhos_P),
                Type = c(rep('Cancer', 195), rep('Adjacent', 195)) %>%
                  factor(levels = c('Cancer', 'Adjacent'), ordered = T))
# add patho and retain only adeno, adenos and NE
Index$Pathology = rep(GBC_Main_Clinical$Pathological_type, 2)
Index = filter(Index, Pathology %in% c('adenocarcinoma','adenosquanmous carcinoma','neuroendocrine carcinoma')) %>%
  mutate(Pathology = case_when(Pathology == 'adenocarcinoma' ~ 'AC',
                               Pathology == 'adenosquanmous carcinoma' ~ 'AS',
                               Pathology == 'neuroendocrine carcinoma' ~ 'NE') %>%
           factor(levels = c('AC','AS','NE'), ordered = T)) 
# add Ident
Index$Ident = ifelse(as.character(Index$Type == 'Cancer'), as.character(Index$Pathology), as.character(Index$Type)) %>%
  factor(levels = c('AC','AS','NE','Adjacent'), ordered = T)

# Have sample retained?
Index$SampleRetained = ifelse(Index$Patient %in% c('GBC_001','GBC_004','GBC_005','GBC_008','GBC_014',
                                                   'GBC_021','GBC_028', 'GBC_031','GBC_032','GBC_035',
                                                   'GBC_052','GBC_056','GBC_067','GBC_081','GBC_117','GBC_122',
                                                   'GBC_154','GBC_185','GBC_191','GBC_193','GBC_199','GBC_202','GBC_206'),
                              'Yes','No')

##### Calc CIN
# add CIN
##### Calc CIN prep
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
Index$CIN = sapply(Index$SampleWES, function(x){as.numeric(filter(Output, Sample == x)$CIN)}) %>% as.numeric()
# Input$CINGroup = ifelse(Input$CIN >= median(Input$CIN, na.rm = T), 'CIN High','CIN Low') %>%
#   factor(levels = c('CIN High','CIN Low'), ordered = T)
#
rm(InputCIN, InputCalcIndex, Output)

##### Calc TMB
TempTMBCalc <- maftools::read.maf(maf = GBC_Main_Mutation)
TempTMB <- tmb(TempTMBCalc)
Index$TMB <- sapply(Index$Patient, function(x){filter(TempTMB, Tumor_Sample_Barcode == x)$total_perMB %>% as.numeric()}) %>% as.numeric()
# rm
rm(TempTMB, TempTMBCalc, Temp)

# Index tumor
Indext <- filter(Index, Type == 'Cancer')

#
Index$Pathology2 = ifelse(Index$Pathology == 'NE', 'NE', 'AC or AS')
Indext$Pathology2 = ifelse(Indext$Pathology == 'NE', 'NE', 'AC or AS')

##### Figure S3A #####
#
Input = filter(Indext, !is.na(SamplePro), !is.na(Pathology))
# adeno：ERBB2
Input$ERBB2_Protein = sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['ERBB2',x]), NA)})
# sqa：KRT6A
Input$KRT6A_RNA = sapply(Input$SampleRNA, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['KRT6A',x]), NA)})
# NE：CHGB,NCAM1
Input$NCAM1_Protein = sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['NCAM1',x]), NA)})
Input$CHGB_RNA = sapply(Input$SampleRNA, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['CHGB',x]), NA)})
# vis
plot_grid(
  # AC:ERBB2
  ggboxplot(Input, 'Pathology', 'ERBB2_Protein', color = 'Pathology', lwd = 1) +
    scale_color_manual(values = ColPatho) + theme(legend.position = 'none') +
    stat_pvalue_manual(compare_means(ERBB2_Protein ~ Pathology, data = Input, method = "t.test") %>% 
                         mutate(y.position = c(5,6,7)), label = 'p.adj'),
  # AS:KRT6A
  ggboxplot(Input, 'Pathology', 'KRT6A_RNA', color = 'Pathology', lwd = 1) +
    scale_color_manual(values = ColPatho) + theme(legend.position = 'none') +
    stat_pvalue_manual(compare_means(KRT6A_RNA ~ Pathology, data = Input, method = "t.test") %>% 
                         mutate(y.position = c(15,17,19)), label = 'p.adj'),
  # NE:NCAM1, CHGB
  ggboxplot(Input, 'Pathology', 'NCAM1_Protein', color = 'Pathology', lwd = 1) +
    scale_color_manual(values = ColPatho) + theme(legend.position = 'none') +
    stat_pvalue_manual(compare_means(NCAM1_Protein ~ Pathology, data = Input, method = "t.test") %>% 
                         mutate(y.position = c(5,6,7)), label = 'p.adj'),
  ggboxplot(Input, 'Pathology', 'CHGB_RNA', color = 'Pathology', lwd = 1) +
    scale_color_manual(values = ColPatho) + theme(legend.position = 'none') +
    stat_pvalue_manual(compare_means(CHGB_RNA ~ Pathology, data = Input, method = "t.test") %>% 
                         mutate(y.position = c(10,11,12)), label = 'p.adj'),
  ncol = 2) 

##### Figure 3A #####
library(factoextra); library(FactoMineR)
library(r.jive)
##### Input prep
Input <- Indext %>% filter(!is.na(SampleWES), !is.na(Pathology)) %>% select(Patient, SampleRNA, SamplePro, SamplePhos, Pathology)
InputRNA <- GBC_Main_RNA_logTPM[,na.omit(Input$SampleRNA)]
InputPro <- GBC_Main_Pro[,na.omit(Input$SamplePro)]
InputPhos <- GBC_Main_Phos_knn[,na.omit(Input$SamplePhos)]

# feature selection
InputRNA <- InputRNA[order(apply(InputRNA, 1, sd),decreasing = T),][1:quantile(1:nrow(InputRNA), 0.3), ] %>% t() %>% as.data.frame()
InputPro <- InputPro[order(apply(InputPro, 1, sd),decreasing = T),][1:quantile(1:nrow(InputPro), 0.3), ] %>% t() %>% as.data.frame()
InputPhos <- InputPhos[order(apply(InputPhos, 1, sd),decreasing = T),][1:quantile(1:nrow(InputPhos), 0.3), ] %>% t() %>% as.data.frame()
InputJIVE <- cbind(InputPro, InputPhos)

##### Get contribution
PCA(as.data.frame(InputRNA), graph = T) # 29.15 9.44
PCA(as.data.frame(InputPro), graph = T) # 19.74 8,8
PCA(as.data.frame(InputPhos), graph = T)# 21.06 14.05
##### PCA
InputRNA <- PCA(as.data.frame(InputRNA), graph = F)
InputPro <- PCA(as.data.frame(InputPro), graph = F)
InputPhos <- PCA(as.data.frame(InputPhos), graph = F)
InputJIVE <- PCA(as.data.frame(InputJIVE), graph = F)
##### obtain coord
InputCoordRNA <- InputRNA$ind$coord %>% as.data.frame() %>% dplyr::select(1,2)
InputCoordRNA$Pathology <- sapply(rownames(InputCoordRNA), function(x){filter(Indext, SampleRNA == x)$Pathology}) %>% as.character()
InputCoordPro <- InputPro$ind$coord %>% as.data.frame() %>% dplyr::select(1,2)
InputCoordPro$Pathology <- sapply(rownames(InputCoordPro), function(x){filter(Indext, SamplePro == x)$Pathology}) %>% as.character()
InputCoordPhos <- InputPhos$ind$coord %>% as.data.frame() %>% dplyr::select(1,2)
InputCoordPhos$Pathology <- sapply(rownames(InputCoordPhos), function(x){filter(Indext, SamplePhos == x)$Pathology}) %>% as.character()
InputCoordJIVE <- InputJIVE$ind$coord %>% as.data.frame() %>% dplyr::select(1,2)
InputCoordJIVE$Pathology <- sapply(rownames(InputCoordJIVE), function(x){filter(Indext, SamplePhos == x)$Pathology}) %>% as.character()

##### vis
plot_grid(ggscatter(InputCoordRNA, x = 'Dim.1', y = 'Dim.2', color = 'Pathology', size = 4, 
                    ellipse = TRUE, ellipse.alpha = 0.2, ellipse.border.remove = T, ellipse.level = 0.95,
                    palette = ColPatho) +
            theme(legend.position = 'none') +
            scale_color_manual(values = ColPatho) +
            labs(x = 'PC1', y = 'PC2', title = 'mRNA'),
          ggscatter(InputCoordPro, x = 'Dim.1', y = 'Dim.2', color = 'Pathology', size = 4, 
                    ellipse = TRUE, ellipse.alpha = 0.2, ellipse.border.remove = T, ellipse.level = 0.95,
                    palette = ColPatho) +
            theme(legend.position = 'none') +
            scale_color_manual(values = ColPatho) +
            labs(x = 'PC1', y = 'PC2', title = 'Protein'),
          ggscatter(InputCoordPhos, x = 'Dim.1', y = 'Dim.2', color = 'Pathology', size = 4, 
                    ellipse = TRUE, ellipse.alpha = 0.2, ellipse.border.remove = T, ellipse.level = 0.95,
                    palette = ColPatho) +
            theme(legend.position = 'none') +
            scale_color_manual(values = ColPatho) +
            labs(x = 'PC1', y = 'PC2', title = 'Phosphosite'), 
          ggscatter(InputCoordJIVE, x = 'Dim.1', y = 'Dim.2', color = 'Pathology', size = 4, 
                    ellipse = TRUE, ellipse.alpha = 0.2, ellipse.border.remove = T, ellipse.level = 0.95,
                    palette = ColPatho) + 
            theme(legend.position = 'none') +
            scale_color_manual(values = ColPatho) +
            labs(x = 'JIVE 1', y = 'JIVE 2', title = 'JIVE'), nrow = 2) 
#
rm(Input, InputCoordJIVE,InputCoordPhos,InputCoordPro,InputCoordRNA,InputJIVE,InputRNA,InputPro,InputPhos,
   OutputCluster, Output, OutputPCs, MutualSample)

##### Figure S3B & S3C#####
##### 1.1 NE score genesets prep
# 10 genes
Ref.10 <- c("SCG3", "CHGA", "CHGB", "CHRNB2", "PCSK1", "ELAVL4", "ENO2", "SCN3A", "SYP", "NKX2-1")
# 50 genes
Ref.50 <- c("BEX1", "ASCL1", "INSM1", "CHGA", "TAGLN3", "KIF5C", "CRMP1", "SCG3", "SYT4", "RTN1", "MYT1",
            "SYP", "KIF1A", "TMSB15A", "SYN1", "SYT11", "RUNDC3A", "TFF3", "CHGB", "FAM57B", "SH3GL2",
            "BSN", "SEZ6", "TMSB15B", "CELF3")
# 70 genes
Ref.70 <- c("ASXL3", "CAND2", "ETV5", "GPX2", "JAKMIP2", "KIAA0408", "SOGA3", "TRIM9", "BRINP1", "C7orf76",
            "GNAO1", "KCNB2", "KCND2", "LRRC16B", "MAP10", "NRSN1", "PCSK1", "PROX1", "RGS7", "SCG3",
            "SEC11C", "SEZ6", "ST8SIA3", "SVOP", "SYT11", "AURKA", "DNMT1", "EZH2", "MYCN")

##### 1.2 index
SeparatedCNV <- function(vec) {
  vec = case_when(vec == -2 ~ 'Deletion',
                  vec == -1 ~ 'ShallowDeletion',
                  vec == 0 ~ 'Diploid', 
                  vec == 1 ~ 'Gain',
                  vec == 2 ~ 'Amplification') %>%
    factor(levels = c('Amplification','Gain', 'Diploid', 'ShallowDeletion', 'Deletion'), ordered = T)
  return(vec)
}
IndexPatho <- rbind(select(GBC_Main_Index, 1,2,4,6) %>% 
                      bioinfoamateur::dfm_change_colnames(c('Sample','WES','RNA','Pro')), 
                    mutate(GBC_Main_Index, Sample = paste0(Sample, '_P')) %>% select(1,3,5,7) %>% 
                      bioinfoamateur::dfm_change_colnames(c('Sample','WES','RNA','Pro')))
IndexPatho$Tumor <- ifelse(!grepl(fixed('_P'),IndexPatho$Sample), 'Yes', NA)
# retain samples having proteomic data
IndexPatho <- filter(IndexPatho, !is.na(Pro))

##### add ERBB3 mutation
IndexPatho$Mut_ERBB3 <- case_when(IndexPatho$Sample %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB3')$Tumor_Sample_Barcode ~ 'Yes',
                                  IndexPatho$Sample %in% GBC_Main_Mutation$Tumor_Sample_Barcode ~ 'No',
                                  !is.na(IndexPatho$Sample) ~ NA) %>%
  factor(levels = c('Yes','No'), ordered = T)

##### add ERBB2 CNV
# ERBB2
IndexPatho$CNV_ERBB2 <- sapply(IndexPatho$WES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV_cate['ERBB2',x]), NA)})
IndexPatho$CNV_ERBB2 <- SeparatedCNV(IndexPatho$CNV_ERBB2)

##### add cli paras
# Age
IndexPatho$Age = sapply(IndexPatho$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Age}) %>% as.numeric()
# Gender
IndexPatho$Gender = sapply(IndexPatho$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Sex}) %>% as.character() %>%
  factor(levels = c('Male','Female'), ordered = T)
# Stage
IndexPatho$Stage = sapply(IndexPatho$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Tumor_stage2}) %>% as.numeric()
IndexPatho$Stage = case_when(IndexPatho$Stage == 1 ~ 'I', IndexPatho$Stage == 2 ~ 'II',
                             IndexPatho$Stage == 3 ~ 'III', IndexPatho$Stage == 4 ~ 'IV') %>%
  factor(levels = c('I','II','III','IV'), ordered = T)
# Patho
IndexPatho$Pathology = sapply(IndexPatho$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Pathological_type}) %>% as.character()
IndexPatho$Pathology = case_when(IndexPatho$Pathology == 'adenocarcinoma' ~ 'AC',
                                 IndexPatho$Pathology == 'adenosquanmous carcinoma' ~ 'AS',
                                 IndexPatho$Pathology == 'neuroendocrine carcinoma' ~ 'NE',
                                 !is.na(IndexPatho$Pathology) ~ 'Others') %>%
  factor(levels = c('AC','AS','NE','Others'), ordered = T)
# LM
IndexPatho$LM <- sapply(IndexPatho$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Regional_invasion}) %>% as.numeric()
IndexPatho$LM = ifelse(IndexPatho$LM == 1, 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)

##### add MYC mRNA and Protein levels and pro-level phos
# MYC
IndexPatho$mRNA_MYC <- sapply(IndexPatho$RNA, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['MYC',x]), NA)}) %>% as.numeric()
IndexPatho$Phos_MYC_S62 <- sapply(IndexPatho$Pro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Phos_NA['MYC:S62',x]), NA)}) %>% as.numeric()
IndexPatho$ProPhos_MYC <- sapply(IndexPatho$Pro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_ProLevelPhos_knn['MYC',x]), NA)}) %>% as.numeric()
# MEIS1
IndexPatho$mRNA_MEIS1 <- sapply(IndexPatho$RNA, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['MEIS1',x]), NA)}) %>% as.numeric()
IndexPatho$Protein_MEIS1 <- sapply(IndexPatho$Pro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['MEIS1',x]), NA)}) %>% as.numeric()

##### 1.3 data arrangement
GBC.rna <- GBC_Main_RNA_logTPM[, na.omit(GBC_Main_Index$RNA_T)]
### index first
# 4 class
Index.gbc <- tibble(Sample = colnames(GBC.rna),
                    Pathology = sapply(colnames(GBC.rna), function(x){filter(IndexPatho, RNA == x, Tumor == "Yes")$Pathology}))
Index.gbc$Pathology = as.numeric(Index.gbc$Pathology)
Index.gbc$Pathology <- case_when(Index.gbc$Pathology == 1 ~ 'AC',
                                 Index.gbc$Pathology == 2 ~ 'AS',
                                 Index.gbc$Pathology == 3 ~ 'NE',
                                 is.na(Index.gbc$Pathology) ~ 'Others') %>%
  factor(levels = c('AC','AS','NE','Others'), ordered = T)
table(Index.gbc$Pathology)
# 2 class
Index.gbc$Pathology2 <- ifelse(Index.gbc$Pathology == 'NE', 'NE','Others')
Index.gbc$Pathology2 <- factor(Index.gbc$Pathology2, levels = c('NE','Others'), ordered = T)
#
Index.gbc <- na.omit(Index.gbc)

##### mtx first
mtx.gbc <- GBC_Main_RNA_logTPM[, Index.gbc$Sample]

### uniformed index and mtx making
# mutual genes
temp.mutual.gene <- intersect(rownames(mtx.gbc), rownames(Pan_RNA))
mtx.gbc <- mtx.gbc[temp.mutual.gene, ]
mtx.pan <- Pan_RNA[temp.mutual.gene, ]
# index
index.pan <- tibble(Sample = Pan_Clinical$sample,
                    Pathology = Pan_Clinical$`cancer type abbreviation`,
                    Pathology2 = Pan_Clinical$`cancer type abbreviation`,)
identical(index.pan$Sample, colnames(mtx.pan))
##### merge
mtx.all <- cbind(mtx.gbc, mtx.pan)
index.all <- rbind(Index.gbc, index.pan)
index.all$Cohort <- c(rep('Fu-GBC', 130), rep('TCGA', 1311)) %>%
  factor(levels = c('Fu-GBC','TCGA'), ordered = T)

##### 1.4 remove batch
library(sva)
# prep
mtx.all.t <- t(mtx.all)
identical(rownames(mtx.all.t), index.all$Sample)
# sva
mtx.all.db <- ComBat(dat = mtx.all, batch = index.all$Cohort, 
                     par.prior = TRUE, prior.plots = FALSE)

##### 1.5 score calc
Ref.ne <- data.frame(NE10 = c(Ref.10, rep(NA, 19)),
                     NE50 = c(Ref.50, rep(NA, 4)),
                     NE70 = c(Ref.70))
# ssgsea 
ssgsea.all <- GSVA::gsva(as.matrix(mtx.all), Ref.ne, method = 'ssgsea')
# add to index
index.all$NE10 <- as.numeric(ssgsea.all['NE10', index.all$Sample])
index.all$NE50 <- as.numeric(ssgsea.all['NE50', index.all$Sample])
index.all$NE70 <- as.numeric(ssgsea.all['NE70', index.all$Sample])

##### 1.6 vis in GBC
vis.gbc <- filter(index.all, Cohort == 'Fu-GBC', Pathology != 'Others')
# Vis
ggboxplot(vis.gbc, x = "Pathology", y = "NE10", color = "Pathology", add = "jitter") +
  scale_color_manual(values = ColPatho2) +
  stat_compare_means(method = "kruskal.test", 
                     label.y = max(vis.gbc$NE10) + 0.2,
                     aes(label = paste0("Kruskal-Wallis p = ", ..p.format..))) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("NE", "AC"), c("NE", "AS"), c("AC", "AS")), 
                     aes(label = paste("p =", sprintf("%.2f", ..p..)))) +
  theme_minimal(base_size = 15) + 
  theme(
    text = element_text(family = "Arial", color = "black"),  
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid = element_blank(), 
    panel.background = element_rect(fill = "white", color = "white"),  
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"), 
    legend.position = "none" 
  ) +
  labs(title = "10-gene NE score", x = NULL, y = "NE Score") # 5*5 NEscore-GBC

##### 1.7 vis in TCGA pan-cancer
vis.pan <- index.all
vis.pan$Pathology2 = factor(vis.pan$Pathology2, levels = c('NE','Others','ACC'))
# Vis
ggboxplot(vis.pan, x = "Pathology2", y = "NE10") +
  stat_compare_means(method = "kruskal.test", 
                     label.y = max(vis.pan$NE10) + 0.2, 
                     aes(label = paste0("Kruskal-Wallis p = ", ..p.format..))) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("NE", "Others"), c("NE", "ACC"), c("NE", "PCPG"),
                                        c('NE','SKCM'),c('NE','THCA')), 
                     aes(label = paste("p =", sprintf("%.2f", ..p..)))) +
  theme_minimal(base_size = 15) + 
  theme(
    text = element_text(family = "Arial", color = "black"),  
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = "white"),  
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"), 
    legend.position = "none" 
  ) +
  labs(title = "10-gene NE score", x = NULL, y = "NE Score") 
# data for Table S4
write.xlsx(vis.pan, file = 'NEscore.xlsx', rowNames = T, colNames = T, overwrite = T)
#
rm(vis.gbc, vis.pan, mtx.all, mtx.gbc, mtx.pan, mtx.all.db, mtx.all.t, 
   ssgsea.all, Index.gbc, index.pan, index.all, GBC.rna, temp.mutual.gene, temp.mutual.sample)

##### Figure S3D #####
##### 1.1 info loading
RefNEfactor <- read.xlsx('./NE receptor-ligand.xlsx', sheet = 2, sep.names = '_') %>% as_tibble()

##### 1.2 selecting
### prep
InputRNA <- GBC_Main_RNA_logTPM[intersect(RefNEfactor$Gene_Symbol, rownames(GBC_Main_RNA_logTPM)), 
                                c(filter(Indext, !is.na(SampleRNA), Pathology2 == 'AC or AS')$SampleRNA,
                                  filter(Indext, !is.na(SampleRNA), Pathology2 == 'NE')$SampleRNA)]
CM = data.table(AC_AS = c(rep(1,125),rep(0,5)),
                NE = c(rep(0,125),rep(1,5)))
rownames(CM) = c(filter(Indext, !is.na(SampleRNA), Pathology2 == 'AC or AS')$SampleRNA,
                 filter(Indext, !is.na(SampleRNA), Pathology2 == 'NE')$SampleRNA)
##### DEA
DEA_NEs <- bioinfoamateur::core_Differential_analysis_continuous(CM, InputRNA, log = T, p.adj = F, show.belong = T)
# add info
DEA_NEs$Gene = rownames(DEA_NEs)
DEA_NEs$NE_Axis = sapply(DEA_NEs$Gene, function(x){filter(RefNEfactor, Gene_Symbol == x)$Axis %>% as.character()})
DEA_NEs$NE_Type = sapply(DEA_NEs$Gene, function(x){filter(RefNEfactor, Gene_Symbol == x)$Type %>% as.character()})

### Gene selection
# RNA
SelectHormone_HPG = c('GNRH1','GNRH2','RLN1','TRH') # RNA
SelectReceptor_NT = filter(DEA_NEs, NE_Axis == 'Neurotransmitters', NE_Type == 'Receptor', belong == 2)$Gene
SelectReceptor_HPG = filter(DEA_NEs, NE_Axis == 'HP-Growth', NE_Type == 'Receptor', belong == 2)$Gene
# Protein
SelectReceptor_Protein = filter(RefNEfactor, Gene_Symbol %in% rownames(GBC_Main_Pro))$Gene_Symbol %>% unique()

##### 1.3 vis prep
# RNA
InputRNAVis <- t(InputRNA) %>% as.data.frame()
InputRNAVis$Pathology = c(rep('AC_AS', 125), rep('NE', 5))
InputRNAVis <- reshape2::melt(InputRNAVis, id.var = c('Pathology'), variable.name = 'Gene', value.name = 'Expression')
# Pro
InputPro <- GBC_Main_Pro[intersect(RefNEfactor$Gene_Symbol, rownames(GBC_Main_Pro)), 
                         c(filter(Indext, !is.na(SamplePro), Pathology2 == 'AC or AS')$SamplePro,
                           filter(Indext, !is.na(SamplePro), Pathology2 == 'NE')$SamplePro)]

InputProVis <- t(InputPro) %>% as.data.frame()
InputProVis$Pathology = c(rep('AC_AS', 179), rep('NE', 6))
InputProVis <- reshape2::melt(InputProVis, id.var = c('Pathology'), variable.name = 'Gene', value.name = 'Expression')

Input = filter(InputRNAVis, Gene %in% SelectHormone_HPG)

##### 1.4 vis
Vis1 <- ggboxplot(Input, "Gene", "Expression", color = "Pathology",
                  palette = ColPatho3) +
  stat_pvalue_manual(compare_means(
    Expression ~ Pathology, data = Input, group.by = 'Gene',
    method = "wilcox") %>% mutate(y.position = 8), label = 'p.signif', x = 'Gene') +
  ggtitle('mRNA Expression\nHormones in HPG Axis') + xlab('')
Vis1

Input = filter(InputRNAVis, Gene %in% SelectReceptor_NT)
# Vis
Vis2 <- ggboxplot(Input, "Gene", "Expression", color = "Pathology",
                  palette = ColPatho3) +
  stat_pvalue_manual(compare_means(
    Expression ~ Pathology, data = Input, group.by = 'Gene',
    method = "wilcox") %>% mutate(y.position = 8), label = 'p.signif', x = 'Gene') +
  ggtitle('mRNA Expression\nNeurotransmitter receptors') + xlab('')
Vis2

#
Input = filter(InputRNAVis, Gene %in% SelectReceptor_HPG)
# Vis
Vis3 <- ggboxplot(Input, "Gene", "Expression", color = "Pathology",
                  palette = ColPatho3) +
  stat_pvalue_manual(compare_means(
    Expression ~ Pathology, data = Input, group.by = 'Gene',
    method = "wilcox") %>% mutate(y.position = 8), label = 'p.signif', x = 'Gene') +
  ggtitle('mRNA Expression\nReceptors in HPG Axis') + xlab('')
Vis3
#
Vis1 # HPG
Vis2 # neurotransmitter
Vis3 # HPG receptor
#
rm(Input, InputPro, InputProVis, InputRNAVis, Vis1, Vis2, Vis3, CM, DEA_NEs,
   SelectHormone_HPG, SelectReceptor_HPG, SelectReceptor_NT, SelectReceptor_Protein)

##### Figure S3E #####
Input <- Indext %>% filter(!is.na(SampleWES), !is.na(Pathology)) %>% select(Pathology, TMB, CIN)
Input <- melt(Input, id.vars = 'Pathology', value.name = 'Value', variable.name = 'Genomic features')
#
ggboxplot(Input, 'Pathology', 'Value', color = 'Pathology', lwd = 01) +
  scale_color_manual(values = ColPatho) + xlab('') +
  stat_pvalue_manual(compare_means(Value ~ Pathology, data = Input, method = "t.test", group.by = 'Genomic features') %>% 
                       mutate(y.position = c(12,14,16,20,22,24)), label = 'p.adj') +
  theme(legend.position = 'right') +
  facet_wrap(~`Genomic features`)
#
rm(Input)

##### Figure 3B & 3C #####
##### 1.1 prep
InputCNV <- GBC_Main_CNV[, c(filter(Indext, !is.na(SampleWES), Pathology2 == 'AC or AS')$SampleWES,
                             filter(Indext, !is.na(SampleWES), Pathology2 == 'NE')$SampleWES)]
CM = data.table(AC_AS = c(rep(1,154),rep(0,6)),
                NE = c(rep(0,154),rep(1,6)))
rownames(CM) = c(filter(Indext, !is.na(SampleWES), Pathology2 == 'AC or AS')$SampleWES,
                 filter(Indext, !is.na(SampleWES), Pathology2 == 'NE')$SampleWES)

##### 1.2 CNV DEA
# data for Table S4
DEA_CNV <- bioinfoamateur::core_Differential_analysis_continuous(CM, InputCNV, log = T, show.belong = F, p.adj = F)
DEA_CNV$Gene = rownames(DEA_CNV)
# add band info
DEA_CNV$Band = sapply(DEA_CNV$Gene, function(x){
  Band = filter(GeneCoord, gene_name == x)$ChrBand
  Band = ifelse(length(Band) == 0, NA, as.character(Band))
}) 
# add COSMIC info
DEA_CNV$COSMIC <- ifelse(DEA_CNV$Gene %in% ref.cosmic$`Gene Symbol`, 'COSMIC', NA)
# add padj
DEA_CNV$padj <- p.adjust(DEA_CNV$`AC_AS-p-value`, method = 'BH')
# 
write.xlsx(DEA_CNV, file = 'DEA cnv.xlsx', rowNames = T, colNames = T, overwrite = T)

#
DEA_CNV <- bioinfoamateur::core_Differential_analysis_continuous(CM, InputCNV, log = T, show.belong = F, p.adj = T)
DEA_CNV$Gene = rownames(DEA_CNV)
# add band info
DEA_CNV$Band = sapply(DEA_CNV$Gene, function(x){
  Band = filter(GeneCoord, gene_name == x)$ChrBand
  Band = ifelse(length(Band) == 0, NA, as.character(Band))
}) 
# add COSMIC info
DEA_CNV$COSMIC <- ifelse(DEA_CNV$Gene %in% ref.cosmic$`Gene Symbol`, 'COSMIC', NA)

##### 1.3 Transform to display form
# calc mean
InputCNVmean <- aggregate(t(InputCNV), by = list(Pathology = c(rep('AC_AS',154), rep('NE',6))), mean)
# transform
InputCNVmean = column_to_rownames(InputCNVmean, var = 'Pathology') %>% t() %>% as_tibble()
InputCNVmean = cbind(InputCNVmean, DEA_CNV[,c(1,5,6,7)])
InputCNVmean = as_tibble(InputCNVmean)
InputCNVmean$Significant = ifelse(InputCNVmean$`AC_AS-p.adj` <= 0.05,'Sig','ns')
InputCNVmean = InputCNVmean[,-3]

##### 1.4 Vis prep
# matrix to long
InputCNVmean_Vis = filter(InputCNVmean, !is.na(Band))
InputCNVmean_Vis$PseudoRank = 1:nrow(InputCNVmean_Vis)
InputCNVmean_Vis = melt(InputCNVmean_Vis, id.vars = c("Gene",'Band','PseudoRank','COSMIC','Significant'), variable.name = 'Pathology', value.name = 'CNV') %>% 
  as_tibble()
InputCNVmean_Vis <- InputCNVmean_Vis %>% mutate(Chr = paste0('chr',str_split_fixed(Band, '[pq]', 2)[,1]) %>%
                                                  factor(levels = levels(Chr$Chr), ordered = T))

##### 1.5 Vis
#
ggscatter(InputCNVmean_Vis, 'PseudoRank', 'CNV', color = 'Pathology', alpha = 'Significant', size = 0.5) +
  scale_color_manual(values = ColPatho3) +
  scale_alpha_manual(values = c('Sig' = 1, 'ns' = 0.2)) +
  xlab('') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(breaks = c(-0.5,0,0.5,1)) +
  geom_vline(xintercept = cumsum(table(InputCNVmean_Vis$Chr)/2),
             linetype = 2, color = 'grey80') 
# maunally mark COSMIC + Sig
filter(InputCNVmean_Vis, COSMIC == 'COSMIC', Significant == 'Sig') %>% View()

##### 1.6 ORA
### chr6q loss
Output6qLoss <- clusterProfiler::enricher(gene = filter(InputCNVmean_Vis,  Significant == 'Sig', Chr == 'chr6')$Gene,
                                          TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.05) %>% mutate(Significant = -log10(p.adjust))
Output6qLoss$ID = factor(Output6qLoss$ID, levels = rev(Output6qLoss$ID), ordered = T)
### chr16q loss
Output16qLoss <- clusterProfiler::enricher(gene = filter(InputCNVmean_Vis,  Significant == 'Sig', Chr == 'chr16')$Gene,
                                           TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.05) %>% mutate(Significant = -log10(p.adjust))
Output16qLoss$ID = factor(Output16qLoss$ID, levels = rev(Output16qLoss$ID), ordered = T)
# data  for Table S4
write.xlsx(rbind(mutate(Output6qLoss, Belong = '6qloss'),
                 mutate(Output16qLoss, Belong = '16qloss')), file = '6q 16q ORA.xlsx', rowNames = T, colNames = T, overwrite = T)

##### Vis
ggbarplot(Output6qLoss, 'ID', 'Significant', color = 'brown', fill = 'brown') +
  coord_flip() + xlab('') + ylab('-LogFDR') +
  scale_y_continuous(expand = c(0,0)) 
ggbarplot(Output16qLoss, 'ID', 'Significant', color = 'brown', fill = 'brown') +
  coord_flip() + xlab('') + ylab('-LogFDR') +
  scale_y_continuous(expand = c(0,0))
#
rm(InputCNVmean, InputCNVmean_Vis, Output16qLoss, Output6qLoss, DEA_CNV, CM, InputCNV, InputRNA)

##### Figure 3D #####
##### 1.1 input prep
InputRNA <- GBC_Main_RNA_logTPM[, c(filter(Indext, !is.na(SampleRNA), Pathology2 == 'AC or AS')$SampleRNA,
                                    filter(Indext, !is.na(SampleRNA), Pathology2 == 'NE')$SampleRNA)]
CMRNA = data.table(AC_AS = c(rep(1,125),rep(0,5)),
                   NE = c(rep(0,125),rep(1,5)))
rownames(CMRNA) = c(filter(Indext, !is.na(SampleRNA), Pathology2 == 'AC or AS')$SampleRNA,
                    filter(Indext, !is.na(SampleRNA), Pathology2 == 'NE')$SampleRNA)
# Pro Phos
InputPro <- GBC_Main_Pro[, c(filter(Indext, !is.na(SamplePro), Pathology2 == 'AC or AS')$SamplePro,
                             filter(Indext, !is.na(SamplePro), Pathology2 == 'NE')$SamplePro)]
InputPhos <- GBC_Main_Phos_knn[, c(filter(Indext, !is.na(SamplePhos), Pathology2 == 'AC or AS')$SamplePhos,
                                   filter(Indext, !is.na(SamplePhos), Pathology2 == 'NE')$SamplePhos)]
CMPP = data.table(AC_AS = c(rep(1,179),rep(0,6)),
                  NE = c(rep(0,179),rep(1,6)))
rownames(CMPP) = c(filter(Indext, !is.na(SamplePro), Pathology2 == 'AC or AS')$SamplePro,
                   filter(Indext, !is.na(SamplePro), Pathology2 == 'NE')$SamplePro)

##### 1.2 DEA
DEA_RNA <- bioinfoamateur::core_Differential_analysis_continuous(CMRNA, InputRNA,log = T,show.belong = F)
DEA_Pro <- bioinfoamateur::core_Differential_analysis_continuous(CMPP, InputPro,log = T,show.belong = F)
DEA_Phos <- bioinfoamateur::core_Differential_analysis_continuous(CMPP, InputPhos,log = T,show.belong = F)

##### 1.3 GSEA
GSEARNA <- clusterProfiler::GSEA(bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = rownames(DEA_RNA),
                                                                                      LogFC = DEA_RNA$`NE-logFC`)),
                                 TERM2GENE = TotalPathway)@result %>% filter(p.adjust <= 0.05) %>% arrange(desc(NES))
GSEAPro <- clusterProfiler::GSEA(bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = rownames(DEA_Pro),
                                                                                      LogFC = DEA_Pro$`NE-logFC`)),
                                 TERM2GENE = TotalPathway)@result %>% filter(p.adjust <= 0.05) %>% arrange(desc(NES))
GSEAPhos <- clusterProfiler::GSEA(bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = str_split_fixed(rownames(DEA_Phos), ':',2)[,1],
                                                                                       LogFC = DEA_Phos$`NE-logFC`)),
                                  TERM2GENE = TotalPathway)@result %>% filter(p.adjust <= 0.05) %>% arrange(desc(NES))

##### 1.4 Pathway selection
# Up
PathwayCC <- c('REACTOME mRNA Splicing','KEGG Spliceosome','REACTOME Transport of Mature Transcript to Cytoplasm',
               'HALLMARK_E2F_TARGETS','REACTOME Cell Cycle','REACTOME DNA Replication',
               'REACTOME DNA Repair','REACTOME Homologous DNA Pairing and Strand Exchange','REACTOME DNA Double-Strand Break Repair')
PathwayOxi <- c('REACTOME Oxidative Stress Induced Senescence','REACTOME DNA Damage/Telomere Stress Induced Senescence')
PathwayMYC <- c('HALLMARK_MYC_TARGETS_V1','HALLMARK_MYC_TARGETS_V2')
PathwayNeuro <- c('KEGG Neuroactive ligand-receptor interaction',
                  'REACTOME Norepinephrine Neurotransmitter Release Cycle',
                  'REACTOME Neurotransmitter release cycle')
# Down
PathwayMetab <- c('HALLMARK_BILE_ACID_METABOLISM','HALLMARK_FATTY_ACID_METABOLISM',
                  'KEGG Metabolism of xenobiotics by cytochrome P450','KEGG Cholesterol metabolism')
PathwayAdhe <- c('KEGG Focal adhesion','KEGG Cell adhesion molecules','HALLMARK_APICAL_JUNCTION',
                 'REACTOME Extracellular matrix organization')
PathwayImm <- c('REACTOME Cytokine Signaling in Immune system','REACTOME Adaptive Immune System',
                'HALLMARK_INFLAMMATORY_RESPONSE')
PathwaySelected <- c(PathwayCC,PathwayOxi,PathwayMYC,PathwayNeuro,
                     PathwayMetab,PathwayAdhe,PathwayImm)

# re-GSEA
GSEARNA <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = rownames(DEA_RNA),
                                                                LogFC = DEA_RNA$`NE-logFC`))
GSEAPro <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = rownames(DEA_Pro),
                                                                LogFC = DEA_Pro$`NE-logFC`))
GSEAPhos <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = str_split_fixed(rownames(DEA_Phos), ':',2)[,1],
                                                                 LogFC = DEA_Phos$`NE-logFC`))
GSEARNA <- clusterProfiler::GSEA(GSEARNA, TERM2GENE = TotalPathway, pvalueCutoff = 1)@result %>% 
  filter(ID %in% PathwaySelected) %>% select(ID, NES, pvalue, p.adjust) %>%
  arrange(desc(NES)) %>% mutate(Omics = 'RNA')
GSEAPro <- clusterProfiler::GSEA(GSEAPro, TERM2GENE = TotalPathway, pvalueCutoff = 1)@result %>% 
  filter(ID %in% PathwaySelected) %>% select(ID, NES, pvalue, p.adjust) %>%
  arrange(desc(NES)) %>% mutate(Omics = 'Protein')
GSEAPhos <- clusterProfiler::GSEA(GSEAPhos, TERM2GENE = TotalPathway, pvalueCutoff = 1)@result %>% 
  filter(ID %in% PathwaySelected) %>% select(ID, NES, pvalue, p.adjust) %>%
  arrange(desc(NES)) %>% mutate(Omics = 'Phosphorylation')

# Transform
GSEAInput <- rbind(GSEARNA, GSEAPro, GSEAPhos) %>% as.tibble()
# data for Table S4
write.xlsx(GSEAInput, file = 'NE-GBC function.xlsx', rowNames = T, colNames = T, overwrite = T)
#
GSEAInput <- GSEAInput[,-3]
GSEAInput$ID = factor(GSEAInput$ID, levels = PathwaySelected, ordered = T)
GSEAInputNES <- dcast(GSEAInput, ID~Omics, value.var = 'NES') %>% as_tibble() %>% column_to_rownames(var = 'ID')
GSEAInputSig <- dcast(GSEAInput, ID~Omics, value.var = 'p.adjust') %>% as_tibble() %>% column_to_rownames(var = 'ID')

##### 1.5 Vis
Heatmap(GSEAInputNES, cluster_rows = F, cluster_columns = F, name = 'NES',
        rect_gp = gpar(col = "white", lwd = 2),
        col = colorRamp2(c(-1.5,0,1.5), c(ColColor$`Low-Indigo`[8], 'white', ColColor$`High-Red`[9])),
        height = unit(15,'cm'),width = unit(4,'cm')) %>%
  draw(heatmap_legend_side = 'left')
#
rm(InputRNA, InputPro, InputPhos, DEA_RNA, DEA_Phos, DEA_Pro, GSEARNA, GSEAPro, GSEAPhos, GSEAInput,
   PathwayAdhe,PathwayCC,PathwayImm,PathwayMYC,PathwayMetab,PathwayNeuro,PathwayOxi,PathwaySelected,
   CMRNA, CMPP)

##### Figure S3F #####
##### 1.1 TF prep
library(decoupleR)
#
load('../OmniPathR_collecTRI.RData') # Inputnet
RefTF = Inputnet
RefTFFamily = fread('../Homo_sapiens_TF.txt') %>% filter(Symbol %in% RefTF$source) %>% as_tibble()

##### 1.2 TF activity inferring
InputTF = GBC_Main_RNA_logTPM[, na.omit(Index$SampleRNA)] %>% as.matrix()
# infer
Output <- run_wmean(mat=InputTF, net=RefTF, .source='source', .target='target', .mor='mor', times = 100, minsize = 5) 
Output
# long2mtx
Outputmtx <- dcast(filter(Output, statistic == 'wmean'), source ~ condition, value.var = "score") %>% 
  column_to_rownames(var = 'source')
# data for Table S4
write.xlsx(Outputmtx, file = 'TF activity.xlsx', rowNames = T, colNames = T, overwrite = T)

#
detach("package:decoupleR", unload = TRUE)
detach("package:OmnipathR", unload = TRUE)
rm(InputTF)

##### 1.3 DEA of TFs
CM <- data.frame(AC_AS = ifelse(colnames(Outputmtx) %in% filter(Index, as.character(Ident) %in% c('AC','AS'))$SampleRNA, 1, 0),
                 NE = ifelse(colnames(Outputmtx) %in% filter(Index, as.character(Ident) == 'NE')$SampleRNA, 1, 0),
                 Adjacent = ifelse(colnames(Outputmtx) %in% filter(Index, as.character(Ident) == 'Adjacent')$SampleRNA, 1, 0))
rownames(CM) = colnames(Outputmtx)
## DEA
DEA_TF <- bioinfoamateur::core_Differential_analysis_continuous(CM, Outputmtx, 
                                                                log = T, p.adj = F, show.belong = F)
DEA_TF_belong <- bioinfoamateur::core_Differential_analysis_continuous(CM, Outputmtx, 
                                                                       log = T, p.adj = F, show.belong = T)
DEA_TF_belong$TF = rownames(DEA_TF_belong)

##### 1.4 TF selection and vis
Input = Outputmtx[c(rownames(filter(DEA_TF_belong, belong == 1) %>% arrange(desc(`AC_AS-logFC`)))[1:8],
                    rownames(filter(DEA_TF_belong, belong == 2) %>% arrange(desc(`NE-logFC`)))[1:4], 'ASCL1','ASCL2','MYC','NEUROG1',
                    rownames(filter(DEA_TF_belong, belong == 3) %>% arrange(desc(`Adjacent-logFC`)))[1:8]), ] %>% t() %>% as.data.frame()

# add patho
Input$Pathology = apply(CM, 1, function(x){colnames(CM)[which(x == 1)]})
# Calc Mean
Input = aggregate(Input[,-ncol(Input)], by = list(Pathology = Input$Pathology), mean) %>% 
  column_to_rownames(var = 'Pathology') %>% scale() %>% t() %>% as.data.frame()
Input = Input[,c(1,3,2)]
# Vis
Heatmap(t(Input), name = 'TF Activity\n(Z-score)', 
        column_split = c(rep(1,8),rep(2,8),rep(3,8)), column_title = NULL,
        row_names_side = 'left', column_names_rot = 30,
        column_names_centered = F, cluster_rows = F, cluster_columns = F,
        colorRamp2(c(-1.5, 0, 1.5), c(ColColor$`Low-LightBlue`[7], "white", ColColor$`High-Red`[7])),
        rect_gp = gpar(col = "black", lwd = 1),
        height = unit(3,'cm'), width = unit(16, 'cm')) 
#
rm(Input, DEA_TF, DEA_TF_belong, CM, Output, Outputmtx)

##### Figure S3G #####
##### 1.1 TF activity
library(decoupleR)
load('../OmniPathR_collecTRI.RData') # Inputnet
RefTF = Inputnet
rm(Inputnet)
RefTFFamily = fread('../Homo_sapiens_TF.txt') %>% filter(Symbol %in% RefTF$source) %>% as_tibble()
#
InputTF = GBC_Main_RNA_logTPM[, na.omit(Index$SampleRNA)] %>% as.matrix()
# infer
Output <- run_wmean(mat=InputTF, net=RefTF, .source='source', .target='target', .mor='mor', times = 100, minsize = 5) 
Output
# long2mtx
Outputmtx <- dcast(filter(Output, statistic == 'wmean'), source ~ condition, value.var = "score") %>% 
  column_to_rownames(var = 'source')
#
detach("package:decoupleR", unload = TRUE)
detach("package:OmnipathR", unload = TRUE)
rm(InputTF)

##### 1.2 TF DEA
CM <- data.frame(AC_AS = ifelse(colnames(Outputmtx) %in% filter(Index, as.character(Ident) %in% c('AC','AS'))$SampleRNA, 1, 0),
                 NE = ifelse(colnames(Outputmtx) %in% filter(Index, as.character(Ident) == 'NE')$SampleRNA, 1, 0),
                 Adjacent = ifelse(colnames(Outputmtx) %in% filter(Index, as.character(Ident) == 'Adjacent')$SampleRNA, 1, 0))
rownames(CM) = colnames(Outputmtx)
# DEA
DEA_TF <- bioinfoamateur::core_Differential_analysis_continuous(CM, Outputmtx, log = T, p.adj = F, show.belong = F)
DEA_TF$TF = rownames(DEA_TF)
DEA_TF_belong <- bioinfoamateur::core_Differential_analysis_continuous(CM, Outputmtx, log = T, p.adj = F, show.belong = T)
DEA_TF_belong$TF = rownames(DEA_TF_belong)

##### 1.3 arrange protein abundance of TFs
InputPro <- GBC_Main_Pro[intersect(DEA_TF$TF, rownames(GBC_Main_Pro)), 
                         c(filter(Index, !is.na(SamplePro), Pathology2 == 'AC or AS', Type == 'Cancer')$SamplePro,
                           filter(Index, !is.na(SamplePro), Pathology2 == 'NE', Type == 'Cancer')$SamplePro,
                           filter(Index, !is.na(SamplePro), Type == 'Adjacent')$SamplePro)]
CM = data.table(AC_AS = c(rep(1,179),rep(0,6),rep(0,130)),
                NE = c(rep(0,179),rep(1,6),rep(0,130)),
                Adjacent = c(rep(0,179),rep(0,6),rep(1,130)))
rownames(CM) = c(filter(Index, !is.na(SamplePro), Pathology2 == 'AC or AS', Type == 'Cancer')$SamplePro,
                 filter(Index, !is.na(SamplePro), Pathology2 == 'NE', Type == 'Cancer')$SamplePro,
                 filter(Index, !is.na(SamplePro), Type == 'Adjacent')$SamplePro)
# DEA
DEA_TF_Pro <- bioinfoamateur::core_Differential_analysis_continuous(CM, InputPro, log = T, p.adj = T, show.belong = F)
DEA_TF_Pro$TF = rownames(DEA_TF_Pro)
DEA_TF_Pro_belong <- bioinfoamateur::core_Differential_analysis_continuous(CM, InputPro, log = T, p.adj = T, show.belong = T)
DEA_TF_Pro_belong$TF = rownames(DEA_TF_Pro_belong)

##### 1.4 Crucial TF identification
# add FC and Pval
Input = tibble(TF = intersect(DEA_TF_Pro$TF, DEA_TF$TF))
Input$LogFC_Activity = sapply(Input$TF, function(x){filter(DEA_TF, TF == x)$`NE-logFC`}) %>% as.numeric()
Input$LogFC_Abundance = sapply(Input$TF, function(x){filter(DEA_TF_Pro, TF == x)$`NE-logFC`}) %>% as.numeric()
Input$Pval_Activity = sapply(Input$TF, function(x){filter(DEA_TF, TF == x)$`NE-p-value`}) %>% as.numeric()
Input$Padj_Abundance = sapply(Input$TF, function(x){filter(DEA_TF_Pro, TF == x)$`NE-p.adj`}) %>% as.numeric()
# judge sig
Input$Sig_Activity = ifelse(Input$Pval_Activity <= 0.05 & Input$LogFC_Activity > 0, '*', 'ns')
Input$Sig_Abundance = ifelse(Input$Padj_Abundance <= 0.05 & Input$LogFC_Abundance > 0, '*', 'ns')
# add marker indicating sig
Input$Marker = case_when(Input$Sig_Activity == '*' & Input$Sig_Abundance == '*' ~ 'Significant both in TF Activity and Abundance',
                         Input$Sig_Activity == 'ns' & Input$Sig_Abundance == 'ns' ~ 'Not Significant',
                         Input$Sig_Activity == '*' & Input$Sig_Abundance == 'ns' ~ 'Significant only in TF Activity', 
                         Input$Sig_Activity == 'ns' & Input$Sig_Abundance == '*' ~ 'Significant only in TF Abundance') %>%
  factor(levels = c('Significant both in TF Activity and Abundance', 'Significant only in TF Activity',
                    'Significant only in TF Abundance','Not Significant'), ordered = T)

##### 2.4 Crucial TF visualization #####
library(ggrepel)
ggscatter(Input, 'LogFC_Activity', 'LogFC_Abundance', color = 'Marker') +
  scale_color_manual(values = c('Significant both in TF Activity and Abundance' = 'red',
                                'Significant only in TF Activity' = ColColor$`High-Orange`[5],
                                'Significant only in TF Abundance' = ColColor$`High-Pink`[5],
                                'Not Significant' = 'grey80')) +
  geom_hline(yintercept = 0, linetype = 2, color = 'black') +
  geom_vline(xintercept = 0, linetype = 2, color = 'black') +
  # geom_hline(yintercept = c(0.5,-0.5), linetype = 2, color = 'grey') +
  # geom_vline(xintercept = c(0.5,-0.5), linetype = 2, color = 'grey') +
  ylab('LogFC (Protein Abundance)') + xlab('LogFC (TF Activity)') +
  theme(legend.position = 'left') +
  geom_text_repel(aes(x=LogFC_Activity, y=LogFC_Abundance, 
                      label=TF), data = filter(Input, Marker == 'Significant both in TF Activity and Abundance'),
                  min.segment.length = 0, seed = 42, box.padding = 0.5)
# 
write.xlsx(Input, file = 'Crucial TF identification.xlsx', rowNames = T, colNames = T, overwrite = T)

##### Figure S3H #####
##### 1.1 SCLC data
# Liu Q, Zhang J, Guo C, et al. Proteogenomic characterization of small cell lung cancer identifies biological insights and subtype-specific therapeutic strategies. Cell. 2024;187(1):184-203.e28. doi:10.1016/j.cell.2023.12.004
ExternalIndex <- openxlsx::read.xlsx('SCLC_Omics.xlsx', sheet = 2, rowNames = F)
ExternalRNA <- openxlsx::read.xlsx('SCLC_Omics.xlsx', sheet = 5, rowNames = T)
ExternalProtein <- openxlsx::read.xlsx('SCLC_Omics.xlsx', sheet = 6, rowNames = T)
ExternalPhos <- openxlsx::read.xlsx('SCLC_Omics.xlsx', sheet = 7, rowNames = T)
#
ExternalIndexNE <- tibble(Sample = rep(ExternalIndex$Sample.ID, 2),
                          OmicID = c(paste0('T',ExternalIndex$Sample.ID),
                                     paste0('N',ExternalIndex$Sample.ID)),
                          Pathology = rep(ExternalIndex$Histologic.type, 2))
ExternalIndexNE$Pathology = ifelse(ExternalIndexNE$Pathology=='small cell lung cancer','SCLC','CSCLC') %>%
  factor(levels = c('CSCLC','SCLC'), ordered = T)
InputExternal <- ExternalIndexNE

##### 1.2 vis
# MEIS1
InputExternal$mRNA_MEIS1 = sapply(InputExternal$OmicID, function(x){ifelse(x %in% colnames(ExternalRNA),as.numeric(ExternalRNA['MEIS1',x]),NA)}) %>% as.numeric()
InputExternal$Protein_MEIS1 = sapply(InputExternal$OmicID, function(x){ifelse(x %in% colnames(ExternalProtein),as.numeric(ExternalProtein['MEIS1',x]),NA)}) %>% as.numeric()
# MYC
InputExternal$mRNA_MYC = sapply(InputExternal$OmicID, function(x){ifelse(x %in% colnames(ExternalRNA),as.numeric(ExternalRNA['MYC',x]),NA)}) %>% as.numeric()
InputExternal$Protein_MYC = sapply(InputExternal$OmicID, function(x){ifelse(x %in% colnames(ExternalProtein),as.numeric(ExternalProtein['MYC',x]),NA)}) %>% as.numeric()
InputExternal$Phos_MYC_S62 = sapply(InputExternal$OmicID, function(x){ifelse(x %in% colnames(ExternalPhos),as.numeric(ExternalPhos['MYC:S62',x]),NA)}) %>% as.numeric()
##### Vis
ggplot(data = InputExternal, aes(Protein_MEIS1, mRNA_MYC, color = Pathology)) +
  geom_point(alpha = 0.5, shape = 16, size = 3) +
  geom_smooth(method = 'lm', se = F) + 
  stat_cor(aes(color = Pathology), method = 'pearson') +
  scale_color_manual(values = c('CSCLC'='#00A08A','SCLC'='#F2AD00')) +
  scale_x_continuous(limits = c(-1.2,1.5)) +
  ylab('mRNA MYC') +
  xlab('Protein MEIS1') +
  theme(legend.position="right") +
  theme_bw() 
# data for Table S4
write.xlsx(InputExternal, file = 'SCLC data.xlsx', rowNames = T, colNames = T, overwrite = T)

##### Figure S3I #####
##### 1.1 prep
Input = t(GBC_Main_Pro[rownames(GBC_Main_Pro) != 'MEIS1',na.omit(Indext)$SamplePro])
apply(Input, 2, function(x){cor(as.numeric(x), as.numeric(GBC_Main_Pro['MEIS1',na.omit(Indext)$SamplePro]))}) %>% as.numeric()
Output = tibble(Gene = colnames(Input),
                Cor = apply(Input, 2, function(x){cor(as.numeric(x), 
                                                      as.numeric(GBC_Main_Pro['MEIS1',na.omit(Indext)$SamplePro]))}) %>% 
                  as.numeric(),
                FDR = apply(Input, 2, function(x){cor.test(as.numeric(x), 
                                                           as.numeric(GBC_Main_Pro['MEIS1',na.omit(Indext)$SamplePro]))$p.value}) %>% 
                  as.numeric())
Output$FDR = p.adjust(Output$FDR, method = 'fdr')
##### gsea input
InputGSEA = bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = Output$Gene, Cor = Output$Cor))
OutputGSEA = clusterProfiler::GSEA(InputGSEA, TERM2GENE = TotalPathway)@result %>%
  filter(NES > 0, p.adjust <= 0.05)

##### 1.2 Calc Pathway
GBC_ssGSEA_Pro <- GSVA::gsva(as.matrix(GBC_Main_Pro), TotalPathwayGSVA, method = 'gsva')
Input = t(GBC_ssGSEA_Pro[,na.omit(Indext)$SamplePro])
Output = tibble(Pathway = colnames(Input),
                Cor = apply(Input, 2, function(x){cor(as.numeric(x), 
                                                      as.numeric(GBC_Main_Pro['MEIS1',na.omit(Indext)$SamplePro]))}) %>% 
                  as.numeric(),
                FDR = apply(Input, 2, function(x){cor.test(as.numeric(x), 
                                                           as.numeric(GBC_Main_Pro['MEIS1',na.omit(Indext)$SamplePro]))$p.value}) %>% 
                  as.numeric())
Output$FDR = p.adjust(Output$FDR, method = 'fdr')

##### pathway select
SelectPathway <- c('REACTOME Nervous system development',
                   'REACTOME Neuronal System',
                   'REACTOME Developmental Biology',
                   'KEGG Axon guidance',
                   'REACTOME Neurotransmitter receptors and postsynaptic signal transmission',
                   'KEGG Neurotrophin signaling pathway')

##### 1.3 Vis
# Order by MEIS1
TempOrder = GBC_Main_Pro['MEIS1', na.omit(Indext$SamplePro)] %>% as.numeric()
names(TempOrder) = na.omit(Indext$SamplePro)
TempOrder = sort(TempOrder, decreasing = T)
#
Input = GBC_ssGSEA_Pro[SelectPathway, names(TempOrder)]
# Anno prep
AnnoExp = columnAnnotation(MEIS1 = TempOrder,
                           NE = ifelse(names(TempOrder) %in% filter(Indext, Pathology2 == 'NE')$SamplePro, 'NE', 'AC_AS'),
                           col = list(
                             MEIS1 = colorRamp2(c(-1.5, 0, 1.5), c(ColColor$`Low-Indigo`[8], "white",
                                                                   ColColor$`High-Red`[8])),
                             NE = c('AC_AS' = ColPatho3[1], 'NE' = ColPatho3[2])
                           ), 
                           gp = gpar(col = "white"),
                           simple_anno_size = unit(1, "cm"))
Heatmap(t(scale(t(Input))), cluster_rows = F, cluster_columns = F, name = 'Pathway activity',
        top_annotation = AnnoExp, show_column_names = F,
        colorRamp2(c(-1.5, 0, 1.5), c(ColColor$`Low-LightBlue`[3], "white",
                                      ColColor$`High-Red`[8])),
        height = unit(6,'cm'), width = unit(14,'cm')) %>%
  draw(heatmap_legend_side = 'left',
       annotation_legend_side = 'left') 
# data for Table S4
Input <- t(Input) %>% as.data.frame()
write.xlsx(Input, file = 'ssGSEA MEIS1.xlsx', rowNames = T, colNames = T, overwrite = T)
#
rm(Input, TempOrder, AnnoExp,SelectPathway,Output,
   InputGSEA,OutputGSEA)


##### Figure S3J #####
Signature = as.data.frame(TotalPathwayGSVA[,'KEGG Axon guidance'])
colnames(Signature) = 'KEGG Axon guidance'
# Calc
Output <- GSVA::gsva(as.matrix(ExternalProtein[,InputExternal$OmicID]), SignatureNE, method = 'gsva')
Output <- tibble(Sample = colnames(Output),
                 NEScore = as.numeric(Output[1,]))
# add metadata
InputExternal$NEScore = sapply(InputExternal$OmicID, function(x){filter(Output, Sample == x)$NEScore}) %>% as.numeric()

#
ggscatter(InputExternal, x = "Protein_MEIS1", y = "NEScore",
          color = "black", size = 3, 
          add = "reg.line",  
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE, 
          cor.coef = TRUE, 
          cor.coeff.args = list(method = "spearman", label.x = 1, label.sep = "\n")) +
  ylab('KEGG Axon guidance') 

##### Figure S3K #####
library(TCseq)
##### 1.1 cell line data read in
# rna
raw.rna <- fread('../counts.txt') %>% as.data.frame()
# pro
raw.pro.fill <- fread('../20240819.GBC.pro.csv') %>% as.data.frame()

##### 1.2 trends identification
load('temp.rda')
temp.trend <- timeclustplot(input.pi.ym.ssgsea.pro, categories="Time", cols = 3)
plot_grid(temp.trend[[6]], temp.trend[[14]], nrow = 1)
# data for Table S4
temp.out <- data.frame(Protein = names(input.pi.ym.pro@cluster),
                       Cluster = unname(input.pi.ym.pro@cluster)) %>%
  filter(Cluster %in% c(6,14)) %>% arrange(Cluster, Protein)
write.xlsx(temp.out, file = 'TCseq trends.xlsx', rowNames = T, colNames = T, overwrite = T)

##### Figure S3L #####
input.ns <- index.pro %>% mutate(`Neuronal system` = as.numeric(ssgsea.pro['REACTOME Neuronal System',index.pro$Sample]),
                                 `Nerve system development` = as.numeric(ssgsea.pro['REACTOME Nervous system development',index.pro$Sample]))
input.ns$Success <- as.factor(input.ns$Success)
filtered_df <- input.ns %>%
  dplyr::filter(Success == "Yes", Condition == "MEIS1") %>%
  dplyr::select(TimeLine, `Neuronal system`, `Nerve system development`)
long_df <- filtered_df %>%
  tidyr::pivot_longer(
    cols = c(`Neuronal system`, `Nerve system development`),
    names_to = "System",
    values_to = "Score"
  )
long_df <- long_df %>%
  group_by(System) %>%
  mutate(Score = scale(Score))
# theme
arial_theme <- theme(
  text = element_text(family = "Arial", color = 'black', size = 12), 
  panel.grid = element_blank(),         
  panel.background = element_rect(fill = "white", color = "white"),
  axis.line = element_line(color = "black"),  
  axis.ticks = element_line(color = "black"),  
  legend.position = "right"  
)
##### Vis
combined_plot <- ggplot(long_df, aes(x = TimeLine, y = Score, color = System, group = System)) +
  scale_color_manual(values = c('Neuronal system' = ColJournal$Nature[1], 'Nerve system development' = ColJournal$Nature[2])) +
  stat_summary(fun = mean, geom = "line", size = 1) +  
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +  
  labs(
    title = "System Score Changes Over Time (Condition: MEIS1, Success: Yes)",
    x = NULL, y = "Normalized pathway score", color = 'Pathway'
  ) + 
  arial_theme 
# data for Table S4
write.xlsx(long_df, file = 'Cell line 2 function.xlsx', rowNames = T, colNames = T, overwrite = T)
#
print(combined_plot)





