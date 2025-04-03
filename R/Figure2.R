##### Figure 2 and Figure S2
rm(list = ls()); try(dev.off(), silent = T)
# packages loading
library(tidyverse); library(data.table)
library(ggpubr); library(cowplot); library(ggthemes); library(ggplot2); library(ggrepel)
library(maftools)
library(ComplexHeatmap); library(circlize)
library(openxlsx)

##### Prep1. Data Loading #####
# 1 omics
load('../DataMain.Rda')
load('../GBC_Benign_Cohort.RData')
# data for Table S3
write.xlsx(GBC_Benign_Pro, file = 'Benign protein.xlsx', rowNames = T, colNames = T, overwrite = T)
write.xlsx(GBC_Benign_Phos_knn, file = 'Benign phos.xlsx', rowNames = T, colNames = T, overwrite = T)
# 2 color
load('../Colors (ggsci).RData')
library(wesanderson)
ColCNA = c(ColColor$`Low-Blue`[7], ColColor$`Low-Blue`[3], ColColor$`Low-BlueGreen`[5],
           ColColor$`High-Red`[3], ColColor$`High-Red`[7])
names(ColCNA) = c('Deletion','ShallowDeletion','Diploid','Gain','Amplification')
ColType2 = c(ColJournal$Nature[3],ColJournal$Nature[7],ColJournal$Nature[2],ColJournal$Nature[5],ColJournal$Nature[8])
names(ColType2) = c('Normal','Adjacent','Benign','ERBB2_others','ERBB2_Amp')
# 3 Genesets
load('../Hallmark_KEGG_Reactome.RData')
# 4 COSMIC Cancer Gene Census
ref.cosmic <- fread('../COSMIC genes.csv')

##### Prep extra. calc ssgsea #####
ssGSEAInput_Pro <- cbind(GBC_Main_Pro, GBC_Benign_Pro)
ssGSEAInput_Phos <- cbind(GBC_Main_ProLevelPhos_knn[intersect(rownames(GBC_Main_ProLevelPhos_knn), rownames(GBC_Benign_ProLevelPhos_knn)),],
                          GBC_Benign_ProLevelPhos_knn[intersect(rownames(GBC_Main_ProLevelPhos_knn), rownames(GBC_Benign_ProLevelPhos_knn)),])
# ssgsea
ssGSEAOutput_Pro <- GSVA::gsva(as.matrix(ssGSEAInput_Pro), TotalPathwayGSVA, method = 'ssgsea') %>% as.data.frame()
ssGSEAOutput_Phos <- GSVA::gsva(as.matrix(ssGSEAInput_Phos), TotalPathwayGSVA, method = 'ssgsea') %>% as.data.frame()

##### Prep2. Index (N = 193 that have proteomic data) #####
IndexERBB <- rbind(select(GBC_Main_Index, 1,2,4,6) %>% 
                     bioinfoamateur::dfm_change_colnames(c('Sample','WES','RNA','Pro')), 
                   mutate(GBC_Main_Index, Sample = paste0(Sample, '_P')) %>% select(1,3,5,7) %>% 
                     bioinfoamateur::dfm_change_colnames(c('Sample','WES','RNA','Pro')))
IndexERBB$Tumor <- ifelse(!grepl(fixed('_P'),IndexERBB$Sample), 'Yes', NA)
# retain samples having proteomic data
IndexERBB <- filter(IndexERBB, !is.na(Pro), !is.na(WES))

##### add ERBB mutation
# EGFR
IndexERBB$Mut_EGFR <- case_when(IndexERBB$Sample %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'EGFR')$Tumor_Sample_Barcode ~ 'Yes',
                                IndexERBB$Sample %in% GBC_Main_Mutation$Tumor_Sample_Barcode ~ NA,
                                !is.na(IndexERBB$Sample) ~ NA)
IndexERBB$Mut_EGFR <- ifelse(IndexERBB$Mut_EGFR == "Yes", filter(GBC_Main_Mutation, Tumor_Sample_Barcode == IndexERBB$Sample)$Variant_Classification, NA)
# ERBB2
IndexERBB$Mut_ERBB2 <- case_when(IndexERBB$Sample %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB2')$Tumor_Sample_Barcode ~ 'Yes',
                                 IndexERBB$Sample %in% GBC_Main_Mutation$Tumor_Sample_Barcode ~ NA,
                                 !is.na(IndexERBB$Sample) ~ NA)
IndexERBB$Mut_ERBB2 <- ifelse(IndexERBB$Mut_ERBB2 == "Yes", filter(GBC_Main_Mutation, Tumor_Sample_Barcode == IndexERBB$Sample)$Variant_Classification, NA)
# ERBB3
IndexERBB$Mut_ERBB3 <- case_when(IndexERBB$Sample %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB3')$Tumor_Sample_Barcode ~ 'Yes',
                                 IndexERBB$Sample %in% GBC_Main_Mutation$Tumor_Sample_Barcode ~ NA,
                                 !is.na(IndexERBB$Sample) ~ NA)
IndexERBB$Mut_ERBB3 <- ifelse(IndexERBB$Mut_ERBB3 == "Yes", filter(GBC_Main_Mutation, Tumor_Sample_Barcode == IndexERBB$Sample)$Variant_Classification, NA)
# ERBB4
IndexERBB$Mut_ERBB4 <- case_when(IndexERBB$Sample %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB4')$Tumor_Sample_Barcode ~ 'Yes',
                                 IndexERBB$Sample %in% GBC_Main_Mutation$Tumor_Sample_Barcode ~ NA,
                                 !is.na(IndexERBB$Sample) ~ NA)
IndexERBB$Mut_ERBB4 <- ifelse(IndexERBB$Mut_ERBB4 == "Yes", filter(GBC_Main_Mutation, Tumor_Sample_Barcode == IndexERBB$Sample)$Variant_Classification, NA)
# ERBB Pathway
IndexERBB$Mut_ERBBPathway <- case_when(IndexERBB$Sample %in% filter(GBC_Main_Mutation, Hugo_Symbol %in% filter(TotalPathway, Pathway == 'KEGG ErbB signaling pathway')$Gene)$Tumor_Sample_Barcode ~ 'Yes',
                                       IndexERBB$Sample %in% GBC_Main_Mutation$Tumor_Sample_Barcode ~ 'No',
                                       !is.na(IndexERBB$Sample) ~ NA) %>%
  factor(levels = c('Yes','No'), ordered = T)

##### add ERBBs CNV
SeparatedCNV <- function(vec) {
  vec = case_when(vec == -2 ~ 'Deletion',
                  vec == -1 ~ 'ShallowDeletion',
                  vec == 0 ~ 'Diploid', 
                  vec == 1 ~ 'Gain',
                  vec == 2 ~ 'Amplification') %>%
    factor(levels = c('Amplification','Gain', 'Diploid', 'ShallowDeletion', 'Deletion'), ordered = T)
  return(vec)
}
# EGFR
IndexERBB$CNV_EGFR <- sapply(IndexERBB$WES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_GBC_Main_GBC_Main_CNV_cate['EGFR',x]), NA)})
IndexERBB$CNV_EGFR <- SeparatedCNV(IndexERBB$CNV_EGFR)
# ERBB2
IndexERBB$CNV_ERBB2 <- sapply(IndexERBB$WES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_GBC_Main_GBC_Main_CNV_cate['ERBB2',x]), NA)})
IndexERBB$CNV_ERBB2 <- SeparatedCNV(IndexERBB$CNV_ERBB2)
# ERBB3
IndexERBB$CNV_ERBB3 <- sapply(IndexERBB$WES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_GBC_Main_GBC_Main_CNV_cate['ERBB3',x]), NA)})
IndexERBB$CNV_ERBB3 <- SeparatedCNV(IndexERBB$CNV_ERBB3)
# ERBB4
IndexERBB$CNV_ERBB4 <- sapply(IndexERBB$WES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_GBC_Main_GBC_Main_CNV_cate['ERBB4',x]), NA)})
IndexERBB$CNV_ERBB4 <- SeparatedCNV(IndexERBB$CNV_ERBB4)

##### add cli paras
# Age
IndexERBB$Age = sapply(IndexERBB$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Age}) %>% as.numeric()
# Gender
IndexERBB$Gender = sapply(IndexERBB$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Sex}) %>% as.character() %>%
  factor(levels = c('Male','Female'), ordered = T)
# Stage
IndexERBB$Stage = sapply(IndexERBB$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Tumor_stage2}) %>% as.numeric()
IndexERBB$Stage = case_when(IndexERBB$Stage == 1 ~ 'I', IndexERBB$Stage == 2 ~ 'II',
                            IndexERBB$Stage == 3 ~ 'III', IndexERBB$Stage == 4 ~ 'IV') %>%
  factor(levels = c('I','II','III','IV'), ordered = T)
# Patho
IndexERBB$Pathology = sapply(IndexERBB$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Pathological_type}) %>% as.character()
IndexERBB$Pathology = case_when(IndexERBB$Pathology == 'adenocarcinoma' ~ 'AC',
                                IndexERBB$Pathology == 'adenosquanmous carcinoma' ~ 'AS',
                                IndexERBB$Pathology == 'neuroendocrine carcinoma' ~ 'NE',
                                !is.na(IndexERBB$Pathology) ~ 'Others') %>%
  factor(levels = c('AC','AS','NE','Others'), ordered = T)
# LM
IndexERBB$LM <- sapply(IndexERBB$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Regional_invasion}) %>% as.numeric()
IndexERBB$LM = ifelse(IndexERBB$LM == 1, 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)

##### add ERBB2 mRNA and Protein levels and pro-level phos
# EGFR
IndexERBB$mRNA_EGFR <- sapply(IndexERBB$RNA, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['EGFR',x]), NA)}) %>% as.numeric()
IndexERBB$Protein_EGFR <- sapply(IndexERBB$Pro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['EGFR',x]), NA)}) %>% as.numeric()
IndexERBB$ProPhos_EGFR <- sapply(IndexERBB$Pro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_ProLevelPhos_knn['EGFR',x]), NA)}) %>% as.numeric()
# ERBB2
IndexERBB$mRNA_ERBB2 <- sapply(IndexERBB$RNA, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['ERBB2',x]), NA)}) %>% as.numeric()
IndexERBB$Protein_ERBB2 <- sapply(IndexERBB$Pro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['ERBB2',x]), NA)}) %>% as.numeric()
IndexERBB$ProPhos_ERBB2 <- sapply(IndexERBB$Pro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_ProLevelPhos_knn['ERBB2',x]), NA)}) %>% as.numeric()
# ERBB3
IndexERBB$mRNA_ERBB3 <- sapply(IndexERBB$RNA, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['ERBB3',x]), NA)}) %>% as.numeric()
IndexERBB$Protein_ERBB3 <- sapply(IndexERBB$Pro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['ERBB3',x]), NA)}) %>% as.numeric()
IndexERBB$ProPhos_ERBB3 <- sapply(IndexERBB$Pro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_ProLevelPhos_knn['ERBB3',x]), NA)}) %>% as.numeric()
# ERBB4
IndexERBB$mRNA_ERBB4 <- sapply(IndexERBB$RNA, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['ERBB4',x]), NA)}) %>% as.numeric()
IndexERBB$Protein_ERBB4 <- sapply(IndexERBB$Pro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['ERBB4',x]), NA)}) %>% as.numeric()
IndexERBB$ProPhos_ERBB4 <- sapply(IndexERBB$Pro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_ProLevelPhos_knn['ERBB4',x]), NA)}) %>% as.numeric()

##### add TMB
GBC_Main_Maf <- maftools::read.maf(maf = GBC_Main_Mutation)
TempTMB <- tmb(GBC_Main_Maf)
IndexERBB$TMB <- sapply(IndexERBB$Sample, function(x){filter(TempTMB, Tumor_Sample_Barcode == x)$total_perMB %>% as.numeric()}) %>% as.numeric()
rm(TempTMB)

##### Figure 2A #####
##### Input prep
InputOp <- filter(IndexERBB[,1:14], Tumor == 'Yes')
InputOp[is.na(InputOp)] = ''
InputOp$CNV_EGFR = as.character(InputOp$CNV_EGFR)
InputOp$CNV_EGFR = ifelse(InputOp$CNV_EGFR == 'Amplification', 'Amplification','')
InputOp$CNV_ERBB2 = as.character(InputOp$CNV_ERBB2) 
InputOp$CNV_ERBB2 = ifelse(InputOp$CNV_ERBB2 == 'Amplification', 'Amplification','')
InputOp$CNV_ERBB3 = as.character(InputOp$CNV_ERBB3)
InputOp$CNV_ERBB3 = ifelse(InputOp$CNV_ERBB3 == 'Amplification', 'Amplification','')
InputOp$CNV_ERBB4 = as.character(InputOp$CNV_ERBB4) 
InputOp$CNV_ERBB4 = ifelse(InputOp$CNV_ERBB4 == 'Amplification', 'Amplification','')

#
InputOp2 = tibble(EGFR = rep('',166), ERBB2 = rep('',166),
                  ERBB3 = rep('',166), ERBB4 = rep('',166))
for (i in 1:nrow(InputOp2)) {
  for (j in 1:4){
    TempMut = InputOp[i,(j+5)]; TempCNV = as.character(InputOp[i,(j+10)]) 
    if (TempMut != '' & TempCNV != ''){ InputOp2[i,j] = paste0(TempCNV,';',TempMut) }
    if (TempMut != '' & TempCNV == ''){ InputOp2[i,j] = TempMut }
    if (TempMut == '' & TempCNV == ''){ InputOp2[i,j] = '' }
    if (TempMut == '' & TempCNV != ''){ InputOp2[i,j] = TempCNV }
  }
}
rm(i, j, TempMut, TempCNV)
# arrange
InputOp2 = as.data.frame(InputOp2)
rownames(InputOp2) = InputOp$Sample
InputOp2 = t(InputOp2) %>% as.matrix()

##### Vis prep
HM_colAnno = columnAnnotation(Stage = filter(IndexERBB, Tumor == 'Yes')$Stage %>% as.character(),
                              Age = filter(IndexERBB, Tumor == 'Yes')$Age,
                              Gender = filter(IndexERBB, Tumor == 'Yes')$Gender,
                              `Liver invasion` = filter(IndexERBB, Tumor == 'Yes')$LM %>% as.character(),
                              TMB = filter(IndexERBB, Tumor == 'Yes')$TMB,
                              `ERBB2 mRNA` = filter(IndexERBB, Tumor == 'Yes')$mRNA_ERBB2,
                              `ERBB2 Protein` = filter(IndexERBB, Tumor == 'Yes')$Protein_ERBB2,
                              col = list(
                                Stage = c('I' = ColColor$`High-Red`[1], 'II' = ColColor$`High-Red`[3],
                                          'III' = ColColor$`High-Red`[5], 'IV' = ColColor$`High-Red`[7]),
                                Age = colorRamp2(c(30,90), c(ColColor$`High-Orange`[1], ColColor$`High-Orange`[8])),
                                Gender = c('Male' = 'black', 'Female' = 'white'),
                                `Liver invasion` = c('No' = 'white', 'Yes' = 'black'),
                                TMB = colorRamp2(c(0, 3, 6), c("white", ColColor$`Low-Indigo`[4],  ColColor$`Low-Indigo`[9])),
                                `ERBB2 mRNA` = colorRamp2(c(2,5,8), c(ColColor$`Low-LightBlue`[8], 'white', ColColor$`High-Red`[8])),
                                `ERBB2 Protein` = colorRamp2(c(-2,0,2), c(ColColor$`Low-LightBlue`[8], 'white', ColColor$`High-Red`[8]))
                              ),
                              simple_anno_size = unit(0.5, "cm"),
                              gp = gpar(col = "white", lwd = 0.2),
                              annotation_name_side = 'left')

##### HM Layout
HM_Col = c(Amplification = ColColor$`High-Red`[10],
           # Gain = ColColor$`High-Red`[1],
           # ShallowDeletion = ColColor$`Low-Blue`[1], 
           # Deletion = ColColor$`Low-Blue`[10],
           Splice_Site = ColColor$`Single-Brown`[9], Missense_Mutation = ColColor$`Low-BlueGreen`[9])
HM_layout = list(
  background = function(x, y, w, h) grid.rect(x, y, w, h, 
                                              gp = gpar(fill = "white", lwd = 0.1)),
  Amplification = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                 gp = gpar(fill = HM_Col["Amplification"], col = NA)),
  # Gain = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
  #                                       gp = gpar(fill = HM_Col["Gain"], col = NA)),
  # ShallowDeletion = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
  #                                                  gp = gpar(fill = HM_Col["ShallowDeletion"], col = NA)),
  # Deletion = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
  #                                           gp = gpar(fill = HM_Col["Deletion"], col = NA)),
  Splice_Site = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.7, 
                                               gp = gpar(fill = HM_Col["Splice_Site"], col = NA)),
  Missense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, 
                                                     gp = gpar(fill = HM_Col["Missense_Mutation"], col = NA))
)

##### Vis
oncoPrint(InputOp2, alter_fun = HM_layout, col = HM_Col, show_column_names = F,
          top_annotation = HM_colAnno, # column_split = ERBB$Cluster,
          left_annotation =  rowAnnotation(
            rbar = anno_oncoprint_barplot(
              axis_param = list(direction = "reverse")
            )), right_annotation = NULL,
          height = unit(4, 'cm'), width = unit(15,'cm')) 
#
rm(HM_Col, HM_colAnno, HM_layout, InputOp, InputOp2, Input, ERBB, GBC_Clusters)

##### Figure 2B #####
Input <- filter(IndexERBB, Tumor == 'Yes') %>%
  select(CNV_ERBB2, 
         mRNA_EGFR, Protein_EGFR, ProPhos_EGFR,
         mRNA_ERBB2, Protein_ERBB2, ProPhos_ERBB2,
         mRNA_ERBB3, Protein_ERBB3, ProPhos_ERBB3,
         mRNA_ERBB4, Protein_ERBB4, ProPhos_ERBB4)
Input = reshape2::melt(Input, id.var = 'CNV_ERBB2', variable.name = 'Omic', value.name = 'Abundance') %>%
  as_tibble()
Input$Omics = str_split_fixed(Input$Omic, '_', 2)[,1] %>% factor(levels = c('mRNA',"Protein",'ProPhos'), ordered = T)
Input$Gene = str_split_fixed(Input$Omic, '_', 2)[,2] %>% factor(levels = c('EGFR','ERBB2','ERBB3','ERBB4'), ordered = T)
Input = Input[,-2] %>% na.omit()

# Vis
ggboxplot(Input, 'Omics', 'Abundance', color = 'CNV_ERBB2') + 
  facet_wrap(~Gene) +
  scale_color_manual(values = ColCNA) +
  stat_pvalue_manual(
    compare_means(
      Abundance ~ CNV_ERBB2, data = mutate(Input, CNV_ERBB2 = as.character(CNV_ERBB2)), group.by = c('Omics','Gene'),
      method = "t.test", ref.group = "Amplification"
    ), x = "Omics", y.position = 11,
    label = "p.signif",
    position = position_dodge(0.8)
  ) +
  xlab('') + 
  theme(legend.position = 'right')
#
rm(Input)

##### Figure 2C & Figure S2A #####
##### 1.1 Index prep
Input <- tibble(Patient = c(rep(GBC_Main_Index$Sample, 2), GBC_Benign_Index$Sample),
                Sample = c(paste0(GBC_Main_Index$Sample,'_T'), paste0(GBC_Main_Index$Sample, '_P'), GBC_Benign_Index$Sample),
                SampleWES = c(GBC_Main_Index$WES_T, GBC_Main_Index$WES_P, rep(NA, 18)),
                SampleRNA = c(GBC_Main_Index$RNA_T, GBC_Main_Index$RNA_P, rep(NA, 18)),
                SamplePro = c(GBC_Main_Index$ProPhos_T, GBC_Main_Index$ProPhos_P, GBC_Benign_Index$Pro),
                SamplePhos = c(GBC_Main_Index$ProPhos_T, GBC_Main_Index$ProPhos_P, GBC_Benign_Index$Phos),
                Type = c(rep('Cancer', 195), rep('Adjacent', 195), GBC_Benign_Index$Type) %>%
                  factor(levels = c('Normal', 'Adjacent', 'Benign', 'Cancer'), ordered = T))
# add cluster
Input$Cluster <- NA

# add position
Input$Position = c(GBC_Main_Clinical$Tumor_location, GBC_Main_Clinical$Tumor_location, rep(NA, 18))
Input$Position = case_when(grepl('1',Input$Position) ~ 'Fundus',
                           grepl('2',Input$Position) ~ 'Body',
                           grepl('3',Input$Position) ~ 'Neck',
                           !is.na(Input$Position) ~ 'Others')

# add regional invasion
Input$LM <- c(GBC_Main_Clinical$Regional_invasion, GBC_Main_Clinical$Regional_invasion, rep(NA, 18)) %>%
  factor(levels = c(0, 1), ordered = T)

# add stage
Input$Stage <- c(GBC_Main_Clinical$Tumor_stage2, GBC_Main_Clinical$Tumor_stage2, rep(NA, 18)) %>%
  factor(levels = c(1,2,3,4), ordered = T)

# add ERBB2/3 mutation
Input$`Mut_ERBB2` <- case_when((Input$Patient %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB2')$Tumor_Sample_Barcode) ~ 'Mut',
                               (Input$Patient %in% GBC_Main_Mutation$Tumor_Sample_Barcode) ~ 'WT')
Input$`Mut_ERBB3` <- case_when((Input$Patient %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB3')$Tumor_Sample_Barcode) ~ 'Mut',
                               (Input$Patient %in% GBC_Main_Mutation$Tumor_Sample_Barcode) ~ 'WT')

## add ERBB2/3 Amp
# separated
Input$`CNV_ERBB2` <- sapply(Input$SampleWES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_GBC_Main_GBC_Main_CNV_cate['ERBB2',x]), NA)})
Input$`CNV_ERBB2` <- SeparatedCNV(Input$`CNV_ERBB2`)
Input$`CNV_ERBB3` <- sapply(Input$SampleWES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_GBC_Main_GBC_Main_CNV_cate['ERBB3',x]), NA)})
Input$`CNV_ERBB3` <- SeparatedCNV(Input$`CNV_ERBB3`)
# continious
Input$`CNV_ERBB2_value` <- sapply(Input$SampleWES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV['ERBB2',x]), NA)})

# add ERBBs RNA
Input$`RNA_EGFR` <- sapply(Input$SampleRNA, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['EGFR',x]), NA)})
Input$`RNA_ERBB2` <- sapply(Input$SampleRNA, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['ERBB2',x]), NA)})
Input$`RNA_ERBB3` <- sapply(Input$SampleRNA, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['ERBB3',x]), NA)})
Input$`RNA_ERBB4` <- sapply(Input$SampleRNA, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['ERBB4',x]), NA)})

# add ERBBs Protein
Input$`Pro_EGFR` <- c(sapply(Input$SamplePro[1:390], function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['EGFR',x]), NA)}),
                      sapply(Input$SamplePro[391:408], function(x){ifelse(!is.na(x), as.numeric(GBC_Benign_Pro['EGFR',x]), NA)}))
Input$`Pro_ERBB2` <- c(sapply(Input$SamplePro[1:390], function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['ERBB2',x]), NA)}),
                       sapply(Input$SamplePro[391:408], function(x){ifelse(!is.na(x), as.numeric(GBC_Benign_Pro['ERBB2',x]), NA)}))

# add ERBB2/3 Phos
Input$`Phos_EGFR` <- c(sapply(Input$SamplePhos[1:390], function(x){ifelse(!is.na(x), as.numeric(GBC_Main_ProLevelPhos_knn['EGFR',x]), NA)}),
                       sapply(Input$SamplePhos[391:408], function(x){ifelse(!is.na(x), as.numeric(GBC_Benign_ProLevelPhos_knn['EGFR',x]), NA)}))
Input$`Phos_ERBB2` <- c(sapply(Input$SamplePhos[1:390], function(x){ifelse(!is.na(x), as.numeric(GBC_Main_ProLevelPhos_knn['ERBB2',x]), NA)}),
                        sapply(Input$SamplePhos[391:408], function(x){ifelse(!is.na(x), as.numeric(GBC_Benign_ProLevelPhos_knn['ERBB2',x]), NA)}))
Input$`Phos_ERBB3` <- c(sapply(Input$SamplePhos[1:390], function(x){ifelse(!is.na(x), as.numeric(GBC_Main_ProLevelPhos_knn['ERBB3',x]), NA)}),
                        sapply(Input$SamplePhos[391:408], function(x){ifelse(!is.na(x), as.numeric(GBC_Benign_ProLevelPhos_knn['ERBB3',x]), NA)}))

## add ERBB pathway
# Pro level
Input$`Pathway_ERBB_Pro` <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Pro['KEGG ErbB signaling pathway',x]), NA)})
Input$`Pathway_ERBB2_Pro` <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Pro['REACTOME Signaling by ERBB2',x]), NA)})
Input$`Pathway_ERBB4_Pro` <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Pro['REACTOME Signaling by ERBB4',x]), NA)})
# phos level
Input$`Pathway_ERBB_Phos` <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Phos['KEGG ErbB signaling pathway',x]), NA)})
Input$`Pathway_ERBB2_Phos` <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Phos['REACTOME Signaling by ERBB2',x]), NA)})
Input$`Pathway_ERBB4_Phos` <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Phos['REACTOME Signaling by ERBB4',x]), NA)})

# add MYC amp
Input$`CNV_MYC` <- sapply(Input$SampleWES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_GBC_Main_GBC_Main_CNV_cate['MYC',x]), NA)})
Input$`CNV_MYC` <- SeparatedCNV(Input$`CNV_MYC`)
Input$`CNV_MYC_value` <- sapply(Input$SampleWES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV['MYC',x]), NA)})

# add MYC RNA
Input$`RNA_MYC` <- sapply(Input$SampleRNA, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['MYC',x]), NA)})

# add MYC Phos
Input$`MYC_S62` <- c(sapply(Input$SamplePhos[1:390], function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Phos_NA['MYC:S62',x]), NA)}),
                     sapply(Input$SamplePhos[391:408], function(x){ifelse(!is.na(x), as.numeric(GBC_Benign_Phos_NA['MYC:S62',x]), NA)}))

# add MYC pathway
Input$`Pathway_MYC_down1` <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Pro['HALLMARK_MYC_TARGETS_V1',x]), NA)})
Input$`Pathway_MYC_down2` <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Pro['HALLMARK_MYC_TARGETS_V2',x]), NA)})

##### 1.2 Vis prep
Input2 <- Input[,c(1, 7, 14, 26:31)]
# add final marker
Input2$Ident <- case_when(Input2$CNV_ERBB2 == 'Amplification' ~ 'ERBB2_Amp',
                          !is.na(Input2$CNV_ERBB2) ~ 'ERBB2_others',
                          as.character(Input2$Type) != 'Cancer' ~ as.character(Input2$Type)) %>%
  factor(levels = c('Normal', 'Adjacent', 'Benign', 'ERBB2_others', 'ERBB2_Amp'), ordered = T)

# add EGFR
Input2$Pathway_EGFR_Pro <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Pro['REACTOME Signaling by EGFR',x]), NA)})
Input2$Pathway_EGFR_Phos <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Phos['REACTOME Signaling by EGFR',x]), NA)})
# add Ras
Input2$Pathway_PI3K_Pro <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Pro['KEGG PI3K-Akt signaling pathway',x]), NA)})
Input2$Pathway_PI3K_Phos <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Phos['KEGG PI3K-Akt signaling pathway',x]), NA)})
# add MAPK
Input2$Pathway_MAPK_Pro <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Pro['KEGG MAPK signaling pathway',x]), NA)})
Input2$Pathway_MAPK_Phos <- sapply(Input$SamplePro, function(x){ifelse(!is.na(x), as.numeric(ssGSEAOutput_Phos['KEGG MAPK signaling pathway',x]), NA)})
#####  data for Table S3
write.xlsx(Input2, file = 'ERBB2 CNV Pathway.xlsx', rowNames = T, colNames = T, overwrite = T)

##### transform
Input2 <- select(Input2, -1, -2, -3)
Input2 = Input2 %>% melt(id.vars = 'Ident', value.name = 'Activity', variable.name = 'Pathway') %>%
  mutate(Omics = ifelse(grepl('Phos',Pathway), 'Phosphoproteomics', 'Proteomics'),
         Signal = case_when(grepl('ERBB2', Pathway) ~ 'ERBB2 Signaling',
                            grepl('ERBB4', Pathway) ~ 'ERBB4 Signaling',
                            grepl('EGFR', Pathway) ~ 'EGFR Signaling',
                            grepl('ERBB', Pathway) ~ 'pan-ERBBs Signaling',
                            grepl('PI3K', Pathway) ~ 'PI3K/AKT Signaling',
                            grepl('MAPK', Pathway) ~ 'MAPK/ERK Signaling') %>% 
           factor(levels = c('pan-ERBBs Signaling', 'EGFR Signaling', 'ERBB2 Signaling', 'ERBB4 Signaling',
                             'PI3K/AKT Signaling', 'MAPK/ERK Signaling'), ordered = T)) %>% na.omit() %>% as_tibble()
##### 1.3 vis
ggboxplot(Input2, 'Ident', 'Activity', color = 'Ident', lwd = 0.75) +
  scale_color_manual(values = ColType2) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black')) +
  facet_wrap(~ Omics + Signal, scales = 'free_y', nrow = 2) +
  stat_pvalue_manual(compare_means(Activity ~ Ident, data = mutate(Input2, Ident = as.character(Ident)), ref.group = 'ERBB2_Amp',
                                   group.by = c("Omics", 'Signal'), method = "t.test") %>%
                       mutate(y.position = rep(c(0.3,0.35,0.4,0.45), 11)),
                     label = "p.signif") +
  xlab('') + ylab('Pathway Activity')

##### Figure S2B #####
##### 1.1 external RNA-Seq data loading
# GSE202479
DataExt1 <- fread('../10.1371:journal.pone.0283770.txt')
DataExt1 <- select(DataExt1, 2, 11:30) %>% distinct(gene_name, .keep_all = T) %>%
  column_to_rownames(var = 'gene_name') %>% as.data.frame()
colnames(DataExt1) = str_split_fixed(colnames(DataExt1), '_', 2)[,1]
# indexing
DataExt1Index <- tibble(Sample = colnames(DataExt1),
                        Type = case_when(colnames(DataExt1) %in% c('Y8','Y12','Y13','Y16') ~ 'Stone',
                                         colnames(DataExt1) %in% c('T11','T24','T30') ~ 'Adenoma',
                                         colnames(DataExt1) %in% c('T5','T12','T13','T18','T31') ~ 'Early Ca',
                                         colnames(DataExt1) %in% c('T1','T19','T22','T27','T32') ~ 'Advanced Ca',
                                         grepl('N',colnames(DataExt1)) ~ 'Normal') %>%
                          factor(levels = c('Normal', 'Stone','Adenoma','Early Ca','Advanced Ca'), ordered = T))

##### 1.2 ssgsea
PathwayNeeded <- c('KEGG ErbB signaling pathway','REACTOME Signaling by EGFR',
                   'REACTOME Signaling by ERBB2','REACTOME Signaling by ERBB4',
                   'KEGG PI3K-Akt signaling pathway','KEGG Ras signaling pathway')
ssGSEA1 <- GSVA::gsva(as.matrix(DataExt1), TotalPathwayGSVA[,PathwayNeeded], method = 'ssgsea') %>% t() %>% as.data.frame()
DataExt1Index = cbind(DataExt1Index, ssGSEA1)
# data for Table S3
write.xlsx(DataExt1Index, file = 'ExternalCohort paradox.xlsx', rowNames = T, colNames = T, overwrite = T)

##### 1.3 vis prep
Input1 <- melt(DataExt1Index[,c(2,3:8)], id.vars = 'Type', value.name = 'Activity', variable.name = 'Pathway') %>%
  mutate(Signal = case_when(grepl('ERBB2', Pathway) ~ 'ERBB2 Signaling',
                            grepl('ERBB4', Pathway) ~ 'ERBB4 Signaling',
                            grepl('EGFR', Pathway) ~ 'EGFR Signaling',
                            grepl('ERBB', Pathway, ignore.case = T) ~ 'pan-ERBBs Signaling',
                            grepl('PI3K', Pathway) ~ 'PI3K/AKT Signaling',
                            grepl('RAS', Pathway, ignore.case = T) ~ 'RAS/MAPK/ERK Signaling') %>% 
           factor(levels = c('pan-ERBBs Signaling', 'EGFR Signaling', 'ERBB2 Signaling', 'ERBB4 Signaling',
                             'PI3K/AKT Signaling', 'RAS/MAPK/ERK Signaling'), ordered = T))

##### 1.4 vis
ColType = c(ColJournal$Nature[3],ColJournal$Nature[7],
            ColJournal$Nature[5],ColJournal$Nature[1],
            ColJournal$Nature[8],ColJournal$Nature[8])
ggboxplot(Input1, 'Type', 'Activity', color = 'Type', add = 'jitter', lwd = 0.75) +
  scale_color_manual(values = ColType) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black')) +
  facet_wrap(~Signal, scales = 'free_y', nrow = 2) +
  stat_pvalue_manual(compare_means(Activity ~ Type, data = Input1 %>% mutate(Type = as.character(Type)),
                                   ref.group = 'Normal', group.by = c('Signal'), method = "t.test") %>%
                       mutate(y.position = rep(c(2.8,2.85,2.9,2.95), 6)),
                     label = "p.signif") 
#
rm(DataExt1, DataExt1Index)

##### Figure S2C #####
##### 1.1 index
IndexTemp <- filter(IndexERBB, Tumor == 'Yes', !is.na(Pro))

##### 1.2 calc ssgsea-gsva
InputGSVAraw <- GSVA::gsva(as.matrix(GBC_Main_ProLevelPhos_knn[,IndexTemp$Pro]), TotalPathwayGSVA, method = 'gsva')

##### 1.3 vis prep pathway
PathwayPanelErbB <- c('ErbBs' = 'KEGG ErbB signaling pathway', 'EGFR' = 'REACTOME Signaling by EGFR',
                      'ERBB2' = 'REACTOME Signaling by ERBB2', 'ERBB4' = 'REACTOME Signaling by ERBB4',
                      'Ras' = 'KEGG Ras signaling pathway', 'PI3K' = 'KEGG PI3K-Akt signaling pathway')
PathwaySelected <- c('REACTOME ERBB2 Activates PTK6 Signaling','REACTOME GRB7 events in ERBB2 signaling',
                     'REACTOME PI3K/AKT activation','REACTOME RET signaling','KEGG FoxO signaling pathway',
                     PathwayPanelErbB, 'REACTOME Signalling to ERKs',
                     'HALLMARK_DNA_REPAIR','HALLMARK_G2M_CHECKPOINT','KEGG Cell cycle','HALLMARK_INTERFERON_GAMMA_RESPONSE',
                     'HALLMARK_MTORC1_SIGNALING','HALLMARK_MYC_TARGETS_V1','HALLMARK_WNT_BETA_CATENIN_SIGNALING',
                     'HALLMARK_PI3K_AKT_MTOR_SIGNALING')
InputHMpathway <- as.data.frame(InputGSVAraw)[PathwaySelected, IndexTemp$Pro]
Inputinfo = arrange(IndexTemp, Mut_ERBBPathway, Mut_ERBB2, Mut_ERBB3)
#
Inputinfo$`Genome ERBB2 status` <- case_when(Inputinfo$`CNV_ERBB2` == 'Amplification' ~ 'Amp',
                                             !is.na(Inputinfo$`Mut_ERBB2`) ~ 'Mut',
                                             Inputinfo$`Mut_ERBB2` == 'WT' ~ 'WT')
Inputinfo$`Genome ERBB3 status` <- case_when(Inputinfo$`CNV_ERBB3` == 'Amplification' ~ 'Amp',
                                             !is.na(Inputinfo$`Mut_ERBB3`) ~ 'Mut',
                                             Inputinfo$`Mut_ERBB3` == 'WT' ~ 'WT')
Inputinfo = arrange(Inputinfo, Mut_ERBBPathway, `Genome ERBB2 status`, `Genome ERBB3 status`)


##### 1.4 Vis
TempPathology = c(sapply(Inputinfo$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Pathological_type}))
AnnoSample = columnAnnotation(ERBB2 = c(Inputinfo$`Genome ERBB2 status`),
                              ERBB3 = c(Inputinfo$`Genome ERBB3 status`),
                              Pro_ERBB2 = sapply(Inputinfo$Pro, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['ERBB2',x]), NA)}),
                              Pathology = case_when(TempPathology == 'adenocarcinoma' ~ 'adenocarcinoma',
                                                    TempPathology == 'adenosquanmous carcinoma' ~ 'adenosquanmous carcinoma',
                                                    TempPathology == 'neuroendocrine carcinoma' ~ 'neuroendocrine carcinoma',
                                                    !is.na(TempPathology) ~ 'Others'),
                              col = list(
                                ERBB2 = c('WT' = 'green', 'Mut' = "orange", 'Amp' = 'red'),
                                ERBB3 = c('WT' = 'green', 'Mut' = "orange", 'Amp' = 'red'),
                                Pro_ERBB2 = colorRamp2(c(-1,0,1), c(ColColor$`Low-Indigo`[5],'white',ColColor$`High-Orange`[5])),
                                Pathology = c('adenocarcinoma' = wes_palette('GrandBudapest2', 4, type = c("discrete"))[1], 
                                              'adenosquanmous carcinoma' = wes_palette('GrandBudapest2', 4, type = c("discrete"))[2],
                                              'neuroendocrine carcinoma' = wes_palette('GrandBudapest2', 4, type = c("discrete"))[3],
                                              'Others' = wes_palette('GrandBudapest2', 4, type = c("discrete"))[4])
                              ),
                              gp = gpar(col = "white"),
                              simple_anno_size = unit(0.5, "cm"))
Heatmap(t(scale(t(InputHMpathway[,c(Inputinfo$Pro)]))), name = 'Scaled features', show_column_names = F,
        column_order = c(Inputinfo$Pro), 
        rect_gp = gpar(col = "white", lwd = 0.5),
        col = colorRamp2(c(-1.5,0,1.5), c(ColColor$`Low-Blue`[5],'white',ColColor$`High-Red`[5])),
        top_annotation = AnnoSample,
        height = unit(20/2.5,'cm'), width = unit(18,'cm'))
# data for Table S3
write.xlsx(InputHMpathway[,c(Inputinfo$Pro)], file = 'Heatmap value.xlsx',
           rowNames = T, colNames = T, overwrite = T)
write.xlsx(Inputinfo, file = 'Heatmap anno.xlsx',
           rowNames = T, colNames = T, overwrite = T)
#
rm(Input1, Input2, InputGSVAraw, InputHMpathway, Inputinfo, IndexTemp,
   AnnoSample, PathwayPanelErbB, PathwayNeeded, PathwaySelected, TempPathology)

##### Figure S2D #####
##### 1.0 known site vis
lollipopPlot(maf = GBC_Main_Maf, gene = 'ERBB2', showDomainLabel = F, 
             AACol = 'AAChange.refGene', showMutationRate = FALSE)    
lollipopPlot(maf = GBC_Main_Maf, gene = 'ERBB3', showDomainLabel = F,
             AACol = 'AAChange.refGene', showMutationRate = FALSE)

##### 1.1 Data loading
# TCGA Pan-cancer cohort (XENA)
load('../TCGA_Pan-cancer.RData')
# ClinVar
load('../ClinVar.RData')
# COSMIC
load('../COSMIC_Mut.RData')

##### 1.2 Vis Prep
library(trackViewer)
library(circlize)
##### Mutation All
##### Protein domain
LolDomain = readRDS(file = system.file('extdata', 'protein_domains.RDs', package = 'maftools'))
LolDomainERBB2 <- LolDomain[HGNC %in% "ERBB2"][refseq.ID == 'NM_004448',]
LolDomainERBB3 <- LolDomain[HGNC %in% "ERBB3"][refseq.ID == 'NM_001982',]
##### Protein Seq
LolSeqERBB2 = c('MELAALCRWGLLLALLPPGAASTQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSFLQDIQEVQGYVLIAHNQVRQVPLQRLRIVRGTQLFEDNYALAVLDNGDPLNNTTPVTGASPGGLRELQLRSLTEILKGGVLIQRNPQLCYQDTILWKDIFHKNNQLALTLIDTNRSRACHPCSPMCKGSRCWGESSEDCQSLTRTVCAGGCARCKGPLPTDCCHEQCAAGCTGPKHSDCLACLHFNHSGICELHCPALVTYNTDTFESMPNPEGRYTFGASCVTACPYNYLSTDVGSCTLVCPLHNQEVTAEDGTQRCEKCSKPCARVCYGLGMEHLREVRAVTSANIQEFAGCKKIFGSLAFLPESFDGDPASNTAPLQPEQLQVFETLEEITGYLYISAWPDSLPDLSVFQNLQVIRGRILHNGAYSLTLQGLGISWLGLRSLRELGSGLALIHHNTHLCFVHTVPWDQLFRNPHQALLHTANRPEDECVGEGLACHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQNGSVTCFGPEADQCVACAHYKDPPFCVARCPSGVKPDLSYMPIWKFPDEEGACQPCPINCTHSCVDLDDKGCPAEQRASPLTSIISAVVGILLVVVLGVVFGILIKRRQQKIRKYTMRRLLQETELVEPLTPSGAMPNQAQMRILKETELRKVKVLGSGAFGTVYKGIWIPDGENVKIPVAIKVLRENTSPKANKEILDEAYVMAGVGSPYVSRLLGICLTSTVQLVTQLMPYGCLLDHVRENRGRLGSQDLLNWCMQIAKGMSYLEDVRLVHRDLAARNVLVKSPNHVKITDFGLARLLDIDETEYHADGGKVPIKWMALESILRRRFTHQSDVWSYGVTVWELMTFGAKPYDGIPAREIPDLLEKGERLPQPPICTIDVYMIMVKCWMIDSECRPRFRELVSEFSRMARDPQRFVVIQNEDLGPASPLDSTFYRSLLEDDDMGDLVDAEEYLVPQQGFFCPDPAPGAGGMVHHRHRSSSTRSGGGDLTLGLEPSEEEAPRSPLAPSEGAGSDVFDGDLGMGAAKGLQSLPTHDPSPLQRYSEDPTVPLPSETDGYVAPLTCSPQPEYVNQPDVRPQPPSPREGPLPAARPAGATLERPKTLSPGKNGVVKDVFAFGGAVENPEYLTPQGGAAPQPHPPPAFSPAFDNLYYWDQDPPERGAPPSTFKGTPTAENPEYLGLDVPV')
LolSeqERBB3 = c('MRANDALQVLGLLFSLARGSEVGNSQAVCPGTLNGLSVTGDAENQYQTLYKLYERCEVVMGNLEIVLTGHNADLSFLQWIREVTGYVLVAMNEFSTLPLPNLRVVRGTQVYDGKFAIFVMLNYNTNSSHALRQLRLTQLTEILSGGVYIEKNDKLCHMDTIDWRDIVRDRDAEIVVKDNGRSCPPCHEVCKGRCWGPGSEDCQTLTKTICAPQCNGHCFGPNPNQCCHDECAGGCSGPQDTDCFACRHFNDSGACVPRCPQPLVYNKLTFQLEPNPHTKYQYGGVCVASCPHNFVVDQTSCVRACPPDKMEVDKNGLKMCEPCGGLCPKACEGTGSGSRFQTVDSSNIDGFVNCTKILGNLDFLITGLNGDPWHKIPALDPEKLNVFRTVREITGYLNIQSWPPHMHNFSVFSNLTTIGGRSLYNRGFSLLIMKNLNVTSLGFRSLKEISAGRIYISANRQLCYHHSLNWTKVLRGPTEERLDIKHNRPRRDCVAEGKVCDPLCSSGGCWGPGPGQCLSCRNYSRGGVCVTHCNFLNGEPREFAHEAECFSCHPECQPMEGTATCNGSGSDTCAQCAHFRDGPHCVSSCPHGVLGAKGPIYKYPDVQNECRPCHENCTQGCKGPELQDCLGQTLVLIGKTHLTMALTVIAGLVVIFMMLGGTFLYWRGRRIQNKRAMRRYLERGESIEPLDPSEKANKVLARIFKETELRKLKVLGSGVFGTVHKGVWIPEGESIKIPVCIKVIEDKSGRQSFQAVTDHMLAIGSLDHAHIVRLLGLCPGSSLQLVTQYLPLGSLLDHVRQHRGALGPQLLLNWGVQIAKGMYYLEEHGMVHRNLAARNVLLKSPSQVQVADFGVADLLPPDDKQLLYSEAKTPIKWMALESIHFGKYTHQSDVWSYGVTVWELMTFGAEPYAGLRLAEVPDLLEKGERLAQPQICTIDVYMVMVKCWMIDENIRPTFKELANEFTRMARDPPRYLVIKRESGPGIAPGPEPHGLTNKKLEEVELEPELDLDLDLEAEEDNLATTTLGSALSLPVGTLNRPRGSQSLLSPSSGYMPMNQGNLGESCQESAVSGSSERCPRPVSLHPMPRGCLASESSEGHVTGSEAELQEKVSMCRSRSRSRSPRPRGDSAYHSQRHSLLTPVTPLSPPGLEEEDVNGYVMPDTHLKGTPSSREGTLSSVGLSSVLGTEEEDEDEEYEYMNRRRRHSPPHPPRPSSLEELGYEYMDVGSDLSASLGSTQSCPLHPVPIMPTAGTTPDEDYEYMNRQRDGGGPGGDYAAMGACPASEQGYEEMRAFQGPGHQAPHVHYARLKTLRSLEATDSAFDNPDYWHSRLFPKANAQRT')
##### replace pattern
LolAAReplace <- c("Ala"="A", "Arg"="R", "Asn"="N", "Asp"="D", "Cys"="C",
                  "Gln"="Q", "Glu"="E", "Gly"="G", "His"="H", "Ile"="I", 
                  "Leu"="L", "Lys"="K", "Met"="M", "Phe"="F", "Pro"="P", 
                  "Ser"="S", "Thr"="T", "Trp"="W", "Tyr"="Y", "Val"="V")

##### 1.3 Mutation extraction ERBB2 and Vis
### [ERBB2]
## all mut info
LolMutERBB2 <- rbind(filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB2') %>% 
                       select(Gene = Hugo_Symbol, aaChange = aaChange) %>%
                       mutate(Cohort = 'GBC', Type = 'GBC'),
                     filter(COSMIC_MutSite, GENE_SYMBOL == 'ERBB2') %>% 
                       select(Gene = GENE_SYMBOL, aaChange = MUTATION_AA) %>%
                       mutate(Cohort = 'COSMIC', Type = 'External'),
                     filter(ClinVar_Mut, `#Symbol` == 'ERBB2', ProteinChange != '-') %>% 
                       select(Gene = `#Symbol`, aaChange = ProteinChange) %>% # 后续调整aa格式
                       mutate(Cohort = 'ClinVar', Type = 'External'),
                     filter(Pan_Mutation, gene == 'ERBB2', effect != 'Silent') %>% 
                       select(Gene = gene, aaChange = Amino_Acid_Change) %>%
                       mutate(Cohort = 'TCGA Pan-cancer', Type = 'External')) %>%
  filter(aaChange != '') %>% filter(aaChange != 'p.?')
LolMutERBB2 = as_tibble(LolMutERBB2)
# uniform the aachange 
for (aa in names(LolAAReplace)) {LolMutERBB2$aaChange <- gsub(aa, LolAAReplace[aa], LolMutERBB2$aaChange)}
rm(aa)
### extract position and freq info
# origin
LolMutERBB2$WT = str_split_fixed(LolMutERBB2$aaChange, fixed('.'), 2)[,2]
LolMutERBB2$WT = str_split_fixed(LolMutERBB2$WT, '[0-9]', 2)[,1]
# position
LolMutERBB2$Position = str_split_fixed(LolMutERBB2$aaChange, '[A-Z*]', 2)[,2]
LolMutERBB2$Position = sub("(\\d+).*", "\\1", LolMutERBB2$Position) %>% as.numeric()
filter(LolMutERBB2, is.na(Position))
LolMutERBB2 = na.omit(LolMutERBB2)
# select the correct WT aa by protein sequences
Output = tibble(Gene = character(), aaChange = character(), Cohort = character(),
                Type = character(), WT = character(), Position = numeric())
for (i in 1:length(table(LolMutERBB2$Position))) {
  TempMutPosition = table(LolMutERBB2$Position)[i] %>% names() %>% as.numeric()
  TempCorrectaa = str_sub(LolSeqERBB2, TempMutPosition, TempMutPosition)
  Tempdf = filter(LolMutERBB2, Position == TempMutPosition, WT == TempCorrectaa)
  Output = rbind(Output, Tempdf)
}
LolMutERBB2 = Output
rm(TempMutPosition, TempCorrectaa, Tempdf, i)
# mutant
LolMutERBB2$Mutant = str_split_fixed(LolMutERBB2$aaChange, fixed('.'), 2)[,2]
LolMutERBB2$Mutant = sapply(1:length(LolMutERBB2$Mutant), function(x){
  str_split_fixed(LolMutERBB2$Mutant[x], as.character(LolMutERBB2$Position[x]), 2)[,2]
})
LolMutERBB2$Mutant = ifelse(LolMutERBB2$Mutant == LolMutERBB2$WT, '=', LolMutERBB2$Mutant)
# remove silent mutation
LolMutERBB2 = filter(LolMutERBB2, Mutant != '=')
# change
LolMutERBB2$Change = case_when(grepl('del', LolMutERBB2$Mutant) ~ 'InDel',
                               grepl('ins', LolMutERBB2$Mutant) ~ 'InDel',
                               grepl('dup', LolMutERBB2$Mutant) ~ 'Dup',
                               grepl('splice', LolMutERBB2$Mutant) ~ 'Splice',
                               grepl('fs', LolMutERBB2$Mutant) ~ 'Frame-shift',
                               grepl('Ter', LolMutERBB2$Mutant) ~ 'Terminated',
                               grepl('[A-Z]', LolMutERBB2$Mutant) ~ 'aaChange',
                               grepl(fixed('*'), LolMutERBB2$Mutant) ~ 'Terminated')
table(LolMutERBB2$Change)
# Freq
LolMutERBB2$Frequency = sapply(1:nrow(LolMutERBB2), function(x){
  Tempdf = filter(LolMutERBB2, Position == LolMutERBB2$Position[x], 
                  Mutant == LolMutERBB2$Mutant[x], Type == LolMutERBB2$Type[x])
  TempN = nrow(Tempdf)
  TempN
})
# final layout
LolMutERBB2$LayOut = ifelse(LolMutERBB2$Change %in% c('aaChange'), 
                            paste0(LolMutERBB2$WT, LolMutERBB2$Position, LolMutERBB2$Mutant),
                            paste0(LolMutERBB2$WT, LolMutERBB2$Position, ' ', LolMutERBB2$Change))
# remove duplicated
LolMutERBB2_clean = distinct(LolMutERBB2, Type, LayOut, .keep_all = T)
# for vis
LolMutERBB2_forvis = na.omit(LolMutERBB2_clean) %>% filter(Frequency >= 3)
# combine same position
Output = tibble(Gene = character(), aaChange = character(), Cohort = character(),
                Type = character(), WT = character(), Position = numeric(),
                Mutant = character(), Change = character(), Frequency = numeric(), 
                LayOut = character(), Multi = character())
for (i in 1:length(table(LolMutERBB2_forvis$Position))) {
  TempMutPosition = table(LolMutERBB2_forvis$Position)[i] %>% names() %>% as.numeric()
  Tempdf = filter(LolMutERBB2_forvis, Position == TempMutPosition)
  if (nrow(Tempdf) == 1){
    Tempdf$Multi = 'No'
    Output = rbind(Output, Tempdf)
  } else {
    Tempdf$Multi = 'Yes'; TempAllFreq = sum(Tempdf$Frequency)
    a = distinct(Tempdf, Type, WT, Position, .keep_all = T)
    a$Frequency = TempAllFreq
    Output = rbind(Output, a)
  }
}
rm(TempMutPosition, Tempdf, i, a)
LolMutERBB2_forvis = filter(Output, Frequency >= 5)
# add GBC
LolMutERBB2_forvis = rbind(LolMutERBB2_forvis,
                           filter(LolMutERBB2_clean, Type == 'GBC') %>% mutate(Multi = 'No'))
# add final layout
LolMutERBB2_forvis$Final_Layout = case_when(LolMutERBB2_forvis$Multi == 'Yes' ~ paste0(LolMutERBB2_forvis$WT,LolMutERBB2_forvis$Position,'X'),
                                            LolMutERBB2_forvis$Multi == 'No' ~ LolMutERBB2_forvis$LayOut)

### ERBB2 Vis Prep
## Mut info
LolVisMut = GRanges("ERBB2", IRanges(LolMutERBB2_forvis$Position, width = 1,
                                     names = LolMutERBB2_forvis$Final_Layout))
# basiccol and alpha
LolVisMut$color <- ifelse(LolMutERBB2_forvis$Type == 'External',
                          ColJournal$COSMICsignature[1], ColJournal$COSMICsignature[3])
LolVisMut$border <- ifelse(LolMutERBB2_forvis$Frequency > 5, 'black', 'grey90')
LolVisMut$alpha = case_when(LolMutERBB2_forvis$Type == 'GBC' ~ 1,
                            LolMutERBB2_forvis$Frequency >= 10 ~ 1)
LolVisMut$alpha[is.na(LolVisMut$alpha)] = 0.6
# anno color
LolVisMut$dashline.col <- ifelse(LolMutERBB2_forvis$Frequency > 100, 'red','black')
LolVisMut$dashline.lwd <- ifelse(LolMutERBB2_forvis$Frequency > 100, 4, 0.5)
# label
LolVisMut$label <- case_when(LolMutERBB2_forvis$Type == 'GBC' ~ LolMutERBB2_forvis$Frequency,
                             LolMutERBB2_forvis$Frequency >= 10 ~ LolMutERBB2_forvis$Frequency)
LolVisMut$label[is.na(LolVisMut$label)] = ''
LolVisMut$label = as.character(LolVisMut$label)
LolVisMut$label.col <- 'black'
LolVisMut$label.cex <- case_when(LolMutERBB2_forvis$Type == 'GBC' ~ 1,
                                 LolMutERBB2_forvis$Frequency >= 7 ~ log10(LolMutERBB2_forvis$Frequency) + 0.5)
LolVisMut$label.parameter.rot <- 90
# label col
LolCol <- list(gpar(col="#ec4646"),gpar(col="black"));names(LolCol) <- c("high","low")
LolCol = LolCol[ifelse(LolMutERBB2_forvis$Frequency > 100, 'high','low')]
LolVisMut$label.parameter.gp = LolCol
# height
LolVisMut$score <- ifelse(LolMutERBB2_forvis$Type == 'GBC', log2(LolMutERBB2_forvis$Frequency) + 0.5,
                          log(LolMutERBB2_forvis$Frequency, base = 4) + 0.1) # 点的高度
# circle large
LolVisMut$cex <- case_when(LolMutERBB2_forvis$Type == 'GBC' ~ 1,
                           LolMutERBB2_forvis$Type != 'GBC' ~ log(LolMutERBB2_forvis$Frequency,base = 30) + 0.1)
# up down
LolVisMut$SNPsideID <- ifelse(LolMutERBB2_forvis$Type == 'GBC', 'bottom', 'top')
# Protein info
LolVisProtein = GRanges("ERBB2",IRanges(c(LolDomainERBB2$Start, 1),
                                        width = c(LolDomainERBB2$End - LolDomainERBB2$Start, 1255),
                                        names = c(LolDomainERBB2$Label, 'Full-length')))
LolVisProtein$fill <- case_when(LolVisProtein@ranges@NAMES == 'Full-length' ~ 'black',
                                LolVisProtein@ranges@NAMES == 'Recep_L_domain' ~ ColJournal$Nature[1],
                                LolVisProtein@ranges@NAMES == 'Furin-like' ~ ColJournal$Nature[7],
                                LolVisProtein@ranges@NAMES == 'FU' ~ ColJournal$Nature[3],
                                LolVisProtein@ranges@NAMES == 'TM_ErbB2' ~ ColJournal$Nature[2],
                                LolVisProtein@ranges@NAMES == 'PTKc_HER2' ~ ColJournal$Nature[4],
                                LolVisProtein@ranges@NAMES == 'Pkinase_Tyr' ~ ColJournal$Nature[9])
LolVisProtein$height <- case_when(LolVisProtein@ranges@NAMES == 'Full-length' ~ 0,
                                  LolVisProtein@ranges@NAMES == 'Recep_L_domain' ~ 0.06,
                                  LolVisProtein@ranges@NAMES == 'Furin-like' ~ 0.04,
                                  LolVisProtein@ranges@NAMES == 'FU' ~ 0.06,
                                  LolVisProtein@ranges@NAMES == 'TM_ErbB2' ~ 0.06,
                                  LolVisProtein@ranges@NAMES == 'PTKc_HER2' ~ 0.04,
                                  LolVisProtein@ranges@NAMES == 'Pkinase_Tyr' ~ 0.06)
### ERBB Vis
lolliplot(LolVisMut, LolVisProtein, jitter="node", type="circle",
          xaxis = c(0, 200, 400, 600, 800, 1000, 1200, 1255))
# data for Table S3
write.xlsx(LolMutERBB2_forvis, file = 'lol ERBB2.xlsx', rowNames = T, colNames = T, overwrite = T)
#
rm(LolMutERBB2_forvis, LolMutERBB2, LolMutERBB2_clean, Output, LolCol)

##### 1.4 Mutation extraction ERBB3 and Vis
### [ERBB3]
## all mut info
LolMutERBB3 <- rbind(filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB3') %>% 
                       select(Gene = Hugo_Symbol, aaChange = aaChange) %>%
                       mutate(Cohort = 'GBC', Type = 'GBC'),
                     filter(COSMIC_MutSite, GENE_SYMBOL == 'ERBB3') %>% 
                       select(Gene = GENE_SYMBOL, aaChange = MUTATION_AA) %>%
                       mutate(Cohort = 'COSMIC', Type = 'External'),
                     filter(ClinVar_Mut, `#Symbol` == 'ERBB3', ProteinChange != '-') %>% 
                       select(Gene = `#Symbol`, aaChange = ProteinChange) %>% # 后续调整aa格式
                       mutate(Cohort = 'ClinVar', Type = 'External'),
                     filter(Pan_Mutation, gene == 'ERBB3', effect != 'Silent') %>% 
                       select(Gene = gene, aaChange = Amino_Acid_Change) %>%
                       mutate(Cohort = 'TCGA Pan-cancer', Type = 'External')) %>%
  filter(aaChange != '') %>% filter(aaChange != 'p.?')
LolMutERBB3 = as_tibble(LolMutERBB3)
# uniform the aachange 
for (aa in names(LolAAReplace)) {LolMutERBB3$aaChange <- gsub(aa, LolAAReplace[aa], LolMutERBB3$aaChange)}
rm(aa)
### extract position and freq info
# origin
LolMutERBB3$WT = str_split_fixed(LolMutERBB3$aaChange, fixed('.'), 2)[,2]
LolMutERBB3$WT = str_split_fixed(LolMutERBB3$WT, '[0-9]', 2)[,1]
# position
LolMutERBB3$Position = str_split_fixed(LolMutERBB3$aaChange, '[A-Z*]', 2)[,2]
LolMutERBB3$Position = sub("(\\d+).*", "\\1", LolMutERBB3$Position) %>% as.numeric()
filter(LolMutERBB3, is.na(Position))
LolMutERBB3 = na.omit(LolMutERBB3)
# select the correct WT aa by protein sequences
Output = tibble(Gene = character(), aaChange = character(), Cohort = character(),
                Type = character(), WT = character(), Position = numeric())
for (i in 1:length(table(LolMutERBB3$Position))) {
  TempMutPosition = table(LolMutERBB3$Position)[i] %>% names() %>% as.numeric()
  TempCorrectaa = str_sub(LolSeqERBB3, TempMutPosition, TempMutPosition)
  Tempdf = filter(LolMutERBB3, Position == TempMutPosition, WT == TempCorrectaa)
  Output = rbind(Output, Tempdf)
}
LolMutERBB3 = Output
rm(TempMutPosition, TempCorrectaa, Tempdf, i)
# mutant
LolMutERBB3$Mutant = str_split_fixed(LolMutERBB3$aaChange, fixed('.'), 2)[,2]
LolMutERBB3$Mutant = sapply(1:length(LolMutERBB3$Mutant), function(x){
  str_split_fixed(LolMutERBB3$Mutant[x], as.character(LolMutERBB3$Position[x]), 2)[,2]
})
LolMutERBB3$Mutant = ifelse(LolMutERBB3$Mutant == LolMutERBB3$WT, '=', LolMutERBB3$Mutant)
LolMutERBB3 = filter(LolMutERBB3, Mutant != '?')
# remove silent mutation
LolMutERBB3 = filter(LolMutERBB3, Mutant != '=')
# change
LolMutERBB3$Change = case_when(grepl('del', LolMutERBB3$Mutant) ~ 'InDel',
                               grepl('ins', LolMutERBB3$Mutant) ~ 'InDel',
                               grepl('dup', LolMutERBB3$Mutant) ~ 'Dup',
                               grepl('splice', LolMutERBB3$Mutant) ~ 'Splice',
                               grepl('fs', LolMutERBB3$Mutant) ~ 'Frame-shift',
                               grepl('Ter', LolMutERBB3$Mutant) ~ 'Terminated',
                               grepl('[A-Z]', LolMutERBB3$Mutant) ~ 'aaChange',
                               grepl(fixed('*'), LolMutERBB3$Mutant) ~ 'Terminated')
table(LolMutERBB3$Change)
# Freq
LolMutERBB3$Frequency = sapply(1:nrow(LolMutERBB3), function(x){
  Tempdf = filter(LolMutERBB3, Position == LolMutERBB3$Position[x], 
                  Mutant == LolMutERBB3$Mutant[x], Type == LolMutERBB3$Type[x])
  TempN = nrow(Tempdf)
  TempN
})
# final layout
LolMutERBB3$LayOut = ifelse(LolMutERBB3$Change %in% c('aaChange'), 
                            paste0(LolMutERBB3$WT, LolMutERBB3$Position, LolMutERBB3$Mutant),
                            paste0(LolMutERBB3$WT, LolMutERBB3$Position, ' ', LolMutERBB3$Change))
# remove duplicated
LolMutERBB3_clean = distinct(LolMutERBB3, Type, LayOut, .keep_all = T)
# for vis
LolMutERBB3_forvis = na.omit(LolMutERBB3_clean) %>% filter(Frequency >= 3, Type != 'GBC')
# combine same position
Output = tibble(Gene = character(), aaChange = character(), Cohort = character(),
                Type = character(), WT = character(), Position = numeric(),
                Mutant = character(), Change = character(), Frequency = numeric(), 
                LayOut = character(), Multi = character())
for (i in 1:length(table(LolMutERBB3_forvis$Position))) {
  TempMutPosition = table(LolMutERBB3_forvis$Position)[i] %>% names() %>% as.numeric()
  Tempdf = filter(LolMutERBB3_forvis, Position == TempMutPosition)
  if (nrow(Tempdf) == 1){
    Tempdf$Multi = 'No'
    Output = rbind(Output, Tempdf)
  } else {
    Tempdf$Multi = 'Yes'; TempAllFreq = sum(Tempdf$Frequency)
    a = distinct(Tempdf, Type, WT, Position, .keep_all = T)
    a$Frequency = TempAllFreq
    Output = rbind(Output, a)
  }
}
rm(TempMutPosition, Tempdf, i, a)
LolMutERBB3_forvis = filter(Output, Frequency >= 5)
# add GBC
LolMutERBB3_forvis = rbind(LolMutERBB3_forvis,
                           filter(LolMutERBB3_clean, Type == 'GBC') %>% mutate(Multi = 'No'))
# add final layout
LolMutERBB3_forvis$Final_Layout = case_when(LolMutERBB3_forvis$Multi == 'Yes' ~ paste0(LolMutERBB3_forvis$WT,LolMutERBB3_forvis$Position,'X'),
                                            LolMutERBB3_forvis$Multi == 'No' ~ LolMutERBB3_forvis$LayOut)

### ERBB3 Vis Prep
library(trackViewer)
## Mut info
LolVisMut = GRanges("ERBB3", IRanges(LolMutERBB3_forvis$Position, width = 1,
                                     names = LolMutERBB3_forvis$Final_Layout))
# basiccol and alpha
LolVisMut$color <- ifelse(LolMutERBB3_forvis$Type == 'External',
                          ColJournal$COSMICsignature[1], ColJournal$COSMICsignature[3])
LolVisMut$border <- ifelse(LolMutERBB3_forvis$Frequency > 10, 'black', 'grey90')
LolVisMut$alpha = case_when(LolMutERBB3_forvis$Type == 'GBC' ~ 1,
                            LolMutERBB3_forvis$Frequency >= 10 ~ 1)
LolVisMut$alpha[is.na(LolVisMut$alpha)] = 0.6
# anno color
LolVisMut$dashline.col <- ifelse(LolMutERBB3_forvis$Frequency > 50, 'red','black')
LolVisMut$dashline.lwd <- ifelse(LolMutERBB3_forvis$Frequency > 50, 4, 0.5)
# label
LolVisMut$label <- case_when(LolMutERBB3_forvis$Type == 'GBC' ~ LolMutERBB3_forvis$Frequency,
                             LolMutERBB3_forvis$Frequency >= 10 ~ LolMutERBB3_forvis$Frequency)
LolVisMut$label[is.na(LolVisMut$label)] = ''
LolVisMut$label = as.character(LolVisMut$label)
LolVisMut$label.col <- 'black'
LolVisMut$label.cex <- case_when(LolMutERBB3_forvis$Type == 'GBC' ~ 1,
                                 LolMutERBB3_forvis$Frequency >= 10 ~ log10(LolMutERBB3_forvis$Frequency) + 0.5)
LolVisMut$label.parameter.rot <- 90
# label col
LolCol <- list(gpar(col="#ec4646"),gpar(col="black"));names(LolCol) <- c("high","low")
LolCol = LolCol[ifelse(LolMutERBB3_forvis$Frequency > 50, 'high','low')]
LolVisMut$label.parameter.gp = LolCol
# height
LolVisMut$score <- ifelse(LolMutERBB3_forvis$Type == 'GBC', log2(LolMutERBB3_forvis$Frequency) + 0.5,
                          log(LolMutERBB3_forvis$Frequency, base = 4) + 0.1) # 点的高度
# circle large
LolVisMut$cex <- case_when(LolMutERBB3_forvis$Type == 'GBC' ~ 1,
                           LolMutERBB3_forvis$Type != 'GBC' ~ log(LolMutERBB3_forvis$Frequency,base = 30) + 0.1)
# up down
LolVisMut$SNPsideID <- ifelse(LolMutERBB3_forvis$Type == 'GBC', 'bottom', 'top')
# Protein info
LolVisProtein = GRanges("ERBB3",IRanges(c(LolDomainERBB3$Start, 1),
                                        width = c(LolDomainERBB3$End - LolDomainERBB3$Start, 1342),
                                        names = c(LolDomainERBB3$Label, 'Full-length')))
LolVisProtein$fill <- case_when(LolVisProtein@ranges@NAMES == 'Full-length' ~ 'black',
                                LolVisProtein@ranges@NAMES == 'Recep_L_domain' ~ ColJournal$Nature[1],
                                LolVisProtein@ranges@NAMES == 'Furin-like' ~ ColJournal$Nature[7],
                                LolVisProtein@ranges@NAMES == 'FU' ~ ColJournal$Nature[3],
                                LolVisProtein@ranges@NAMES == 'TM_ErbB3' ~ ColJournal$Nature[2],
                                LolVisProtein@ranges@NAMES == 'PTK_HER3' ~ ColJournal$Nature[4],
                                LolVisProtein@ranges@NAMES == 'Pkinase_Tyr' ~ ColJournal$Nature[9])
LolVisProtein$height <- case_when(LolVisProtein@ranges@NAMES == 'Full-length' ~ 0,
                                  LolVisProtein@ranges@NAMES == 'Recep_L_domain' ~ 0.06,
                                  LolVisProtein@ranges@NAMES == 'Furin-like' ~ 0.04,
                                  LolVisProtein@ranges@NAMES == 'FU' ~ 0.06,
                                  LolVisProtein@ranges@NAMES == 'TM_ErbB3' ~ 0.06,
                                  LolVisProtein@ranges@NAMES == 'PTK_HER3' ~ 0.04,
                                  LolVisProtein@ranges@NAMES == 'Pkinase_Tyr' ~ 0.06)
### ERBB Vis
lolliplot(LolVisMut, LolVisProtein, jitter="node", type="circle",
          xaxis = c(0, 200, 400, 600, 800, 1000, 1200, 1342))
# data for Table S3
write.xlsx(LolMutERBB3_forvis, file = 'lol ERBB3.xlsx', rowNames = T, colNames = T, overwrite = T)
#
rm(LolCol, LolDomain, LolDomainERBB2, LolDomainERBB3, LolAAReplace, LolSeqERBB2, LolSeqERBB3,
   LolMutERBB3, LolMutERBB3_clean, LolMutERBB3_forvis, maf_aachange, TempAllFreq,
   LolVisMut, LolVisProtein, Output, COSMIC_Genes, COSMIC_MutSite, ClinVar_Disease, ClinVar_Mut)

##### Figure 2E #####
##### 1.1 index prep
Input <- select(GBC_Main_Clinical, 1, 3, 4, 7:13, 16, 18, 23:24, 31:32) %>% as_tibble()
Input$Cluster <- NA
# Mutation
Input$`Mut EGFR` <- case_when((Input$Patient_ID %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'EGFR')$Tumor_Sample_Barcode) ~ 'Mut',
                              (Input$Patient_ID %in% GBC_Main_Mutation$Tumor_Sample_Barcode) ~ 'WT')
Input$`Mut ERBB2` <- case_when((Input$Patient_ID %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB2')$Tumor_Sample_Barcode) ~ 'Mut',
                               (Input$Patient_ID %in% GBC_Main_Mutation$Tumor_Sample_Barcode) ~ 'WT')
Input$`Mut ERBB3` <- case_when((Input$Patient_ID %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB3')$Tumor_Sample_Barcode) ~ 'Mut',
                               (Input$Patient_ID %in% GBC_Main_Mutation$Tumor_Sample_Barcode) ~ 'WT')
Input$`Mut ERBB4` <- case_when((Input$Patient_ID %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB4')$Tumor_Sample_Barcode) ~ 'Mut',
                               (Input$Patient_ID %in% GBC_Main_Mutation$Tumor_Sample_Barcode) ~ 'WT')
Input$`Mut ERBBs` <- case_when((Input$Patient_ID %in% filter(GBC_Main_Mutation, 
                                                             Hugo_Symbol %in% c('EGFR','ERBB2','ERBB4','ERBB3'))$Tumor_Sample_Barcode) ~ 'Mut',
                               (Input$Patient_ID %in% GBC_Main_Mutation$Tumor_Sample_Barcode) ~ 'WT')
## Pathway Mutation
Input$`Mut ERBB Pathway` <- case_when((Input$Patient_ID %in% filter(GBC_Main_Mutation, 
                                                                    Hugo_Symbol %in% filter(TotalPathway, Pathway == 'KEGG ErbB signaling pathway')$Gene,
                                                                    Hugo_Symbol %in% ref.cosmic$`Gene Symbol`)$Tumor_Sample_Barcode) ~ 'Mut',
                                      (Input$Patient_ID %in% GBC_Main_Mutation$Tumor_Sample_Barcode) ~ 'WT')
table(Input$`Mut ERBB Pathway`)

# add fusion info
Input$`Fusion ERBB2` <- NA

# CNV
# calc
Input$`CNV s EGFR` <- sapply(GBC_Main_Index$WES_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV_cate['EGFR',x]), NA)})
Input$`CNV s EGFR` <- SeparatedCNV(Input$`CNV s EGFR`)
Input$`CNV c EGFR` <- sapply(GBC_Main_Index$WES_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV['EGFR',x]), NA)})
Input$`CNV s ERBB2` <- sapply(GBC_Main_Index$WES_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV_cate['ERBB2',x]), NA)})
Input$`CNV s ERBB2` <- SeparatedCNV(Input$`CNV s ERBB2`)
Input$`CNV c ERBB2` <- sapply(GBC_Main_Index$WES_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV['ERBB2',x]), NA)})
Input$`CNV s ERBB3` <- sapply(GBC_Main_Index$WES_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV_cate['ERBB3',x]), NA)})
Input$`CNV s ERBB3` <- SeparatedCNV(Input$`CNV s ERBB3`)
Input$`CNV c ERBB3` <- sapply(GBC_Main_Index$WES_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV['ERBB3',x]), NA)})
Input$`CNV s ERBB4` <- sapply(GBC_Main_Index$WES_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV_cate['ERBB4',x]), NA)})
Input$`CNV s ERBB4` <- SeparatedCNV(Input$`CNV s ERBB4`)
Input$`CNV c ERBB4` <- sapply(GBC_Main_Index$WES_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV['ERBB4',x]), NA)})

##### RNA 
Input$`RNA EGFR` <- sapply(GBC_Main_Index$RNA_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['EGFR',x]), NA)})
Input$`RNA ERBB2` <- sapply(GBC_Main_Index$RNA_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['ERBB2',x]), NA)})
Input$`RNA ERBB3` <- sapply(GBC_Main_Index$RNA_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['ERBB3',x]), NA)})
Input$`RNA ERBB4` <- sapply(GBC_Main_Index$RNA_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_RNA_logTPM['ERBB4',x]), NA)})

##### Protein
Input$`Pro EGFR` <- sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['EGFR',x]), NA)})
Input$`Pro ERBB2` <- sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['ERBB2',x]), NA)})
Input$`Pro ERBB3` <- sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['ERBB3',x]), NA)})
Input$`Pro ERBB4` <- sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_Pro['ERBB4',x]), NA)})

##### Phos
Input$`Phos EGFR` <- sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_ProLevelPhos_knn['EGFR',x]), NA)})
Input$`Phos ERBB2` <- sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_ProLevelPhos_knn['ERBB2',x]), NA)})
Input$`Phos ERBB3` <- sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_ProLevelPhos_knn['ERBB3',x]), NA)})
Input$`Phos ERBB4` <- sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_ProLevelPhos_knn['ERBB4',x]), NA)})

# remove columns with all NA
Input = Input[,apply(Input, 2, function(x){sum(is.na(x)) != 195})]

##### factorize
Input$Sex = factor(Input$Sex)
Input$Pathological_type = factor(Input$Pathological_type, levels = c('adenocarcinoma','adenosquanmous carcinoma','neuroendocrine carcinoma',
                                                                     'adenocarcinoma; carcinosarcoma','carcinosarcoma','low grade intraepithelial neoplasia',
                                                                     'clear cell carcinoma','squamous cell carcinoma'))
Input$Differentiation = factor(Input$Differentiation, levels = c(1, 2, 3), ordered = T) 
Input$Regional_invasion = factor(Input$Regional_invasion)
Input$Regional_lymph_node_metastasis = factor(Input$Regional_lymph_node_metastasis)
Input$Distal_metastasis = factor(Input$Distal_metastasis)
Input$Vascular_invasion = factor(Input$Vascular_invasion)
Input$Perineural_invasion = factor(Input$Perineural_invasion)
Input$Tumor_stage2 = factor(Input$Tumor_stage2, levels = c(1,2,3,4), ordered = T)
Input$cholelithiasis = factor(Input$cholelithiasis)

# add domain
InputDomain <- mutate(Input,
                      `Domain Mut ERBB2` = case_when(Patient_ID %in% filter(GBC_Main_Mutation, grepl('S656F',aaChange))$Tumor_Sample_Barcode  ~ 'TMD',
                                                     `Mut ERBB2` == 'Mut' ~ 'KD'),
                      `Domain Mut ERBB3` = case_when(Patient_ID %in% c('GBC_079','GBC_164')  ~ 'RD',
                                                     Patient_ID %in% c('GBC_174','GBC_191','GBC_128','GBC_152',
                                                                       'GBC_068','GBC_100','GBC_192','GBC_237')  ~ 'Furin',
                                                     Patient_ID %in% c('GBC_094','GBC_112','GBC_278')  ~ 'TMD',
                                                     Patient_ID %in% c('GBC_186')  ~ 'Signal transduction'))
# conclusion
InputDomain$`Genome ERBB2` <- case_when(InputDomain$`CNV s ERBB2` == 'Amplification' ~ 'Amp',
                                        InputDomain$`Domain Mut ERBB2` == 'KD' ~ 'KD',
                                        InputDomain$`Domain Mut ERBB2` == 'TMD' ~ 'TMD',
                                        InputDomain$`Mut ERBB2` == 'WT' ~ 'WT')
InputDomain$`Genome ERBB3` <- case_when(InputDomain$`CNV s ERBB3` == 'Amplification' ~ 'Amp',
                                        InputDomain$`Domain Mut ERBB3` == 'Furin' ~ 'Furin',
                                        InputDomain$`Domain Mut ERBB3` == 'RD' ~ 'RD',
                                        InputDomain$`Domain Mut ERBB3` == 'Signal transduction' ~ 'Signal transduction',
                                        InputDomain$`Domain Mut ERBB3` == 'TMD' ~ 'TMD',
                                        InputDomain$`Mut ERBB3` == 'WT' ~ 'WT')
InputDomain$`Genome ERBB2 status` <- case_when(InputDomain$`CNV s ERBB2` == 'Amplification' ~ 'Amp',
                                               InputDomain$`Mut ERBB2` == 'Mut' ~ 'Mut',
                                               InputDomain$`Mut ERBB2` == 'WT' ~ 'WT')
InputDomain$`Genome ERBB3 status` <- case_when(InputDomain$`CNV s ERBB3` == 'Amplification' ~ 'Amp',
                                               InputDomain$`Mut ERBB3` == 'Mut' ~ 'Mut',
                                               InputDomain$`Mut ERBB3` == 'WT' ~ 'WT')

##### 1.2 calc ssgsea
PathwayPanelErbB <- c('ErbBs' = 'KEGG ErbB signaling pathway', 'EGFR' = 'REACTOME Signaling by EGFR',
                      'ERBB2' = 'REACTOME Signaling by ERBB2', 'ERBB4' = 'REACTOME Signaling by ERBB4',
                      'Ras' = 'KEGG Ras signaling pathway', 'PI3K' = 'KEGG PI3K-Akt signaling pathway')
InputGSVAErbB <- GSVA::gsva(as.matrix(GBC_Main_Pro), TotalPathwayGSVA[, PathwayPanelErbB], method = 'ssgsea') %>%
  t() %>% scale() %>% t() %>% as.data.frame()
# add to index
InputDomain$`Pathway ErbBs` = sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(InputGSVAErbB['KEGG ErbB signaling pathway',x]), NA)})
InputDomain$`Pathway EGFR` = sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(InputGSVAErbB['REACTOME Signaling by EGFR',x]), NA)})
InputDomain$`Pathway ERBB2` = sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(InputGSVAErbB['REACTOME Signaling by ERBB2',x]), NA)})
InputDomain$`Pathway ERBB4` = sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(InputGSVAErbB['REACTOME Signaling by ERBB4',x]), NA)})
InputDomain$`Pathway Ras` = sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(InputGSVAErbB['KEGG Ras signaling pathway',x]), NA)})
InputDomain$`Pathway PI3K` = sapply(GBC_Main_Index$ProPhos_T, function(x){ifelse(!is.na(x), as.numeric(InputGSVAErbB['KEGG PI3K-Akt signaling pathway',x]), NA)})

##### 1.3 vis prep
Input2 = tibble(`Genome ERBB2` = rep('TAC', 137),
                `Genome ERBB3` = rep('TAC', 137),
                `Genome ERBB2 status` = rep('TAC', 137),
                `Genome ERBB3 status` = rep('TAC', 137),
                `Pathway ErbBs` = as.numeric(InputGSVAErbB[1,grepl('P',colnames(InputGSVAErbB))]),
                `Pathway EGFR` = as.numeric(InputGSVAErbB[2,grepl('P',colnames(InputGSVAErbB))]),
                `Pathway ERBB2` = as.numeric(InputGSVAErbB[3,grepl('P',colnames(InputGSVAErbB))]),
                `Pathway ERBB4` = as.numeric(InputGSVAErbB[4,grepl('P',colnames(InputGSVAErbB))]),
                `Pathway Ras` = as.numeric(InputGSVAErbB[5,grepl('P',colnames(InputGSVAErbB))]),
                `Pathway PI3K` = as.numeric(InputGSVAErbB[6,grepl('P',colnames(InputGSVAErbB))]))
InputVis = rbind(InputDomain[,42:51], Input2)
# factorize
InputVis$`Genome ERBB2` = factor(InputVis$`Genome ERBB2`, levels = c('TAC','WT','RD','Furin','TMD','KD','Signal transduction','Amp'), ordered = T)
InputVis$`Genome ERBB2 status` = factor(InputVis$`Genome ERBB2 status`, levels = c('TAC','WT','Mut','Amp'), ordered = T)
InputVis$`Genome ERBB3` = factor(InputVis$`Genome ERBB3`, levels = c('TAC','WT','RD','Furin','TMD','KD','Signal transduction','Amp'), ordered = T)
InputVis$`Genome ERBB3 status` = factor(InputVis$`Genome ERBB3 status`, levels = c('TAC','WT','Mut','Amp'), ordered = T)
table(InputVis$`Genome ERBB2`)

##### 1.4 vis
InputCol = c(ColColor$`Low-Blue`[5],'green', wes_palette('Moonrise3', 5, type = c("discrete")),ColColor$`High-Orange`[5] ,'red')
names(InputCol) = c('TAC','WT','RD','Furin','TMD','KD','Signal transduction','Mut','Amp')
## domain
InputVisDomain = reshape2::melt(InputVis[,c(1,2,5:10)], id.var = c('Genome ERBB2','Genome ERBB3'),
                                variable.name = 'Pathway', value.name = 'Pathway score') %>% tibble()
InputVisDomain = reshape2::melt(InputVisDomain, id.var = c('Pathway','Pathway score'),
                                variable.name = 'Gene', value.name = 'Genome status')
InputVisDomain$Gene = factor(InputVisDomain$Gene, levels = c('Genome ERBB2','Genome ERBB3'), ordered = T)
InputVisDomain$`Genome status` = factor(InputVisDomain$`Genome status`, levels = c('TAC','WT','RD','Furin','TMD','KD','Signal transduction','Amp'), ordered = T)
#
InputVisDomain.mean.mutant.erbb3 <- filter(InputVisDomain, Gene == 'Genome ERBB3', !(`Genome status` %in% c('WT','TAC','Amp')))
InputVisDomain.mean.mutant.erbb3 <- filter(InputVisDomain.mean.mutant.erbb3, !is.na(`Genome status`))
InputVisDomain.mean.mutant.erbb3$`Genome status` = 'Mutant'
#integ
InputVisDomain.erbb3 <- rbind(filter(InputVisDomain, Gene == 'Genome ERBB3', !is.na(`Genome status`)),
                              InputVisDomain.mean.mutant.erbb3)
#factorize
InputVisDomain.erbb3$`Genome status` = factor(InputVisDomain.erbb3$`Genome status`,
                                              levels = c('TAC','WT','RD','Furin','TMD','KD','Signal transduction','Mutant','Amp'), ordered = T)
# new color
InputCol2 = c(ColColor$`Low-Blue`[5],'green', wes_palette('Moonrise3', 5, type = c("discrete")),ColColor$`High-Orange`[5] ,ColColor$`Single-BlueGrey`[9],'red')
names(InputCol2) = c('TAC','WT','RD','Furin','TMD','KD','Signal transduction','Mutant','Amp')
# 
ggboxplot(na.omit(InputVisDomain.erbb3), 'Genome status', 'Pathway score', 
          # color = 'Genome status', 
          #add = 'jitter', add.params = list(size = 0.2)
          ) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "WT", na.rm = T) +
  scale_color_manual(values = InputCol) +
  facet_wrap(~Gene+Pathway, scales = 'free_y', nrow = 2) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('') 
# data for Table S3
write.xlsx(InputVisDomain.erbb3, file = 'ERBB3 domain pathway.xlsx', rowNames = T, colNames = T, overwrite = T)
#
rm(InputVis, InputVisDomain, InputVisNoDomain, InputGSVAErbB, InputGSVAErbB_TAC, Input2, a, ssGSEA1, ssGSEAInput_Phos, ssGSEAInput_Pro, ssGSEAOutput_Phos, ssGSEAInput_Pro,
   InputDomain, InputVisDomain.erbb3, InputVisDomain.mean.mutant.erbb3)

##### Figure 2K #####
##### 1.1 data loading
load('../20240303.GBC.ERBB2.RData')
GBCphosclean = GBCphosclean[, colnames(GBCproclean)]
# median filtering, imputation
GBCphosclean_halfNA <- GBCphosclean[apply(GBCphosclean, 1, function(x){sum(is.na(x)) <= 15}),]
GBCphosclean_halfNA <- GBCphosclean_halfNA[,apply(GBCphosclean_halfNA, 2, function(x){sum(is.na(x)) <= 15496/2})]
GBCphosclean_knn <- DreamAI(GBCphosclean_halfNA, method = 'KNN')$Ensemble %>% as.data.frame()

##### 1.2 index
identical(colnames(GBCphosclean), colnames(GBCproclean))
colnames(GBCphosclean)
Index <- tibble(Sample = colnames(GBCphosclean),
                Type = c(rep('EV',3),rep('ERBB2OE',3),rep('ERBB3',3),rep('PGAP3Fusion',3),
                         rep('ERBB2_S310F',3),rep('ERBB2_S656F',3),
                         rep('ERBB3_V104L',3),rep('ERBB3_V653E',3),rep('ERBB3_R667L',3),rep('ERBB3_R679Q',3)))
Index$Pro = Index$Sample
Index$Phos = ifelse(Index$Sample %in% colnames(GBCphosclean_knn), Index$Sample, NA)
# remove pgap3 fusion
Index = filter(Index, Type != 'PGAP3Fusion')
# factorize
table(Index$Type)
Index$Type = factor(Index$Type, levels = c('EV','ERBB2OE','ERBB2_S310F','ERBB2_S656F',
                                           'ERBB3','ERBB3_V104L','ERBB3_V653E','ERBB3_R667L','ERBB3_R679Q'), ordered = T)

##### normalization data frame before heatmap
GBCphosclean <- cbind(GBCphosclean, 
                      data.frame(E3_104_rep3 = rep(NA, nrow(GBCphosclean)),
                                 E3_104_rep2 = rep(NA, nrow(GBCphosclean)),
                                 E3_104_rep3 = rep(NA, nrow(GBCphosclean)),
                                 E3_667_rep3 = rep(NA, nrow(GBCphosclean))))
GBCphosclean = GBCphosclean[, colnames(GBCproclean)]
# data for Table S3
GBCphosclean.rec <- GBCphosclean_knn[,na.omit(Index$Phos)]
colnames(GBCphosclean.rec) <- sapply(colnames(GBCphosclean.rec), function(x){filter(Index, Phos == x)$Type %>% as.character()})
write.xlsx(GBCphosclean.rec, file = 'ERBB cellline proteomics.xlsx', rowNames = T, colNames = T, overwrite = T)

##### 1.3 Kinase activity calc
library(KSEAapp)
ref.kinase = KSData
# Index
Index.phos = filter(Index, !is.na(Phos))
# func
Enrich_create_GSVA_ref_list <- function(gmt){
  require(dplyr)
  colnames(gmt) = c('term', 'gene')
  term <- gmt$term[!duplicated(gmt[,1])] %>% as.vector()
  hallmark_gsva <- as.data.frame(c(rep(1, length(term)))) %>% t() %>% as.data.frame()
  colnames(hallmark_gsva) <- term
  for(i in 1:length(term)){
    a <- gmt[gmt[,1] == term[i], 2] %>% as.data.frame()
    for(j in 1:(dim(a)[1])){
      hallmark_gsva[j, i] = a[j, 1]
    }
  }
  
  return(hallmark_gsva)
}
kinase.set = filter(ref.kinase, GENE %in% c('AKT1','AKT2','MAPK3','MAPK1','MAPK10','MAPK9','MAPK8',
                                            'MAPK7','MAPK15','MAPK14'),
                    SUB_ORGANISM == 'human')
kinase.set = kinase.set[,c(3,8,10)]
kinase.set = tibble(Kinase = kinase.set$GENE, 
                    Site = paste0(kinase.set$SUB_GENE,':',kinase.set$SUB_MOD_RSD))
kinase.set.enrich = Enrich_create_GSVA_ref_list(kinase.set)
# GSVA
library(GSVA)
kinase.enrich.output <- gsva(as.matrix(GBCphosclean_knn), kinase.set.enrich, method = 'ssgsea') %>%
  t() %>% scale() %>% t() %>% as.data.frame()
# order kinase
kinase.enrich.output = kinase.enrich.output[c('AKT1','AKT2','MAPK3','MAPK1','MAPK10','MAPK9','MAPK8',
                                              'MAPK7','MAPK15','MAPK14'), Index.phos$Phos]

##### 1.4 vis
Heatmap(kinase.enrich.output, cluster_rows = F, cluster_columns = F, 
        na_col = 'grey90', name = 'Kinase Activity',
        column_split = Index.phos$Type, column_title_rot = 45,
        show_column_names = F,
        height = unit(5*1.1,'cm'), width = unit(11.5*1.1,'cm'),
        row_split = c('AKTs','AKTs','ERK1-2','ERK1-2',
                      'JNKs','JNKs','JNKs','ERK5','ERK7-8','p38'))
# data for Table S3
colnames(kinase.enrich.output) = colnames(GBCphosclean.rec)
write.xlsx(kinase.enrich.output, file = 'Kinase activity.xlsx', rowNames = T, colNames = T, overwrite = T)







