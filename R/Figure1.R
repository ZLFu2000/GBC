##### Figure 1 and Figure S1
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
# 2 color
load('../Colors (ggsci).RData')
# 3 Genesets
load('../Hallmark_KEGG_Reactome.RData')
# 4 COSMIC Cancer Gene Census
ref.cosmic <- fread('../COSMIC genes.csv')

##### Prep2. Index #####
index.raw <- tibble(Patient = c(rep(GBC_Main_Index$Sample, 2)),
                Sample = c(paste0(GBC_Main_Index$Sample,'_T'), paste0(GBC_Main_Index$Sample, '_P')),
                SampleWES = c(GBC_Main_Index$WES_T, GBC_Main_Index$WES_P),
                SampleRNA = c(GBC_Main_Index$RNA_T, GBC_Main_Index$RNA_P),
                SamplePro = c(GBC_Main_Index$ProPhos_T, GBC_Main_Index$ProPhos_P),
                SamplePhos = c(GBC_Main_Index$ProPhos_T, GBC_Main_Index$ProPhos_P),
                Type = c(rep('Cancer', 195), rep('Adjacent', 195)) %>%
                  factor(levels = c('Adjacent', 'Cancer'), ordered = T))
# Calc TMB
TempTMB <- tmb(GBC_Main_Maf)
index.raw$TMB <- sapply(index.raw$Patient, function(x){filter(TempTMB, Tumor_Sample_Barcode == x)$total_perMB %>% as.numeric()}) %>% as.numeric()
# rm
rm(TempTMB, TempTMBCalc, Temp)

# add pathotype
index.raw$Pathology = rep(GBC_Main_Clinical$Pathological_type, 2)
index.raw = index.raw %>%
  mutate(Pathology = case_when(Pathology == 'adenocarcinoma' ~ 'AC',
                               Pathology == 'adenosquanmous carcinoma' ~ 'AS',
                               Pathology == 'neuroendocrine carcinoma' ~ 'NE',
                               !is.na(Pathology) ~ 'Others') %>%
           factor(levels = c('AC','AS','NE','Others'), ordered = T)) 
# add regional invasion
index.raw$LM <- c(GBC_Main_Clinical$Regional_invasion, GBC_Main_Clinical$Regional_invasion) %>%
  factor(levels = c(0, 1), ordered = T)
# add stage
index.raw$Stage <- c(GBC_Main_Clinical$Tumor_stage2, GBC_Main_Clinical$Tumor_stage2) %>%
  factor(levels = c(1,2,3,4), ordered = T)

# add ERBB2/3 and ERBB pathway mutation
index.raw$`Mut_ERBB2` <- case_when((index.raw$Patient %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB2')$Tumor_Sample_Barcode) & index.raw$Type == 'Cancer' ~ 'Mut',
                               (index.raw$Patient %in% GBC_Main_Mutation$Tumor_Sample_Barcode) & index.raw$Type == 'Cancer'~ 'WT')
index.raw$`Mut_ERBB3` <- case_when((index.raw$Patient %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ERBB3')$Tumor_Sample_Barcode) & index.raw$Type == 'Cancer' ~ 'Mut',
                               (index.raw$Patient %in% GBC_Main_Mutation$Tumor_Sample_Barcode) & index.raw$Type == 'Cancer' ~ 'WT')
index.raw$`Mut_ERBB_Pathway` <- case_when((index.raw$Patient %in% filter(GBC_Main_Mutation, Hugo_Symbol %in% filter(TotalPathway, Pathway == 'KEGG ErbB signaling pathway')$Gene)$Tumor_Sample_Barcode) & index.raw$Type == 'Cancer' ~ 'Mut',
                                      (index.raw$Patient %in% GBC_Main_Mutation$Tumor_Sample_Barcode) & index.raw$Type == 'Cancer' ~ 'WT')

## add ERBB2/3 CNV
# separated
SeparatedCNV <- function(vec) {
  vec = case_when(vec == -2 ~ 'Deletion',
                  vec == -1 ~ 'Shallow Deletion',
                  vec == 0 ~ 'Diploid', 
                  vec == 1 ~ 'Gain',
                  vec == 2 ~ 'Amplification') %>%
    factor(levels = c('Amplification','Gain', 'Diploid', 'Shallow Deletion', 'Deletion'), ordered = T)
  return(vec)
}
index.raw$`CNV_ERBB2` <- sapply(index.raw$SampleWES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV_cate['ERBB2',x]), NA)})
index.raw$`CNV_ERBB2` <- SeparatedCNV(index.raw$`CNV_ERBB2`)
index.raw$`CNV_ERBB3` <- sapply(index.raw$SampleWES, function(x){ifelse(!is.na(x), as.numeric(GBC_Main_CNV_cate['ERBB3',x]), NA)})
index.raw$`CNV_ERBB3` <- SeparatedCNV(index.raw$`CNV_ERBB3`)

# Add survival info
index.raw$OStime = sapply(index.raw$Patient, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$OS_days/30 %>% as.numeric()})
index.raw$OS = sapply(index.raw$Patient, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$OS_event %>% as.numeric()})

# tumor only
index.tumor <- filter(index.raw, Type == 'Cancer')

##### Figure 1A was manually created using Adobe Illustrator #####
##### Figure S1A were manually created using Adobe Illustrator #####
##### Values for the bar plot in Figure S1B were directly obtained from the row counts of the original data matrix #####

##### Figure S1C #####
library(FactoMineR)
library(factoextra)
# input prep
index.pca <- tibble(Sample = colnames(GBC_Main_Pro),
                    Type = ifelse(grepl('T', colnames(GBC_Main_Pro)), 'Tumor', 'Adjacent') %>%
                      factor(levels = c('Tumor','Adjacent'), ordered = T))
# Add batch
index.pca$Batch <- table.s1$Batch # could be found in Table S1
# check sample order
identical(colnames(GBC_Main_Pro), index.pca$Sample)
# PCA
pca.result <- PCA(t(GBC_Main_Pro), graph = FALSE)
# Vis.df
pca.df <- data.frame(
  Sample = colnames(GBC_Main_Pro),
  PC1 = pca.result$ind$coord[, 1],  
  PC2 = pca.result$ind$coord[, 2],  
  Batch = factor(index.pca$Batch),  
  Type = index.pca$Type
)
# Vis
ggplot(pca.df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +  
  stat_ellipse(aes(group = Batch), type = "norm", level = 0.95, linetype = 2) +  
  theme_minimal(base_size = 15) +  
  theme(
    text = element_text(family = "Arial", color = "black"),  
    panel.grid = element_blank(), 
    panel.background = element_rect(fill = "white", color = "white"),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    legend.position = "right", 
    axis.text = element_text(color = "black")  
  ) +
  labs(
    title = "PCA of Protein Expression Matrix",
    x = paste0("PC1 (", round(pca.result$eig[, 2][1], 2), "%)"),
    y = paste0("PC2 (", round(pca.result$eig[, 2][2], 2), "%)"),
    color = "TMT Batch" 
  )  # 5*6 TMT Batch
#
rm(index.pca, pca.result, pca.df)

##### Figuer S1D & Figuer S1E #####
##### 1.1 index prep
index.t <- filter(select(index.raw, Patient, Sample, SampleWES, SampleRNA, SamplePro, SamplePhos, Type), 
                  Type == 'Cancer') %>% na.omit() # N = 129
index.n <- filter(select(index.raw, Patient, Sample, SampleWES, SampleRNA, SamplePro, SamplePhos, Type),
                  Type == 'Adjacent') %>% na.omit() # N = 63
##### mtx prep
# rna
input.t.rna <- GBC_Main_RNA_logTPM[, index.t$SampleRNA]
input.n.rna <- GBC_Main_RNA_logTPM[, index.n$SampleRNA]
# pro
input.t.pro <- GBC_Main_Pro[, index.t$SamplePro]
input.n.pro <- GBC_Main_Pro[, index.n$SamplePro]

##### mutual feature prep
mutual.gene <- intersect(rownames(input.t.rna), rownames(input.t.pro))
input.t.rna <- input.t.rna[mutual.gene, ]; input.t.pro <- input.t.pro[mutual.gene, ]
input.n.rna <- input.n.rna[mutual.gene, ]; input.n.pro <- input.n.pro[mutual.gene, ]

##### 1.2 calc cor
# calc cor
calc.cor <- function(mtx.rna, mtx.pro) {
  output <- tibble(Gene = character(), Cor = numeric(), Pval = numeric(), FDR = numeric())
  # calc
  for (i in 1:length(mutual.gene)) {
    temp.out <- tibble(Gene = mutual.gene[i], 
                       Cor = cor(as.numeric(mtx.rna[mutual.gene[i],]),
                                 as.numeric(mtx.pro[mutual.gene[i],]),
                                 method = 'spearman'),
                       Pval = cor.test(as.numeric(mtx.rna[mutual.gene[i],]),
                                       as.numeric(mtx.pro[mutual.gene[i],]),
                                       method = 'spearman')$p.value, 
                       FDR = NA)
    # rbind
    output <- rbind(output, temp.out)
  }
  # fdr
  output$FDR <- p.adjust(output$Pval, method = 'fdr')
  # return
  return(output)
}
cor.t <- calc.cor(input.t.rna, input.t.pro)
cor.n <- calc.cor(input.n.rna, input.n.pro)
##### data for Table S1
library(openxlsx)
write.xlsx(arrange(cor.t, Gene), file = 'RNA-Pro Cor tumor.xlsx', rowNames = T, colNames = T, overwrite = T)
write.xlsx(arrange(cor.n, Gene), file = 'RNA-Pro Cor NAT.xlsx', rowNames = T, colNames = T, overwrite = T)

# vis cor
cor.vis <- tibble(Gene = cor.t$Gene,
                  CorTumor = cor.t$Cor, CorNAT = cor.n$Cor)
cor.vis <- melt(cor.vis, id.vars = 'Gene', variable.name = 'Type', value.name = 'Cor') %>% as_tibble()

##### 1.3 vis
plot_grid(gghistogram(cor.t, x = "Cor", fill = ColJournal$Nature[1],
                      add = "median", bins = 50) +
            scale_x_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.5), labels = seq(-1,1,0.5)) +
            xlab(NULL) + ylab('Probability density')
          , 
          gghistogram(cor.n, x = "Cor", fill = ColJournal$Nature[2],
                      add = "median", bins = 50) +
            scale_x_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.5), labels = seq(-1,1,0.5)) +
            xlab('Spearman’s correlation') + ylab(''),
          nrow = 2) # 5*7.5 RNA-protein cor
# stat
median(cor.t$Cor)
median(cor.n$Cor)
table(cor.t$Cor > 0)[2] / length(cor.t$Cor)
table(cor.t$FDR <= 0.01 & cor.t$Cor > 0)[2] / length(cor.t$Cor)
mean(cor.t$Cor)
mean(cor.n$Cor)
table(cor.n$Cor > 0)[2] / length(cor.n$Cor)
table(cor.n$FDR <= 0.01 & cor.n$Cor > 0)[2] / length(cor.n$Cor)

##### 1.4 cor-based gsea
gsea.input.t <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = cor.t$Gene, FC = cor.t$Cor))
gsea.input.n <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = cor.n$Gene, FC = cor.n$Cor))
# gsea
gsea.t <- clusterProfiler::GSEA(gsea.input.t, TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.05) %>% arrange(desc(NES))
gsea.n <- clusterProfiler::GSEA(gsea.input.n, TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.05) %>% arrange(desc(NES))
##### data for Table S1
write.xlsx(gsea.t, file = 'Cor GSEA.xlsx', rowNames = T, colNames = T, overwrite = T)
# Vis
# tumor
gsea.t.vis <- gsea.t[c('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
                       'REACTOME Collagen biosynthesis and modifying enzymes',
                       'REACTOME Extracellular matrix organization',
                       'HALLMARK_INFLAMMATORY_RESPONSE',
                       'HALLMARK_HYPOXIA',
                       'KEGG Basal transcription factors',
                       'REACTOME Mitochondrial tRNA aminoacylation',
                       'KEGG Protein export',
                       'KEGG Aminoacyl-tRNA biosynthesis'),]
gsea.t.vis$Type <- ifelse(gsea.t.vis$NES > 0, 'GBC', 'NAT')
ggplot(gsea.t.vis, aes(x = reorder(ID, NES), y = NES)) +
  geom_bar(stat = "identity", aes(fill = Type)) +
  coord_flip() +
  theme_minimal() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c('GBC' = ColJournal$Nature[1], 'NAT' = ColJournal$Nature[2])) +
  theme(
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black"),  
    text = element_text(family = "Arial", color = "black")
  ) +
  labs(title = "Spearman’s correlation in GBC tissues",
       x = "Pathway",
       y = "-log10(p.adjust)")
#
rm(gsea.n, gsea.t, gsea.t.vis, index.t, index.n, gsea.input.n, gsea.input.t, mutual.gene,
   input.n.pro, input.n.rna, input.t.pro, input.t.rna,
   cor.n, cor.t, cor.vis)

##### Figure S1F (oncoprint) #####
length(table(GBC_Main_Maf@data$Tumor_Sample_Barcode)) #  = 167
##### 1.1 gene selection and Maf prep
# In cosmic
maf.filtered <-  GBC_Main_Maf@data %>% filter(Hugo_Symbol %in% ref.cosmic$`Gene Symbol`)
# top 20 freq
temp.top20 <- maf.filtered %>%
  count(Hugo_Symbol, name = "Mutation_Count") %>%
  arrange(desc(Mutation_Count)) %>%
  slice_head(n = 20)
maf.filtered <- maf.filtered %>% filter(Hugo_Symbol %in% temp.top20$Hugo_Symbol)
# before writing: retain all samples 
temp.samples <- levels(maf.filtered$Tumor_Sample_Barcode)
# dataframe 2 maf
write.table(maf.filtered, file = "filtered_top20_genes.maf", sep = "\t", quote = FALSE, row.names = FALSE)
maf.filtered <- maftools::read.maf(maf = 'filtered_top20_genes.maf')
# check samples
print(getSampleSummary(maf.filtered))
# samples with mutation
temp.samples_with_mutations <- unique(maf.filtered@data$Tumor_Sample_Barcode)
# samples without mutation
temp.no_mutation_samples <- setdiff(temp.samples, temp.samples_with_mutations)
##### join
# empty df
no_mutation_df <- maf.filtered@data[0, ]
for(i in 1:length(temp.no_mutation_samples)) {
  temp.row <- maf.filtered@data[1, ] 
  temp.row[!is.na(temp.row)] = NA
  temp.row$Tumor_Sample_Barcode <- temp.no_mutation_samples[i]
  no_mutation_df <- rbind(no_mutation_df, temp.row)
  rm(i, temp.row)
}
# join maf
maf.filtered.na <- bind_rows(maf.filtered@data, no_mutation_df)
maf.filtered.na$Tumor_Sample_Barcode <- factor(maf.filtered.na$Tumor_Sample_Barcode, 
                                               levels = levels(getSampleSummary(GBC_Main_Maf)$Tumor_Sample_Barcode), ordered = T)
maf.filtered.na <- arrange(maf.filtered.na, Tumor_Sample_Barcode)
# final maf for Figure S1F
write.table(maf.filtered.na, file = "filtered_top20_genes_retainNA.maf", sep = "\t", quote = FALSE, row.names = FALSE)
maf.filtered.na <- maftools::read.maf(maf = 'filtered_top20_genes_retainNA.maf')

##### 1.2 Maf visualization
oncoplot(
  maf = maf.filtered.na,
  top = 20, 
  colors = c(
    "Frame_Shift_Del" = "#4CAF50",  
    "Frame_Shift_Ins" = "#B3E5FC",  
    "Missense_Mutation" = "#1565C0", 
    "Nonsense_Mutation" = "#F44336",  
    "In_Frame_Del" = "#FFCDD2",  
    "In_Frame_Ins" = "#FFEB3B", 
    "Splice_Site" = "#FF9800",  
    "Multi_Hit" = "#000000", 
    'Translation_Start_Site' = '#8BC34A'
  )
  ,  
  drawRowBar = TRUE,  
  drawColBar = TRUE,  
  annotationFontSize = 0.8,  
  fontSize = 0.7,  
  legendFontSize = 1.0,  
  bgCol = "grey94",
  removeNonMutated = F,
  showTumorSampleBarcodes = F  
) 
#
rm(maf.filtered, maf.filtered.na, no_mutation_df, temp.top20, 
   temp.no_mutation_samples, temp.samples_with_mutations, temp.samples)

##### Figure S1G #####
##### Genes
InputGeneMut <- c('TP53','ELF3','EGFR','ERBB2','ERBB3','SMAD4','PIK3CA','CTNNB1','KRAS')
##### Cohorts
# 1. CCR: Giraldo cohort
InputMafCCR <- fread('../data_mutations.txt')
# 2. JOH: Nepal cohort
InputMafJOH <- read.xlsx('../WES.xlsx',
                         startRow = 2, rowNames = F, colNames = T)
# 3. NG: Li cohort (Liu YB)
InputMafNGWES <- read.xlsx('../WES and Panel.xlsx',
                           sheet = 1, startRow = 2, rowNames = F, colNames = T)
InputMafNGTar <- read.xlsx('../WES and Panel.xlsx',
                           sheet = 2, startRow = 2, rowNames = F, colNames = T)
length(table(InputMafNGWES$Mutated_sample.6)) # 32
length(table(InputMafNGTar$Mutated_sample.6)) # 49
# 4. Gut: Li cohort (Liu YB)
InputMafGut <- read.xlsx('..mutation.xlsx',
                         startRow = 2, rowNames = F, colNames = T)
InputMafGut <- tibble(Gene = str_split_fixed(InputMafGut$Hugo_Symbol, ' ', 2)[,1],
                      Sample = paste0(InputMafGut$Tumor_Sa.mple_Barc.ode, InputMafGut$X19))

##### Maf Gene freq calc manually
Input <- tibble(Gene = InputGeneMut,
                `FU-GBC` = 0, CCR = 0, JOH = 0, NG = 0, Gut = 0)
# FU-GBC
Input$`FU-GBC`[1] = 80/167; Input$`FU-GBC`[2] = 15/167; Input$`FU-GBC`[3] = 1/167; Input$`FU-GBC`[4] = 3/167
Input$`FU-GBC`[5] = 14/167; Input$`FU-GBC`[6] = 14/167; Input$`FU-GBC`[7] = 8/167; Input$`FU-GBC`[8] = 8/167
Input$`FU-GBC`[9] = 7/167
# CCR
Input$CCR[1] = 154/244; Input$CCR[2] = 21/244; Input$CCR[3] = 3/244; Input$CCR[4] = 19/244
Input$CCR[5] = 16/244; Input$CCR[6] = 52/244; Input$CCR[7] = 26/244; Input$CCR[8] = 15/244
Input$CCR[9] = 18/244
# JOH
Input$JOH[1] = 27/92; Input$JOH[2] = 6/92; Input$JOH[3] = 0/92; Input$JOH[4] = 5/92
Input$JOH[5] = 3/92; Input$JOH[6] = 1/92; Input$JOH[7] = 3/92; Input$JOH[8] = 3/92
Input$JOH[9] = 1/92
# NG
Input$NG[1] = 0.471; Input$NG[2] = 2/57; Input$NG[3] = 3/57; Input$NG[4] = 7/57
Input$NG[5] = 7/57; Input$NG[6] = 3/57; Input$NG[7] = 3/57; Input$NG[8] = 3/57
Input$NG[9] = 0.078
# Gut
Input$Gut[1] = 43/157; Input$Gut[2] = 6/157; Input$Gut[3] = 8/157; Input$Gut[4] = 12/157
Input$Gut[5] = 13/157; Input$Gut[6] = 18/157; Input$Gut[7] = 8/157; Input$Gut[8] = 5/157
Input$Gut[9] = 8/157
##### data for Table S2
write.xlsx(Input, file = 'MutFreq of external cohorts.xlsx',
           rowNames = T, colNames = T, overwrite = T)

##### arrange
Input <- melt(Input, id.vars = 'Gene', variable.name = 'Cohort', value.name = 'Frequency') %>% as_tibble()
Input$Gene = factor(Input$Gene, levels = rev(InputGeneMut), ordered = T)
Input$Cohort = factor(Input$Cohort, levels = c('FU-GBC','CCR','JOH','NG','Gut'), ordered = T)

##### Vis
ggbarplot(Input, "Gene", "Frequency", fill = "Cohort", color = "Cohort", 
          label = F, position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0,0)) + xlab('') + ylab('Mutation Frequency') +
  scale_color_manual(values = c('FU-GBC'='red','CCR'=ColJournal$Nature[5],'JOH'=ColJournal$Nature[4],
                                'NG'=ColJournal$Nature[2],'Gut'='grey')) +
  scale_fill_manual(values = c('FU-GBC'='red','CCR'=ColJournal$Nature[5],'JOH'=ColJournal$Nature[4],
                               'NG'=ColJournal$Nature[2],'Gut'='grey')) +
  coord_flip() + theme(legend.position = 'left') 

#
rm(Input, InputMafCCR, InputMafGut, InputMafJOH, InputMafNGTar, InputMafNGWES, InputGeneCNV, InputGeneMut)

##### Figure 1B #####
##### 1.1 calc 
Input <- somaticInteractions(maf = GBC_Main_Maf, top = 25, pvalue = c(0.05), genes = names(sort(table(GBC_Main_Mutation$Hugo_Symbol), decreasing = T))[1:100])
Input <- filter(Input, pValue <= 0.05)
##### data for Table S2
write.xlsx(Input, file = 'Co-mutation(sig).xlsx', rowNames = T, colNames = T, overwrite = T)
##### 1.2 Visualization 
##### 1. SMAD4 + ARID1A 3.67e-4
##### 2. SMAD4 + ATM    1.45e-2
##### 3. ERBB3 + KMT2C  2.50e-2
##### 4. TP53 + KEAP1   2.70e-2
VisCoMut <- function(gene1, gene2){
  Output <- tibble(Sample = names(table(GBC_Main_Mutation$Tumor_Sample_Barcode)),
                   Gene1 = rep('Mut', 167), Gene2 = rep("Mut", 167))
  Output$Gene1 = ifelse(Output$Sample %in% filter(GBC_Main_Mutation, Hugo_Symbol == gene1)$Tumor_Sample_Barcode, 'Mut', 'WT')
  Output$Gene2 = ifelse(Output$Sample %in% filter(GBC_Main_Mutation, Hugo_Symbol == gene2)$Tumor_Sample_Barcode, 'Mut', 'WT')
  Output = arrange(Output, Gene1, Gene2)
  colnames(Output)[2:3] = c(gene1, gene2)
  # return(Output)
  #### Vis
  TempCol = c(ColColor$`Single-Brown`[10],'grey95'); names(TempCol) = c('Mut','WT')
  HmInput <- as.matrix(Output[,2:3]) %>% t()
  Heatmap(HmInput, col = TempCol, name = 'Mutation status', 
          height = unit(2,'cm'),width = unit(6,'cm'))
}
draw(VisCoMut('ARID1A','SMAD4') %v%
       VisCoMut('SMAD4','ATM') %v%
       VisCoMut('ERBB3','KMT2C') %v%
       VisCoMut('TP53','KEAP1') %v%
       VisCoMut('TP53','KRAS') %v%
       VisCoMut('TP53','ATM'))      
#
rm(VisCoMut, Input)

##### Figure 1C #####
##### 1.1 manually pick genes
Input.KeyMutation <- c('TP53','ELF3','EPHA2','ARID1A','KMT2C','STK11','ERBB3','SMAD4')
Input.Target <- c('TP53','CDK1','NCAM1','PCNA',     # TP53:CDK1,PCNA
                  'ELF3','CCL21','IDH1','FCGBP',  # ELF3:S100A7
                  'EPHA2','IDO1','EGFR','IGF2BP1',  # EPHA2:EGFR
                  'ARID1A','BAAT','SULT1A3','HMGB2', # ARID1A
                  'KMT2C','BCKDK','ITGA3','TMF1', # KMT2C
                  'STK11','S100P','TRIM29','HMGA2', # STK11
                  'ERBB3','PRG2','ITGA3','PRMT1', # ERBB3
                  'SMAD4','HNF4A','ITGB4','MMP28') # SMAD4:HNF4A
##### calc prep
OutputGene <- tibble(expand.grid(MutantGene = Input.KeyMutation,
                                 AlteredGene = Input.Target,
                                 Omics = c('mRNA','Protein'),
                                 LogFC = NA,
                                 Pval = NA))
OutputGene$MutantGene = as.character(OutputGene$MutantGene)
OutputGene$AlteredGene = as.character(OutputGene$AlteredGene)
OutputGene$Omics = as.character(OutputGene$Omics)

##### 1.2 Sample and mtx prep
# prep
SamplesWES = table(GBC_Main_Mutation$Tumor_Sample_Barcode) %>% names()
# sample name prophos
SamplePP = filter(GBC_Main_Index, Sample %in% SamplesWES)$ProPhos_T %>% na.omit() %>% as.character()
names(SamplePP) = filter(GBC_Main_Index, ProPhos_T %in% SamplePP)$Sample
# sample name rna
SampleRNA = filter(GBC_Main_Index, Sample %in% SamplesWES)$RNA_T %>% na.omit() %>% as.character()
names(SampleRNA) = filter(GBC_Main_Index, RNA_T %in% SampleRNA)$Sample
# 
InputRNA = GBC_Main_RNA_logTPM[,SampleRNA]
InputPro = GBC_Main_Pro[,SamplePP]

##### 1.3 Start calc
##### calc
### Gene
for (i in 1:nrow(OutputGene)) {
  # exist in omics?
  TempCurrentOmic = OutputGene$Omics[i]
  TempMutGene = OutputGene$MutantGene[i]; TempTargetGene = OutputGene$AlteredGene[i]
  #
  if (TempCurrentOmic == 'mRNA') {
    if (!TempTargetGene %in% rownames(GBC_Main_RNA_logTPM)) {
      OutputGene$LogFC[i] = NA; OutputGene$Pval[i] = NA
      next
    }
  } else {
    if (!TempTargetGene %in% rownames(GBC_Main_Pro)) {
      OutputGene$LogFC[i] = NA; OutputGene$Pval[i] = NA
      next
    }
  }
  ##### start calc
  if (TempCurrentOmic == 'mRNA') {
    TempMutSample = SampleRNA[names(SampleRNA) %in% filter(GBC_Main_Mutation, Hugo_Symbol == TempMutGene)$Tumor_Sample_Barcode] %>% na.omit() %>% as.character()
    TempWTSample = SampleRNA[!names(SampleRNA) %in% filter(GBC_Main_Mutation, Hugo_Symbol == TempMutGene)$Tumor_Sample_Barcode] %>% na.omit() %>% as.character()
    OutputGene$LogFC[i] = try(mean(as.numeric(InputRNA[TempTargetGene,TempMutSample])) - mean(as.numeric(InputRNA[TempTargetGene,TempWTSample])))
    OutputGene$Pval[i] = try(wilcox.test(as.numeric(InputRNA[TempTargetGene,TempMutSample]), 
                                         as.numeric(InputRNA[TempTargetGene,TempWTSample]))$p.value)
  }
  if (TempCurrentOmic == 'Protein') {
    TempMutSample = SamplePP[names(SamplePP) %in% filter(GBC_Main_Mutation, Hugo_Symbol == TempMutGene)$Tumor_Sample_Barcode] %>% na.omit() %>% as.character()
    TempWTSample = SamplePP[!names(SamplePP) %in% filter(GBC_Main_Mutation, Hugo_Symbol == TempMutGene)$Tumor_Sample_Barcode] %>% na.omit() %>% as.character()
    OutputGene$LogFC[i] = try(mean(as.numeric(InputPro[TempTargetGene,TempMutSample])) - mean(as.numeric(InputPro[TempTargetGene,TempWTSample])))
    OutputGene$Pval[i] = try(wilcox.test(as.numeric(InputPro[TempTargetGene,TempMutSample]), 
                                         as.numeric(InputPro[TempTargetGene,TempWTSample]))$p.value)
  }
  rm(i, TempMutSample, TempWTSample)
}
OutputGene$LogFC = as.numeric(OutputGene$LogFC); OutputGene$Pval = as.numeric(OutputGene$Pval)
# FDR
OutputGene$FDR <- p.adjust(OutputGene$Pval, method = 'fdr')
# rm empty rows
OutputGene <- na.omit(OutputGene)
##### data for Table S2
write.xlsx(OutputGene, file = 'CrucialGene Mut-RNA&Pro FC.xlsx',
           rowNames = T, colNames = T, overwrite = T)

##### 1.4 vis
OutputGene$Significant <- ifelse(OutputGene$FDR < 0.05, "Significant", "Not Significant")
OutputGene <- OutputGene %>%
  arrange(LogFC) %>%
  mutate(
    Rank = row_number(), 
    Significant = ifelse(FDR <= 0.05, TRUE, FALSE),  
    Size = -log10(FDR), 
    Label = ifelse(Significant, paste(MutantGene, AlteredGene, sep = "-"), NA)  
  )
ggplot(OutputGene, aes(x = Rank, y = LogFC, size = Size, color = Significant, alpha = Significant)) +
  geom_point() +
  geom_text_repel(aes(label = Label), size = 3.5, family = "Arial", color = "black", na.rm = TRUE) +  
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey")) + 
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4)) + 
  theme_minimal() +
  theme(
    text = element_text(family = "Arial", color = "black"), 
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black")
  ) +
  labs(
    title = "LogFC Ranking Plot",
    x = "Ranking (LogFC from smallest to largest)",
    y = "LogFC",
    size = expression(-log[10](FDR)),
    color = "Significant"
  )

# 
rm(Input.KeyMutation, Input.Target, TempCurrentOmic, TempMutGene, TempTargetGene,
   SamplePP, SampleRNA, SamplesWES, OutputGene, InputRNA, InputPro)

##### Figure 1D #####
##### 1.1 prep
# gene id prep
# Integrin α subunit
Genes.ITGA <- c("ITGA1", "ITGA2", "ITGA2B", "ITGA3", "ITGA4", "ITGA5", "ITGA6", "ITGA7", "ITGA8",
                "ITGA9", "ITGA10", "ITGA11", "ITGAD", "ITGAE", "ITGAL", "ITGAM", "ITGAV", "ITGAX")
# Integrin β subunit
Genes.ITGB <- c("ITGB1", "ITGB2", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8")
##### TP53 mutation and EPHA2 mutation conclusion
Index.comut <- tibble(Sample = names(table(GBC_Main_Maf@data$Tumor_Sample_Barcode)))
Index.comut$TP53 <- ifelse(Index.comut$Sample %in% filter(GBC_Main_Maf@data, Hugo_Symbol == 'TP53')$Tumor_Sample_Barcode, 'Mut', 'WT')
Index.comut$EPHA2 <- ifelse(Index.comut$Sample %in% filter(GBC_Main_Maf@data, Hugo_Symbol == 'EPHA2')$Tumor_Sample_Barcode, 'Mut', 'WT')
Index.comut$Condition <- paste0(Index.comut$TP53, '-', Index.comut$EPHA2)
Index.comut$Condition <- factor(Index.comut$Condition, levels = c('WT-WT','Mut-WT','WT-Mut','Mut-Mut'), ordered = T)
##### mtx arrangement
Mtx.rna <- InputRNA[c('TP53', 'EPHA2', Genes.ITGA, Genes.ITGB),]
Mtx.pro <- InputPro[c('TP53', 'EPHA2', Genes.ITGA, Genes.ITGB),]
rownames(Mtx.pro) <- c('TP53', 'EPHA2', Genes.ITGA, Genes.ITGB)
# change colnames
colnames(Mtx.rna) <- ifelse(
  !is.na(match(colnames(Mtx.rna), GBC_Main_Index$RNA_T)), 
  GBC_Main_Index$Sample[match(colnames(Mtx.rna), GBC_Main_Index$RNA_T)], 
  colnames(Mtx.rna)
)
colnames(Mtx.pro) <- ifelse(
  !is.na(match(colnames(Mtx.pro), GBC_Main_Index$ProPhos_T)), 
  GBC_Main_Index$Sample[match(colnames(Mtx.pro), GBC_Main_Index$ProPhos_T)], 
  colnames(Mtx.pro)
)
# transform
Mtx.rna <- t(Mtx.rna); Mtx.pro <- t(Mtx.pro)
# add condition
Mtx.rna <- aggregate(Mtx.rna, by = list(Condition = filter(Index.comut, Sample %in% rownames(Mtx.rna))$Condition),
                     mean, na.rm = T) %>% column_to_rownames(var = 'Condition') %>% t() %>% as.data.frame()
Mtx.pro <- aggregate(Mtx.pro, by = list(Condition = filter(Index.comut, Sample %in% rownames(Mtx.pro))$Condition),
                     mean, na.rm = T) %>% column_to_rownames(var = 'Condition') %>% t() %>% as.data.frame()
# scale
Mtx.rna <- t(scale(t(Mtx.rna))) %>% as.data.frame()
Mtx.pro <- t(scale(t(Mtx.pro))) %>% as.data.frame()

##### 1.2 Vis
# table to long
Mtx.long.rna <- Mtx.rna %>%
  rownames_to_column(var = "Gene") %>%  
  pivot_longer(cols = -Gene, names_to = "Group", values_to = "Expression") %>%
  mutate(Type = "RNA")
Mtx.long.pro <- Mtx.pro %>%
  rownames_to_column(var = "Gene") %>%  
  pivot_longer(cols = -Gene, names_to = "Group", values_to = "Expression") %>%
  mutate(Type = "Protein")
Mtx.long <- bind_rows(Mtx.long.rna, Mtx.long.pro)
# add index
Mtx.long <- Mtx.long %>%
  mutate(
    Gene_numeric = as.numeric(factor(Gene, levels = c('TP53', 'EPHA2', Genes.ITGA, Genes.ITGB),
                                     ordered = T)),  
    Gene_adjusted = ifelse(Type == "RNA", Gene_numeric - 0.25, Gene_numeric + 0.25)  
  )
##### data for Table S2
write.xlsx(Mtx.long, file = 'TP53-EPHA2 ITGN.xlsx',
           rowNames = T, colNames = T, overwrite = T)

##### Vis
ggplot(Mtx.long, aes(x = Gene_adjusted, y = Group, fill = Expression, shape = Type)) +
  geom_point(
    size = 4, stroke = 1,  
    color = "black" 
  ) +
  scale_fill_gradient2(low = ColColor$`Low-Blue`[5], mid = "white", high = ColColor$`High-Red`[9], midpoint = 0, na.value = "grey100") +  
  scale_shape_manual(values = c("RNA" = 24, "Protein" = 21)) +  
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") + 
  scale_x_continuous(
    breaks = unique(Mtx.long$Gene_numeric),  
    labels = unique(Mtx.long$Gene),  
    expand = c(0, 0.4)  
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    text = element_text(family = "Arial", color = "black"),  
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    panel.grid.major = element_line(color = "grey95", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.line = element_blank() 
  ) +
  labs(
    title = "Scatterplot of Mutant vs Altered Genes",
    x = "Altered Genes in mRNA and Protein Levels",
    y = "Mutant Genes",
    color = "LogFC"
  ) 
#
rm(Mtx.long, Mtx.long.pro, Mtx.long.rna, Mtx.pro, Mtx.rna, Genes.ITGA, Genes.ITGB, Index.comut)

##### Figure 1E #####
# Key Mutation: TP53, ELF3 
##### 1.1 DEA 
## prep
SamplesWES = table(GBC_Main_Mutation$Tumor_Sample_Barcode) %>% names()
# sample name prophos
SamplePP = filter(GBC_Main_Index, Sample %in% SamplesWES)$ProPhos_T %>% na.omit() %>% as.character()
names(SamplePP) = filter(GBC_Main_Index, ProPhos_T %in% SamplePP)$Sample
# sample name rna
SampleRNA = filter(GBC_Main_Index, Sample %in% SamplesWES)$RNA_T %>% na.omit() %>% as.character()
names(SampleRNA) = filter(GBC_Main_Index, RNA_T %in% SampleRNA)$Sample
# 
InputRNA = GBC_Main_RNA_logTPM[,SampleRNA]
InputPro = GBC_Main_Pro[,SamplePP]
InputPhos = GBC_Main_ProLevelPhos_knn[,SamplePP]
#
MutualGene = intersect(rownames(InputRNA), rownames(InputPro))

##### 1.2 DEA for pathway selection
### TP53
TempMutSampleRNA = SampleRNA[names(SampleRNA) %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'TP53')$Tumor_Sample_Barcode] %>% na.omit() %>% as.character()
TempWTSampleRNA = SampleRNA[!names(SampleRNA) %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'TP53')$Tumor_Sample_Barcode] %>% na.omit() %>% as.character()
TempMutSample = SamplePP[names(SamplePP) %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'TP53')$Tumor_Sample_Barcode] %>% na.omit() %>% as.character()
TempWTSample = SamplePP[!names(SamplePP) %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'TP53')$Tumor_Sample_Barcode] %>% na.omit() %>% as.character()
# 
OutputTP53 <- tibble(Gene = character(),
                     RNA_LogFC = numeric(), RNA_Padj = numeric(),
                     Pro_LogFC = numeric(), Pro_Padj = numeric(),
                     Phos_LogFC = numeric(), Phos_Padj = numeric())
for (i in 1:length(MutualGene)) {
  a <- tibble(Gene = MutualGene[i],
              RNA_LogFC = mean(as.numeric(InputRNA[MutualGene[i],TempMutSampleRNA])) - mean(as.numeric(InputRNA[MutualGene[i],TempWTSampleRNA])), 
              RNA_Padj = wilcox.test(as.numeric(InputRNA[MutualGene[i],TempMutSampleRNA]),
                                     as.numeric(InputRNA[MutualGene[i],TempWTSampleRNA]))$p.value,
              Pro_LogFC = mean(as.numeric(InputPro[MutualGene[i],TempMutSample])) - mean(as.numeric(InputPro[MutualGene[i],TempWTSample])), 
              Pro_Padj = wilcox.test(as.numeric(InputPro[MutualGene[i],TempMutSample]),
                                     as.numeric(InputPro[MutualGene[i],TempWTSample]))$p.value,
              Phos_LogFC = try(mean(as.numeric(InputPhos[MutualGene[i],TempMutSample])) - mean(as.numeric(InputPhos[MutualGene[i],TempWTSample]))), 
              Phos_Padj = try(wilcox.test(as.numeric(InputPhos[MutualGene[i],TempMutSample]),
                                          as.numeric(InputPhos[MutualGene[i],TempWTSample]))$p.value))
  OutputTP53 <- rbind(OutputTP53, a)
}
OutputTP53$Phos_Padj = as.numeric(OutputTP53$Phos_Padj)

### ELF3
TempMutSampleRNA = SampleRNA[names(SampleRNA) %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ELF3')$Tumor_Sample_Barcode] %>% na.omit() %>% as.character()
TempWTSampleRNA = SampleRNA[!names(SampleRNA) %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ELF3')$Tumor_Sample_Barcode] %>% na.omit() %>% as.character()
TempMutSample = SamplePP[names(SamplePP) %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ELF3')$Tumor_Sample_Barcode] %>% na.omit() %>% as.character()
TempWTSample = SamplePP[!names(SamplePP) %in% filter(GBC_Main_Mutation, Hugo_Symbol == 'ELF3')$Tumor_Sample_Barcode] %>% na.omit() %>% as.character()
# 
OutputELF3 <- tibble(Gene = character(),
                     RNA_LogFC = numeric(), RNA_Padj = numeric(),
                     Pro_LogFC = numeric(), Pro_Padj = numeric(),
                     Phos_LogFC = numeric(), Phos_Padj = numeric())
for (i in 1:length(MutualGene)) {
  a <- tibble(Gene = MutualGene[i],
              RNA_LogFC = mean(as.numeric(InputRNA[MutualGene[i],TempMutSampleRNA])) - mean(as.numeric(InputRNA[MutualGene[i],TempWTSampleRNA])), 
              RNA_Padj = wilcox.test(as.numeric(InputRNA[MutualGene[i],TempMutSampleRNA]),
                                     as.numeric(InputRNA[MutualGene[i],TempWTSampleRNA]))$p.value,
              Pro_LogFC = mean(as.numeric(InputPro[MutualGene[i],TempMutSample])) - mean(as.numeric(InputPro[MutualGene[i],TempWTSample])), 
              Pro_Padj = wilcox.test(as.numeric(InputPro[MutualGene[i],TempMutSample]),
                                     as.numeric(InputPro[MutualGene[i],TempWTSample]))$p.value,
              Phos_LogFC = try(mean(as.numeric(InputPhos[MutualGene[i],TempMutSample])) - mean(as.numeric(InputPhos[MutualGene[i],TempWTSample]))), 
              Phos_Padj = try(wilcox.test(as.numeric(InputPhos[MutualGene[i],TempMutSample]),
                                          as.numeric(InputPhos[MutualGene[i],TempWTSample]))$p.value))
  OutputELF3 <- rbind(OutputELF3, a)
}
OutputELF3$Phos_Padj = as.numeric(OutputELF3$Phos_Padj)

##### 1.3 ORA for pathway selection
##### Up in TP53
# HALLMARK_INTERFERON_GAMMA_RESPONSE
# KEGG DNA replication
# REACTOME Antigen Presentation
# KEGG Phagosome
# HALLMARK_MYC_TARGETS_V2
ORAresultTP53_Up <- clusterProfiler::enricher(filter(OutputTP53, Pro_LogFC > 0, Pro_Padj <= 0.05)$Gene,
                                              TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.05)
##### Down in TP53
# KEGG ErbB signaling pathway
# HALLMARK_PI3K_AKT_MTOR_SIGNALING
# HALLMARK_APICAL_JUNCTION
# REACTOME Nucleotide Excision Repair
ORAresultTP53_Down <- clusterProfiler::enricher(filter(OutputTP53, Pro_LogFC < 0, Pro_Padj <= 0.05)$Gene,
                                                TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.05)

##### Up in ELF3
# HALLMARK_OXIDATIVE_PHOSPHORYLATION
# HALLMARK_FATTY_ACID_METABOLISM
# HALLMARK_GLYCOLYSIS
ORAresultELF3_Up <- clusterProfiler::enricher(filter(OutputELF3, Pro_LogFC > 0, Pro_Padj <= 0.05)$Gene,
                                              TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.05)

##### 1.4 Heatmap for schematic plot
HMInputGene <- rbind(filter(OutputTP53, Gene %in% filter(TotalPathway, Pathway == 'HALLMARK_INTERFERON_GAMMA_RESPONSE')$Gene,
                            Pro_Padj <= 0.05) %>%
                       mutate(Belong = 'TP53 Up', Pathway = 'IFNy Signaling'),
                     filter(OutputTP53, Gene %in% filter(TotalPathway, Pathway == 'KEGG DNA replication')$Gene,
                            Pro_Padj <= 0.05) %>%
                       mutate(Belong = 'TP53 Up', Pathway = 'Cell Cycle'),
                     filter(OutputTP53, Gene %in% filter(TotalPathway, Pathway == 'KEGG Phagosome')$Gene,
                            Pro_Padj <= 0.05) %>%
                       mutate(Belong = 'TP53 Up', Pathway = 'Phagosome'),
                     filter(OutputTP53, Gene %in% filter(TotalPathway, Pathway == 'HALLMARK_MYC_TARGETS_V2')$Gene,
                            Pro_Padj <= 0.05) %>%
                       mutate(Belong = 'TP53 Up', Pathway = 'MYC Targets'),
                     filter(OutputTP53, Gene %in% filter(TotalPathway, Pathway == 'KEGG ErbB signaling pathway')$Gene,
                            Pro_Padj <= 0.05) %>%
                       mutate(Belong = 'TP53 Down', Pathway = 'ERBB Signaling'),
                     filter(OutputTP53, Gene %in% filter(TotalPathway, Pathway == 'HALLMARK_PI3K_AKT_MTOR_SIGNALING')$Gene,
                            Pro_Padj <= 0.05) %>%
                       mutate(Belong = 'TP53 Down', Pathway = 'PI3K-AKT Signaling'),
                     filter(OutputTP53, Gene %in% filter(TotalPathway, Pathway == 'HALLMARK_APICAL_JUNCTION')$Gene,
                            Pro_Padj <= 0.05) %>%
                       mutate(Belong = 'TP53 Down', Pathway = 'Cell Adhasion'),
                     filter(OutputTP53, Gene %in% filter(TotalPathway, Pathway == 'REACTOME Nucleotide Excision Repair')$Gene,
                            Pro_Padj <= 0.05) %>%
                       mutate(Belong = 'TP53 Down', Pathway = 'DNA Repair'),
                     filter(OutputELF3, Gene %in% filter(TotalPathway, Pathway == 'HALLMARK_OXIDATIVE_PHOSPHORYLATION')$Gene,
                            Pro_Padj <= 0.05) %>%
                       mutate(Belong = 'ELF3 Up', Pathway = 'OXPHOS'),
                     filter(OutputELF3, Gene %in% filter(TotalPathway, Pathway == 'HALLMARK_FATTY_ACID_METABOLISM')$Gene,
                            Pro_Padj <= 0.05) %>%
                       mutate(Belong = 'ELF3 Up', Pathway = 'FFA Metabolism'),
                     filter(OutputELF3, Gene %in% filter(TotalPathway, Pathway == 'HALLMARK_GLYCOLYSIS')$Gene,
                            Pro_Padj <= 0.05) %>%
                       mutate(Belong = 'ELF3 Up', Pathway = 'Glycolysis')
) %>% as_tibble()
##### data for Table S2
write.xlsx(HMInputGene, file = 'Schematic data.xlsx',
           rowNames = T, colNames = T, overwrite = T)

##### 1.5 Vis
HMInput = HMInputGene[,c(1,2,4,6)] %>% distinct(Gene, .keep_all = T) %>% column_to_rownames(var = 'Gene')
# 
Heatmap(HMInput, cluster_rows = F, cluster_columns = F,
        col = colorRamp2(c(-1.5,0,1.5), c(ColColor$`Low-Indigo`[8], 'white', ColColor$`High-Red`[9])),
        height = unit(60,'cm'), width = unit(4,'cm')) # 

#
rm(HMInput, HMInputGene, a,
   InputPhos, InputPro, InputRNA, 
   ORAresult, ORAresultELF3_Down, ORAresultELF3_Up, ORAresultTP53_Down, ORAresultTP53_Up,
   Output, OutputELF3, OutputTP53, i,
   InputCNVPeak, InputGene, MutualGene, SamplePP, SampleRNA, SamplesWES, TempMutSample, TempMutSampleRNA,
   TempWTSample, TempWTSampleRNA)

##### Figure 1F was generated by GISTIC2 #####

##### Figure 1G #####
##### 1.0 prep
# CNV peak
InputCNVPeak = fread('/Users/fuzile/Desktop/たいようけい/我参与的项目/主要参与/2023.08 GBC MO/1 原始数据/原始数据整理/0905/WES_CNV/amp_genes.conf_90.txt')
InputCNVPeak = InputCNVPeak[-c(1:4),-1]
InputCNVPeak = as.matrix(InputCNVPeak) %>% as.character()
InputCNVPeak = InputCNVPeak[InputCNVPeak != ''] %>% na.omit() %>% as.character()
InputCNVPeak
# mutual sample
MutualSample = filter(GBC_Main_Index, !is.na(WES_T), !is.na(RNA_T), !is.na(ProPhos_T))

##### 1.1 Gene selection
### Step 1. Mutual: 6299
InputGene = intersect(rownames(GBC_Main_CNV), rownames(GBC_Main_RNA_logTPM)) %>%
  intersect(rownames(GBC_Main_Pro))
### Step 2. In Amp peak: 198
InputGene = intersect(InputGene, InputCNVPeak)
### Step 3. calc CNV-RNA, CNV-Protein
Output = tibble(Gene = character(),
                CorRNA = numeric(),
                CorRNAPval = numeric(),
                CorPro = numeric(),
                CorProPval = numeric()
)
for (i in 1:length(InputGene)){
  a = tibble(Gene = InputGene[i],
             CorRNA = cor(as.numeric(GBC_Main_CNV[InputGene[i], MutualSample$WES_T]),
                          as.numeric(GBC_Main_RNA_logTPM[InputGene[i],MutualSample$RNA_T])),
             CorRNAPval = cor.test(as.numeric(GBC_Main_CNV[InputGene[i], MutualSample$WES_T]),
                                   as.numeric(GBC_Main_RNA_logTPM[InputGene[i],MutualSample$RNA_T]), method = 'spearman')$p.value,
             CorPro = cor(as.numeric(GBC_Main_CNV[InputGene[i], MutualSample$WES_T]),
                          as.numeric(GBC_Main_Pro[InputGene[i],MutualSample$ProPhos_T])),
             CorProPval = cor.test(as.numeric(GBC_Main_CNV[InputGene[i], MutualSample$WES_T]),
                                   as.numeric(GBC_Main_Pro[InputGene[i],MutualSample$ProPhos_T]), method = 'spearman')$p.value
  )
  Output = rbind(Output, a)
}
### Step 4. judge
## direction
Output$Direction <- ifelse(Output$CorRNA * Output$CorPro > 0, 'Yes', 'No')
table(Output$Direction) # 174
## significancy
Output$SigRNA = ifelse(Output$CorRNAPval <= 0.05, 'Yes', 'No') # 155
table(Output$SigRNA)
Output$SigPro = ifelse(Output$CorProPval <= 0.05, 'Yes', 'No') # 110
table(Output$SigPro)
Output$SigBoth = ifelse(Output$SigRNA == 'Yes' & Output$SigPro == 'Yes', 'Yes', 'No')
table(Output$SigBoth) # 106
## direaction and significancy
Output$DicAndSig = ifelse(Output$Direction == 'Yes' & Output$SigBoth == 'Yes', 'Yes', 'No')
table(Output$DicAndSig) # 106

##### 1.2 enrichment
ORAresult = clusterProfiler::enricher(filter(Output, DicAndSig == 'Yes')$Gene, TERM2GENE = TotalPathway)@result %>%
  # select(ID, pvalue, GeneRatio) %>%
  filter(pvalue <= 0.05) %>% mutate(LogFDR = -log10(pvalue)) %>%
  arrange(desc(LogFDR))
# data for Table S2
write.xlsx(ORAresult, file = 'CNV cis ORA.xlsx',
           rowNames = T, colNames = T, overwrite = T)

# vis prep
ORAresult = clusterProfiler::enricher(filter(Output, DicAndSig == 'Yes')$Gene, TERM2GENE = TotalPathway)@result %>%
  select(ID, pvalue, GeneRatio) %>%
  filter(pvalue <= 0.05) %>% mutate(LogFDR = -log10(pvalue)) %>%
  arrange(desc(LogFDR))
ORAresultVis = ORAresult[c('REACTOME Signaling by FGFR1 in disease',
                        'KEGG Central carbon metabolism in cancer',
                        'REACTOME MET activates RAS signaling',
                        'REACTOME PI3K/AKT Signaling in Cancer',
                        'REACTOME Signaling by ERBB2',
                        'HALLMARK_MYC_TARGETS_V2',
                        'KEGG Adherens junction',
                        'HALLMARK_E2F_TARGETS') %>% rev(),]

##### 1.3 Vis
ggbarplot(ORAresultVis, 'ID', 'LogFDR', color = 'brown', fill = 'brown',
          position = position_dodge()) +
  coord_flip() + xlab('') + ylab('-LogPval') +
  scale_color_manual(values = c('Upregulated' = ColColor$`High-OrangeRed`[5], 'Downregulated' = ColColor$`Low-Indigo`[5])) +
  scale_fill_manual(values = c('Upregulated' = ColColor$`High-OrangeRed`[5], 'Downregulated' = ColColor$`Low-Indigo`[5])) +
  scale_y_continuous(expand = c(0,0)) 
#
rm(InputCNVPeak, MutualSample, ORAresult, ORAresultVis, i, InputGene)

##### Figure S1H #####
##### 1.1 DEA
# Index
IndexTvsN.rna <- filter(index.raw, Type != 'Benign') %>% select(1,4,7) %>%
  filter(!is.na(SampleRNA))
IndexTvsN.pro <- filter(index.raw, Type != 'Benign') %>% select(1,5,7) %>%
  filter(!is.na(SamplePro))
# Cluster matrix
CM.rna <- data.frame(Tumor = ifelse(IndexTvsN.rna$Type == 'Cancer', 1, 0),
                     NAT = ifelse(IndexTvsN.rna$Type == 'Adjacent', 1, 0))
CM.pro <- data.frame(Tumor = ifelse(IndexTvsN.pro$Type == 'Cancer', 1, 0),
                     NAT = ifelse(IndexTvsN.pro$Type == 'Adjacent', 1, 0))
# DEA
DEA_RNA <- core_Differential_analysis_continuous(CM.rna, GBC_Main_RNA_logTPM[,IndexTvsN.rna$SampleRNA],
                                                 log = T, p.adj = T, method = 'Wilcox', show.belong = T)
DEA_Pro <- core_Differential_analysis_continuous(CM.pro, GBC_Main_Pro[,IndexTvsN.pro$SamplePro],
                                                 log = T, p.adj = T, method = 'Wilcox', show.belong = T)
DEA_Phos <- core_Differential_analysis_continuous(CM.pro, GBC_Main_ProLevelPhos_knn[,IndexTvsN.pro$SamplePro],
                                                  log = T, p.adj = T, method = 'Wilcox', show.belong = T)
table(DEA_RNA$belong); table(DEA_Pro$belong); table(DEA_Phos$belong)

##### 1.2 ORA (LogFC > 0.5)
# Data for Table S2
# Pro
ORA_up.pro <- clusterProfiler::enricher(filter(DEA_Pro, belong == 1, `Tumor-logFC` >= 0.5) %>% rownames(), TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.01) %>% arrange(p.adjust) %>% mutate(Belong = 'Tumor', Omic = 'Protein')
ORA_down.pro <- clusterProfiler::enricher(filter(DEA_Pro, belong == 2, `NAT-logFC` >= 0.5) %>% rownames(), TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.01) %>% arrange(p.adjust) %>% mutate(Belong = 'NAT', Omic = 'Protein')
ORA.pro <- rbind(ORA_up.pro[1:6, ], ORA_down.pro[1:6,] %>% arrange(desc(p.adjust)))
# RNA
ORA_up.rna <- clusterProfiler::enricher(filter(DEA_RNA, belong == 1, `Tumor-logFC` >= 0.5) %>% rownames(), TERM2GENE = TotalPathway, pvalueCutoff = 1)@result %>%
  filter(ID %in% ORA_up.pro[1:6, ]$ID) %>% arrange(p.adjust) %>% mutate(Belong = 'Tumor', Omic = 'mRNA')
ORA_down.rna <- clusterProfiler::enricher(filter(DEA_RNA, belong == 2, `NAT-logFC` >= 0.5) %>% rownames(), TERM2GENE = TotalPathway, pvalueCutoff = 1)@result %>%
  filter(ID %in% ORA_down.pro[1:6,]$ID) %>% arrange(p.adjust) %>% mutate(Belong = 'NAT', Omic = 'mRNA')
ORA.rna <- rbind(ORA_up.rna, ORA_down.rna %>% arrange(desc(p.adjust)))
# Phos
ORA_up.phos <- clusterProfiler::enricher(filter(DEA_Phos, belong == 1, `Tumor-logFC` >= 0.5) %>% rownames(), TERM2GENE = TotalPathway)@result %>%
  filter(ID %in% ORA_up.pro[1:6, ]$ID) %>% arrange(p.adjust) %>% mutate(Belong = 'Tumor', Omic = 'Protein phosphorylation')
ORA_down.phos <- clusterProfiler::enricher(filter(DEA_Phos, belong == 2, `NAT-logFC` >= 0.5) %>% rownames(), TERM2GENE = TotalPathway)@result %>%
  filter(ID %in% ORA_down.pro[1:6,]$ID) %>% arrange(p.adjust) %>% mutate(Belong = 'NAT', Omic = 'Protein phosphorylation')
ORA.phos <- rbind(ORA_up.phos, ORA_down.phos %>% arrange(desc(p.adjust)))
## Integ
ORA <- rbind(ORA.rna, ORA.pro, ORA.phos)
# 
write.xlsx(ORA, file = 'ORA TvsN.xlsx', rowNames = T, colNames = T, overwrite = T)

# Pro
ORA_up.pro <- clusterProfiler::enricher(filter(DEA_Pro, belong == 1, `Tumor-logFC` >= 0.5) %>% rownames(), TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.01) %>% select(ID, p.adjust) %>% arrange(p.adjust) %>% mutate(Belong = 'Tumor', Omic = 'Protein')
ORA_down.pro <- clusterProfiler::enricher(filter(DEA_Pro, belong == 2, `NAT-logFC` >= 0.5) %>% rownames(), TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.01) %>% select(ID, p.adjust) %>% arrange(p.adjust) %>% mutate(Belong = 'NAT', Omic = 'Protein')
ORA.pro <- rbind(ORA_up.pro[1:6, ], ORA_down.pro[1:6,] %>% arrange(desc(p.adjust)))
# RNA
ORA_up.rna <- clusterProfiler::enricher(filter(DEA_RNA, belong == 1, `Tumor-logFC` >= 0.5) %>% rownames(), TERM2GENE = TotalPathway, pvalueCutoff = 1)@result %>%
  filter(ID %in% ORA_up.pro[1:6, ]$ID) %>% select(ID, p.adjust) %>% arrange(p.adjust) %>% mutate(Belong = 'Tumor', Omic = 'mRNA')
ORA_down.rna <- clusterProfiler::enricher(filter(DEA_RNA, belong == 2, `NAT-logFC` >= 0.5) %>% rownames(), TERM2GENE = TotalPathway, pvalueCutoff = 1)@result %>%
  filter(ID %in% ORA_down.pro[1:6,]$ID) %>% select(ID, p.adjust) %>% arrange(p.adjust) %>% mutate(Belong = 'NAT', Omic = 'mRNA')
ORA.rna <- rbind(ORA_up.rna, ORA_down.rna %>% arrange(desc(p.adjust)))
# Phos
ORA_up.phos <- clusterProfiler::enricher(filter(DEA_Phos, belong == 1, `Tumor-logFC` >= 0.5) %>% rownames(), TERM2GENE = TotalPathway)@result %>%
  filter(ID %in% ORA_up.pro[1:6, ]$ID) %>% select(ID, p.adjust) %>% arrange(p.adjust) %>% mutate(Belong = 'Tumor', Omic = 'Protein phosphorylation')
ORA_down.phos <- clusterProfiler::enricher(filter(DEA_Phos, belong == 2, `NAT-logFC` >= 0.5) %>% rownames(), TERM2GENE = TotalPathway)@result %>%
  filter(ID %in% ORA_down.pro[1:6,]$ID) %>% select(ID, p.adjust) %>% arrange(p.adjust) %>% mutate(Belong = 'NAT', Omic = 'Protein phosphorylation')
ORA.phos <- rbind(ORA_up.phos, ORA_down.phos %>% arrange(desc(p.adjust)))
## Integ
ORA <- rbind(ORA.rna, ORA.pro, ORA.phos)
# retain pathway
Pathway.select <- ORA.pro$ID
#
rm(ORA_up.rna, ORA_down.rna, ORA_up.pro, ORA_down.pro, ORA_up.phos, ORA_down.phos)

##### 1.3 Vis multi-omic functions
# prep
ORA <- mutate(ORA, LogFDR = ifelse(Belong == 'Tumor', -log10(p.adjust), log10(p.adjust)))
ORA$Omic = factor(ORA$Omic, levels = c('mRNA','Protein','Protein phosphorylation'), ordered = T)
ORA$ID = factor(ORA.pro$ID, levels = Pathway.select, ordered = T)
# vis
ggbarplot(ORA, 'ID', 'LogFDR', color = 'Omic', fill = 'Omic', order = Pathway.select,
          position = position_dodge()) +
  coord_flip() + xlab('') + ylab('-LogFDR') +
  scale_color_manual(values = c('mRNA' = '#F1A42A', 'Protein' = '#4E185F', 'Protein phosphorylation' = '#20706C')) +
  scale_fill_manual(values = c('mRNA' = '#F1A42A', 'Protein' = '#4E185F', 'Protein phosphorylation' = '#20706C')) +
  scale_y_continuous(expand = c(0,0)) 
#
rm(ORA, OutputKSEA, InputKSEA, Pathway.select)

##### Figure S1I #####
##### 1.1 HCC Input
# Ref: Gao Q, Zhu H, Dong L, et al. Integrated Proteogenomic Characterization of HBV-Related Hepatocellular Carcinoma [published correction appears in Cell. 2019 Nov 14;179(5):1240. doi: 10.1016/j.cell.2019.10.038.]. Cell. 2019;179(2):561-577.e22. doi:10.1016/j.cell.2019.08.052
load('/Users/fuzile/Desktop/たいようけい/基本文件/常用外部队列数据/CPTAC HCC/CPTAC组学数据.RData')
rm(CPTAC_HCC_cnv_matrix, CPTAC_HCC_mutation_long,CPTAC_HCC_mutation_matrix,CPTAC_HCC_phos,
   CPTAC_HCC_proteomics,CPTAC_HCC_RNA_count,CPTAC_HCC_RNA_TPM,CPTAC_HCC_survival)

##### 1.2 GBC Input
InputGBC <- GBC_Main_RNA_logTPM[,c(filter(index.raw, Type == 'Cancer')$SampleRNA %>% na.omit(),
                                   filter(index.raw, Type != 'Cancer')$SampleRNA %>% na.omit())]
CM <- data.frame(Tumor = ifelse(colnames(InputGBC) %in% filter(index.raw, Type == 'Cancer')$SampleRNA %>% na.omit(),1,0),
                 NAT = ifelse(colnames(InputGBC) %in% filter(index.raw, Type != 'Cancer')$SampleRNA %>% na.omit(),1,0))
OutputGBC <- bioinfoamateur::core_Differential_analysis_continuous(CM, InputGBC, log = T, method = 't.test', show.belong = F) %>%
  arrange(desc(`Tumor-logFC`))
OutputGBC <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = rownames(OutputGBC),
                                                                  LogFC = OutputGBC$`Tumor-logFC`))

##### 1.3 CCA Input
# Ref: Deng M, Ran P, Chen L, et al. Proteogenomic characterization of cholangiocarcinoma. Hepatology. 2023;77(2):411-429. doi:10.1002/hep.32624
# iCCA
load('/Users/fuzile/Desktop/たいようけい/基本文件/常用外部队列数据/Hepatology CCA/HepCCA.RData')
rm(HepCCA_Phos, HepCCA_Pro, HepCCA_ProLevelPhos)
# eCCA
InputeCCA <- HepCCA_logRNA[,c(filter(HepCCA_Index, Type == 'eCCA')$RNA_T %>% na.omit(),
                              filter(HepCCA_Index, Type == 'eCCA')$RNA_P %>% na.omit())]
CM <- data.frame(Tumor = ifelse(colnames(InputeCCA) %in% filter(HepCCA_Index, Type == 'eCCA')$RNA_T %>% na.omit(),1,0),
                 NAT = ifelse(colnames(InputeCCA) %in% filter(HepCCA_Index, Type == 'eCCA')$RNA_P %>% na.omit(),1,0))
OutputeCCA <- bioinfoamateur::core_Differential_analysis_continuous(CM, InputeCCA, log = T, method = 't.test', show.belong = F) %>%
  arrange(desc(`Tumor-logFC`))
OutputeCCA <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = rownames(OutputeCCA),
                                                                   LogFC = OutputeCCA$`Tumor-logFC`))
# iCCA
InputiCCA <- HepCCA_logRNA[,c(filter(HepCCA_Index, Type == 'iCCA')$RNA_T %>% na.omit(),
                              filter(HepCCA_Index, Type == 'iCCA')$RNA_P %>% na.omit())]
CM <- data.frame(Tumor = ifelse(colnames(InputiCCA) %in% filter(HepCCA_Index, Type == 'iCCA')$RNA_T %>% na.omit(),1,0),
                 NAT = ifelse(colnames(InputiCCA) %in% filter(HepCCA_Index, Type == 'iCCA')$RNA_P %>% na.omit(),1,0))
OutputiCCA <- bioinfoamateur::core_Differential_analysis_continuous(CM, InputiCCA, log = T, method = 't.test', show.belong = F) %>%
  arrange(desc(`Tumor-logFC`))
OutputiCCA <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = rownames(OutputiCCA),
                                                                   LogFC = OutputiCCA$`Tumor-logFC`))

##### 1.4 HCC Input
InputHCC <- CPTAC_HCC_RNA_logTPM
CM <- data.frame(Tumor = ifelse(grepl('T',colnames(CPTAC_HCC_RNA_logTPM)), 1, 0),
                 NAT = ifelse(grepl('N',colnames(CPTAC_HCC_RNA_logTPM)), 1, 0))
OutputHCC <- bioinfoamateur::core_Differential_analysis_continuous(CM, InputHCC, log = T, method = 't.test', show.belong = F) %>%
  arrange(desc(`Tumor-logFC`))
OutputHCC <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = rownames(OutputHCC),
                                                                  LogFC = OutputHCC$`Tumor-logFC`))

##### 1.5 GSEA and pathway selection
GSEAGBC <- clusterProfiler::GSEA(OutputGBC, TERM2GENE = filter(TotalPathway, grepl('KEGG',Pathway)))@result %>%
  select(ID, NES, p.adjust) %>% filter(NES > 0, p.adjust <= 0.01) %>% arrange(desc(NES))
GSEAeCCA <- clusterProfiler::GSEA(OutputeCCA, TERM2GENE = filter(TotalPathway, grepl('KEGG',Pathway)))@result %>%
  select(ID, NES, p.adjust) %>% filter(NES > 0, p.adjust <= 0.01) %>% arrange(desc(NES))
GSEAiCCA <- clusterProfiler::GSEA(OutputiCCA, TERM2GENE = filter(TotalPathway, grepl('KEGG',Pathway)))@result %>%
  select(ID, NES, p.adjust) %>% filter(NES > 0, p.adjust <= 0.01) %>% arrange(desc(NES))
GSEAHCC <- clusterProfiler::GSEA(OutputHCC, TERM2GENE = filter(TotalPathway, grepl('KEGG',Pathway)))@result %>%
  select(ID, NES, p.adjust) %>% filter(NES > 0, p.adjust <= 0.01) %>% arrange(desc(NES))
# Pathway selection
PathwaySelected <- c('KEGG Cell cycle','KEGG Spliceosome','KEGG Proteasome','KEGG Base excision repair',
                     'KEGG Ribosome','KEGG Mucin type O-glycan biosynthesis','KEGG Oxidative phosphorylation','KEGG Endocytosis',
                     'KEGG ECM-receptor interaction','KEGG PI3K-Akt signaling pathway','KEGG Hippo signaling pathway','KEGG Cellular senescence')

##### 1.6 re-GSEA and visualization
# data for Table S2
GSEAGBC <- clusterProfiler::GSEA(OutputGBC, TERM2GENE = filter(TotalPathway, grepl('KEGG',Pathway)), pvalueCutoff = 1)@result %>%
  filter(ID %in% PathwaySelected) %>% arrange(desc(NES)) %>% mutate(Belong = 'GBC')
GSEAeCCA <- clusterProfiler::GSEA(OutputeCCA, TERM2GENE = filter(TotalPathway, grepl('KEGG',Pathway)), pvalueCutoff = 1)@result %>%
  filter(ID %in% PathwaySelected) %>% arrange(desc(NES)) %>% mutate(Belong = 'eCCA')
GSEAiCCA <- clusterProfiler::GSEA(OutputiCCA, TERM2GENE = filter(TotalPathway, grepl('KEGG',Pathway)), pvalueCutoff = 1)@result %>%
  filter(ID %in% PathwaySelected) %>% arrange(desc(NES)) %>% mutate(Belong = 'iCCA')
GSEAHCC <- clusterProfiler::GSEA(OutputHCC, TERM2GENE = filter(TotalPathway, grepl('KEGG',Pathway)), pvalueCutoff = 1)@result %>%
  filter(ID %in% PathwaySelected) %>% arrange(desc(NES)) %>% mutate(Belong = 'HCC')
Input <- rbind(GSEAGBC, GSEAeCCA, GSEAiCCA, GSEAHCC) %>% as_tibble() 
write.xlsx(Input, file = 'Hep-Bil tumor func.xlsx',
           rowNames = T, colNames = T, overwrite = T)
#
GSEAGBC <- clusterProfiler::GSEA(OutputGBC, TERM2GENE = filter(TotalPathway, grepl('KEGG',Pathway)), pvalueCutoff = 1)@result %>%
  select(ID, NES, p.adjust) %>% filter(ID %in% PathwaySelected) %>% arrange(desc(NES)) %>% mutate(Belong = 'GBC')
GSEAeCCA <- clusterProfiler::GSEA(OutputeCCA, TERM2GENE = filter(TotalPathway, grepl('KEGG',Pathway)), pvalueCutoff = 1)@result %>%
  select(ID, NES, p.adjust) %>% filter(ID %in% PathwaySelected) %>% arrange(desc(NES)) %>% mutate(Belong = 'eCCA')
GSEAiCCA <- clusterProfiler::GSEA(OutputiCCA, TERM2GENE = filter(TotalPathway, grepl('KEGG',Pathway)), pvalueCutoff = 1)@result %>%
  select(ID, NES, p.adjust) %>% filter(ID %in% PathwaySelected) %>% arrange(desc(NES)) %>% mutate(Belong = 'iCCA')
GSEAHCC <- clusterProfiler::GSEA(OutputHCC, TERM2GENE = filter(TotalPathway, grepl('KEGG',Pathway)), pvalueCutoff = 1)@result %>%
  select(ID, NES, p.adjust) %>% filter(ID %in% PathwaySelected) %>% arrange(desc(NES)) %>% mutate(Belong = 'HCC')
Input <- rbind(GSEAGBC, GSEAeCCA, GSEAiCCA, GSEAHCC) %>% as_tibble() 
# factorize
Input$Belong = factor(Input$Belong, levels = c('GBC','eCCA','iCCA','HCC'), ordered = T)
Input$ID = factor(Input$ID, levels = PathwaySelected, ordered = T)
Input$Significance <- -log10(Input$p.adjust)
# Vis
InputVal <- dcast(Input, ID ~ Belong, value.var = 'NES') %>% column_to_rownames(var = 'ID') %>%
  t() %>% scale() %>% t()
InputSig <- dcast(Input, ID ~ Belong, value.var = 'Significance') %>% column_to_rownames(var = 'ID')

Heatmap(InputVal, cluster_rows = F, cluster_columns = F, name = 'NES',
        row_split = c(rep(1,4),rep(2,4),rep(3,4)), row_title = NULL, column_names_rot = 30, column_names_centered = T,
        col = circlize::colorRamp2(c(-1, 0, 1), c(ColColor$`Low-Indigo`[5], "white", ColColor$`High-Red`[5])),
        rect_gp = gpar(col = "white", lwd = 2),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(InputSig[i, j] > -log10(0.05))
            grid.text('*', x, y, gp = gpar(fontsize = 15))
        },
        height = unit(12,'cm'), width = unit(4, 'cm')) %>% draw(heatmap_legend_side = 'left') # 10*10 SF1G TvsN功能
#
rm(CM, CPTAC_HCC_clinical, CPTAC_HCC_RNA_logTPM, GSEAeCCA, GSEAGBC, GSEAHCC, GSEAiCCA,
   HepCCA_Index, HepCCA_Clinical, HepCCA_RNA, HepCCA_logRNA,
   Input, InputeCCA, InputGBC, InputHCC, InputiCCA, InputSig, InputVal,
   Output, OutputeCCA, OutputiCCA, OutputGBC, OutputHCC, PathwaySelected)




