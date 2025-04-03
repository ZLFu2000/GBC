##### Figure 4 and Figure S4
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
load('./GBC_LM_exploration_cohort.RData')
load('./GBC_LM_validation_cohort.RData')
# data for Table S5
write.xlsx(GBC_LMv_Clinical, file = 'LMv clinical.xlsx', rowNames = T, colNames = T, overwrite = T)
write.xlsx(GBC_LMv_Pro, file = 'LMv pro.xlsx', rowNames = T, colNames = T, overwrite = T)
write.xlsx(GBC_LMv_Phos_knn, file = 'LMv phos.xlsx', rowNames = T, colNames = T, overwrite = T)


#### 2 color
load('../Colors (ggsci).RData')
library(wesanderson)

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
Ref.genome <- rtracklayer::import('./Homo_sapiens.GRCh38.108.gtf') %>%
  as.data.frame() %>% select(1,2,3,7,12) %>% distinct(seqnames, gene_name, .keep_all = T) %>% na.omit()
Ref.genome <- mutate(Ref.genome, seqnames = paste0('chr', seqnames) %>% factor(levels = levels(Chr$Chr), ordered = T)) %>%
  na.omit() %>% arrange(seqnames, start)
Ref.genome = filter(Ref.genome, )
Ref.genome$ChrBand = apply(Ref.genome, 1, function(x){
  CurrentBand = filter(Chr, Chr == as.character(x[1]), Start <= as.numeric(x[2]), 
                       End >= as.numeric(x[3]))$ChrBand 
  CurrentBand = ifelse(length(CurrentBand) == 0, NA, as.character(CurrentBand))
}) %>%
  unlist() %>% as.character()
Ref.genome = filter(Ref.genome, ChrBand != '')
Ref.genome = as_tibble(Ref.genome)


##### prep2. index #####
index <- tibble(Patient = rep('A',21),
                Condition = rep('A',21))
index$Patient <- rep(GBC_LMe_Index$Sample[1:7],3)
index$Condition <- c(rep('NAT',7),rep('Cancer Insitu',7),rep('Liver invasion',7))
index$Sample <- c(GBC_LMe_Index$ProPhos_P, GBC_LMe_Index$ProPhos_T, GBC_LMe_Index$ProPhos_LM)
index = na.omit(index)
### primary
index.prim <- rbind(select(GBC_Main_Index, 1,2,4,6) %>% 
                      bioinfoamateur::dfm_change_colnames(c('Sample','WES','RNA','Pro')), 
                    mutate(GBC_Main_Index, Sample = paste0(Sample, '_P')) %>% select(1,3,5,7) %>% 
                      bioinfoamateur::dfm_change_colnames(c('Sample','WES','RNA','Pro')))
index.prim$Tumor <- ifelse(!grepl(fixed('_P'),index.prim$Sample), 'Yes', NA)
# retain samples having proteomic data
index.prim <- filter(index.prim, !is.na(Pro), Tumor == 'Yes')
index.prim$LM <- sapply(index.prim$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Regional_invasion}) %>% as.numeric()
index.prim$LM = ifelse(index.prim$LM == 1, 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)

##### Figure 4A #####
library(ggsurvfit); library(survival)
#
survfit2(Surv(OS_days, OS_event) ~ Regional_invasion, data = GBC_Main_Clinical) |> # Input
  ggsurvfit(linewidth = 1) + ggtitle('Invasion') +
  scale_color_manual(values = c('#82AC7C', 'red')) +
  add_quantile(y_value = 0.6, color = "gray50", linewidth = 0.75) + 
  add_pvalue(caption = "Log-rank {p.value}", location  = "annotation",
             x = 200, y = 0.25)

##### Figure 4B #####
InputMaf_LI <- maftools::read.maf(maf = filter(GBC_Main_Mutation, Tumor_Sample_Barcode %in% filter(GBC_Main_Clinical,Regional_invasion==1)$Patient_ID))
InputMaf_noLI <- maftools::read.maf(maf = filter(GBC_Main_Mutation, Tumor_Sample_Barcode %in% filter(GBC_Main_Clinical,Regional_invasion==0)$Patient_ID))

##### Calc
Output <- mafCompare(m1 = InputMaf_LI, m2 = InputMaf_noLI, 
                     m1Name = 'LI', m2Name = 'non-LI', 
                     minMut = 5)
print(Output)
forestPlot(mafCompareRes = Output, pVal = 0.05) 
# Data for Table S5
write.xlsx(Output$results, file = 'Diff Mut.xlsx', rowNames = T, colNames = T, overwrite = T)
#
rm(InputMaf_LI, InputMaf_noLI, Output)

##### Figure S4A #####
##### 1.1 TCGA data read in
tcgapan.meta <- fread('../TCGA pan metadata') %>% as.data.frame()
tcgapan.mut <- fread('../TCGA pan Mut matrix')
tcgapan.mut <- column_to_rownames(na.omit(tcgapan.mut), var = 'sample') %>% as.data.frame()
# mutual sample 
temp.mutual.sample <- intersect(colnames(tcgapan.mut), tcgapan.meta$sample)
tcgapan.meta <- filter(tcgapan.meta, sample %in% temp.mutual.sample)
tcgapan.mut <- tcgapan.mut[, temp.mutual.sample]
apply(tcgapan.mut[c('MAGEC1','MACF1','PXDN','ADGRV1','NF1'),], 1, function(x){table(x)})
# check sample
tcgapan.mut <- tcgapan.mut[, tcgapan.meta$sample]
identical(colnames(tcgapan.mut), tcgapan.meta$sample)
# change Stage
tcgapan.meta <- tcgapan.meta %>%
  mutate(Stage = case_when(
    ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA", "Stage IB", "I/II NOS") ~ "Stage I",
    ajcc_pathologic_tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC") ~ "Stage II",
    ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "Stage III",
    ajcc_pathologic_tumor_stage %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC") ~ "Stage IV",
    TRUE ~ NA  
  ))
##### 1.2 Func prep
rm(geneX_status, geneX, data_combined, percentage_table, contingency_table, fisher_test, p_value,
   count_table, fisher_test, heatmap_title)
GeneStage <- function(geneX) {
  # 
  geneX_status <- tibble(
    sample = colnames(tcgapan.mut),
    MutationStatus = ifelse(as.numeric(tcgapan.mut[geneX, ]) == 0, "Non-Mutated", "Mutated")
  )
  
  # 
  data_combined <- tcgapan.meta %>%
    left_join(geneX_status, by = "sample") %>%
    filter(!is.na(Stage)) %>%  # remove NA
    mutate(StageIV_Status = ifelse(Stage == "Stage IV", "Stage IV", "Non-Stage IV"))  
  # mutate(StageIV_Status = ifelse(Stage %in% c("Stage III","Stage IV"), "Stage III/IV", "Non-Stage III/IV"))  
  # mutate(StageIV_Status = ifelse(Stage == "Stage III", "Stage III", "Non-Stage III")) 
  
  # 
  contingency_table <- table(data_combined$StageIV_Status, data_combined$MutationStatus)
  fisher_test <- fisher.test(contingency_table)  # Fisher 
  p_value <- fisher_test$p.value
  
  # 
  percentage_table <- data_combined %>%
    count(StageIV_Status, MutationStatus) %>%
    group_by(StageIV_Status) %>%
    mutate(Percentage = n / sum(n))  
  
  #vis
  ggplot(percentage_table, aes(x = StageIV_Status, y = Percentage, fill = MutationStatus)) +
    geom_bar(stat = "identity", position = "fill", color = "black") +  
    scale_y_continuous(labels = scales::percent, expand = c(0,0)) +  
    scale_fill_manual(values = c("Non-Mutated" = "#B3B3B3", "Mutated" = ColJournal$Nature[1])) +  
    theme_minimal() +
    theme(
      text = element_text(family = "Arial", color = "black"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    labs(
      title = paste0("Association between ", geneX, " Mutation and Stage IV"),
      # title = paste0("Association between ", geneX, " Mutation and Stage III/IV"),
      # title = paste0("Association between ", geneX, " Mutation and Stage III"),
      x = "Stage Group",
      y = "Proportion",
      fill = "Mutation Status"
    ) +
    geom_text(aes(label = scales::percent(Percentage, accuracy = 1)), 
              position = position_fill(vjust = 0.5), color = "white", size = 4) + 
    annotate("text", x = 1.5, y = 1.05, label = paste("Fisher's p-value:", format(p_value, digits = 3)),
             hjust = 0.5, size = 5, color = "black")  # p
  
}
plot_grid(GeneStage('MAGEC1'), GeneStage('MACF1'), GeneStage('PXDN'), GeneStage('NF1'),
          nrow = 2)
GeneStage('MACF1') 


##### Figure 4C #####
load('CNVIndex.RData')
Input <- filter(Input, Type == 'Cancer')
##### add MYC
Input <- select(Input, LM, CNV_ERBB2_value, CNV_MYC_value)
Input$LM <- ifelse(Input$LM == 1, 'Liver invasion', 'Cancer Insitu') %>%
  factor(levels = c('Cancer Insitu','Liver invasion'), ordered = T)
##### Vis
plot_grid(ggboxplot(Input, 'LM', 'CNV_ERBB2_value', color = 'LM', add = 'jitter') +
            scale_color_manual(values = c('Cancer Insitu' = ColColor$`Low-Indigo`[5],
                                          'Liver invasion' = ColColor$`High-OrangeRed`[5])) +
            xlab('') + ylab('ERBB2 CNV') +
            theme(legend.position = 'none') +
            stat_compare_means(method = 'wilcox'), # 4*4 F6E NE
          ggboxplot(Input, 'LM', 'CNV_MYC_value', color = 'LM', add = 'jitter') +
            scale_color_manual(values = c('Cancer Insitu' = ColColor$`Low-Indigo`[5],
                                          'Liver invasion' = ColColor$`High-OrangeRed`[5])) +
            xlab('') + ylab('MYC CNV') +
            theme(legend.position = 'none') +
            stat_compare_means(method = 'wilcox'),
          nrow = 1)
#
rm(Input)

##### Figure S4B #####
index.cnv <- filter(index.prim, WES %in% colnames(GBC_Main_CNV))
input.cnv <- GBC_Main_CNV[, index.cnv$WES]
identical(index.cnv$WES, colnames(input.cnv))
# CM
cm.cnv <- data.frame(Yes = ifelse(index.cnv$LM == 'Yes', 1, 0),
                     No = ifelse(index.cnv$LM == 'Yes', 0, 1))
rownames(cm.cnv) <- index.cnv$Sample
# DEA
dea.cnv <- bioinfoamateur::core_Differential_analysis_continuous(cm.cnv, input.cnv, log = T, method = 't.test', show.belong = F, p.adj = F)
dea.cnv <- as.tibble(dea.cnv) %>% mutate(Gene = rownames(dea.cnv))
dea.cnv$FDR <- p.adjust(dea.cnv$`Yes-p-value`, method = 'fdr')
# add coord
dea.cnv$Coord <- sapply(dea.cnv$Gene, function(x){
  temp.band <- filter(Ref.genome, gene_name == x)$ChrBand})
dea.cnv$Coord <- as.character(dea.cnv$Coord)
dea.cnv$Coord <- ifelse(dea.cnv$Coord != 'character(0)', dea.cnv$Coord, NA)
dea.cnv$Coord = factor(dea.cnv$Coord, levels = levels(Chr$ChrBand), ordered = T)
##### Vis
dea.cnv.vis <- dea.cnv %>%
  arrange(`Yes-logFC`) %>%
  mutate(
    Rank = row_number(),  # rank
    Size = -log10(`Yes-p-value`),  
    ColorGroup = case_when(
      `Yes-p-value` > 0.05 ~ "Not Significant",  
      Coord %in% c("7p22.3", "7p22.2", "7p22.1") ~ "Region 7p22 (Significant)",  
      Coord %in% c("8q24.21", "8q24.22", "8q24.23", '8q24.3') ~ "Region 8q24 (Significant)", 
      TRUE ~ "Other Significant"
    )
  )
dea.cnv.vis$Transparency <- case_when(dea.cnv.vis$`Yes-p-value` < 0.05 & dea.cnv.vis$ColorGroup %in% c("Region 7p22 (Significant)", "Region 8q24 (Significant)") ~ 1,
                                      dea.cnv.vis$`Yes-p-value` <= 0.05 ~ 0.1,
                                      dea.cnv.vis$`Yes-p-value` > 0.05 ~ 0.02)

# Vis
ggplot(dea.cnv.vis, aes(x = Rank, y = `Yes-logFC`, size = Size, color = ColorGroup, alpha = Transparency)) +
  geom_point() +
  scale_color_manual(
    values = c("Not Significant" = "grey90", "Region 7p22 (Significant)" = "blue", 
               "Region 8q24 (Significant)" = "red", "Other Significant" = "black")
  ) +
  scale_size_continuous(range = c(1, 10)) +  
  theme_minimal() +
  theme(
    text = element_text(family = "Arial", color = "black"), 
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black")  
  ) +
  labs(
    title = "Differential CNV Analysis: Yes vs No",
    x = "LogFC Ranking",
    y = "Yes vs No LogFC",
    size = expression(-log[10]("Yes-p-value")),
    color = "Significance Category"
  ) +
  guides(alpha = "none")  
# data for Table S5
write.xlsx(dea.cnv.vis, file = 'DEA CNV.xlsx', rowNames = T, colNames = T, overwrite = T)

#
rm(dea.cnv.vis)

##### Figure S4C #####
# data for Table S5
ORA.7p <- clusterProfiler::enricher(Ref.genome %>%
                                      filter(ChrBand %in% c("7p22.3", "7p22.2", "7p22.1")) %>%
                                      pull(gene_name), TERM2GENE = TotalPathway)
ORA.8q <- clusterProfiler::enricher(Ref.genome %>%
                                      filter(ChrBand %in% c("8q24.21", "8q24.22", "8q24.23", "8q24.3")) %>%
                                      pull(gene_name), TERM2GENE = TotalPathway)
ORA.7p <- as.data.frame(ORA.7p@result) %>% filter(pvalue <= 0.05) %>% arrange(pvalue) %>% mutate(Chr = '7p22')
ORA.8q <- as.data.frame(ORA.8q@result) %>% filter(pvalue <= 0.05) %>% arrange(pvalue) %>% mutate(Chr = '8q24')
write.xlsx(rbind(ORA.7p, ORA.8q), file = 'Chr7-8 func.xlsx', rowNames = T, colNames = T, overwrite = T)
##
ORA.7p <- clusterProfiler::enricher(Ref.genome %>%
                                      filter(ChrBand %in% c("7p22.3", "7p22.2", "7p22.1")) %>%
                                      pull(gene_name), TERM2GENE = TotalPathway)
ORA.8q <- clusterProfiler::enricher(Ref.genome %>%
                                      filter(ChrBand %in% c("8q24.21", "8q24.22", "8q24.23", "8q24.3")) %>%
                                      pull(gene_name), TERM2GENE = TotalPathway)
# pathway selection
ORA.7p <- as.data.frame(ORA.7p@result) %>% filter(pvalue <= 0.05) %>% arrange(pvalue) %>% head(10)
ORA.8q <- as.data.frame(ORA.8q@result) %>% filter(pvalue <= 0.05) %>% arrange(pvalue) %>% head(10)

# Vis
ggplot(ORA.7p, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  scale_y_continuous(expand = c(0,0)) +
  theme(
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black"), 
    text = element_text(family = "Arial", color = "black") 
  ) +
  labs(title = "Top 10 Enriched Pathways for 7p22 Region",
       x = "Pathway",
       y = "-log10(p.adjust)")
ggplot(ORA.8q, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  scale_y_continuous(expand = c(0,0)) +
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    text = element_text(family = "Arial", color = "black") 
  ) +
  labs(title = "Top 10 Enriched Pathways for 8q Region",
       x = "Pathway",
       y = "-log10(p.adjust)") 
#
rm(ORA, ORA.7p, ORA.8q, index.cnv, input.cnv)

##### Figure 4D #####
##### 1.1 prep
# primary
index.prim
# exploration
index.expl <- index
index.expl$Condition = as.character(index.expl$Condition)
index.expl$Condition = ifelse(index.expl$Condition == 'Cancer Insitu', 'Primary cancer', index.expl$Condition)
index.expl$Condition = factor(index.expl$Condition, levels = c('NAT','Primary cancer','Liver invasion'), ordered = T)
# validation
index.vali <- tibble(Patient = rep(GBC_LMv_Index$Sample, 3),
                     Condition = c(rep('NAT',11), rep('Primary cancer',11), rep('Liver invasion',11)) %>%
                       factor(levels = c('NAT','Primary cancer','Liver invasion'), ordered = T),
                     Sample = c(GBC_LMv_Index$ProPhos_P, GBC_LMv_Index$ProPhos_T, GBC_LMv_Index$ProPhos_LM))
index.vali <- na.omit(index.vali)

##### 1.2 Gene sets prep
# Genes
set.emt <- c("FN1", "VIM", "S100A4", "CDH2", "CTNNB1", "SNAI1", 
             "SNAI2", "TWIST1", "ZEB2", "DDR2", "ZEB1")
table(filter(TotalPathway, Gene %in% set.emt)$Pathway) %>% names()
set.mmp <- c("MMP1", "MMP2", "MMP3", "MMP7", "MMP8", "MMP9", "MMP10", "MMP11", "MMP12", "MMP13", "MMP14")
table(filter(TotalPathway, Gene %in% set.emt)$Pathway) %>% names()
# selected marker pathways
union(table(filter(TotalPathway, Gene %in% set.emt)$Pathway) %>% names(),
      table(filter(TotalPathway, Gene %in% set.emt)$Pathway) %>% names())
set.vis <- c('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', # EMT
             'REACTOME Cell-cell junction organization', 'KEGG ECM-receptor interaction',
             'KEGG Pathways in cancer')
# invasion pathway/GOs
set.inv.pathways <- c('GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX',
                      'KEGG ECM-receptor interaction',
                      'KEGG Cell adhesion molecules',
                      'KEGG Focal adhesion')

##### 1.3 DEA
# DEA
index.prim
cm.prim <- data.frame(Yes = ifelse(index.prim$LM == 'Yes', 1, 0),
                      No = ifelse(index.prim$LM == 'No', 1, 0))
dea.prim <- bioinfoamateur::core_Differential_analysis_continuous(cm.prim, GBC_Main_Pro[,index.prim$Pro],
                                                                  p.adj = F, show.belong = F, log = T)
dea.prim$Gene = rownames(dea.prim)
filter(dea.prim, Gene %in% c(set.emt, set.mmp), `Yes-p-value` <= 0.05)
filter(dea.prim, Gene %in% filter(TotalPathway, Pathway %in% set.inv.pathways)$Gene,
       `Yes-p-value` <= 0.05)
set.select <- c(filter(dea.prim, Gene %in% c(set.emt, set.mmp), `Yes-p-value` <= 0.05)$Gene,
                filter(dea.prim, Gene %in% filter(TotalPathway, Pathway %in% set.inv.pathways)$Gene,
                       `Yes-p-value` <= 0.05)$Gene) %>% unique()

##### 1.4 Vis upper
##### df prep
dea.prim$`-log10(Yes-p-value)` <- -log10(dea.prim$`Yes-p-value`)
dea.prim.vis <- filter(dea.prim, Gene %in% set.select)
#
p1 <- ggplot(dea.prim.vis, aes(x = reorder(Gene, -`Yes-logFC`), y = `-log10(Yes-p-value)`, 
                               size = abs(`Yes-logFC`), color = `Yes-logFC`)) +
  geom_point(alpha = 1) +
  scale_size_continuous(range = c(4,10)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(
    text = element_text(color = 'black'),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),  
    axis.title.x = element_blank(),  # Hide x-axis label for better alignment
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)  # Adjust margins for tight fit
  ) +
  labs(
    title = NULL,
    y = "-log10(Yes-p-value)",
    color = "Log Fold Change (LogFC)",
    size = "Absolute LogFC"
  )
p1

##### 1.5 Vis lower
# Create a binary matrix indicating the presence of each gene in pathways
TotalPathway_full <- data.frame(Gene = set.select) %>%
  left_join(filter(TotalPathway, Pathway %in% set.vis), by = "Gene")
# Create a binary matrix indicating the presence of each gene in pathways
gene_pathway_matrix <- table(TotalPathway_full$Gene, TotalPathway_full$Pathway)
# Convert to data frame for plotting
gene_pathway_matrix <- as.data.frame(as.table(gene_pathway_matrix))
colnames(gene_pathway_matrix) <- c("Gene", "Pathway", "Count")
gene_pathway_matrix <- gene_pathway_matrix %>%
  mutate(
    Gene = factor(Gene, levels = arrange(dea.prim.vis, desc(`Yes-logFC`))$Gene),  # Ensure Gene is ordered as in set.select
    Count = ifelse(Count > 0, 1, NA)  # Use NA to represent no pathway association
  )

# Create the dot plot for pathway-gene relationships
p2 <- ggplot(gene_pathway_matrix, aes(x = Gene, y = Pathway)) +
  geom_point(aes(size = Count, color = Pathway), alpha = 0.8) +
  # scale_size_manual(range = c(3, 6), na.value = 0) + 
  scale_color_manual(values = ColJournal$Nature[1:3]) +  # Assign different colors to each pathway
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_blank(),  # Remove axis lines
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)  # Adjust margins for tight fit
  ) +
  labs(title = NULL, x = NULL, y = NULL)
print(p2)
# combine
library(patchwork)
combined_plot <- p1 / p2 + plot_layout(heights = c(2, 1))
print(combined_plot)

# data for Table S5
write.xlsx(dea.prim.vis, file = 'Up-lower up.xlsx', rowNames = T, colNames = T, overwrite = T)
write.xlsx(gene_pathway_matrix, file = 'Up-lower lower.xlsx', rowNames = T, colNames = T, overwrite = T)

# rm
rm(TotalPathway_filtered, gene_pathway_matrix, p1, p2, combined_plot, dea.prim.vis)

##### Figure S4D #####
##### 1.0 index
Index <- tibble(Patient = rep('A',21),
                Condition = rep('A',21))
Index$Patient <- rep(GBC_LMe_Index$Sample[1:7],3)
Index$Condition <- c(rep('NAT',7),rep('Cancer Insitu',7),rep('Liver invasion',7))
Index$Sample <- c(GBC_LMe_Index$ProPhos_P, GBC_LMe_Index$ProPhos_T, GBC_LMe_Index$ProPhos_LM)
Index = na.omit(Index)

##### 1.1 signature
SignaturepEMT <- read.xlsx('./pEMT_genes.xlsx', colNames = T)

##### 1.2 ssgsea
Output <- GSVA::gsva(as.matrix(GBC_LMe_Pro[,Index$Sample]), SignaturepEMT, method = 'ssgsea') %>% as.data.frame()
ssGSEA_LMe <- GSVA::gsva(as.matrix(GBC_LMe_Pro[,Index$Sample]), TotalPathwayGSVA, method = 'ssgsea') %>% as.data.frame()

##### vis prep
Index$pEMT <- Output[1,] %>% as.numeric()
Index$cEMT <- as.numeric(ssGSEA_LMe['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',])
##### Vis
ggboxplot(Index, 'Condition', 'pEMT', color = 'Condition', add = 'jitter') +
  scale_color_manual(values = c('NAT' = ColColor$`Low-Blue`[5],
                                'Cancer Insitu' = ColColor$`Low-Indigo`[5],
                                'Liver invasion' = ColColor$`High-OrangeRed`[5])) +
  xlab('') + ylab('pEMT Score') +
  theme(legend.position = 'none') +
  stat_compare_means(method = 'kruskal') +
  stat_pvalue_manual(compare_means(
    pEMT ~ Condition, data = Index,
    method = "t.test") %>% mutate(y.position = c(0.6,0.7,0.8)), label = "p") # 4*4 F6E pEMT
ggboxplot(Index, 'Condition', 'cEMT', color = 'Condition', add = 'jitter') +
  scale_color_manual(values = c('NAT' = ColColor$`Low-Blue`[5],
                                'Cancer Insitu' = ColColor$`Low-Indigo`[5],
                                'Liver invasion' = ColColor$`High-OrangeRed`[5])) +
  xlab('') + ylab('cEMT Score') +
  theme(legend.position = 'none') +
  stat_compare_means(method = 'kruskal') +
  stat_pvalue_manual(compare_means(
    cEMT ~ Condition, data = Index,
    method = "t.test") %>% mutate(y.position = c(0.6,0.7,0.8)), label = "p") # 4*4 F6E cEMT
# data for Table S5
write.xlsx(Index, file = 'EMT score.xlsx', rowNames = T, colNames = T, overwrite = T)

##### Figure 4E #####
##### 1.1 index
IndexPrim <- rbind(select(GBC_Main_Index, 1,2,4,6) %>% 
                     bioinfoamateur::dfm_change_colnames(c('Sample','WES','RNA','Pro')), 
                   mutate(GBC_Main_Index, Sample = paste0(Sample, '_P')) %>% select(1,3,5,7) %>% 
                     bioinfoamateur::dfm_change_colnames(c('Sample','WES','RNA','Pro')))
IndexPrim$Tumor <- ifelse(!grepl(fixed('_P'),IndexPrim$Sample), 'Yes', NA)
# retain samples having proteomic data
IndexPrim <- filter(IndexPrim, !is.na(Pro), Tumor == 'Yes')
# LI
# LM
IndexPrim$LM <- sapply(IndexPrim$Sample, function(x){filter(GBC_Main_Clinical, Patient_ID == x)$Regional_invasion}) %>% as.numeric()
IndexPrim$LM = ifelse(IndexPrim$LM == 1, 'Yes', 'No') %>%
  factor(levels = c('Yes','No'), ordered = T)

##### 1.2 dea
CM <- data.frame(noLI = ifelse(IndexPrim$LM == 'No', 1, 0),
                 LI = ifelse(IndexPrim$LM == 'Yes', 1, 0))
DEA_Pro <- bioinfoamateur::core_Differential_analysis_continuous(CM, GBC_Main_Pro[,IndexPrim$Pro],
                                                 log = T, p.adj = T)
DEA_Pro$Protein = rownames(DEA_Pro)
DEA_Phos <- bioinfoamateur::core_Differential_analysis_continuous(CM, GBC_Main_ProLevelPhos_knn[,IndexPrim$Pro], 
                                                  log = T, p.adj = T)
DEA_Phos$Protein = rownames(DEA_Phos)

##### 1.3 GSEA 
GSEAPro <- bioinfoamateur::Enrich_create_GSEA_object(data.frame(Gene = DEA_Pro$Protein, FC = DEA_Pro$`LI-logFC`))
#
GSEAOutputPro <- clusterProfiler::GSEA(GSEAPro, TERM2GENE = TotalPathway, minGSSize = 5)@result %>% 
  filter(p.adjust <= 0.05) %>% arrange(desc(NES)) %>% mutate(Belong = 'Pro') 
# data for Table S5
write.xlsx(GSEAOutputPro, file = 'Primary (LI vs noLI).xlsx', rowNames = T, colNames = T, overwrite = T)

#
SelectPathway <- c('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','REACTOME Extracellular matrix organization',
                   'REACTOME Collagen degradation','REACTOME Degradation of the extracellular matrix',
                   'KEGG TGF-beta signaling pathway','HALLMARK_INFLAMMATORY_RESPONSE',
                   'REACTOME Signaling by PDGF','KEGG Linoleic acid metabolism', # LI
                   'HALLMARK_HEME_METABOLISM','REACTOME Apoptotic execution phase',
                   'KEGG Oxidative phosphorylation','REACTOME Common Pathway of Fibrin Clot Formation')
#
GSEAOutputPro <- clusterProfiler::GSEA(GSEAPro, TERM2GENE = TotalPathway, minGSSize = 5, pvalueCutoff = 1)@result %>% 
  filter(ID %in% SelectPathway) %>% arrange(desc(NES)) %>% select(ID, NES, p.adjust, core_enrichment)
##### Vis
GSEAOutputPro$LogFDR = -log10(GSEAOutputPro$p.adjust)
GSEAOutputPro$Belong = ifelse(GSEAOutputPro$NES >= 0, 'Upregulated','Downregulated')
ggbarplot(GSEAOutputPro, 'ID', 'NES', color = 'Belong', fill = 'Belong',
          position = position_dodge()) +
  coord_flip() + xlab('') + ylab('NES') +
  scale_color_manual(values = c('Upregulated' = ColColor$`High-OrangeRed`[5], 'Downregulated' = ColColor$`Low-Indigo`[5])) +
  scale_fill_manual(values = c('Upregulated' = ColColor$`High-OrangeRed`[5], 'Downregulated' = ColColor$`Low-Indigo`[5])) +
  scale_y_continuous(expand = c(0,0)) 

##### Figure 4F #####
##### 1.1 protein selection
SelectGeneUp = c('INPP5K','COA6','OCC1')
SelectGeneDown = c('A1BG', 'MASP1', 'FBXO21')
##### 1.2 vis prep
InputLine = GBC_LMe_Pro[c(SelectGeneUp, SelectGeneDown),
                        c(GBC_LMe_Index$ProPhos_P[1:6], GBC_LMe_Index$ProPhos_T[1:6], GBC_LMe_Index$ProPhos_LM[1:6])] %>%
  t() %>% as.tibble()
# add patient
InputLine$Patient <- rep(GBC_LMe_Index$Sample[1:6],3)
# add condition
InputLine$Condition <- c(rep('NAT',6),rep('Primary cancer',6),rep('Liver invasion',6))
# melt
InputLine = reshape2::melt(InputLine, id.var = c('Patient','Condition'), variable.name = 'Protein',
                           value.name = 'Abundance')
InputLine = as_tibble(InputLine)
# factorize
InputLine$Condition = factor(InputLine$Condition, levels = c('NAT','Primary cancer','Liver invasion'), ordered = T)

##### 1.3 Vis
ggplot(data = InputLine, 
       aes(x = Condition, y = Abundance, color = Condition)) +
  # stat_friedman_test(group.by = 'Condition', wid = 'Patient') +
  stat_compare_means(method = 'kruskal') +
  geom_line(aes(group = Patient), color = 'black') +
  scale_color_manual(values = c('NAT' = ColColor$`Low-Blue`[5],
                                'Cancer Insitu' = ColColor$`Low-Indigo`[5],
                                'Liver invasion' = ColColor$`High-OrangeRed`[5])) + 
  geom_point() + xlab('') + ylab('Protein Abundance') +
  theme(legend.position = 'none') +
  facet_wrap(~Protein, scales = 'free_y', nrow = 2) +
  theme_half_open() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')

##### Figure 4G #####
##### 1.1 load mFuzz results
load('mFuzz clusters.RData')
OutputGene = cbind(Output$cluster, Output$membership) %>% as.data.frame()
colnames(OutputGene)[1] = 'Cluster'
# 
InputUp <- filter(OutputGene, Cluster == 5) # the cluster where proteins were always upregulated in invasave lesions
OutputORAUp <- clusterProfiler::enricher(rownames(InputUp), TERM2GENE = TotalPathway)@result %>%
  select(ID, p.adjust, GeneRatio) %>% filter(p.adjust <= 0.05) %>%
  mutate(LogFDR = -log10(p.adjust),
         Belong = 'Upregulated')
##### Down always:C2
InputDown <- filter(OutputGene, Cluster == 2)
OutputORADown <- clusterProfiler::enricher(rownames(InputDown), TERM2GENE = TotalPathway)@result %>%
  select(ID, p.adjust, GeneRatio) %>% filter(p.adjust <= 0.05) %>%
  mutate(LogFDR = log10(p.adjust),
         Belong = 'Downregulated')
# data for Table S5
write.xlsx(InputUp, file = 'Always up protein.xlsx', rowNames = T, colNames = T, overwrite = T)
write.xlsx(InputDown, file = 'Always down protein.xlsx', rowNames = T, colNames = T, overwrite = T)
# data for Table S5
OutputORAUp <- clusterProfiler::enricher(rownames(InputUp), TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.05) %>%
  mutate(LogFDR = -log10(p.adjust),
         Belong = 'Upregulated')
write.xlsx(OutputORAUp, file = 'Always up protein ORA.xlsx', rowNames = T, colNames = T, overwrite = T)
OutputORADown <- clusterProfiler::enricher(rownames(InputDown), TERM2GENE = TotalPathway)@result %>%
  filter(p.adjust <= 0.05) %>%
  mutate(LogFDR = log10(p.adjust),
         Belong = 'Downregulated')
write.xlsx(OutputORADown, file = 'Always down protein ORA.xlsx', rowNames = T, colNames = T, overwrite = T)

# combine
OutputORA <- rbind(OutputORAUp[c('REACTOME Eukaryotic Translation Elongation',
                                 'REACTOME mRNA Splicing',
                                 'REACTOME Mitochondrial translation',
                                 'REACTOME RHOC GTPase cycle',
                                 'KEGG Steroid biosynthesis',
                                 'HALLMARK_MYC_TARGETS_V1',
                                 'REACTOME Cholesterol biosynthesis'),],
                   OutputORADown[c('KEGG Complement and coagulation cascades',
                                   'HALLMARK_INTERFERON_ALPHA_RESPONSE',
                                   'HALLMARK_COAGULATION',
                                   'REACTOME Extracellular matrix organization',
                                   'REACTOME Platelet degranulation',
                                   'REACTOME Neutrophil degranulation',
                                   'KEGG Cell adhesion molecules'),])
OutputORA = arrange(OutputORA, desc(LogFDR))
##### Vis
ggbarplot(OutputORA, 'ID', 'LogFDR', color = 'Belong', fill = 'Belong',
          position = position_dodge()) +
  coord_flip() + xlab('') + ylab('NES') +
  scale_color_manual(values = c('Upregulated' = ColColor$`High-OrangeRed`[5], 'Downregulated' = ColColor$`Low-Indigo`[5])) +
  scale_fill_manual(values = c('Upregulated' = ColColor$`High-OrangeRed`[5], 'Downregulated' = ColColor$`Low-Indigo`[5])) +
  scale_y_continuous(expand = c(0,0))

##### Figure S4E #####
colnames(ref.go) <- c('Pathway','Gene')
#
gene_list <- bioinfoamateur::Enrich_create_GSEA_object(data_frame(Gene = dea.prim$Gene, 
                                                                  logFC = dea.prim$`Yes-logFC`))
# func
gsea_results <- list()
for (pathway_name in set.inv.pathways) {
  # Extract genes for the specific pathway from TotalPathway
  pathway_genes <- rbind(TotalPathway, ref.go) %>%
    filter(Pathway == pathway_name) %>%
    pull(Gene)
  
  # Perform GSEA using clusterProfiler
  gsea_res <- clusterProfiler::GSEA(
    geneList = gene_list,
    TERM2GENE = data.frame(TERM = pathway_name, GENE = pathway_genes),
    pvalueCutoff = 1,
    verbose = FALSE
  )
  
  # Store results for plotting
  gsea_results[[pathway_name]] <- gsea_res
}
temp <- clusterProfiler::GSEA(gene_list, TERM2GENE = rbind(TotalPathway, ref.go), pvalueCutoff = 1)@result
##### Vis
# Plot GSEA results
p1 <- enrichplot::gseaplot2(gsea_results$GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX,
                            geneSetID = 1, title = 'GO Collagen', subplots = 1:2) +
  annotate("text", x = Inf, y = Inf, label = paste("NES:", round(gsea_results$GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX@result$NES, 2)), 
           hjust = 2, vjust = 2, size = 5, color = "black")
p2 <- enrichplot::gseaplot2(gsea_results$`KEGG ECM-receptor interaction`,
                            geneSetID = 1, title = 'KEGG ECM-receptor interaction', subplots = 1:2) +
  annotate("text", x = Inf, y = Inf, label = paste("NES:", round(gsea_results$`KEGG ECM-receptor interaction`@result$NES, 2)), 
           hjust = 2, vjust = 2, size = 5, color = "black")
p3 <- enrichplot::gseaplot2(gsea_results$`KEGG Cell adhesion molecules`,
                            geneSetID = 1, title = 'KEGG Cell adhesion', subplots = 1:2) +
  annotate("text", x = Inf, y = Inf, label = paste("NES:", round(gsea_results$`KEGG Cell adhesion molecules`@result$NES, 2)), 
           hjust = 2, vjust = 2, size = 5, color = "black")
p4 <- enrichplot::gseaplot2(gsea_results$`KEGG Focal adhesion`,
                            geneSetID = 1, title = 'KEGG Focal adhesion', subplots = 1:2) +
  annotate("text", x = Inf, y = Inf, label = paste("NES:", round(gsea_results$`KEGG Focal adhesion`@result$NES, 2)), 
           hjust = 2, vjust = 2, size = 5, color = "black")
plot_grid(p1,p2,p3,p4, nrow = 2) 
#
rm(p1,p2,p3,p4,temp,gsea_results,gene_list)

##### Figure S4F #####
##### 1.1 index making
index.expl.phos <- filter(index.expl, Condition != "NAT")
##### 1.2 KSEA
cm.prim.phos <- data.frame(LI = ifelse(index.prim$LM == 'Yes', 1, 0),
                           No = ifelse(index.prim$LM == 'No', 1, 0))
cm.expl.phos <- data.frame(LI = ifelse(index.expl.phos$Condition == 'Liver invasion', 1, 0),
                           PC = ifelse(index.expl.phos$Condition != 'Liver invasion', 1, 0))
##### 1.3 dea
dea.prim.phos <- bioinfoamateur::core_Differential_analysis_continuous(cm.prim.phos, GBC_Main_Phos_knn[, index.prim$Pro],
                                                                       p.adj = T, show.belong = F, log = T)
dea.expl.phos <- bioinfoamateur::core_Differential_analysis_continuous(cm.expl.phos, GBC_LMe_Phos_knn[, index.expl.phos$Sample],
                                                                       p.adj = T, show.belong = F, log = T)
##### 1.4 KSEA
library(KSEAapp)
# input
ksea.input.prim <- tibble(Protein = rownames(dea.prim.phos),
                          Gene = str_split_fixed(rownames(dea.prim.phos), ':', 2)[, 1],
                          Peptide = rownames(dea.prim.phos),
                          Residue.Both = str_split_fixed(rownames(dea.prim.phos), ':', 2)[, 2],
                          p = dea.prim.phos$`LI-p.adj`, 
                          FC = dea.prim.phos$`LI-logFC`)
ksea.input.expl <- tibble(Protein = rownames(dea.expl.phos),
                          Gene = str_split_fixed(rownames(dea.expl.phos), ':', 2)[, 1],
                          Peptide = rownames(dea.expl.phos),
                          Residue.Both = str_split_fixed(rownames(dea.expl.phos), ':', 2)[, 2],
                          p = dea.expl.phos$`LI-p.adj`, 
                          FC = dea.expl.phos$`LI-logFC`)
### KSEA running
# prim
KSEA.Complete(KSData,
              ksea.input.prim,
              NetworKIN = T,
              NetworKIN.cutoff = 5,
              m.cutoff = 5,
              p.cutoff = 0.05)
ksea.output.prim <- fread('KSEA Kinase Scores.csv')
# expl
KSEA.Complete(KSData,
              ksea.input.expl,
              NetworKIN = T,
              NetworKIN.cutoff = 5,
              m.cutoff = 5,
              p.cutoff = 0.05)
ksea.output.expl <- fread('KSEA Kinase Scores.csv')

##### Vis prep
ksea.combined <- merge(ksea.output.expl, ksea.output.prim, by = "Kinase.Gene", suffixes = c(".expl", ".prim"))
# data for Table S5
write.xlsx(ksea.combined, file = 'KSEA in two cohort.xlsx', rowNames = T, colNames = T, overwrite = T)

# Filter to keep only significant kinases in either comparison (p-value <= 0.05)
ksea.filtered <- ksea.combined %>%
  filter(p.value.expl <= 0.05 | p.value.prim <= 0.05) %>%
  arrange(desc(z.score.expl))  # Sort by z-score from ksea.out.expl in descending order

# Mark significance (assuming p-value <= 0.05 is significant)
ksea.filtered <- ksea.filtered %>%
  mutate(
    sig.expl = ifelse(p.value.expl <= 0.05, TRUE, FALSE),
    sig.prim = ifelse(p.value.prim <= 0.05, TRUE, FALSE),
    color.expl = ifelse(sig.expl, "blue", "grey"),
    color.prim = ifelse(sig.prim, "red", "grey")
  )

# Reshape data for ggplot
ksea.melted <- ksea.filtered %>%
  select(Kinase.Gene, z.score.expl, z.score.prim, color.expl, color.prim) %>%
  pivot_longer(
    cols = starts_with("z.score"),
    names_to = "Condition",
    values_to = "z.score"
  ) %>%
  mutate(
    Condition = ifelse(Condition == "z.score.expl", "Liver Invasion vs Primary", "Primary (Liver Invasion) vs Primary (No Invasion)"),
    color = ifelse(Condition == "Liver Invasion vs Primary", color.expl, color.prim)
  )
ksea.melted$Kinase.Gene = factor(ksea.melted$Kinase.Gene, 
                                 levels = ksea.filtered$Kinase.Gene, ordered = T)

##### Vis
ggplot(ksea.melted, aes(x = Kinase.Gene, y = z.score, fill = color)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  scale_fill_manual(values = c("blue", "red", "grey")) +
  labs(
    title = NULL,
    x = NULL,
    y = "Relative kinase activity",
    fill = "Condition"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) # 

#
rm(ksea.combined, ksea.filtered, ksea.input.expl, ksea.input.prim, ksea.melted,
   ksea.output.expl, ksea.output.prim)

##### Figure 4H #####
##### 1.1 dea
CMe <- data.frame(Insitu = c(rep(1,7), rep(0,7)),
                  LI = c(rep(0,7), rep(1,7)))
CMv <- data.frame(Insitu = c(rep(1,11), rep(0,11)),
                  LI = c(rep(0,11), rep(1,11)))
# prep
InputProe <- GBC_LMe_Pro[, c(GBC_LMe_Index$ProPhos_T,GBC_LMe_Index$ProPhos_LM)]
InputProv <- GBC_LMv_Pro[, c(GBC_LMv_Index$ProPhos_T,GBC_LMv_Index$ProPhos_LM)]
InputPhose <- GBC_LMe_Phos_knn[, c(GBC_LMe_Index$ProPhos_T,GBC_LMe_Index$ProPhos_LM)]
InputPhosv <- GBC_LMv_Phos_knn[, c(GBC_LMv_Index$ProPhos_T,GBC_LMv_Index$ProPhos_LM)]
##### DEA
DEA_Proe <- bioinfoamateur::core_Differential_analysis_continuous(CMe, InputProe, log = T, p.adj = F, show.belong = F)
DEA_Prov <- bioinfoamateur::core_Differential_analysis_continuous(CMv, InputProv, log = T, p.adj = F, show.belong = F)
DEA_Phose <- bioinfoamateur::core_Differential_analysis_continuous(CMe, InputPhose, log = T, p.adj = F, show.belong = F)
DEA_Phosv <- bioinfoamateur::core_Differential_analysis_continuous(CMv, InputPhosv, log = T, p.adj = F, show.belong = F)
DEA_Proe$Protein = rownames(DEA_Proe); DEA_Prov$Protein = rownames(DEA_Prov)
DEA_Phose$PhosphoSite = rownames(DEA_Phose)
DEA_Phosv$PhosphoSite = rownames(DEA_Phosv)

# 
DEA_Proe = DEA_Proe[,-c(1:2)]; colnames(DEA_Proe) = c('LI_pval_e','LI_LogFC_e','Protein')
DEA_Prov = DEA_Prov[,-c(1:2)]; colnames(DEA_Prov) = c('LI_pval_v','LI_LogFC_v','Protein')
DEA_Phose = DEA_Phose[,-c(1:2)]; colnames(DEA_Phose) = c('LI_pval_e','LI_LogFC_e','PhosphoSite')
DEA_Phosv = DEA_Phosv[,-c(1:2)]; colnames(DEA_Phosv) = c('LI_pval_v','LI_LogFC_v','PhosphoSite')
# 
DEA_Pro <- full_join(DEA_Proe, DEA_Prov, by = 'Protein')
DEA_Phos <- full_join(DEA_Phose, DEA_Phosv, by = 'PhosphoSite')
DEA_Phos$Protein = str_split_fixed(DEA_Phos$PhosphoSite, ':', 2)[,1]

##### 1.2 sig cross-validation
##### add meta
DEA_Pro$Sig_E <- case_when(DEA_Pro$LI_pval_e <= 0.05 & DEA_Pro$LI_LogFC_e > 0 ~ 'Sig Up',
                           DEA_Pro$LI_pval_e <= 0.05 & DEA_Pro$LI_LogFC_e < 0 ~ 'Sig Down',)
DEA_Pro$Sig_V <- case_when(DEA_Pro$LI_pval_v <= 0.05 & DEA_Pro$LI_LogFC_v > 0 ~ 'Sig Up',
                           DEA_Pro$LI_pval_v <= 0.05 & DEA_Pro$LI_LogFC_v < 0 ~ 'Sig Down',)
DEA_Pro$Sig <- case_when(DEA_Pro$Sig_E == 'Sig Up' & DEA_Pro$Sig_V == 'Sig Up' ~ 'Sig Up',
                         DEA_Pro$Sig_E == 'Sig Down' & DEA_Pro$Sig_V == 'Sig Down' ~ 'Sig Down',
                         DEA_Pro$Sig_E == 'Sig Up' & DEA_Pro$Sig_V != 'Sig Down' ~ 'Only Up in E',
                         DEA_Pro$Sig_E != 'Sig Down' & DEA_Pro$Sig_V == 'Sig Up' ~ 'Only Up in V',
                         !is.na(DEA_Pro$Protein) ~ 'Others')
DEA_Pro$Druggable <- ifelse(DEA_Pro$Protein %in% RefDrug$gene_name, 'Druggable', NA)

##### 1.3 Vis
library(ggrepel)
ggscatter(DEA_Pro, 'LI_LogFC_e', 'LI_LogFC_v', color = 'Sig') +
  scale_color_manual(values = c('Sig Up' = 'red',
                                'Only Up in E' = ColColor$`High-Orange`[5],
                                'Only Up in V' = ColColor$`High-Orange`[5],
                                'Sig Down' = 'blue',
                                'Others' = 'grey90')) +
  geom_hline(yintercept = 0, linetype = 2, color = 'black') +
  geom_vline(xintercept = 0, linetype = 2, color = 'black') +
  # geom_hline(yintercept = c(0.5,-0.5), linetype = 2, color = 'grey') +
  # geom_vline(xintercept = c(0.5,-0.5), linetype = 2, color = 'grey') +
  ylab('LogFC (Validation)') + xlab('LogFC (Exploration)') +
  theme(legend.position = 'left') +
  geom_text_repel(aes(x=LI_LogFC_e, y=LI_LogFC_v, 
                      label=Protein), data = filter(DEA_Pro, Sig %in% c('Sig Up','Sig Down')),
                  min.segment.length = 0, seed = 42, box.padding = 0.5) 
# data for Table S5
DEA_Pro$padj_e <- p.adjust(DEA_Pro$LI_pval_e, method = 'BH')
DEA_Pro$padj_v <- p.adjust(DEA_Pro$LI_pval_v, method = 'BH')
write.xlsx(DEA_Pro, file = 'Target fc.xlsx', rowNames = T, colNames = T, overwrite = T)

##### Figure S4G & S4H #####
index.expl <- index.expl %>%
  mutate(PHGDH = as.numeric(GBC_LMe_Pro['PHGDH', index.expl$Sample]),
         ACAT1 = as.numeric(GBC_LMe_Pro['ACAT1', index.expl$Sample]))
index.vali <- index.vali %>%
  mutate(PHGDH = as.numeric(GBC_LMv_Pro['PHGDH', index.vali$Sample]),
         ACAT1 = as.numeric(GBC_LMv_Pro['ACAT1', index.vali$Sample]))
##### Vis
plot_grid(ggboxplot(index.expl, 'Condition', 'PHGDH', color = 'Condition') +
            scale_color_manual(values = c('NAT' = ColColor$`Low-Blue`[5],
                                          'Primary cancer' = ColColor$`Low-Indigo`[5],
                                          'Liver invasion' = ColColor$`High-OrangeRed`[5])) +
            stat_compare_means(method = "kruskal.test") + 
            stat_compare_means(
              method = "wilcox.test",
              comparisons = list(
                c("NAT", "Primary cancer"),
                c("NAT", "Liver invasion"),
                c("Primary cancer", "Liver invasion")
              ),
              label = "p.signif",
              label.y = c(max(index.expl$PHGDH) + 0.4, max(index.expl$PHGDH) + 0.8, max(index.expl$PHGDH) + 1.2) 
            ) + xlab(NULL) + theme(legend.position = 'none'),
          ggboxplot(index.expl, 'Condition', 'ACAT1', color = 'Condition') +
            scale_color_manual(values = c('NAT' = ColColor$`Low-Blue`[5],
                                          'Primary cancer' = ColColor$`Low-Indigo`[5],
                                          'Liver invasion' = ColColor$`High-OrangeRed`[5])) +
            stat_compare_means(method = "kruskal.test") + 
            stat_compare_means(
              method = "wilcox.test",
              comparisons = list(
                c("NAT", "Primary cancer"),
                c("NAT", "Liver invasion"),
                c("Primary cancer", "Liver invasion")
              ),
              label = "p.signif",
              label.y = c(max(index.expl$ACAT1) + 0.4, max(index.expl$ACAT1) + 0.8, max(index.expl$ACAT1) + 1.2)  
            ) + xlab(NULL) + theme(legend.position = 'none'),
          ggboxplot(index.vali, 'Condition', 'PHGDH', color = 'Condition') +
            scale_color_manual(values = c('NAT' = ColColor$`Low-Blue`[5],
                                          'Primary cancer' = ColColor$`Low-Indigo`[5],
                                          'Liver invasion' = ColColor$`High-OrangeRed`[5])) +
            stat_compare_means(method = "kruskal.test") + 
            stat_compare_means(
              method = "wilcox.test",
              comparisons = list(
                c("NAT", "Primary cancer"),
                c("NAT", "Liver invasion"),
                c("Primary cancer", "Liver invasion")
              ),
              label = "p.signif",
              label.y = c(max(index.vali$PHGDH) + 0.4, max(index.vali$PHGDH) + 0.8, max(index.vali$PHGDH) + 1.2)  
            ) + xlab(NULL) + theme(legend.position = 'none'),
          ggboxplot(index.vali, 'Condition', 'ACAT1', color = 'Condition') +
            scale_color_manual(values = c('NAT' = ColColor$`Low-Blue`[5],
                                          'Primary cancer' = ColColor$`Low-Indigo`[5],
                                          'Liver invasion' = ColColor$`High-OrangeRed`[5])) +
            stat_compare_means(method = "kruskal.test") +  
            stat_compare_means(
              method = "wilcox.test",
              comparisons = list(
                c("NAT", "Primary cancer"),
                c("NAT", "Liver invasion"),
                c("Primary cancer", "Liver invasion")
              ),
              label = "p.signif",
              label.y = c(max(index.vali$ACAT1) + 0.3, max(index.vali$ACAT1) + 0.6, max(index.vali$ACAT1) + 0.9)  
            ) + xlab(NULL) + theme(legend.position = 'none'),
          nrow = 2) 



