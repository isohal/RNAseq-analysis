# RNAseq-analysis
This page describes the code used to conduct bioinformatics analysis in this publication. To read the full publication, please go to this link - 

## A first-in-class fully modified version of miR-34a with outstanding stability, activity, and anti-tumor efficacy

Ahmed M. Abdelaal<sup>1</sup>, Ikjot S. Sohal<sup>1*</sup>, Shreyas Iyer<sup>1</sup>, Kasireddy Sudarshan<sup>2</sup>, Nadia A. Lanman<sup>3,4</sup>, Harish Kothandaraman<sup>3</sup>, Philip S. Low<sup>2,3</sup>, and Andrea L. Kasinski<sup>1,3*</sup>

<sup>1</sup>Department of Biological Sciences, <sup>2</sup>Department of Chemistry, <sup>3</sup>Purdue Institute for Cancer Research, <sup>4</sup>Department of Comparative Pathobiology, Purdue University West Lafayette, Indiana, 47907, USA

\*To whom correspondence should be addressed:

Andrea L. Kasinski :Tel: 765-496-1658; Email: [akasinski\@purdue.edu](mailto:akasinski@purdue.edu)

Ikjot S. Sohal: Email [isohal\@purdue.edu](mailto:isohal@purdue.edu)

## R script for global transcriptomics analysis of PM-miR-34a and FM-miR-34a activity

The following script describes the code to normalize raw read count matrix and generate expression dataset using *DESeq2* package. Downstream analysis of the expression dataset, which includes Volcano plot, overlap of downregulated geneset with predicted targets of miR-34a, ontology analysis for miRNA prediction, and gene ontology analysis is shown for FM34a vs NC comparison.

### Generating normalized count matrix

``` r
library(readr)
library(readxl)
library('biomaRt')

#Read raw counts file
CountMatrix <- read_csv("<filepath>/Quantification_Batch_Correction.csv")
CountMatrix.df = as.data.frame(CountMatrix)
#mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
geneIDs = CountMatrix$Gene_ID

rownames(CountMatrix.df) = geneIDs
CountMatrix.df = CountMatrix.df[,-1]
CountMatrix.df[is.na(CountMatrix.df)] = 0

##Determining normalizing factor
x = colSums(CountMatrix.df)
x = x/x["NC2"]

##Generating normalized matrix
nFM34a1 = as.integer(CountMatrix.df$FM34a1/x["FM34a1"])
nFM34a2 = as.integer(CountMatrix.df$FM34a2/x["FM34a2"])
nFM34a3 = as.integer(CountMatrix.df$FM34a3/x["FM34a3"])
nNC1 = as.integer(CountMatrix.df$NC1/x["NC1"])
nNC2 = as.integer(CountMatrix.df$NC2/x["NC2"])
nNC3 = as.integer(CountMatrix.df$NC3/x["NC3"])
nPM34a1 = as.integer(CountMatrix.df$PM34a1/x["PM34a1"])
nPM34a2 = as.integer(CountMatrix.df$PM34a2/x["PM34a2"])
nPM34a3 = as.integer(CountMatrix.df$PM34a3/x["PM34a3"])
nUn1 = as.integer(CountMatrix.df$Un1/x["Un1"])
nUn2 = as.integer(CountMatrix.df$Un2/x["Un2"])
nUn3 = as.integer(CountMatrix.df$Un3/x["Un3"])

Norm.CountMatrix.df = cbind.data.frame(nFM34a1, nFM34a2, nFM34a3, nNC1, nNC2, nNC3, nPM34a1, nPM34a2, nPM34a3, nUn1, nUn2, nUn3)
rownames(Norm.CountMatrix.df) = CountMatrix$Gene_ID

##Saving normalized matrix
write.csv(Norm.CountMatrix.df, file = "<filepath>/NormalizedCountMatrix_allgenes.csv")

##Creating condition labels for samples in the normalized matrix
deseqData = cbind.data.frame(condition = c("FM34a", "FM34a", "FM34a", "NC", "NC", "NC", "PM34a", "PM34a", "PM34a", "Un", "Un", "Un"))
rownames(deseqData) = colnames(Norm.CountMatrix.df)

##Confirming that colnames of normalized matrix is equal to
##rownames of "condition" vector
all(colnames(Norm.CountMatrix.df) %in% rownames(deseqData))
```

### RNA sequencing analysis using *DESeq2* package

``` r
library(DESeq2)

##Generating DESeq expression dataset from normalized count matrix
dds.normalized = DESeqDataSetFromMatrix(Norm.CountMatrix.df, colData = deseqData, design = ~condition)
colSums(counts(dds.normalized[,1:12]))

##Run DESeq
dds = DESeq(dds.normalized)

##Extracting results by contrasting two conditions (FM34a vs NC is shown here)
results_FM34avsNC = results(dds, contrast = c("condition", "FM34a", "NC"))
results_FM34avsNC_df = as.data.frame(results_FM34avsNC)
results_FM34avsNC_df = results_FM34avsNC_df[results_FM34avsNC_df$pvalue < 0.05,]
write.csv(results_FM34avsNC_df, file = "<filepath>/results_FM34avsNC.csv")

##Top downregulated genes - FM34a vs NC
D_FM34avsNC = results_FM34avsNC[results_FM34avsNC$pvalue < 0.05 & results_FM34avsNC$log2FoldChange < 0,]
allD_FM34avsNC = rownames(D_FM34avsNC[order(D_FM34avsNC$log2FoldChange, decreasing = F),]) #All downregulated genes - FM34a vs NC
topD_FM34avsNC = c(rownames(D_FM34avsNC[order(D_FM34avsNC$pvalue, decreasing = F),])[1:6], rownames(D_FM34avsNC[order(D_FM34avsNC$log2FoldChange, decreasing = F),])[1:2])

##Top upregulated genes - FM34a vs NC
U_FM34avsNC = results_FM34avsNC[results_FM34avsNC$pvalue < 0.05 & results_FM34avsNC$log2FoldChange > 0,]
allU_FM34avsNC = rownames(U_FM34avsNC[order(U_FM34avsNC$log2FoldChange, decreasing = T),]) #All upregulated genes - FM34a vs NC
topU_FM34avsNC = c(rownames(U_FM34avsNC[order(U_FM34avsNC$pvalue, decreasing = F),])[1:6], rownames(U_FM34avsNC[order(U_FM34avsNC$log2FoldChange, decreasing = T),])[1:2])
```

#### Volcano plot for various comparisons (script for "FM34a vs NC" shown here)

``` r
##Setting custom colors for upregulated, downregulated and non-significant genes - FM34a vs NC
keyvals.color_FM34avsNC = ifelse(
  rownames(results_FM34avsNC) %in% allU_FM34avsNC, '#000000',
  ifelse(rownames(results_FM34avsNC) %in% allD_FM34avsNC, '#da4040',
         'grey'))
keyvals.color_FM34avsNC[is.na(keyvals.color_FM34avsNC)] <- 'grey'
names(keyvals.color_FM34avsNC)[keyvals.color_FM34avsNC == 'grey'] = 'Not significant'
names(keyvals.color_FM34avsNC)[keyvals.color_FM34avsNC == '#000000'] = 'Upregulated'
names(keyvals.color_FM34avsNC)[keyvals.color_FM34avsNC == '#da4040'] = 'Downregulated'

##Generate volcano plot
library(EnhancedVolcano)

##Saving tiff file for FM34a vs NC comparison:
tiff("<filepath>/FM34a_vs_NC_20221029.tiff", units = "in", width = 6.5, height = 9, res = 300)

EnhancedVolcano(results_FM34avsNC,
    lab = rownames(results_FM34avsNC),
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = c(topD_FM34avsNC, topU_FM34avsNC),
    pCutoff = 5e-02,
    title = "FM34a vs NC",
    subtitle = NULL,
    border = 'full',
    borderWidth = 1,
    ylim = c(0,55),
    #xlim = c(-8,8),
    pointSize = 3,
    labSize = 6,
    axisLabSize = 24,
    boxedLabels = F,
    drawConnectors = T,
    widthConnectors = 0.8,
    arrowheads = F,
    lengthConnectors = unit(0.01, 'npc'),
    typeConnectors = 'open',
    max.overlaps = 200,
    caption = NULL,
    FCcutoff = 0,
    #col = c("grey30","grey30", "grey30", "#ffad73"),
    colCustom = keyvals.color_FM34avsNC,
    colAlpha = 1/3)
```
Volcano plot of up-regulated and down-regulated genes between cells transfected with FM-miR-34a and NC.

<img src="https://github.com/isohal/RNAseq-analysis/assets/42177931/4e759823-a711-45e0-9414-7c1b5645fed5" width="450"/>

#### Overlap of downregulated genes and predicted miR-34a targets (script for "FM34a vs NC" shown here)

``` r
#Getting list of predicted targets of miR-34a from miRdb database
PredictedTargets_miR34a = read_excel("<filepath>/PredictedTargets.miR34a.xlsx")
rownames(PredictedTargets_miR34a) = PredictedTargets_miR34a$`Gene Symbol`

#OVerlap with FM34a targets
FM34a_pTargets = intersect(PredictedTargets_miR34a$`Gene Symbol`, allD_FM34avsNC)
Rank.dist.FM34a = PredictedTargets_miR34a[FM34a_pTargets,]$`Target Rank`
Target.Score.FM34a = PredictedTargets_miR34a[FM34a_pTargets,]$`Target Score`
FM34a.target.overlap = cbind.data.frame(GeneName = FM34a_pTargets, TargetRank = Rank.dist.FM34a, TargetScore = Target.Score.FM34a)
write.csv(FM34a.target.overlap, file = "<filepath>/FM34a.target.overlap.results.csv")
```

#### Ontology analysis for miRNA prediction

Using the downregulated geneset of FM34a vs NC, miRNA ontology was performed using the online version of *g:Profiler2* package - the parameters are described in the code below. Final data was exported in an excel file and plotted in R using *ggplot* package.

``` r
##Ontology analysis for miRNA (performed using online version of g:Profiler2). The parameters used are mentioned below:
Down.FM34avsNC_GO_miRNA = gost(query = allD_FM34avsNC,
                    organism = "hsapiens",
                    ordered_query = F,
                    multi_query = F,
                    significant = T,
                    exclude_iea = F,
                    measure_underrepresentation = F,
                    evcodes = F,
                    user_threshold = 0.05,
                    correction_method = "g_SCS",
                    domain_scope = "annotated",
                    custom_bg = NULL,
                    sources = c("MIRNA"))

##Import excel file saved from g:Profiler2 miRNA ontology analysis 
FM_miR_34a <- read_excel("<filepath>/FM-miR-34a.xlsx")

##GOplot for miRNA category
ggplot(data = FM_miR_34a) +
  geom_point(aes(x = 1,y = -log10(adj_pval)),
             shape = 19, 
             color = red3, 
             size = 20*FM_miR_34a$'intersection_size'/max(FM_miR_34a$'intersection_size'), 
             position = "jitter") + 
  scale_y_continuous(limits = c(0,30)) +
  scale_x_continuous(limits = c(0.6,1.4)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold", size = 24, colour = "black"),
        axis.ticks.length = unit(5,"points"),
        axis.line = element_line(colour = "black", size = 1))
Down.All_GOplot_FM_miRNA = publish_gostplot(gostplot(Down.FM34avsNC_GO_miRNA, capped = F, interactive = F, pal = c(MIRNA = red3)), highlight_terms = c("MIRNA:hsa-miR-34a-5p", "MIRNA:hsa-miR-34c-5p", "MIRNA:hsa-miR-449b-5p"))
ggsave("Down.All_GOplot_FM_miRNA.tiff",
       plot = last_plot(),
       device = "tiff",
       scale = 1.2,
       width = 4.5,
       height = 6,
       units = "in",
       dpi = 300)
```
MiRNA ontology based on genes downregulated in the FM-miR-34a vs NC comparison.
<img src="https://github.com/isohal/RNAseq-analysis/assets/42177931/ae10bd67-c253-4aa4-97f1-41bfc914f0be" width="450"/>


#### Gene ontology analysis for downregulated and upregulated genes

``` r
##Gene ontology analysis for downregulated genes
library(gprofiler2)
red3 = rgb(205,0,0,115, maxColorValue = 255)
gold2 = rgb(240, 213, 127, 100, maxColorValue = 255)
##Ontology analysis for all categories (for both PM and FM34a)
Down.All_GO = gost(query = list("PM34avsNC" = allD_PM34avsNC,
                           "FM34avsNC" = allD_FM34avsNC),
                    organism = "hsapiens",
                    ordered_query = F,
                    multi_query = T,
                    significant = T,
                    exclude_iea = F,
                    measure_underrepresentation = F,
                    evcodes = F,
                    user_threshold = 0.05,
                    correction_method = "g_SCS",
                    domain_scope = "annotated",
                    custom_bg = NULL,
                    sources = c("GO:BP", "KEGG", "REAC"))

##GOplot for all categories, with specific terms highlighted, for both PM34a vs NC and FM34a vs NC comparisons
Down.All_GOplot = publish_gostplot(gostplot(Down.All_GO, capped = F, interactive = F), highlight_terms = c("GO:0000278", "GO:0008283", "GO:0016477", "GO:0048870", "GO:0012501", "KEGG:04218", "KEGG:04110", "KEGG:04218", "REAC:R-HSA-453279", "REAC:R-HSA-72764", "REAC:R-HSA-927802"))
ggsave("Down.All_GOplot.tiff",
       plot = last_plot(),
       device = "tiff",
       scale = 1.5,
       width = 6,
       height = 9,
       units = "in",
       dpi = 300)
```
Gene ontology of biological processes, KEGG pathways and Reactome pathways based on genes downregulated in the FM-miR-34a vs NC and PM-miR-34a vs NC comparison.
<img src="https://github.com/isohal/RNAseq-analysis/assets/42177931/66122640-64f2-499f-808d-ff43e936582d" width="450"/>

``` {.r .R}
##Gene ontology analysis for upregulated genes
Up.All_GO = gost(query = list("PM34avsNC" = allU_PM34avsNC,
                           "FM34avsNC" = allU_FM34avsNC),
                    organism = "hsapiens",
                    ordered_query = F,
                    multi_query = T,
                    significant = T,
                    exclude_iea = T,
                    measure_underrepresentation = F,
                    evcodes = F,
                    user_threshold = 0.05,
                    correction_method = "g_SCS",
                    domain_scope = "annotated",
                    custom_bg = NULL,
                    sources = c("GO:BP", "KEGG", "REAC"))

##GOplot for all categories, with specific terms highlighted, for both PM34a vs NC and FM34a vs NC comparisons
Up.All_GOplot = publish_gostplot(gostplot(Up.All_GO, capped = F, interactive = F), highlight_terms = c("GO:0030030", "GO:0044782", "GO:0007018", "KEGG:00970", "REAC:R-HSA-5617833", "REAC:R-HSA-379716"))
ggsave("Up.All_GOplot.tiff",
       plot = last_plot(),
       device = "tiff",
       scale = 1.5,
       width = 6,
       height = 9,
       units = "in",
       dpi = 300)
```
Gene ontology of biological processes, KEGG pathways and Reactome pathways based on genes downregulated in the FM-miR-34a vs NC and PM-miR-34a vs NC comparison.
<img src="https://github.com/isohal/RNAseq-analysis/assets/42177931/6fa93027-f520-409e-8a38-0dd1c70da7e6" width="450"/>
