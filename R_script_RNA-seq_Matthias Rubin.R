# Loading DESEq2 library
library(DESeq2)


# Read count data
countData <- read.table("featurecounts_paired.txt", header = TRUE, row.names = 1)
# List of columns to exclude
columns_to_exclude <- c("Chr", "Start", "End", "Strand", "Length")
# Subset to keep only the relevant columns for DESeq2
countData_subset <- countData[, !(colnames(countData) %in% columns_to_exclude)]
# Create a vector of new row names
new_column_names  <- c("HER21", "HER22", "HER23", "NonTNBC1", "NonTNBC2", "NonTNBC3","Normal1","Normal2","Normal3","TNBC1","TNBC2", "TNBC3")
# Assign the new row names to the countData
colnames(countData_subset) <- new_column_names 

# read in sample info
colData<- read.table("sample_info.txt", header = TRUE, row.names = 1)
#making sure the row names in colData matches the column names in countData_subset
all(colnames(countData_subset) %in% rownames(colData))
# are they in the same order?
all(colnames(countData_subset) == rownames(colData))

# generating the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData_subset,
                              colData = colData,
                              design = ~ origin)


dds

# setring factor level
dds$origin <- relevel(dds$origin, ref ="Normal")

# running DESeq
dds <- DESeq(dds)
res <- results(dds)
res

## Visualize the data and quality control
# Extracting transformed values
vsd <- vst(dds, blind=TRUE)
# Principal component analysis of the samples
plotPCA(vsd, intgroup=c("origin"))

### Contrasting out the comparisons of the 3 tumor subgroups
## TNBC_VS_HER2
res_TNBC_vs_HER2 <- results(dds, contrast=c("origin","TNBC","HER2"), alpha= 0.05)
## Adding gene names
res_TNBC_vs_HER2$ensembl <- sapply( strsplit( rownames(res_TNBC_vs_HER2), split="\\+" ), "[", 1 )
library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res_TNBC_vs_HER2$ensembl,
                  mart = ensembl )
idx <- match( res_TNBC_vs_HER2$ensembl, genemap$ensembl_gene_id )
res_TNBC_vs_HER2$entrez <- genemap$entrezgene[ idx ]
res_TNBC_vs_HER2$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(res_TNBC_vs_HER2,4)

# Numbers of genes up or downregulated
summary(res_TNBC_vs_HER2)
# filtering out the genes with adjusted p-value < 0.05
sum( res_TNBC_vs_HER2$padj < 0.05, na.rm=TRUE )
resSig_TNBC_vs_HER2 <- res_TNBC_vs_HER2[ which(res_TNBC_vs_HER2$padj < 0.05 ), ]
summary(resSig_TNBC_vs_HER2)
# significant genes with the strongest down-regulation
head( resSig_TNBC_vs_HER2[ order( resSig_TNBC_vs_HER2$log2FoldChange ), ] )
# with the strongest upregulation
tail( resSig_TNBC_vs_HER2[ order( resSig_TNBC_vs_HER2$log2FoldChange ), ] )

## NonTNBC_VS_HER2
res_NonTNBC_vs_HER2 <- results(dds, contrast=c("origin","NonTNBC","HER2"), alpha= 0.05)
## Adding gene names
res_NonTNBC_vs_HER2$ensembl <- sapply( strsplit( rownames(res_NonTNBC_vs_HER2), split="\\+" ), "[", 1 )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res_NonTNBC_vs_HER2$ensembl,
                  mart = ensembl )
idx <- match( res_NonTNBC_vs_HER2$ensembl, genemap$ensembl_gene_id )
res_NonTNBC_vs_HER2$entrez <- genemap$entrezgene[ idx ]
res_NonTNBC_vs_HER2$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(res_NonTNBC_vs_HER2,4)

# Numbers of genes up or downregulated
summary(res_NonTNBC_vs_HER2)
# filtering out the genes with adjusted p-value < 0.05
sum( res_NonTNBC_vs_HER2$padj < 0.05, na.rm=TRUE )
resSig_NonTNBC_vs_HER2 <- res_NonTNBC_vs_HER2[ which(res_NonTNBC_vs_HER2$padj < 0.05 ), ]
summary(resSig_NonTNBC_vs_HER2)
# significant genes with the strongest down-regulation
head( resSig_NonTNBC_vs_HER2[ order( resSig_NonTNBC_vs_HER2$log2FoldChange ), ] )
#  with the strongest upregulation
tail( resSig_NonTNBC_vs_HER2[ order( resSig_NonTNBC_vs_HER2$log2FoldChange ), ] )

## NonTNBC_VS_TNBC
res_NonTNBC_vs_TNBC <- results(dds, contrast=c("origin","NonTNBC","TNBC"), alpha= 0.05)
## Adding gene names
res_NonTNBC_vs_TNBC$ensembl <- sapply( strsplit( rownames(res_NonTNBC_vs_TNBC), split="\\+" ), "[", 1 )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res_NonTNBC_vs_TNBC$ensembl,
                  mart = ensembl )
idx <- match( res_NonTNBC_vs_TNBC$ensembl, genemap$ensembl_gene_id )
res_NonTNBC_vs_TNBC$entrez <- genemap$entrezgene[ idx ]
res_NonTNBC_vs_TNBC$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(res_NonTNBC_vs_TNBC,4)

# Numbers of genes up or downregulated
summary(res_NonTNBC_vs_TNBC)
# filtering out the genes with adjusted p-value < 0.05
sum( res_NonTNBC_vs_TNBC$padj < 0.05, na.rm=TRUE )
resSig_NonTNBC_vs_TNBC <- res_NonTNBC_vs_TNBC[ which(res_NonTNBC_vs_TNBC$padj < 0.05 ), ]
summary(resSig_NonTNBC_vs_TNBC)
# significant genes with the strongest down-regulation
head( resSig_NonTNBC_vs_TNBC[ order( resSig_NonTNBC_vs_TNBC$log2FoldChange ), ] )
#  with the strongest upregulation
tail( resSig_NonTNBC_vs_TNBC[ order( resSig_NonTNBC_vs_TNBC$log2FoldChange ), ] )

## TNBC_VS_Normal (this was done to just have an idea how the difference between cancer and normal could look; the data was not used in the report)
res_TNBC_vs_Normal <- results(dds, contrast=c("origin","TNBC","Normal"), alpha= 0.05)
## Adding gene names
res_TNBC_vs_Normal$ensembl <- sapply( strsplit( rownames(res_TNBC_vs_Normal), split="\\+" ), "[", 1 )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res_TNBC_vs_Normal$ensembl,
                  mart = ensembl )
idx <- match( res_TNBC_vs_Normal$ensembl, genemap$ensembl_gene_id )
res_TNBC_vs_Normal$entrez <- genemap$entrezgene[ idx ]
res_TNBC_vs_Normal$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(res_TNBC_vs_Normal,4)

# Numbers of genes up or downregulated
summary(res_TNBC_vs_Normal)
# filtering out the genes with adjusted p-value < 0.05
sum( res_TNBC_vs_Normal$padj < 0.05, na.rm=TRUE )
resSig_TNBC_vs_Normal <- res_TNBC_vs_Normal[ which(res_TNBC_vs_Normal$padj < 0.05 ), ]
summary(resSig_TNBC_vs_Normal)
# significant genes with the strongest down-regulation
head( resSig_TNBC_vs_Normal[ order( resSig_TNBC_vs_Normal$log2FoldChange ), ] )
#  with the strongest upregulation
tail( resSig_TNBC_vs_Normal[ order( resSig_TNBC_vs_Normal$log2FoldChange ), ] )


# Export the results of the contrast and filtering
write.csv( as.data.frame(resSig_TNBC_vs_HER2), file="results_TNBC_vs_HER2_padj_0.05.csv" )
write.csv( as.data.frame(resSig_NonTNBC_vs_HER2), file="results_NonTNBC_vs_HER2_padj_0.05.csv" )
write.csv( as.data.frame(resSig_NonTNBC_vs_TNBC), file="results_NonTNBC_vs_TNBC_padj_0.05.csv" )
write.csv( as.data.frame(resSig_TNBC_vs_Normal), file="results_TNBC_vs_Normal_padj_0.05.csv" )


## Enhanced volcano to show up- and down-regulated genes
# TNBC_VS_HER2
library(EnhancedVolcano)
EnhancedVolcano(resSig_TNBC_vs_HER2,
                lab = resSig_TNBC_vs_HER2$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'TNBC vs HER2',
                pointSize = 3.0,
                labSize = 4.0,
                colAlpha = 2,
                legendPosition = 'None',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)

# Comparing certain genes between groups
SNORA53 <- plotCounts(dds, "ENSG00000212443", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = SNORA53, main = "Expression of SNORA53")

GRB7 <- plotCounts(dds, "ENSG00000141738", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = GRB7, main = "Expression of GRB7")

SNORA22 <- plotCounts(dds, "ENSG00000206634", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = SNORA22, main = "Expression of SNORA22")

CDKN2A <- plotCounts(dds, "ENSG00000147889", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = CDKN2A, main = "Expression of CDKN2A")

ERBB2 <- plotCounts(dds, "ENSG00000141736", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = ERBB2, main = "Expression of ERBB2")

# NonTNBC_VS_HER2
EnhancedVolcano(resSig_NonTNBC_vs_HER2,
                lab = resSig_NonTNBC_vs_HER2$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'NonTNBC vs HER2',
                pointSize = 3.0,
                labSize = 4.0,
                colAlpha = 2,
                legendPosition = 'None',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)

# Comparing certain genes between groups
AGR3 <- plotCounts(dds, "ENSG00000173467", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = AGR3, main = "Expression of AGR3")

ESR1 <- plotCounts(dds, "ENSG00000091831", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = ESR1, main = "Expression of ESR1")

CYP2B7P1 <- plotCounts(dds, "ENSG00000291083", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = CYP2B7P1, main = "Expression of Pseudogen CYP2B7P1")

NKAIN1 <- plotCounts(dds, "ENSG00000084628", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = NKAIN1, main = "Expression of NKAIN1")

# NonTNBC_VS_TNBC
EnhancedVolcano(resSig_NonTNBC_vs_TNBC,
                lab = resSig_NonTNBC_vs_TNBC$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'NonTNBC vs TNBC',
                pointSize = 3.0,
                labSize = 4.0,
                colAlpha = 2,
                legendPosition = 'None',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)

# Comparing certain genes between groups
ELF5 <- plotCounts(dds, "ENSG00000135374", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = ELF5, main = "Expression of ELF5")

AGR2 <- plotCounts(dds, "ENSG00000106541", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = AGR2, main = "Expression of AGR2")

TFF1 <- plotCounts(dds, "ENSG00000160182", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = TFF1, main = "Expression of TFF1")

PGR <- plotCounts(dds, "ENSG00000082175", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = PGR, main = "Expression of PGR")

RAB21 <- plotCounts(dds, "ENSG00000080371", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = RAB21, main = "Expression of RAB21")

# TNBC_VS_Normal (this was done to just have an idea how the difference between cancer and normal could look; the data was not used in the report)
EnhancedVolcano(resSig_TNBC_vs_Normal,
                lab = resSig_TNBC_vs_Normal$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'TNBC vs Normal',
                pointSize = 3.0,
                labSize = 4.0,
                colAlpha = 2,
                legendPosition = 'None',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)

# Comparing certain genes between groups
RN7SLP1  <- plotCounts(dds, "ENSG00000263740", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = RN7SLP1 , main = "Expression of Pseudogen RN7SLP1")

SCARNA7 <- plotCounts(dds, "ENSG00000238741", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = SCARNA7, main = "Expression of SCARNA7")

RN7SL3 <- plotCounts(dds, "ENSG00000278771", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = RN7SL3, main = "Expression of RN7SL3")

UBC <- plotCounts(dds, "ENSG00000150991", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = UBC, main = "Expression of UBC")

PPP1R15A <- plotCounts(dds, "ENSG00000087074", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = PPP1R15A, main = "Expression of PPP1R15A")

RN7SL1 <- plotCounts(dds, "ENSG00000276168", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = RN7SL1, main = "Expression of RN7SL1")

RMRP <- plotCounts(dds, "ENSG00000277027", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = RMRP, main = "Expression of RMRP")

ATF3 <- plotCounts(dds, "ENSG00000162772", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = ATF3, main = "Expression of ATF3")

TMEM219 <- plotCounts(dds, "ENSG00000149932", intgroup = c("origin"), returnData = TRUE)
boxplot(count ~ origin, data = TMEM219, main = "Expression of TMEM219")


## GO overrepresentation analysis
library(clusterProfiler)
library(org.Hs.eg.db)

# TNBC_VS_HER2
geneListAll_TNBC_VS_HER2 <- rownames(res_TNBC_vs_HER2)
geneList_TNBC_VS_HER2 <- rownames(resSig_TNBC_vs_HER2)

ego <- enrichGO(gene         = geneList_TNBC_VS_HER2,
                universe      = geneListAll_TNBC_VS_HER2,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego, 3)  
dotplot(ego) + ggtitle("GO terms TNBC VS HER2")

write.csv( as.data.frame(ego), file="results_enrichGO_TNBC_VS_HER2.csv" )


# NonTNBC_VS_HER2
geneListAll_NonTNBC_VS_HER2 <- rownames(res_NonTNBC_vs_HER2)
geneList_NonTNBC_VS_HER2 <- rownames(resSig_NonTNBC_vs_HER2)

ego <- enrichGO(gene         = geneList_NonTNBC_VS_HER2,
                universe      = geneListAll_NonTNBC_VS_HER2,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENSEMBL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego, 3)  
dotplot(ego) + ggtitle("GO terms NonTNBC VS HER2")

write.csv( as.data.frame(ego), file="results_enrichGO_NonTNBC_VS_HER2.csv" )



# NonTNBC_VS_TNBC
geneListAll_NonTNBC_VS_TNBC <- rownames(res_NonTNBC_vs_TNBC)
geneList_NonTNBC_VS_TNBC <- rownames(resSig_NonTNBC_vs_TNBC)


ego <- enrichGO(gene         = geneList_NonTNBC_VS_TNBC,
                universe      = geneListAll_NonTNBC_VS_TNBC,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego, 3)  
dotplot(ego) + ggtitle("GO terms NonTNBC VS TNBC")

write.csv( as.data.frame(ego), file="results_enrichGO_NonTNBC_VS_TNBC.csv" )

# TNBC_VS_Normal (this was done to just have an idea how the difference between cancer and normal could look; the data was not used in the report)
geneListAll_TNBC_VS_Normal <- rownames(res_TNBC_vs_Normal)
geneList_TNBC_VS_Normal <- rownames(resSig_TNBC_vs_Normal)

ego <- enrichGO(gene         = geneList_TNBC_VS_Normal,
                universe      = geneListAll_TNBC_VS_Normal,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego, 3)  
dotplot(ego) + ggtitle("GO terms TNBC VS Normal")

write.csv( as.data.frame(ego), file="results_enrichGO_TNBC_VS_Normal.csv" )

