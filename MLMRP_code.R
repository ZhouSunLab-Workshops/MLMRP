#1.lncRNA DE expression
load("1.Rdata")
library(DESeq2)

lncRNA_dds=DESeqDataSetFromMatrix(countData = train_count, 
                                     colData = train_label, 
                                     design = ~ patient_ID + tissue_type)
lncRNA_dds=DESeq(lncRNA_dds)
lncRNA_res=results(lncRNA_dds, pAdjustMethod = "fdr")

DE_lncRNA=lncRNA_res[c(which(lncRNA_res$log2FoldChange<(-1)),which(lncRNA_res$log2FoldChange>1)),]
DE_lncRNA=subset(DE_lncRNA,padj<0.05)
DE_lncRNA=data.frame(DE_lncRNA)

write.table(DE_lncRNA, "DE_lncRNA.txt", sep = "\t")

## Volcano plot
logFC=lncRNA_res$log2FoldChange
padj=lncRNA_res$padj

data=data.frame(logFC = logFC, padj = padj, row.names = rownames(lncRNA_res))
data$sig[(data$padj >= 0.05)|((data$logFC <= 1) & data$logFC >= -1)]="no"
data$sig[data$padj <= 0.05 & data$logFC >= 1]="up"
data$sig[data$padj <= 0.05 & data$logFC <= -1]="down"
data$label=c()
data[which(rownames(data) %in% logi.model$coefnames), "label"]=rownames(data)[which(rownames(data) %in% logi.model$coefnames)]
data[which(!is.na(data$label)), "label"]=c("ADAMTS9-AS1", "LINC03016", "MIR503HG", "AP003548.1", "PGM5-AS1", "LINC01082")

pdf(file = "pheatmap.pdf", width = 3, height = 4)
theme_set(theme_bw())
p=ggplot(data, aes(logFC, -1*log10(padj), color = sig)) + 
  geom_point() +  
  #xlim(-1.5, 2.5)+
  labs(x = "log2(FoldChange)",y = "-log10(FDR)") + 
  scale_color_manual(values = c("#0072B5","grey","#BC3C28")) + 
  geom_hline(yintercept = -log10(0.05), linetype=4) + 
  geom_vline(xintercept = c(-1, 1),linetype=4) +
  #geom_text_repel(data = data, aes(x = logFC, y = -1*log10(padj), label = label),
  # size = 3, box.padding = unit(0.5, "lines"), segment.color = "black",
  # point.padding = unit(0.4, "lines")) + 
  geom_text(aes(label = label), size = 3) +
  theme(panel.grid = element_blank(), axis.text=element_text(size = 12),
        axis.title = element_text(size=12), legend.text = element_text(size=12)) 
p
dev.off()

## Clustering
library(pheatmap)

train_DE_TPM=lncRNA_TPM[rownames(DE_lncRNA), colnames(train_count)]

annotation_col=data.frame(Tissue_type = ifelse(train_label$tissue_type == "N", "ANT", "ESCC"))
rownames(annotation_col)=rownames(train_label)
annotation_col$Cluster = ifelse(cluster == 1, "Cluster1", "Cluster2")[rownames(annotation_col)]

ann_colors=list(Tissue_type = c(ANT = "#d5c0d6", ESCC = "#7875a3"),
                Cluster = c(Cluster1 = "#f2c4ad", Cluster2 = "#d0866f"))
bk=unique(c(seq(-2,2, length=100)))
pdf(file = "pheatmap.pdf", width=8, height=8)
p=pheatmap(train_DE_TPM,
           breaks = bk, 
           scale = 'row', 
           show_colnames = F, 
           show_rownames = F,
           annotation_col = annotation_col, 
           border=FALSE, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           fontsize = 12, 
           annotation_colors = ann_colors,
           clustering_distance_rows = "canberra", 
           clustering_distance_cols = "canberra",
           clustering_method = "ward.D")
cluster=cutree(p$tree_col, 2)
table(cluster, annotation_col$tissue_type)
dev.off()

# 2.Feature selection
load("2.Rdata")
library(caret)

## 2.1 Remove Zero- and Near Zero-Variance Predictors
DE_lncRNA_TPM=lncRNA_TPM[, rownames(DE_lncRNA)]

nzv=nearZeroVar(DE_lncRNA_TPM, saveMetrics= TRUE)
filtered_DE_lncRNA_TPM=DE_lncRNA_TPM[, !nzv$nzv]

## 2.2 Remove high correlative lncRNAs
cor_filteredlncR=cor(filtered_DE_lncRNA_TPM)

highlyCorlncR=findCorrelation(cor_filteredlncR, cutoff = .9)
filtered_DE_lncRNA_TPM=filtered_DE_lncRNA_TPM[, -highlyCorlncR]

## 2.3 Data normalizasion
preProcValues=preProcess(train[, colnames(filtered_DE_lncRNA_TPM)], method = c("center", "scale"))

train_zscore=predict(preProcValues, train[, colnames(filtered_DE_lncRNA_TPM)])
train_zscore$tissue_type=train_label[rownames(train_zscore), "tissue_type"]

## 2.4 Recursive Feature Elimination
rfeCtrl=rfeControl(functions = rfFuncs,
                   rerank = FALSE,
                   method = "repeatedcv",
                   saveDetails = TRUE,
                   number = 10,
                   repeats = 5,
                   verbose = TRUE,
                   returnResamp = "all",
                   p = 0.8,
                   index = NULL,
                   indexOut = NULL,
                   timingSamps = 0,
                   seeds = NA,
                   allowParallel = TRUE)

set.seed(123)
rf_rfe = rfe(x = train_zscore[, -ncol(train_zscore)],  
             y = train_zscore[, ncol(train_zscore)], 
             sizes = c(1:50),
             rfeControl = rfeCtrl)

plot(rf_rfe, type = c('g','o'), xlim = c(1,50))

best=pickSizeBest(x = rf_rfe$results, metric = "Accuracy", maximize = T)
within1Pct=pickSizeTolerance(x = rf_rfe$results, metric = "Accuracy", maximize = T, tol = 0.5)

signature=rf_rfe$optVariables[1:within1Pct]

## 2.5 Cosistent expression validation
GSE53622=data.frame(t(GSE53622_expr[signature, ]))
GSE53622$tissue_type=GSE53622_clinical[rownames(GSE53622), "tissue_type"]

tend=apply(GSE53622[, 1:7], 2, FUN = function(x){t.test(x[which(GSE53622$tissue_type == "T")], 
                                            x[which(GSE53622$tissue_type == "N")])$statistic})
rm_sig=names(c(which(tend > 0 & DE_lncRNA[signature, "log2FoldChange"] < 0), 
         which(tend < 0 & DE_lncRNA[signature, "log2FoldChange"] > 0)))

signature=setdiff(signature, rm_sig)

## 2.6 Box plot
library(ggplot2)
library(gridExtra)
library(ggpubr)

plot_column_boxplot = function(data, column) {
  ggplot(data = data, aes(x = tissue_type, y = get(column))) +
    theme_classic() +
    geom_boxplot(color="black") + 
    geom_point(aes(colour = tissue_type, shape = tissue_type), size = 2)  + 
    geom_line(aes(group = patient), size = 0.8, colour = "#9C9C9C") + 
    scale_color_manual(values = c(N = "#0000EE", T = "#EE0000")) +
    scale_shape_manual(values = c(N = 16, T = 16)) +
    stat_compare_means(method = "t.test", paired = T) +
    labs(y = column) +
    theme(legend.position = "none") + 
    theme(panel.grid = element_blank(), axis.text=element_text(size = 12),
          axis.title = element_text(size = 12), legend.text = element_text(size = 12))
}

train_signature_zscore=train_zscore[, signature]
train_signature_zscore$tissue_type=train_label[rownames(train_signature_zscore), "tissue_type"]
train_signature_zscore$patient=train_label[rownames(train_signature_zscore), "patient_ID"]

myplots_boxplot=lapply(signature, plot_column_boxplot, data = train_signature_zscore)
grid.arrange(myplots_boxplot[[1]], myplots_boxplot[[2]], myplots_boxplot[[3]], 
             myplots_boxplot[[4]], myplots_boxplot[[5]], myplots_boxplot[[6]],  ncol = 6)

## 2.7 Functional enrichment
library(Hmisc)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

data=t(GSE53622_mRNA_expr[, rownames(GSE53622_clinical)[GSE53622_clinical$tissue_type == "T"]])

M=rcorr(as.matrix(data), type = "pearson")  # 计算相关系数和p值

a1=apply(M$r[7:nrow(M$r), 1:6], 2, FUN = function(x){names(x)[which(x > 0.7)]})
b1=unlist(a1)
b1=unique(b1)
b1_ID=bitr(b1, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

a2=apply(M$r[7:nrow(M$r), 1:6], 2, FUN = function(x){names(x)[which(x < (-0.7))]})
b2=unlist(a2)
b2=unique(b2)
b2_ID=bitr(b2, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

b1_kegg=enrichKEGG(b1_ID$ENTREZID, organism = 'hsa', 
                   keyType = 'kegg', pvalueCutoff = 0.01, qvalueCutoff = 0.01, 
                   minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'fdr',
                   use_internal_data = FALSE)
b2_kegg=enrichKEGG(b2_ID$ENTREZID, organism = 'hsa', 
                   keyType = 'kegg', pvalueCutoff = 0.01, qvalueCutoff = 0.01, 
                   minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'fdr',
                   use_internal_data = FALSE)

b1_keggx=setReadable(b1_kegg, 'org.Hs.eg.db', 'ENTREZID')
p1=heatplot(b1_keggx, showCategory = 10)
b2_keggx=setReadable(b2_kegg, 'org.Hs.eg.db', 'ENTREZID')
p2=heatplot(b2_keggx, showCategory = 10)

cowplot::plot_grid(p1, p2, ncol=1, labels = LETTERS[1:2])

b1_keggx2=pairwise_termsim(b1_keggx, showCategory = 10)
p1=treeplot(b1_keggx2, showCategory = 10, hilight = FALSE, offset_tiplab = 1.5)
b2_keggx2=pairwise_termsim(b2_keggx)
p2=treeplot(b2_keggx2)

aplot::plot_list(p1, p2, tag_levels = 'A', ncol = 1)

# 3. Establishment and verification of a multi-lncRNA diagnostic signature
load("3.Rdata")
## (1) Establishment
train_signature_TDM=data.frame(train_signature_TDM)
train_signature_TDM$tissue_type=train_label[rownames(train_signature_TDM), "tissue_type"]

logi.model_TDM=train(form = tissue_type ~ .,
                     data = train_signature_TDM[, c(signature, "tissue_type")],
                     trControl = trainControl(method = "cv", number = 5),
                     control=list(maxit=100),
                     method = "glm",
                     family = "binomial")

## (2) validation: SCH validation cohort
test_signature_TDM=data.frame(test_signature_TDM)
test_signature_TDM$tissue_type=test_label[rownames(test_signature_TDM), "tissue_type"]

test.score_TDM=predict(logi.model_TDM$finalModel, 
                       newdata = data.frame(test_signature_TDM), 
                       type = "link")
test_signature_TDM$score2=(test.score_TDM - min(test.score_TDM))/(max(test.score_TDM) - min(test.score_TDM))

roc=roc(response = test_signature_TDM$tissue_type, 
        predictor = test_signature_TDM$score2,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 

confusionMatrix(factor(ifelse(test_signature_TDM$score2 >= 0.5, "T", "N")), 
                factor(test_signature_TDM$tissue_type),
                positive = "T")

data=test_signature_TDM
data$predict=ifelse(data$score2 >= 0.5, "T", "N")

data$score3=floor(10 * data$score2) + 1
data$score3[data$score3 > 10]=10
data$score3=as.character(data$score3)

annotation_col=data.frame(Real = data$tissue_type, 
                          Predict = data$predict,
                          Score = data$score3)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(N = "#d5c0d6", T = "#7875a3"),
                Predict = c(N = "#f2c4ad", T = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-1, 1, length=100)))

p=pheatmap(t(data[order(data$score2), c(1, 5, 2, 3, 6, 4)]),
           breaks = bk, 
           scale = 'row', 
           show_colnames = F, 
           show_rownames = T,
           annotation_col = annotation_col, 
           border=FALSE, 
           cellheight = 30,
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           fontsize = 12, 
           annotation_colors = ann_colors,
           cluster_cols = F,
           cluster_rows = F)

## (3) validation: CAMS cohort
matrix.score_TDM=predict(logi.model_TDM$finalModel, 
                         newdata = data.frame(matrix_TDM), 
                         type = "link")
matrix.predict_TDM=predict(logi.model_TDM, newdata = matrix_TDM)
matrix_TDM$score2=(matrix.score_TDM - min(matrix.score_TDM))/(max(matrix.score_TDM) - min(matrix.score_TDM))

roc=roc(response = matrix_TDM$tissue_type, 
        predictor = matrix_TDM$score2,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 

confusionMatrix(factor(ifelse(matrix_TDM$score2 >= 0.5, "T", "N")), 
                factor(matrix_TDM$tissue_type),
                positive = "T")

data=matrix_TDM
data$predict=ifelse(data$score2 >= 0.5, "T", "N")

data$score3=floor(10 * data$score2) + 1
data$score3[data$score3 > 10]=10
data$score3=as.character(data$score3)

annotation_col=data.frame(Real = data$tissue_type, 
                          Predict = data$predict,
                          Score = data$score3)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(N = "#d5c0d6", T = "#7875a3"),
                Predict = c(N = "#f2c4ad", T = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-1, 1, length=100)))

p=pheatmap(t(data[order(data$score2), c(1, 5, 2, 3, 6, 4)]),
           breaks = bk, 
           scale = 'row', 
           show_colnames = F, 
           show_rownames = T,
           annotation_col = annotation_col, 
           border=FALSE, 
           cellheight = 30,
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           fontsize = 12, 
           annotation_colors = ann_colors,
           cluster_cols = F,
           cluster_rows = F)

# 4. Independent validation of the MLMRPscore
load("3.Rdata")
## (1) You cohort
GSE130078.score_TDM=predict(logi.model_TDM$finalModel, 
                            newdata = data.frame(GSE130078_signature_TDM), 
                            type = "link")
GSE130078.predict_TDM=predict(logi.model_TDM, newdata = GSE130078_signature_TDM)
GSE130078_signature_TDM$score2=(GSE130078.score_TDM - min(GSE130078.score_TDM))/(max(GSE130078.score_TDM) - min(GSE130078.score_TDM))

roc=roc(response = GSE130078_signature_TDM$tissue_type, 
        predictor = GSE130078_signature_TDM$score2,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 

confusionMatrix(factor(ifelse(GSE130078_signature_TDM$score2 >= 0.5, "T", "N")), 
                factor(GSE130078_signature_TDM$tissue_type),
                positive = "T")

data=GSE130078_signature_TDM
data$predict=ifelse(data$score2 >= 0.5, "T", "N")

data$score3=floor(10 * data$score2) + 1
data$score3[data$score3 > 10]=10
data$score3=as.character(data$score3)

annotation_col=data.frame(Real = data$tissue_type, 
                          Predict = data$predict,
                          Score = data$score3)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(N = "#d5c0d6", T = "#7875a3"),
                Predict = c(N = "#f2c4ad", T = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-1, 1, length=100)))

p=pheatmap(t(data[order(data$score2), c(1, 5, 2, 3, 6, 4)]),
           breaks = bk, 
           scale = 'row', 
           show_colnames = F, 
           show_rownames = T,
           annotation_col = annotation_col, 
           border=FALSE, 
           cellheight = 30,
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           fontsize = 12, 
           annotation_colors = ann_colors,
           cluster_cols = F,
           cluster_rows = F)

## (2) TCGA-GTEx cohort
TCGA_GTEx.score_TDM=predict(logi.model_TDM$finalModel, 
                            newdata = data.frame(TCGA_GTEx_signature_TDM), 
                            type = "link")
TCGA_GTEx.predict_TDM=predict(logi.model_TDM, newdata = TCGA_GTEx_signature_TDM)
TCGA_GTEx_signature_TDM$score2=(TCGA_GTEx.score_TDM - min(TCGA_GTEx.score_TDM))/(max(TCGA_GTEx.score_TDM) - min(TCGA_GTEx.score_TDM))

roc=roc(response = TCGA_GTEx_signature_TDM$tissue_type, 
        predictor = TCGA_GTEx_signature_TDM$score2,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 

confusionMatrix(factor(ifelse(TCGA_GTEx_signature_TDM$score2 >= 0.5, "T", "N")), 
                factor(TCGA_GTEx_signature_TDM$tissue_type),
                positive = "T")

data=TCGA_GTEx_signature_TDM
data$predict=ifelse(data$score2 >= 0.5, "T", "N")

data$score3=floor(10 * data$score2) + 1
data$score3[data$score3 > 10]=10
data$score3=as.character(data$score3)

annotation_col=data.frame(Real = data$tissue_type, 
                          Predict = data$predict,
                          Score = data$score3)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(N = "#d5c0d6", T = "#7875a3"),
                Predict = c(N = "#f2c4ad", T = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-1, 1, length=100)))

p=pheatmap(t(data[order(data$score2), c(1, 5, 2, 3, 6, 4)]),
           breaks = bk, 
           scale = 'row', 
           show_colnames = F, 
           show_rownames = T,
           annotation_col = annotation_col, 
           border=FALSE, 
           cellheight = 30,
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           fontsize = 12, 
           annotation_colors = ann_colors,
           cluster_cols = F,
           cluster_rows = F)

## (3) Li cohort-1
GSE53624.score=predict(logi.model_TDM$finalModel, 
                       newdata = data.frame(GSE53624_signature), 
                       type = "link")
GSE53624.predict=predict(logi.model_TDM, newdata = GSE53624_signature)
GSE53624_signature$score2=(GSE53624.score - min(GSE53624.score))/(max(GSE53624.score) - min(GSE53624.score))

roc=roc(response = GSE53624_signature$tissue_type, 
        predictor = GSE53624_signature$score2,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 

confusionMatrix(factor(ifelse(GSE53624_signature$score2 >= 0.5, "T", "N")), 
                factor(GSE53624_signature$tissue_type),
                positive = "T")

data=GSE53624_signature
data$predict=ifelse(data$score2 >= 0.5, "T", "N")

data$score3=floor(10 * data$score2) + 1
data$score3[data$score3 > 10]=10
data$score3=as.character(data$score3)

annotation_col=data.frame(Real = data$tissue_type, 
                          Predict = data$predict,
                          Score = data$score3)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(N = "#d5c0d6", T = "#7875a3"),
                Predict = c(N = "#f2c4ad", T = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-1, 1, length=100)))

p=pheatmap(t(data[order(data$score2), c(1, 5, 2, 3, 6, 4)]),
           breaks = bk, 
           scale = 'row', 
           show_colnames = F, 
           show_rownames = T,
           annotation_col = annotation_col, 
           border=FALSE, 
           cellheight = 30,
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           fontsize = 12, 
           annotation_colors = ann_colors,
           cluster_cols = F,
           cluster_rows = F)

## (4) Li cohort-2
GSE53622.score=predict(logi.model_TDM$finalModel, 
                       newdata = data.frame(GSE53622_signature), 
                       type = "link")
GSE53622.predict=predict(logi.model_TDM, newdata = GSE53622_signature)
GSE53622_signature$score2=(GSE53622.score - min(GSE53622.score))/(max(GSE53622.score) - min(GSE53622.score))

roc=roc(response = GSE53622_signature$tissue_type, 
        predictor = GSE53622_signature$score2,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 

confusionMatrix(factor(ifelse(GSE53622_signature$score2 >= 0.5, "T", "N")), 
                factor(GSE53622_signature$tissue_type),
                positive = "T")

data=GSE53622_signature
data$predict=ifelse(data$score2 >= 0.5, "T", "N")

data$score3=floor(10 * data$score2) + 1
data$score3[data$score3 > 10]=10
data$score3=as.character(data$score3)

annotation_col=data.frame(Real = data$tissue_type, 
                          Predict = data$predict,
                          Score = data$score3)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(N = "#d5c0d6", T = "#7875a3"),
                Predict = c(N = "#f2c4ad", T = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-1, 1, length=100)))

p=pheatmap(t(data[order(data$score2), c(1, 5, 2, 3, 6, 4)]),
           breaks = bk, 
           scale = 'row', 
           show_colnames = F, 
           show_rownames = T,
           annotation_col = annotation_col, 
           border=FALSE, 
           cellheight = 30,
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           fontsize = 12, 
           annotation_colors = ann_colors,
           cluster_cols = F,
           cluster_rows = F)

# 5. Early diagnosis
load("3.Rdata")
library(OptimalCutpoints)

## (1) SCH cohort
train_signature_TDM=data.frame(train_signature_TDM)
train_signature_TDM$stage=train_label[rownames(train_signature_TDM), "stage"]

train.score_TDM=predict(logi.model_TDM$finalModel, 
                        newdata = data.frame(train_signature_TDM), 
                        type = "link")
train_signature_TDM$score2=(train.score_TDM - min(train.score_TDM))/(max(train.score_TDM) - min(train.score_TDM))
train_signature_TDM=train_signature_TDM[train_signature_TDM$stage %in% c("I", "II", "N"), ]
train_signature_TDM$stage=ifelse(train_signature_TDM$stage == "N", "N", "I_II")

optimal.I.II.result=optimal.cutpoints(X = "score2", 
                                      status = "stage", 
                                      tag.healthy = "N",
                                      methods = "Youden", 
                                      data = train_signature_TDM,
                                      control = control.cutpoints(),
                                      ci.fit = FALSE,
                                      conf.level = 0.95,
                                      trace = FALSE)

cutoff=optimal.I.II.result$Youden$Global$optimal.cutoff$cutoff  #0.551

test_signature_TDM=data.frame(test_signature_TDM)
test_signature_TDM$stage=test_label[rownames(test_signature_TDM), "stage"]

test.score_TDM=predict(logi.model_TDM$finalModel, 
                       newdata = data.frame(test_signature_TDM), 
                       type = "link")
test_signature_TDM$score2=(test.score_TDM - min(test.score_TDM))/(max(test.score_TDM) - min(test.score_TDM))
test_signature_TDM=test_signature_TDM[test_signature_TDM$stage %in% c("I", "II", "N"), ]
test_signature_TDM$stage=ifelse(test_signature_TDM$stage == "N", "N", "I_II")

SCH=rbind(train_signature_TDM, test_signature_TDM)

roc=roc(response = SCH$stage, 
        predictor = SCH$score2,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 

confusionMatrix(factor(ifelse(SCH$score2 >= cutoff, "I_II", "N")), 
                factor(SCH$stage),
                positive = "I_II")

data=SCH
data$predict=ifelse(data$score2 >= cutoff, "I_II", "N")

data$score3=floor(10 * data$score2) + 1
data$score3[data$score3 > 10]=10
data$score3=as.character(data$score3)

annotation_col=data.frame(Real = data$stage, 
                          Predict = data$predict,
                          Score = data$score3)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(N = "#d5c0d6", I_II = "#7875a3"),
                Predict = c(N = "#f2c4ad", I_II = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-1, 1, length=100)))

p=pheatmap(t(data[order(data$score2), c(1, 5, 2, 3, 6, 4)]),
           breaks = bk, 
           scale = 'row', 
           show_colnames = F, 
           show_rownames = T,
           annotation_col = annotation_col, 
           border=FALSE, 
           cellheight = 30,
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           fontsize = 12, 
           annotation_colors = ann_colors,
           cluster_cols = F,
           cluster_rows = F)

## (2) CAMS validation cohort
matrix.score_TDM=predict(logi.model_TDM$finalModel, 
                         newdata = data.frame(matrix_TDM), 
                         type = "link")
matrix.predict_TDM=predict(logi.model_TDM, newdata = matrix_TDM)
matrix_TDM$score2=(matrix.score_TDM - min(matrix.score_TDM))/(max(matrix.score_TDM) - min(matrix.score_TDM))

stageI_II_sample=intersect(as.numeric(as.character(matrix_TDM[matrix_TDM$tissue_type == "T", "patient"]))/2,
                           RTPCR_clinical[RTPCR_clinical$stage %in% c("I", "II", "I_II"), "patient"])
stageI_II_sample=which(matrix_TDM$patient %in% as.character(stageI_II_sample*2))
matrix_TDM=matrix_TDM[c(1:15, stageI_II_sample), ]
matrix_TDM$stage=ifelse(matrix_TDM$tissue_type == "N", "N", "I_II")

roc=roc(response = matrix_TDM$stage, 
        predictor = matrix_TDM$score2,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 

confusionMatrix(factor(ifelse(matrix_TDM$score2 >= cutoff, "I_II", "N")), 
                factor(matrix_TDM$stage),
                positive = "I_II")

data=matrix_TDM
data$predict=ifelse(data$score2 >= cutoff, "I_II", "N")

data$score3=floor(10 * data$score2) + 1
data$score3[data$score3 > 10]=10
data$score3=as.character(data$score3)

annotation_col=data.frame(Real = data$stage, 
                          Predict = data$predict,
                          Score = data$score3)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(N = "#d5c0d6", I_II = "#7875a3"),
                Predict = c(N = "#f2c4ad", I_II = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-1, 1, length=100)))

p=pheatmap(t(data[order(data$score2), c(1, 5, 2, 3, 6, 4)]),
           breaks = bk, 
           scale = 'row', 
           show_colnames = F, 
           show_rownames = T,
           annotation_col = annotation_col, 
           border=FALSE, 
           cellheight = 30,
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           fontsize = 12, 
           annotation_colors = ann_colors,
           cluster_cols = F,
           cluster_rows = F)

## (3) You cohort
GSE130078.score_TDM=predict(logi.model_TDM$finalModel, 
                            newdata = data.frame(GSE130078_signature_TDM), 
                            type = "link")
GSE130078.predict_TDM=predict(logi.model_TDM, newdata = GSE130078_signature_TDM)
GSE130078_signature_TDM$score2=(GSE130078.score_TDM - min(GSE130078.score_TDM))/(max(GSE130078.score_TDM) - min(GSE130078.score_TDM))

GSE130078_signature_TDM=GSE130078_signature_TDM[GSE130078_signature_TDM$stage %in% c("I", "II", "N"), ]
GSE130078_signature_TDM$stage=ifelse(GSE130078_signature_TDM$stage == "N", "N", "I_II")

roc=roc(response = GSE130078_signature_TDM$stage, 
        predictor = GSE130078_signature_TDM$score2,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 

confusionMatrix(factor(ifelse(GSE130078_signature_TDM$score2 >= cutoff, "I_II", "N")), 
                factor(GSE130078_signature_TDM$stage),
                positive = "I_II")

data=GSE130078_signature_TDM
data$predict=ifelse(data$score2 >= cutoff, "I_II", "N")

data$score3=floor(10 * data$score2) + 1
data$score3[data$score3 > 10]=10
data$score3=as.character(data$score3)

annotation_col=data.frame(Real = data$stage, 
                          Predict = data$predict,
                          Score = data$score3)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(N = "#d5c0d6", I_II = "#7875a3"),
                Predict = c(N = "#f2c4ad", I_II = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-1, 1, length=100)))

p=pheatmap(t(data[order(data$score2), c(1, 5, 2, 3, 6, 4)]),
           breaks = bk, 
           scale = 'row', 
           show_colnames = F, 
           show_rownames = T,
           annotation_col = annotation_col, 
           border=FALSE, 
           cellheight = 30,
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           fontsize = 12, 
           annotation_colors = ann_colors,
           cluster_cols = F,
           cluster_rows = F)

## (4) TCGA-GTEx xohort
TCGA_GTEx.score_TDM=predict(logi.model_TDM$finalModel, 
                            newdata = data.frame(TCGA_GTEx_signature_TDM), 
                            type = "link")
TCGA_GTEx.predict_TDM=predict(logi.model_TDM, newdata = TCGA_GTEx_signature_TDM)
TCGA_GTEx_signature_TDM$score2=(TCGA_GTEx.score_TDM - min(TCGA_GTEx.score_TDM))/(max(TCGA_GTEx.score_TDM) - min(TCGA_GTEx.score_TDM))

TCGA_GTEx_signature_TDM=TCGA_GTEx_signature_TDM[TCGA_GTEx_signature_TDM$stage %in% c("I", "II", "N"), ]
TCGA_GTEx_signature_TDM$stage=ifelse(TCGA_GTEx_signature_TDM$stage == "N", "N", "I_II")

roc=roc(response = TCGA_GTEx_signature_TDM$stage, 
        predictor = TCGA_GTEx_signature_TDM$score2,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 

confusionMatrix(factor(ifelse(TCGA_GTEx_signature_TDM$score2 >= 0.551, "I_II", "N")), 
                factor(TCGA_GTEx_signature_TDM$stage),
                positive = "I_II")

data=TCGA_GTEx_signature_TDM
data$predict=ifelse(data$score2 >= 0.551, "I_II", "N")

data$score3=floor(10 * data$score2) + 1
data$score3[data$score3 > 10]=10
data$score3=as.character(data$score3)

annotation_col=data.frame(Real = data$stage, 
                          Predict = data$predict,
                          Score = data$score3)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(N = "#d5c0d6", I_II = "#7875a3"),
                Predict = c(N = "#f2c4ad", I_II = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-1, 1, length=100)))

p=pheatmap(t(data[order(data$score2), c(1, 5, 2, 3, 6, 4)]),
           breaks = bk, 
           scale = 'row', 
           show_colnames = F, 
           show_rownames = T,
           annotation_col = annotation_col, 
           border=FALSE, 
           cellheight = 30,
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           fontsize = 12, 
           annotation_colors = ann_colors,
           cluster_cols = F,
           cluster_rows = F)

## (5) Li cohort-1
GSE53624.score=predict(logi.model_TDM$finalModel, 
                       newdata = data.frame(GSE53624_signature), 
                       type = "link")
GSE53624.predict=predict(logi.model_TDM, newdata = GSE53624_signature)
GSE53624_signature$score2=(GSE53624.score - min(GSE53624.score))/(max(GSE53624.score) - min(GSE53624.score))

GSE53624_signature=GSE53624_signature[GSE53624_signature$stage %in% c("I", "II", "N"), ]
GSE53624_signature$stage=ifelse(GSE53624_signature$stage == "N", "N", "I_II")

roc=roc(response = GSE53624_signature$stage, 
        predictor = GSE53624_signature$score2,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 

confusionMatrix(factor(ifelse(GSE53624_signature$score2 >= cutoff, "I_II", "N")), 
                factor(GSE53624_signature$stage),
                positive = "I_II")

data=GSE53624_signature
data$predict=ifelse(data$score2 >= cutoff, "I_II", "N")

data$score3=floor(10 * data$score2) + 1
data$score3[data$score3 > 10]=10
data$score3=as.character(data$score3)

annotation_col=data.frame(Real = data$stage, 
                          Predict = data$predict,
                          Score = data$score3)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(N = "#d5c0d6", I_II = "#7875a3"),
                Predict = c(N = "#f2c4ad", I_II = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-1, 1, length=100)))

p=pheatmap(t(data[order(data$score2), c(1, 5, 2, 3, 6, 4)]),
           breaks = bk, 
           scale = 'row', 
           show_colnames = F, 
           show_rownames = T,
           annotation_col = annotation_col, 
           border=FALSE, 
           cellheight = 30,
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           fontsize = 12, 
           annotation_colors = ann_colors,
           cluster_cols = F,
           cluster_rows = F)

## (6) Li cohort-2
GSE53622.score=predict(logi.model_TDM$finalModel, 
                       newdata = data.frame(GSE53622_signature), 
                       type = "link")
GSE53622.predict=predict(logi.model_TDM, newdata = GSE53622_signature)
GSE53622_signature$score2=(GSE53622.score - min(GSE53622.score))/(max(GSE53622.score) - min(GSE53622.score))

GSE53622_signature=GSE53622_signature[GSE53622_signature$stage %in% c("I", "II", "N"), ]
GSE53622_signature$stage=ifelse(GSE53622_signature$stage == "N", "N", "I_II")

roc=roc(response = GSE53622_signature$stage, 
        predictor = GSE53622_signature$score2,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 

confusionMatrix(factor(ifelse(GSE53622_signature$score2 >= cutoff, "I_II", "N")), 
                factor(GSE53622_signature$stage),
                positive = "I_II")

data=GSE53622_signature
data$predict=ifelse(data$score2 >= cutoff, "I_II", "N")

data$score3=floor(10 * data$score2) + 1
data$score3[data$score3 > 10]=10
data$score3=as.character(data$score3)

annotation_col=data.frame(Real = data$stage, 
                          Predict = data$predict,
                          Score = data$score3)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(N = "#d5c0d6", I_II = "#7875a3"),
                Predict = c(N = "#f2c4ad", I_II = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-1, 1, length=100)))

p=pheatmap(t(data[order(data$score2), c(1, 5, 2, 3, 6, 4)]),
           breaks = bk, 
           scale = 'row', 
           show_colnames = F, 
           show_rownames = T,
           annotation_col = annotation_col, 
           border=FALSE, 
           cellheight = 30,
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           fontsize = 12, 
           annotation_colors = ann_colors,
           cluster_cols = F,
           cluster_rows = F)
