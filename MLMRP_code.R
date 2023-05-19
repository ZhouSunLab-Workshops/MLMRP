#Figure 2
##A
library(readxl)
data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 1, skip = 1))
data$padj=as.numeric(data$padj)
data[data$label == "NA", "label"]=NA

pdf(file = "Figure2A.pdf", width = 3, height = 4)
theme_set(theme_bw())
p=ggplot(data, aes(logFC, -1*log10(padj), color = sig)) + 
  geom_point() +  
  #xlim(-1.5, 2.5)+
  labs(x = "log2(FoldChange)",y = "-log10(FDR)") + 
  scale_color_manual(values = c(down="#0072B5", no="grey", up="#BC3C28")) + 
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

##B
library(pheatmap)

train_DE_TPM=data.frame(read_excel(path ="./source_data.xlsx", sheet = 2, skip = 1), row.names = 1)

annotation_col=data.frame(Tissue_type = ifelse(train_DE_TPM$tissue_type == "N", "ANT", "ESCC"),
                          Cluster = train_DE_TPM$Cluster)
rownames(annotation_col)=rownames(train_DE_TPM)

ann_colors=list(Tissue_type = c(ANT = "#d5c0d6", ESCC = "#7875a3"),
                Cluster = c(Cluster1 = "#f2c4ad", Cluster2 = "#d0866f"))
bk=unique(c(seq(-2,2, length=100)))
pdf(file = "Figure2B.pdf", width=8, height=8)
p=pheatmap(t(train_DE_TPM[, 1:(ncol(train_DE_TPM)-2)]),
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
dev.off()

##C
library(ggplot2)
library(gridExtra)
library(ggpubr)

train_signature_zscore=data.frame(read_excel(path ="./source_data.xlsx", sheet = 3, skip = 1), row.names = 1, check.names = F)

plot_column_boxplot = function(data, column) {
  ggplot(data = data, aes(x = tissue_type, y = get(column))) +
    theme_bw() +
    geom_boxplot(color="black") + 
    geom_point(aes(colour = tissue_type, shape = tissue_type), size = 2)  + 
    geom_line(aes(group = patient), size = 0.8, colour = "#9C9C9C") + 
    scale_color_manual(values = c(ANT = "#0000EE", ESCC = "#EE0000")) +
    scale_shape_manual(values = c(ANT = 16, ESCC = 16)) +
    stat_compare_means(method = "t.test", paired = T) +
    labs(title = column, y = "Expression value") +
    theme(legend.position = "none") + 
    theme(panel.grid = element_blank(), axis.text=element_text(size = 12),
          axis.title = element_text(size = 12), legend.text = element_text(size = 12))
}

pdf(file="Figure2C.pdf", width=13, height=3)
myplots_boxplot=lapply(colnames(train_signature_zscore)[1:6], plot_column_boxplot, data = train_signature_zscore)
grid.arrange(myplots_boxplot[[1]], myplots_boxplot[[2]], myplots_boxplot[[3]], 
             myplots_boxplot[[4]], myplots_boxplot[[5]], myplots_boxplot[[6]],  ncol = 6)
dev.off()

p.value=apply(train_signature_zscore[, 1:6], 2, FUN = function(x){t.test(x[which(train_signature_zscore$tissue_type == "ESCC")], x[which(train_signature_zscore$tissue_type == "ANT")], paired = T)$p.value})

##D
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

func_gene=data.frame(read_excel(path ="./source_data.xlsx", sheet = 4, skip = 1))

b1_kegg=enrichKEGG(func_gene$Left_ENTREZID, organism = 'hsa', 
                   keyType = 'kegg', pvalueCutoff = 0.01, qvalueCutoff = 0.01, 
                   minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'fdr',
                   use_internal_data = FALSE)
b2_kegg=enrichKEGG(func_gene$Right_ENTREZID, organism = 'hsa', 
                   keyType = 'kegg', pvalueCutoff = 0.01, qvalueCutoff = 0.01, 
                   minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'fdr',
                   use_internal_data = FALSE)

b1_keggx=setReadable(b1_kegg, 'org.Hs.eg.db', 'ENTREZID')
p1=heatplot(b1_keggx, showCategory = 10)
b2_keggx=setReadable(b2_kegg, 'org.Hs.eg.db', 'ENTREZID')
p2=heatplot(b2_keggx, showCategory = 10)

pdf(file="Figure2D.pdf", width=10, height=8)
cowplot::plot_grid(p1, p2, ncol=1, labels = LETTERS[1:2])
dev.off()

#Figure 3
##A
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 5, skip = 1), row.names = 1)

df=data.frame(value = as.numeric(table(data$Real)), group = names(table(data$Real)))
p1=ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = NA) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#2172AB", "#B14234")) +
  theme_void()

data_m=data.frame(Real = factor(c(rep("ANT", 2), rep("ESCC", 2)), levels = c("ANT", "ESCC")),
                  Predict = factor(rep(c("ANT", "ESCC"), 2), levels = c("ESCC", "ANT")),
                  value = c(as.numeric(table(data$Predict, data$Real))))
data_m$color=paste(data_m$Real, data_m$Predict, sep = "_")
p2=ggplot(data_m, aes(x = Real, y = Predict)) +  
  theme_minimal() + #设置系统自带主题
  theme(panel.grid.major = element_blank()) +  #设置主项网格
  theme(legend.key = element_blank()) + #去掉背景颜色
  geom_tile(aes(fill = color)) +  #设置填充的值
  scale_fill_manual(values = c(ANT_ANT = "#2172AB", ANT_ESCC = "white", ESCC_ANT = "white", ESCC_ESCC = "#B14234")) +
  #scale_fill_gradient(low = "white", high = "red") + #设置颜色梯度
  geom_text(aes(label=value)) + #设置文本显示数值
  geom_vline(xintercept = c(0.5, 1.5, 2.5), size = 0.5) +
  geom_hline(yintercept = c(0.5, 1.5, 2.5), size = 0.5)
p2 #输出图片 

pdf(file="Figure3A.pdf", width=8, height=2.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

##B
library(pheatmap)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 6, skip = 1), row.names = 1, check.names = F)

annotation_col=data.frame(Real = data$Real, 
                          Predict = data$Predict,
                          Score = data$Risk_probability)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(ANT = "#d5c0d6", ESCC = "#7875a3"),
                Predict = c(ANT = "#f2c4ad", ESCC = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-2, 2, length=100)))

p=pheatmap(t(data[order(data$score), 1:6]),
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

pdf(file="Figure3B.pdf", width=8, height=4)
p
dev.off()

##C
library(pROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 7, skip = 1), row.names = 1)

pdf(file="Figure3C.pdf", width=4, height=4)
roc=roc(response = data$Real, 
        predictor = data$score,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 
dev.off()

##D
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 8, skip = 1), row.names = 1)

df=data.frame(value = as.numeric(table(data$Real)), group = names(table(data$Real)))
p1=ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = NA) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#2172AB", "#B14234")) +
  theme_void()

data_m=data.frame(Real = factor(c(rep("ANT", 2), rep("ESCC", 2)), levels = c("ANT", "ESCC")),
                  Predict = factor(rep(c("ANT", "ESCC"), 2), levels = c("ESCC", "ANT")),
                  value = c(as.numeric(table(data$Predict, data$Real))))
data_m$color=paste(data_m$Real, data_m$Predict, sep = "_")
p2=ggplot(data_m, aes(x = Real, y = Predict)) +  
  theme_minimal() + #设置系统自带主题
  theme(panel.grid.major = element_blank()) +  #设置主项网格
  theme(legend.key = element_blank()) + #去掉背景颜色
  geom_tile(aes(fill = color)) +  #设置填充的值
  scale_fill_manual(values = c(ANT_ANT = "#2172AB", ANT_ESCC = "white", ESCC_ANT = "white", ESCC_ESCC = "#B14234")) +
  #scale_fill_gradient(low = "white", high = "red") + #设置颜色梯度
  geom_text(aes(label=value)) + #设置文本显示数值
  geom_vline(xintercept = c(0.5, 1.5, 2.5), size = 0.5) +
  geom_hline(yintercept = c(0.5, 1.5, 2.5), size = 0.5)
p2 #输出图片 

pdf(file="Figure3D.pdf", width=8, height=2.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

##E
library(pheatmap)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 9, skip = 1), row.names = 1, check.names = F)

annotation_col=data.frame(Real = data$Real, 
                          Predict = data$Predict,
                          Score = data$Risk_probability)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(ANT = "#d5c0d6", ESCC = "#7875a3"),
                Predict = c(ANT = "#f2c4ad", ESCC = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-2, 2, length=100)))

p=pheatmap(t(data[order(data$score), 1:6]),
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

pdf(file="Figure3E.pdf", width=8, height=4)
p
dev.off()

##F
library(pROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 10, skip = 1), row.names = 1)

pdf(file="Figure3F.pdf", width=4, height=4)
roc=roc(response = data$Real, 
        predictor = data$score,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 
dev.off()

#Figure 4
##A
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 11, skip = 1), row.names = 1)

df=data.frame(value = as.numeric(table(data$Real)), group = names(table(data$Real)))
p1=ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = NA) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#2172AB", "#B14234")) +
  theme_void()

data_m=data.frame(Real = factor(c(rep("ANT", 2), rep("ESCC", 2)), levels = c("ANT", "ESCC")),
                  Predict = factor(rep(c("ANT", "ESCC"), 2), levels = c("ESCC", "ANT")),
                  value = c(as.numeric(table(data$Predict, data$Real))))
data_m$color=paste(data_m$Real, data_m$Predict, sep = "_")
p2=ggplot(data_m, aes(x = Real, y = Predict)) +  
  theme_minimal() + #设置系统自带主题
  theme(panel.grid.major = element_blank()) +  #设置主项网格
  theme(legend.key = element_blank()) + #去掉背景颜色
  geom_tile(aes(fill = color)) +  #设置填充的值
  scale_fill_manual(values = c(ANT_ANT = "#2172AB", ANT_ESCC = "white", ESCC_ANT = "white", ESCC_ESCC = "#B14234")) +
  #scale_fill_gradient(low = "white", high = "red") + #设置颜色梯度
  geom_text(aes(label=value)) + #设置文本显示数值
  geom_vline(xintercept = c(0.5, 1.5, 2.5), size = 0.5) +
  geom_hline(yintercept = c(0.5, 1.5, 2.5), size = 0.5)
p2 #输出图片 

pdf(file="Figure4A.pdf", width=8, height=2.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

##B
library(pheatmap)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 12, skip = 1), row.names = 1, check.names = F)

annotation_col=data.frame(Real = data$Real, 
                          Predict = data$Predict,
                          Score = data$Risk_probability)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(ANT = "#d5c0d6", ESCC = "#7875a3"),
                Predict = c(ANT = "#f2c4ad", ESCC = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-2, 2, length=100)))

p=pheatmap(t(data[order(data$score), 1:6]),
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

pdf(file="Figure4B.pdf", width=8, height=4)
p
dev.off()

##C
library(pROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 13, skip = 1), row.names = 1)

pdf(file="Figure4C.pdf", width=4, height=4)
roc=roc(response = data$Real, 
        predictor = data$score,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 
dev.off()

##D
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 14, skip = 1), row.names = 1)
data$Real=factor(data$Real, levels = c("HD", "ESCC"))
data$Predict=factor(data$Predict, levels = c("HD", "ESCC"))

df=data.frame(value = as.numeric(table(data$Real)), group = names(table(data$Real)))
df$group=factor(df$group, levels = c("HD", "ESCC"))
p1=ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = NA) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#2172AB", "#B14234")) +
  theme_void()

data_m=data.frame(Real = factor(c(rep("HD", 2), rep("ESCC", 2)), levels = c("HD", "ESCC")),
                  Predict = factor(rep(c("HD", "ESCC"), 2), levels = c("ESCC", "HD")),
                  value = c(as.numeric(table(data$Predict, data$Real))))
data_m$color=paste(data_m$Real, data_m$Predict, sep = "_")
p2=ggplot(data_m, aes(x = Real, y = Predict)) +  
  theme_minimal() + #设置系统自带主题
  theme(panel.grid.major = element_blank()) +  #设置主项网格
  theme(legend.key = element_blank()) + #去掉背景颜色
  geom_tile(aes(fill = color)) +  #设置填充的值
  scale_fill_manual(values = c(HD_HD = "#2172AB", HD_ESCC = "white", ESCC_HD = "white", ESCC_ESCC = "#B14234")) +
  #scale_fill_gradient(low = "white", high = "red") + #设置颜色梯度
  geom_text(aes(label=value)) + #设置文本显示数值
  geom_vline(xintercept = c(0.5, 1.5, 2.5), size = 0.5) +
  geom_hline(yintercept = c(0.5, 1.5, 2.5), size = 0.5)
p2 #输出图片 

pdf(file="Figure4D.pdf", width=8, height=2.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

##E
library(pheatmap)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 15, skip = 1), row.names = 1, check.names = F)

annotation_col=data.frame(Real = data$Real, 
                          Predict = data$Predict,
                          Score = data$Risk_probability)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(HD = "#d5c0d6", ESCC = "#7875a3"),
                Predict = c(HD = "#f2c4ad", ESCC = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-2, 2, length=100)))

p=pheatmap(t(data[order(data$score), 1:6]),
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

pdf(file="Figure4E.pdf", width=8, height=4)
p
dev.off()

##F
library(pROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 16, skip = 1), row.names = 1)

pdf(file="Figure3F.pdf", width=4, height=4)
roc=roc(response = data$Real, 
        predictor = data$score,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 
dev.off()

##G
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 17, skip = 1), row.names = 1)

df=data.frame(value = as.numeric(table(data$Real)), group = names(table(data$Real)))
p1=ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = NA) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#2172AB", "#B14234")) +
  theme_void()

data_m=data.frame(Real = factor(c(rep("ANT", 2), rep("ESCC", 2)), levels = c("ANT", "ESCC")),
                  Predict = factor(rep(c("ANT", "ESCC"), 2), levels = c("ESCC", "ANT")),
                  value = c(as.numeric(table(data$Predict, data$Real))))
data_m$color=paste(data_m$Real, data_m$Predict, sep = "_")
p2=ggplot(data_m, aes(x = Real, y = Predict)) +  
  theme_minimal() + #设置系统自带主题
  theme(panel.grid.major = element_blank()) +  #设置主项网格
  theme(legend.key = element_blank()) + #去掉背景颜色
  geom_tile(aes(fill = color)) +  #设置填充的值
  scale_fill_manual(values = c(ANT_ANT = "#2172AB", ANT_ESCC = "white", ESCC_ANT = "white", ESCC_ESCC = "#B14234")) +
  #scale_fill_gradient(low = "white", high = "red") + #设置颜色梯度
  geom_text(aes(label=value)) + #设置文本显示数值
  geom_vline(xintercept = c(0.5, 1.5, 2.5), size = 0.5) +
  geom_hline(yintercept = c(0.5, 1.5, 2.5), size = 0.5)
p2 #输出图片 

pdf(file="Figure4G.pdf", width=8, height=2.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

##H
library(pheatmap)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 18, skip = 1), row.names = 1, check.names = F)

annotation_col=data.frame(Real = data$Real, 
                          Predict = data$Predict,
                          Score = data$Risk_probability)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(ANT = "#d5c0d6", ESCC = "#7875a3"),
                Predict = c(ANT = "#f2c4ad", ESCC = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-2, 2, length=100)))

p=pheatmap(t(data[order(data$score), 1:6]),
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

pdf(file="Figure4H.pdf", width=8, height=4)
p
dev.off()

##I
library(pROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 19, skip = 1), row.names = 1)

pdf(file="Figure4I.pdf", width=4, height=4)
roc=roc(response = data$Real, 
        predictor = data$score,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 
dev.off()

##J
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 20, skip = 1), row.names = 1)

df=data.frame(value = as.numeric(table(data$Real)), group = names(table(data$Real)))
p1=ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = NA) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#2172AB", "#B14234")) +
  theme_void()

data_m=data.frame(Real = factor(c(rep("ANT", 2), rep("ESCC", 2)), levels = c("ANT", "ESCC")),
                  Predict = factor(rep(c("ANT", "ESCC"), 2), levels = c("ESCC", "ANT")),
                  value = c(as.numeric(table(data$Predict, data$Real))))
data_m$color=paste(data_m$Real, data_m$Predict, sep = "_")
p2=ggplot(data_m, aes(x = Real, y = Predict)) +  
  theme_minimal() + #设置系统自带主题
  theme(panel.grid.major = element_blank()) +  #设置主项网格
  theme(legend.key = element_blank()) + #去掉背景颜色
  geom_tile(aes(fill = color)) +  #设置填充的值
  scale_fill_manual(values = c(ANT_ANT = "#2172AB", ANT_ESCC = "white", ESCC_ANT = "white", ESCC_ESCC = "#B14234")) +
  #scale_fill_gradient(low = "white", high = "red") + #设置颜色梯度
  geom_text(aes(label=value)) + #设置文本显示数值
  geom_vline(xintercept = c(0.5, 1.5, 2.5), size = 0.5) +
  geom_hline(yintercept = c(0.5, 1.5, 2.5), size = 0.5)
p2 #输出图片 

pdf(file="Figure4J.pdf", width=8, height=2.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

##K
library(pheatmap)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 21, skip = 1), row.names = 1, check.names = F)

annotation_col=data.frame(Real = data$Real, 
                          Predict = data$Predict,
                          Score = data$Risk_probability)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(ANT = "#d5c0d6", ESCC = "#7875a3"),
                Predict = c(ANT = "#f2c4ad", ESCC = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-2, 2, length=100)))

p=pheatmap(t(data[order(data$score), 1:6]),
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

pdf(file="Figure4K.pdf", width=8, height=4)
p
dev.off()

##L
library(pROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 22, skip = 1), row.names = 1)

pdf(file="Figure4L.pdf", width=4, height=4)
roc=roc(response = data$Real, 
        predictor = data$score,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 
dev.off()

#Figure 5
##A
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 23, skip = 1), row.names = 1, check.names = F)

df=data.frame(value = as.numeric(table(data$Real)), group = names(table(data$Real)))
p1=ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = NA) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#2172AB", "#B14234")) +
  theme_void()

data_m=data.frame(Real = factor(c(rep("ANT", 2), rep("I & II", 2)), levels = c("ANT", "I & II")),
                  Predict = factor(rep(c("ANT", "I & II"), 2), levels = c("I & II", "ANT")),
                  value = c(as.numeric(table(data$Predict, data$Real))))
data_m$color=paste(data_m$Real, data_m$Predict, sep = "_")
p2=ggplot(data_m, aes(x = Real, y = Predict)) +  
  theme_minimal() + #设置系统自带主题
  theme(panel.grid.major = element_blank()) +  #设置主项网格
  theme(legend.key = element_blank()) + #去掉背景颜色
  geom_tile(aes(fill = color)) +  #设置填充的值
  scale_fill_manual(values = c(ANT_ANT = "#2172AB", `ANT_I & II` = "white", `I & II_ANT` = "white", `I & II_I & II` = "#B14234")) +
  #scale_fill_gradient(low = "white", high = "red") + #设置颜色梯度
  geom_text(aes(label=value)) + #设置文本显示数值
  geom_vline(xintercept = c(0.5, 1.5, 2.5), size = 0.5) +
  geom_hline(yintercept = c(0.5, 1.5, 2.5), size = 0.5)
p2 #输出图片 

pdf(file="Figure5A.1.pdf", width=8, height=2.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

annotation_col=data.frame(Real = data$Real, 
                          Predict = data$Predict,
                          Score = data$Risk_probability)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(ANT = "#d5c0d6", `I & II` = "#7875a3"),
                Predict = c(ANT = "#f2c4ad", `I & II` = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-2, 2, length=100)))

p3=pheatmap(t(data[order(data$score), 1:6]),
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

pdf(file="Figure5A.2.pdf", width=8, height=4)
p3
dev.off()

pdf(file="Figure5A.3.pdf", width=4, height=4)
roc=roc(response = data$Real, 
        predictor = data$score,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 
dev.off()

##B
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 24, skip = 1), row.names = 1, check.names = F)

df=data.frame(value = as.numeric(table(data$Real)), group = names(table(data$Real)))
p1=ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = NA) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#2172AB", "#B14234")) +
  theme_void()

data_m=data.frame(Real = factor(c(rep("ANT", 2), rep("I & II", 2)), levels = c("ANT", "I & II")),
                  Predict = factor(rep(c("ANT", "I & II"), 2), levels = c("I & II", "ANT")),
                  value = c(as.numeric(table(data$Predict, data$Real))))
data_m$color=paste(data_m$Real, data_m$Predict, sep = "_")
p2=ggplot(data_m, aes(x = Real, y = Predict)) +  
  theme_minimal() + #设置系统自带主题
  theme(panel.grid.major = element_blank()) +  #设置主项网格
  theme(legend.key = element_blank()) + #去掉背景颜色
  geom_tile(aes(fill = color)) +  #设置填充的值
  scale_fill_manual(values = c(ANT_ANT = "#2172AB", `ANT_I & II` = "white", `I & II_ANT` = "white", `I & II_I & II` = "#B14234")) +
  #scale_fill_gradient(low = "white", high = "red") + #设置颜色梯度
  geom_text(aes(label=value)) + #设置文本显示数值
  geom_vline(xintercept = c(0.5, 1.5, 2.5), size = 0.5) +
  geom_hline(yintercept = c(0.5, 1.5, 2.5), size = 0.5)
p2 #输出图片 

pdf(file="Figure5B.1.pdf", width=8, height=2.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

annotation_col=data.frame(Real = data$Real, 
                          Predict = data$Predict,
                          Score = data$Risk_probability)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(ANT = "#d5c0d6", `I & II` = "#7875a3"),
                Predict = c(ANT = "#f2c4ad", `I & II` = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-2, 2, length=100)))

p3=pheatmap(t(data[order(data$score), 1:6]),
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

pdf(file="Figure5B.2.pdf", width=8, height=4)
p3
dev.off()

pdf(file="Figure5B.3.pdf", width=4, height=4)
roc=roc(response = data$Real, 
        predictor = data$score,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 
dev.off()

##C
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 25, skip = 1), row.names = 1, check.names = F)

df=data.frame(value = as.numeric(table(data$Real)), group = names(table(data$Real)))
p1=ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = NA) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#2172AB", "#B14234")) +
  theme_void()

data_m=data.frame(Real = factor(c(rep("ANT", 2), rep("I & II", 2)), levels = c("ANT", "I & II")),
                  Predict = factor(rep(c("ANT", "I & II"), 2), levels = c("I & II", "ANT")),
                  value = c(as.numeric(table(data$Predict, data$Real))))
data_m$color=paste(data_m$Real, data_m$Predict, sep = "_")
p2=ggplot(data_m, aes(x = Real, y = Predict)) +  
  theme_minimal() + #设置系统自带主题
  theme(panel.grid.major = element_blank()) +  #设置主项网格
  theme(legend.key = element_blank()) + #去掉背景颜色
  geom_tile(aes(fill = color)) +  #设置填充的值
  scale_fill_manual(values = c(ANT_ANT = "#2172AB", `ANT_I & II` = "white", `I & II_ANT` = "white", `I & II_I & II` = "#B14234")) +
  #scale_fill_gradient(low = "white", high = "red") + #设置颜色梯度
  geom_text(aes(label=value)) + #设置文本显示数值
  geom_vline(xintercept = c(0.5, 1.5, 2.5), size = 0.5) +
  geom_hline(yintercept = c(0.5, 1.5, 2.5), size = 0.5)
p2 #输出图片 

pdf(file="Figure5C.1.pdf", width=8, height=2.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

annotation_col=data.frame(Real = data$Real, 
                          Predict = data$Predict,
                          Score = data$Risk_probability)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(ANT = "#d5c0d6", `I & II` = "#7875a3"),
                Predict = c(ANT = "#f2c4ad", `I & II` = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-2, 2, length=100)))

p3=pheatmap(t(data[order(data$score), 1:6]),
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

pdf(file="Figure5C.2.pdf", width=8, height=4)
p3
dev.off()

pdf(file="Figure5C.3.pdf", width=4, height=4)
roc=roc(response = data$Real, 
        predictor = data$score,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 
dev.off()

##D
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 26, skip = 1), row.names = 1, check.names = F)

df=data.frame(value = as.numeric(table(data$Real)), group = names(table(data$Real)))
p1=ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = NA) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#2172AB", "#B14234")) +
  theme_void()

data_m=data.frame(Real = factor(c(rep("ANT", 2), rep("I & II", 2)), levels = c("ANT", "I & II")),
                  Predict = factor(rep(c("ANT", "I & II"), 2), levels = c("I & II", "ANT")),
                  value = c(as.numeric(table(data$Predict, data$Real))))
data_m$color=paste(data_m$Real, data_m$Predict, sep = "_")
p2=ggplot(data_m, aes(x = Real, y = Predict)) +  
  theme_minimal() + #设置系统自带主题
  theme(panel.grid.major = element_blank()) +  #设置主项网格
  theme(legend.key = element_blank()) + #去掉背景颜色
  geom_tile(aes(fill = color)) +  #设置填充的值
  scale_fill_manual(values = c(ANT_ANT = "#2172AB", `ANT_I & II` = "white", `I & II_ANT` = "white", `I & II_I & II` = "#B14234")) +
  #scale_fill_gradient(low = "white", high = "red") + #设置颜色梯度
  geom_text(aes(label=value)) + #设置文本显示数值
  geom_vline(xintercept = c(0.5, 1.5, 2.5), size = 0.5) +
  geom_hline(yintercept = c(0.5, 1.5, 2.5), size = 0.5)
p2 #输出图片 

pdf(file="Figure5D.1.pdf", width=8, height=2.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

annotation_col=data.frame(Real = data$Real, 
                          Predict = data$Predict,
                          Score = data$Risk_probability)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(ANT = "#d5c0d6", `I & II` = "#7875a3"),
                Predict = c(ANT = "#f2c4ad", `I & II` = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-2, 2, length=100)))

p3=pheatmap(t(data[order(data$score), 1:6]),
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

pdf(file="Figure5D.2.pdf", width=8, height=4)
p3
dev.off()

pdf(file="Figure5D.3.pdf", width=4, height=4)
roc=roc(response = data$Real, 
        predictor = data$score,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 
dev.off()

##E
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 27, skip = 1), row.names = 1, check.names = F)

df=data.frame(value = as.numeric(table(data$Real)), group = names(table(data$Real)))
p1=ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = NA) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#2172AB", "#B14234")) +
  theme_void()

data_m=data.frame(Real = factor(c(rep("ANT", 2), rep("I & II", 2)), levels = c("ANT", "I & II")),
                  Predict = factor(rep(c("ANT", "I & II"), 2), levels = c("I & II", "ANT")),
                  value = c(as.numeric(table(data$Predict, data$Real))))
data_m$color=paste(data_m$Real, data_m$Predict, sep = "_")
p2=ggplot(data_m, aes(x = Real, y = Predict)) +  
  theme_minimal() + #设置系统自带主题
  theme(panel.grid.major = element_blank()) +  #设置主项网格
  theme(legend.key = element_blank()) + #去掉背景颜色
  geom_tile(aes(fill = color)) +  #设置填充的值
  scale_fill_manual(values = c(ANT_ANT = "#2172AB", `ANT_I & II` = "white", `I & II_ANT` = "white", `I & II_I & II` = "#B14234")) +
  #scale_fill_gradient(low = "white", high = "red") + #设置颜色梯度
  geom_text(aes(label=value)) + #设置文本显示数值
  geom_vline(xintercept = c(0.5, 1.5, 2.5), size = 0.5) +
  geom_hline(yintercept = c(0.5, 1.5, 2.5), size = 0.5)
p2 #输出图片 

pdf(file="Figure5E.1.pdf", width=8, height=2.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

annotation_col=data.frame(Real = data$Real, 
                          Predict = data$Predict,
                          Score = data$Risk_probability)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(ANT = "#d5c0d6", `I & II` = "#7875a3"),
                Predict = c(ANT = "#f2c4ad", `I & II` = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-2, 2, length=100)))

p3=pheatmap(t(data[order(data$score), 1:6]),
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

pdf(file="Figure5E.2.pdf", width=8, height=4)
p3
dev.off()

pdf(file="Figure5E.3.pdf", width=4, height=4)
roc=roc(response = data$Real, 
        predictor = data$score,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 
dev.off()

##F
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 28, skip = 1), row.names = 1, check.names = F)

df=data.frame(value = as.numeric(table(data$Real)), group = names(table(data$Real)))
p1=ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = NA) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#2172AB", "#B14234")) +
  theme_void()

data_m=data.frame(Real = factor(c(rep("ANT", 2), rep("I & II", 2)), levels = c("ANT", "I & II")),
                  Predict = factor(rep(c("ANT", "I & II"), 2), levels = c("I & II", "ANT")),
                  value = c(as.numeric(table(data$Predict, data$Real))))
data_m$color=paste(data_m$Real, data_m$Predict, sep = "_")
p2=ggplot(data_m, aes(x = Real, y = Predict)) +  
  theme_minimal() + #设置系统自带主题
  theme(panel.grid.major = element_blank()) +  #设置主项网格
  theme(legend.key = element_blank()) + #去掉背景颜色
  geom_tile(aes(fill = color)) +  #设置填充的值
  scale_fill_manual(values = c(ANT_ANT = "#2172AB", `ANT_I & II` = "white", `I & II_ANT` = "white", `I & II_I & II` = "#B14234")) +
  #scale_fill_gradient(low = "white", high = "red") + #设置颜色梯度
  geom_text(aes(label=value)) + #设置文本显示数值
  geom_vline(xintercept = c(0.5, 1.5, 2.5), size = 0.5) +
  geom_hline(yintercept = c(0.5, 1.5, 2.5), size = 0.5)
p2 #输出图片 

pdf(file="Figure5F.1.pdf", width=8, height=2.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

annotation_col=data.frame(Real = data$Real, 
                          Predict = data$Predict,
                          Score = data$Risk_probability)
rownames(annotation_col)=rownames(data)
ann_colors=list(Real = c(ANT = "#d5c0d6", `I & II` = "#7875a3"),
                Predict = c(ANT = "#f2c4ad", `I & II` = "#d0866f"),
                Score = c(`1` = "#006837", `2` = "#1A9850", `3` = "#66BD63",
                          `4` = "#A6D96A", `5` = "#D9EF8B", `6` = "#FEE08B",
                          `7` = "#FDAE61", `8` = "#F46D43", `9` = "#D73027",
                          `10` = "#A50026"))
bk=unique(c(seq(-2, 2, length=100)))

p3=pheatmap(t(data[order(data$score), 1:6]),
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

pdf(file="Figure5F.2.pdf", width=8, height=4)
p3
dev.off()

pdf(file="Figure5F.3.pdf", width=4, height=4)
roc=roc(response = data$Real, 
        predictor = data$score,
        ci=TRUE, ci.alpha=0.05, print.auc=TRUE, plot=TRUE)
roc_ci=ci.sp(roc, sensitivities = seq(0, 1, .001), boot.n = 1000) 
plot(roc_ci, type = "shape", col = "lightblue", alpha = 0.5) 
dev.off()

#Figure 6
##A
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 29, skip = 1), row.names = 1, check.names = F)
data$label=factor(data$label, levels = c("HD", "EIN", "ESCC"))

plot_column_boxplot = function(data, column) {
  ylim=c(boxplot.stats(data[data$label == "HD", column])$stats[c(1,5)],
         boxplot.stats(data[data$label == "EIN", column])$stats[c(1,5)],
         boxplot.stats(data[data$label == "ESCC", column])$stats[c(1,5)])
  ggplot(data = data, aes(x = label, y = get(column))) +
    geom_boxplot(width=0.7) +
    geom_point(aes(colour = label), size = 1.5, position = position_jitterdodge(jitter.width = 1))  + ##在箱线图的基础上叠加散点图
    scale_color_manual(values = c(HD = "#258279", ESCC = "#B80200", EIN ="#E7AF30")) +
    scale_fill_manual(values = c(HD = "#258279", ESCC = "#B80200", EIN ="#E7AF30")) +
    coord_cartesian(ylim = c(min(ylim), max(ylim))*1.05)+  
    stat_compare_means(method = "wilcox.test", ref.group = "HD") +
    labs(title = column, y = "Expression value") +
    theme_bw() +
    theme(legend.position="none", ##隐藏图例
          #axis.text.y = element_blank(), 
          panel.grid = element_blank(), axis.text=element_text(size = 10),
          axis.title = element_text(size=10), legend.text = element_text(size=10))
}

pdf(file="Figure6A.pdf", width=10, height=3)
myplots_boxplot=lapply(colnames(data)[1:5], plot_column_boxplot, data = data)
grid.arrange(myplots_boxplot[[1]], myplots_boxplot[[2]], myplots_boxplot[[3]], 
             myplots_boxplot[[4]], myplots_boxplot[[5]], ncol = 5)
dev.off()

##B
library(pROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 30, skip = 1), row.names = 1)
data$label=factor(data$label, levels = c("HD", "ESCC"))

pdf(file="Figure6B.pdf", width=8, height=8)
r1=plot.roc(label ~ AP003548.1, data = data, col = "#C0342C", direction = ">", lwd = 2, legacy.axes = T)
r2=lines.roc(label ~ PGM5.AS1, data = data, col = "#191C3A", direction = ">", lwd = 2)
r3=lines.roc(label ~ ADAMTS9.AS1, data = data, col = "#6F863F", direction = ">", lwd = 2)
r4=lines.roc(label ~ LINC01082, data = data, col = "#ECA55B", direction = ">", lwd = 2)
r5=lines.roc(label ~ LINC03016, data = data, col = "#316CA1", direction = ">", lwd = 2)
r6=lines.roc(label ~ SCC, data = data, col = "#805D9D", direction = "<", lwd = 1.25, lty = 2)
r7=lines.roc(label ~ CEA, data = data, col = "#489238", direction = "<", lwd = 1.25, lty = 2)
r8=lines.roc(label ~ CYFRA21.1, data = data, col = "#7B7A7C", direction = "<", lwd = 1.25, lty = 2)

legend('bottomright', 
       c(paste('AP003548.1', "AUC =", round(as.numeric(r1$auc), 3), "(", round(ci.auc(r1)[1], 3), "-", round(ci.auc(r1)[3], 3), ")"),
         paste('PGM5-AS1', "AUC =", round(as.numeric(r2$auc),3), "(", round(ci.auc(r2)[1], 3), "-", round(ci.auc(r2)[3], 3), ")"),
         paste('ADAMTS9-AS1', "AUC =", round(as.numeric(r3$auc),3), "(", round(ci.auc(r3)[1], 3), "-", round(ci.auc(r3)[3], 3), ")"),
         paste('LINC01082', "AUC =", round(as.numeric(r4$auc),3), "(", round(ci.auc(r4)[1], 3), "-", round(ci.auc(r4)[3], 3), ")"),
         paste('LINC03016', "AUC =", round(as.numeric(r5$auc),3), "(", round(ci.auc(r5)[1], 3), "-", round(ci.auc(r5)[3], 3), ")"),
         paste('SCC-Ag', "AUC =", round(as.numeric(r6$auc),3), "(", round(ci.auc(r6)[1], 3), "-", round(ci.auc(r6)[3], 3), ")"),
         paste('CEA', "AUC =", round(as.numeric(r7$auc),3), "(", round(ci.auc(r7)[1], 3), "-", round(ci.auc(r7)[3], 3), ")"),
         paste('CYFRA21-1', "AUC =", round(as.numeric(r8$auc),3), "(", round(ci.auc(r8)[1], 3), "-", round(ci.auc(r8)[3], 3), ")")), 
       lty = c(rep(1, 5), rep(2, 3)), text.font = 1, 
       col=c("#C0342C", "#191C3A", "#6F863F", "#ECA55B", "#316CA1", "#805D9D", "#489238", "#7B7A7C"))
dev.off()

##C
library(pROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 31, skip = 1), row.names = 1)
data$label=factor(data$label, levels = c("HD", "EIN"))

pdf(file="Figure6B.pdf", width=8, height=8)
r1=plot.roc(label ~ AP003548.1, data = data, col = "#C0342C", direction = ">", lwd = 2, legacy.axes = T)
r2=lines.roc(label ~ PGM5.AS1, data = data, col = "#191C3A", direction = ">", lwd = 2)
r3=lines.roc(label ~ ADAMTS9.AS1, data = data, col = "#6F863F", direction = ">", lwd = 2)
r4=lines.roc(label ~ LINC01082, data = data, col = "#ECA55B", direction = ">", lwd = 2)
r5=lines.roc(label ~ LINC03016, data = data, col = "#316CA1", direction = ">", lwd = 2)
r6=lines.roc(label ~ SCC, data = data, col = "#805D9D", direction = "<", lwd = 1.25, lty = 2)
r7=lines.roc(label ~ CEA, data = data, col = "#489238", direction = "<", lwd = 1.25, lty = 2)
r8=lines.roc(label ~ CYFRA21.1, data = data, col = "#7B7A7C", direction = "<", lwd = 1.25, lty = 2)

legend('bottomright', 
       c(paste('AP003548.1', "AUC =", round(as.numeric(r1$auc), 3), "(", round(ci.auc(r1)[1], 3), "-", round(ci.auc(r1)[3], 3), ")"),
         paste('PGM5-AS1', "AUC =", round(as.numeric(r2$auc),3), "(", round(ci.auc(r2)[1], 3), "-", round(ci.auc(r2)[3], 3), ")"),
         paste('ADAMTS9-AS1', "AUC =", round(as.numeric(r3$auc),3), "(", round(ci.auc(r3)[1], 3), "-", round(ci.auc(r3)[3], 3), ")"),
         paste('LINC01082', "AUC =", round(as.numeric(r4$auc),3), "(", round(ci.auc(r4)[1], 3), "-", round(ci.auc(r4)[3], 3), ")"),
         paste('LINC03016', "AUC =", round(as.numeric(r5$auc),3), "(", round(ci.auc(r5)[1], 3), "-", round(ci.auc(r5)[3], 3), ")"),
         paste('SCC-Ag', "AUC =", round(as.numeric(r6$auc),3), "(", round(ci.auc(r6)[1], 3), "-", round(ci.auc(r6)[3], 3), ")"),
         paste('CEA', "AUC =", round(as.numeric(r7$auc),3), "(", round(ci.auc(r7)[1], 3), "-", round(ci.auc(r7)[3], 3), ")"),
         paste('CYFRA21-1', "AUC =", round(as.numeric(r8$auc),3), "(", round(ci.auc(r8)[1], 3), "-", round(ci.auc(r8)[3], 3), ")")), 
       lty = c(rep(1, 5), rep(2, 3)), text.font = 1, 
       col=c("#C0342C", "#191C3A", "#6F863F", "#ECA55B", "#316CA1", "#805D9D", "#489238", "#7B7A7C"))
dev.off()

#Figure 7
library(pROC)
library(customLayout)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 32, skip = 1), row.names = 1)

dca <- function(data, outcome, predictors, xstart=0.01, xstop=0.99, xby=0.01, 
                ymin=-0.05, probability=NULL, harm=NULL,graph=TRUE, intervention=FALSE, title = "DCA",
                interventionper=100, smooth=FALSE,loess.span=0.10) {
  
  # LOADING REQUIRED LIBRARIES
  require(stats)
  
  # data MUST BE A DATA FRAME
  if (class(data)!="data.frame") {
    stop("Input data must be class data.frame")
  }
  
  #ONLY KEEPING COMPLETE CASES
  data=data[complete.cases(data[append(outcome,predictors)]),append(outcome,predictors)]
  
  # outcome MUST BE CODED AS 0 AND 1
  if (max(data[[outcome]])>1 | min(data[[outcome]])<0) {
    stop("outcome cannot be less than 0 or greater than 1")
  }
  # xstart IS BETWEEN 0 AND 1
  if (xstart<0 | xstart>1) {
    stop("xstart must lie between 0 and 1")
  }
  
  # xstop IS BETWEEN 0 AND 1
  if (xstop<0 | xstop>1) {
    stop("xstop must lie between 0 and 1")
  }
  
  # xby IS BETWEEN 0 AND 1
  if (xby<=0 | xby>=1) {
    stop("xby must lie between 0 and 1")
  }
  
  # xstart IS BEFORE xstop
  if (xstart>=xstop) {
    stop("xstop must be larger than xstart")
  }
  
  #STORING THE NUMBER OF PREDICTORS SPECIFIED
  pred.n=length(predictors)
  
  #IF probability SPECIFIED ENSURING THAT EACH PREDICTOR IS INDICATED AS A YES OR NO
  if (length(probability)>0 & pred.n!=length(probability)) {
    stop("Number of probabilities specified must be the same as the number of predictors being checked.")
  }
  
  #IF harm SPECIFIED ENSURING THAT EACH PREDICTOR HAS A SPECIFIED HARM
  if (length(harm)>0 & pred.n!=length(harm)) {
    stop("Number of harms specified must be the same as the number of predictors being checked.")
  }
  
  #INITIALIZING DEFAULT VALUES FOR PROBABILITES AND HARMS IF NOT SPECIFIED
  if (length(harm)==0) {
    harm=rep(0,pred.n)
  }
  if (length(probability)==0) {
    probability=rep(TRUE,pred.n)
  }
  
  
  #CHECKING THAT EACH probability ELEMENT IS EQUAL TO YES OR NO, 
  #AND CHECKING THAT PROBABILITIES ARE BETWEEN 0 and 1
  #IF NOT A PROB THEN CONVERTING WITH A LOGISTIC REGRESSION
  for(m in 1:pred.n) { 
    if (probability[m]!=TRUE & probability[m]!=FALSE) {
      stop("Each element of probability vector must be TRUE or FALSE")
    }
    if (probability[m]==TRUE & (max(data[predictors[m]])>1 | min(data[predictors[m]])<0)) {
      stop(paste(predictors[m],"must be between 0 and 1 OR sepcified as a non-probability in the probability option",sep=" "))  
    }
    if(probability[m]==FALSE) {
      model=NULL
      pred=NULL
      model=glm(data.matrix(data[outcome]) ~ data.matrix(data[predictors[m]]), family=binomial("logit"))
      pred=data.frame(model$fitted.values)
      pred=data.frame(pred)
      names(pred)=predictors[m]
      data=cbind(data[names(data)!=predictors[m]],pred)
      print(paste(predictors[m],"converted to a probability with logistic regression. Due to linearity assumption, miscalibration may occur.",sep=" "))
    }
  }
  
  # THE PREDICTOR NAMES CANNOT BE EQUAL TO all OR none.
  if (length(predictors[predictors=="all" | predictors=="none"])) {
    stop("Prediction names cannot be equal to all or none.")
  }  
  
  #########  CALCULATING NET BENEFIT   #########
  N=dim(data)[1]
  event.rate=colMeans(data[outcome])
  
  # CREATING DATAFRAME THAT IS ONE LINE PER THRESHOLD PER all AND none STRATEGY
  nb=data.frame(seq(from=xstart, to=xstop, by=xby))
  names(nb)="threshold"
  interv=nb
  
  nb["all"]=event.rate - (1-event.rate)*nb$threshold/(1-nb$threshold)
  nb["none"]=0
  
  # CYCLING THROUGH EACH PREDICTOR AND CALCULATING NET BENEFIT
  for(m in 1:pred.n){
    for(t in 1:length(nb$threshold)){
      # COUNTING TRUE POSITIVES AT EACH THRESHOLD
      tp=mean(data[data[[predictors[m]]]>=nb$threshold[t],outcome])*sum(data[[predictors[m]]]>=nb$threshold[t])
      # COUNTING FALSE POSITIVES AT EACH THRESHOLD
      fp=(1-mean(data[data[[predictors[m]]]>=nb$threshold[t],outcome]))*sum(data[[predictors[m]]]>=nb$threshold[t])
      #setting TP and FP to 0 if no observations meet threshold prob.
      if (sum(data[[predictors[m]]]>=nb$threshold[t])==0) {
        tp=0
        fp=0
      }
      
      # CALCULATING NET BENEFIT
      nb[t,predictors[m]]=tp/N - fp/N*(nb$threshold[t]/(1-nb$threshold[t])) - harm[m]
    }
    interv[predictors[m]]=(nb[predictors[m]] - nb["all"])*interventionper/(interv$threshold/(1-interv$threshold))
  }
  
  # CYCLING THROUGH EACH PREDICTOR AND SMOOTH NET BENEFIT AND INTERVENTIONS AVOIDED 
  for(m in 1:pred.n) {
    if (smooth==TRUE){
      lws=loess(data.matrix(nb[!is.na(nb[[predictors[m]]]),predictors[m]]) ~ data.matrix(nb[!is.na(nb[[predictors[m]]]),"threshold"]),span=loess.span)
      nb[!is.na(nb[[predictors[m]]]),paste(predictors[m],"_sm",sep="")]=lws$fitted
      
      lws=loess(data.matrix(interv[!is.na(nb[[predictors[m]]]),predictors[m]]) ~ data.matrix(interv[!is.na(nb[[predictors[m]]]),"threshold"]),span=loess.span)
      interv[!is.na(nb[[predictors[m]]]),paste(predictors[m],"_sm",sep="")]=lws$fitted
    }
  }
  
  # PLOTTING GRAPH IF REQUESTED
  if (graph==TRUE) {
    require(graphics)
    
    # PLOTTING INTERVENTIONS AVOIDED IF REQUESTED
    if(intervention==TRUE) {
      # initialize the legend label, color, and width using the standard specs of the none and all lines
      legendlabel <- NULL
      legendcolor <- NULL
      legendwidth <- NULL
      legendpattern <- NULL
      
      #getting maximum number of avoided interventions
      ymax=max(interv[predictors],na.rm = TRUE)
      
      #INITIALIZING EMPTY PLOT WITH LABELS
      plot(x=nb$threshold, y=nb$all,type="n",xlim=c(xstart, xstop), ylim=c(ymin, ymax), xlab="Threshold probability", ylab=paste("Net reduction in interventions per",interventionper,"patients"), main=title)
      
      #PLOTTING INTERVENTIONS AVOIDED FOR EACH PREDICTOR
      for(m in 1:pred.n) {
        if (smooth==TRUE){
          lines(interv$threshold,data.matrix(interv[paste(predictors[m],"_sm",sep="")]),col=m,lty=2)
        } else {
          lines(interv$threshold,data.matrix(interv[predictors[m]]),col=m,lty=2)
        }
        
        # adding each model to the legend
        legendlabel <- c(legendlabel, predictors[m])
        legendcolor <- c(legendcolor, m)
        legendwidth <- c(legendwidth, 1)
        legendpattern <- c(legendpattern, 2)
      }
    } else {
      # PLOTTING NET BENEFIT IF REQUESTED
      
      # initialize the legend label, color, and width using the standard specs of the none and all lines
      legendlabel <- c("None", "All")
      #legendcolor <- c(17, 8)
      legendcolor <- c('black', 'grey')
      legendwidth <- c(3, 3)
      legendpattern <- c(1, 1)
      
      #getting maximum net benefit
      ymax=max(nb[names(nb)!="threshold"],na.rm = TRUE)
      
      # inializing new benfit plot with treat all option
      plot(x=nb$threshold, y=nb$all, type="l",col="grey", lwd=3 ,xlim=c(xstart, xstop), ylim=c(ymin, ymax), xlab="Threshold probability", ylab="Net benefit", main=title)
      # adding treat none option
      lines(x=nb$threshold, y=nb$none,lwd=3,col="black")
      #PLOTTING net benefit FOR EACH PREDICTOR
      for(m in 1:pred.n) {
        if (smooth==TRUE){
          lines(nb$threshold,data.matrix(nb[paste(predictors[m],"_sm",sep="")]),col='red',lty=1,lwd = 3) 
        } else {
          lines(nb$threshold,data.matrix(nb[predictors[m]]),col='red',lty=1,lwd = 3)
        }
        # adding each model to the legend
        legendlabel <- c(legendlabel, predictors[m])
        #legendcolor <- c(legendcolor, m)
        legendcolor <- c(legendcolor, 'red')
        legendwidth <- c(3, 3)
        #legendpattern <- c(legendpattern, 2)
        legendpattern <- c(legendpattern, 1)
      }
    }
    # then add the legend
    legend("topright", legendlabel, cex=0.8, col=legendcolor, lwd=legendwidth, lty=legendpattern)
    
  }
  
  #RETURNING RESULTS
  results=list() 
  results$N=N
  results$predictors=data.frame(cbind(predictors,harm,probability))
  names(results$predictors)=c("predictor","harm.applied","probability")
  results$interventions.avoided.per=interventionper
  results$net.benefit=nb
  results$interventions.avoided=interv
  
  return(results)
  
}  

pdf(file = "Figure7.pdf", width = 12, height = 7)
lay=lay_new(matrix(1:8, ncol = 4))
lay_set(lay)
a1=dca(data = subset(data, cohort == "train"), outcome = "label", predictors = "MLMRPscore", title = "SCH training cohort", smooth = "TRUE")
a2=dca(data = subset(data, cohort == "test"), outcome = "label", predictors = "MLMRPscore", title = "SCH validation cohort", smooth = "TRUE")
a3=dca(data = subset(data, cohort == "RTPCR"), outcome = "label", predictors = "MLMRPscore", title = "CAMS cohort", smooth = "TRUE")
a4=dca(data = subset(data, cohort == "GSE130078"), outcome = "label", predictors = "MLMRPscore", title = "You cohort", smooth = "TRUE")
a5=dca(data = subset(data, cohort == "TCGA_GTEx"), outcome = "label", predictors = "MLMRPscore", title = "TCGA-GTEx cohort", smooth = "TRUE")
a6=dca(data = subset(data, cohort == "GSE53624"), outcome = "label", predictors = "MLMRPscore", title = "Li cohort-1", smooth = "TRUE")
a7=dca(data = subset(data, cohort == "GSE53622"), outcome = "label", predictors = "MLMRPscore", title = "Li cohort-2", smooth = "TRUE")
dev.off()

#Supplementary figure S1
##A
library(PRROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 33, skip = 1), row.names = 1)

pdf(file = "Supplementary figure S1A.pdf", width = 5, height = 4.5)
pr=pr.curve(scores.class0 = data[data$label == "ESCC", "MLMRPscore"], 
            scores.class1 = data[data$label == "ANT", "MLMRPscore"], curve=T)
plot(pr)
dev.off()

##B
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 34, skip = 1), row.names = 1, check.names = F)

plot_column_boxplot = function(data, column) {
  ggplot(data = data, aes(x = tissue_type, y = get(column))) +
    theme_bw() +
    geom_boxplot(color="black") + 
    geom_point(aes(colour = tissue_type, shape = tissue_type), size = 2)  + 
    geom_line(aes(group = patient), size = 0.8, colour = "#9C9C9C") + 
    scale_color_manual(values = c(ANT = "#0000EE", ESCC = "#EE0000")) +
    scale_shape_manual(values = c(ANT = 16, ESCC = 16)) +
    stat_compare_means(method = "t.test", paired = T) +
    labs(title = column, y = "Expression value") +
    theme(legend.position = "none") + 
    theme(panel.grid = element_blank(), axis.text=element_text(size = 12),
          axis.title = element_text(size = 12), legend.text = element_text(size = 12))
}

pdf(file="Supplementary figure S1B.pdf", width=7, height=7)
myplots_boxplot=lapply(colnames(data)[1:6], plot_column_boxplot, data = data)
grid.arrange(myplots_boxplot[[1]], myplots_boxplot[[2]], myplots_boxplot[[3]], 
             myplots_boxplot[[4]], myplots_boxplot[[5]], myplots_boxplot[[6]],  ncol = 3)
dev.off()

#Supplementary figure S2
##A
library(PRROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 35, skip = 1), row.names = 1)

pdf(file = "Supplementary figure S2A.pdf", width = 5, height = 4.5)
pr=pr.curve(scores.class0 = data[data$label == "ESCC", "MLMRPscore"], 
            scores.class1 = data[data$label == "ANT", "MLMRPscore"], curve=T)
plot(pr)
dev.off()

##B
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 36, skip = 1), row.names = 1, check.names = F)

plot_column_boxplot = function(data, column) {
  ggplot(data = data, aes(x = tissue_type, y = get(column))) +
    theme_bw() +
    geom_boxplot(color="black") + 
    geom_point(aes(colour = tissue_type, shape = tissue_type), size = 2)  + 
    geom_line(aes(group = patient), size = 0.8, colour = "#9C9C9C") + 
    scale_color_manual(values = c(ANT = "#0000EE", ESCC = "#EE0000")) +
    scale_shape_manual(values = c(ANT = 16, ESCC = 16)) +
    stat_compare_means(method = "t.test", paired = T) +
    labs(title = column, y = "Expression value") +
    theme(legend.position = "none") + 
    theme(panel.grid = element_blank(), axis.text=element_text(size = 12),
          axis.title = element_text(size = 12), legend.text = element_text(size = 12))
}

pdf(file="Supplementary figure S2B.pdf", width=7, height=7)
myplots_boxplot=lapply(colnames(data)[1:6], plot_column_boxplot, data = data)
grid.arrange(myplots_boxplot[[1]], myplots_boxplot[[2]], myplots_boxplot[[3]], 
             myplots_boxplot[[4]], myplots_boxplot[[5]], myplots_boxplot[[6]],  ncol = 3)
dev.off()

#Supplementary figure S3
##A
library(PRROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 37, skip = 1), row.names = 1)

pdf(file = "Supplementary figure S3A.pdf", width = 5, height = 4.5)
pr=pr.curve(scores.class0 = data[data$label == "ESCC", "MLMRPscore"], 
            scores.class1 = data[data$label == "ANT", "MLMRPscore"], curve=T)
plot(pr)
dev.off()

##B
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 38, skip = 1), row.names = 1, check.names = F)

plot_column_boxplot = function(data, column) {
  ggplot(data = data, aes(x = tissue_type, y = get(column))) +
    theme_bw() +
    geom_boxplot(color="black") + 
    geom_point(aes(colour = tissue_type, shape = tissue_type), size = 2)  + 
    geom_line(aes(group = patient), size = 0.8, colour = "#9C9C9C") + 
    scale_color_manual(values = c(ANT = "#0000EE", ESCC = "#EE0000")) +
    scale_shape_manual(values = c(ANT = 16, ESCC = 16)) +
    stat_compare_means(method = "t.test", paired = T) +
    labs(title = column, y = "Expression value") +
    theme(legend.position = "none") + 
    theme(panel.grid = element_blank(), axis.text=element_text(size = 12),
          axis.title = element_text(size = 12), legend.text = element_text(size = 12))
}

pdf(file="Supplementary figure S3B.pdf", width=7, height=7)
myplots_boxplot=lapply(colnames(data)[1:6], plot_column_boxplot, data = data)
grid.arrange(myplots_boxplot[[1]], myplots_boxplot[[2]], myplots_boxplot[[3]], 
             myplots_boxplot[[4]], myplots_boxplot[[5]], myplots_boxplot[[6]],  ncol = 3)
dev.off()

##C
library(PRROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 39, skip = 1), row.names = 1)

pdf(file = "Supplementary figure S3C.pdf", width = 5, height = 4.5)
pr=pr.curve(scores.class0 = data[data$label == "ESCC", "MLMRPscore"], 
            scores.class1 = data[data$label == "ANT", "MLMRPscore"], curve=T)
plot(pr)
dev.off()

##D
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 40, skip = 1), row.names = 1, check.names = F)
data$tissue_type=factor(data$tissue_type, levels = c("HD", "ESCC"))

plot_column_boxplot = function(data, column) {
  ylim=c(boxplot.stats(data[data$tissue_type == "HD", column])$stats[c(1,5)],
         boxplot.stats(data[data$tissue_type == "ESCC", column])$stats[c(1,5)])
  ggplot(data = data, aes(x = tissue_type, y = get(column))) +
    theme_bw() +
    geom_boxplot(color="black") + 
    geom_point(aes(colour = tissue_type, shape = tissue_type), size = 2, position = position_jitterdodge())  + ##在箱线图的基础上叠加散点图
    scale_color_manual(values = c(HD = "#0000EE", ESCC = "#EE0000")) +
    scale_shape_manual(values = c(HD = 16, ESCC = 16)) +
    coord_cartesian(ylim = c(min(ylim), max(ylim))*1.05)+  
    stat_compare_means(method = "wilcox.test") +
    labs(title = column, y = "Expression value") +
    theme(legend.position = "none") + 
    theme(panel.grid = element_blank(), axis.text=element_text(size = 12),
          axis.title = element_text(size = 12), legend.text = element_text(size = 12))
}

pdf(file="Supplementary figure S3D.pdf", width=7, height=7)
myplots_boxplot=lapply(colnames(data)[1:6], plot_column_boxplot, data = data)
grid.arrange(myplots_boxplot[[1]], myplots_boxplot[[2]], myplots_boxplot[[3]], 
             myplots_boxplot[[4]], myplots_boxplot[[5]], myplots_boxplot[[6]],  ncol = 3)
dev.off()

p.value=apply(data[, 1:6], 2, FUN = function(x){wilcox.test(x[which(data$tissue_type == "ESCC")], x[which(data$tissue_type == "HD")])$p.value})

#Supplementary figure S4
##A
library(PRROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 41, skip = 1), row.names = 1)

pdf(file = "Supplementary figure S4A.pdf", width = 5, height = 4.5)
pr=pr.curve(scores.class0 = data[data$label == "ESCC", "MLMRPscore"], 
            scores.class1 = data[data$label == "ANT", "MLMRPscore"], curve=T)
plot(pr)
dev.off()

##B
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 42, skip = 1), row.names = 1, check.names = F)

plot_column_boxplot = function(data, column) {
  ggplot(data = data, aes(x = tissue_type, y = get(column))) +
    theme_bw() +
    geom_boxplot(color="black") + 
    geom_point(aes(colour = tissue_type, shape = tissue_type), size = 2)  + 
    geom_line(aes(group = patient), size = 0.8, colour = "#9C9C9C") + 
    scale_color_manual(values = c(ANT = "#0000EE", ESCC = "#EE0000")) +
    scale_shape_manual(values = c(ANT = 16, ESCC = 16)) +
    stat_compare_means(method = "t.test", paired = T) +
    labs(title = column, y = "Expression value") +
    theme(legend.position = "none") + 
    theme(panel.grid = element_blank(), axis.text=element_text(size = 12),
          axis.title = element_text(size = 12), legend.text = element_text(size = 12))
}

pdf(file="Supplementary figure S4B.pdf", width=7, height=7)
myplots_boxplot=lapply(colnames(data)[1:6], plot_column_boxplot, data = data)
grid.arrange(myplots_boxplot[[1]], myplots_boxplot[[2]], myplots_boxplot[[3]], 
             myplots_boxplot[[4]], myplots_boxplot[[5]], myplots_boxplot[[6]],  ncol = 3)
dev.off()

p.value=apply(data[, 1:6], 2, FUN = function(x){t.test(x[which(data$tissue_type == "ESCC")], x[which(data$tissue_type == "ANT")], paired = T)$p.value})

##C
library(PRROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 43, skip = 1), row.names = 1)

pdf(file = "Supplementary figure S4C.pdf", width = 5, height = 4.5)
pr=pr.curve(scores.class0 = data[data$label == "ESCC", "MLMRPscore"], 
            scores.class1 = data[data$label == "ANT", "MLMRPscore"], curve=T)
plot(pr)
dev.off()

##D
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 44, skip = 1), row.names = 1, check.names = F)

plot_column_boxplot = function(data, column) {
  ggplot(data = data, aes(x = tissue_type, y = get(column))) +
    theme_bw() +
    geom_boxplot(color="black") + 
    geom_point(aes(colour = tissue_type, shape = tissue_type), size = 2)  + 
    geom_line(aes(group = patient), size = 0.8, colour = "#9C9C9C") + 
    scale_color_manual(values = c(ANT = "#0000EE", ESCC = "#EE0000")) +
    scale_shape_manual(values = c(ANT = 16, ESCC = 16)) +
    stat_compare_means(method = "t.test", paired = T) +
    labs(title = column, y = "Expression value") +
    theme(legend.position = "none") + 
    theme(panel.grid = element_blank(), axis.text=element_text(size = 12),
          axis.title = element_text(size = 12), legend.text = element_text(size = 12))
}

pdf(file="Supplementary figure S4D.pdf", width=7, height=7)
myplots_boxplot=lapply(colnames(data)[1:6], plot_column_boxplot, data = data)
grid.arrange(myplots_boxplot[[1]], myplots_boxplot[[2]], myplots_boxplot[[3]], 
             myplots_boxplot[[4]], myplots_boxplot[[5]], myplots_boxplot[[6]],  ncol = 3)
dev.off()

p.value=apply(data[, 1:6], 2, FUN = function(x){t.test(x[which(data$tissue_type == "ESCC")], x[which(data$tissue_type == "ANT")], paired = T)$p.value})

#Supplementary figure S5
library(ggplot2)
library(gridExtra)
library(ggpubr)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 45, skip = 1), row.names = 1)

dodge=position_dodge(width = 0.95)
plot_column_violin = function(data, column) {
  ggplot(data = data, aes(x = label, y = MLMRPscore, fill = get(column))) +
    geom_violin(position = dodge) +
    geom_boxplot(width=0.1, position = dodge) +
    scale_color_manual(values = c(Cue = "#F08961", Free = "#64BA9F", No = "#AEC2D1", Yes = "#B6A6B8", Female = "#3E8BA7", Male = "#F0C028")) +
    scale_fill_manual(values = c(Cue = "#F08961", Free = "#64BA9F", No = "#AEC2D1", Yes = "#B6A6B8", Female = "#3E8BA7", Male = "#F0C028")) +
    stat_compare_means(method = "wilcox.test") +
    labs(title = unique(data$cohort), y = "Risk Probability") +
    theme_bw() + 
    theme(legend.position="none",
          panel.grid = element_blank(), axis.text=element_text(size = 12),
          axis.title = element_text(size=12), legend.text = element_text(size=12))
}

pdf(file="Supplementary figure S5.pdf", width=18, height=9)
train_violin=lapply(colnames(data)[4:6], plot_column_violin, data = subset(data, cohort == "SCH discovery cohort"))
test_violin=lapply(colnames(data)[4:6], plot_column_violin, data = subset(data, cohort == "SCH validation cohort"))
GSE53624_violin=lapply(colnames(data)[4:6], plot_column_violin, data = subset(data, cohort == "Li cohort-1"))
GSE53622_violin=lapply(colnames(data)[4:6], plot_column_violin, data = subset(data, cohort == "Li cohort-2"))
TCGA_violin=lapply(colnames(data)[4:6], plot_column_violin, data = subset(data, cohort == "TCGA"))
GSE130078_violin=lapply(colnames(data)[4:5], plot_column_violin, data = subset(data, cohort == "You cohort"))
grid.arrange(train_violin[[1]], test_violin[[1]], GSE53624_violin[[1]], GSE53622_violin[[1]], TCGA_violin[[1]], GSE130078_violin[[1]],  
             train_violin[[2]], test_violin[[2]], GSE53624_violin[[2]], GSE53622_violin[[2]], TCGA_violin[[2]], GSE130078_violin[[2]], 
             train_violin[[3]], test_violin[[3]], GSE53624_violin[[3]], GSE53622_violin[[3]], TCGA_violin[[3]], ncol = 6)
dev.off()

#Supplementary figure S6
##A
library(PRROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 46, skip = 1), row.names = 1)

pdf(file = "Supplementary figure S6A.pdf", width = 5, height = 4.5)
pr=pr.curve(scores.class0 = data[data$Real == "I & II", "score"], 
            scores.class1 = data[data$Real == "ANT", "score"], curve=T)
plot(pr)
dev.off()

##B
library(PRROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 47, skip = 1), row.names = 1)

pdf(file = "Supplementary figure S6B.pdf", width = 5, height = 4.5)
pr=pr.curve(scores.class0 = data[data$Real == "I & II", "score"], 
            scores.class1 = data[data$Real == "ANT", "score"], curve=T)
plot(pr)
dev.off()

##C
library(PRROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 48, skip = 1), row.names = 1)

pdf(file = "Supplementary figure S6C.pdf", width = 5, height = 4.5)
pr=pr.curve(scores.class0 = data[data$Real == "I & II", "score"], 
            scores.class1 = data[data$Real == "ANT", "score"], curve=T)
plot(pr)
dev.off()

##D
library(PRROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 49, skip = 1), row.names = 1)

pdf(file = "Supplementary figure S6D.pdf", width = 5, height = 4.5)
pr=pr.curve(scores.class0 = data[data$Real == "I & II", "score"], 
            scores.class1 = data[data$Real == "ANT", "score"], curve=T)
plot(pr)
dev.off()

##E
library(PRROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 50, skip = 1), row.names = 1)

pdf(file = "Supplementary figure S6E.pdf", width = 5, height = 4.5)
pr=pr.curve(scores.class0 = data[data$Real == "I & II", "score"], 
            scores.class1 = data[data$Real == "ANT", "score"], curve=T)
plot(pr)
dev.off()

##F
library(PRROC)

data=data.frame(read_excel(path ="./source_data.xlsx", sheet = 51, skip = 1), row.names = 1)

pdf(file = "Supplementary figure S6F.pdf", width = 5, height = 4.5)
pr=pr.curve(scores.class0 = data[data$Real == "I & II", "score"], 
            scores.class1 = data[data$Real == "ANT", "score"], curve=T)
plot(pr)
dev.off()