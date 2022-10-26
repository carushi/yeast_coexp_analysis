library(ComplexHeatmap)
library(ggplot2)
library(matrixStats)

ANN_DIR="../../../ann/"

data <- read.table('yeast_count_data.tsv', sep="\t", header=T)
dds <- read.table(file.path(ANN_DIR, 'design_matrix.tsv'), sep="\t", header=T)

rownames(data) = data[,1]
data <- data[,-c(1)]
data <- apply(data, c(1,2), function(x) {log10(x+1)})
annot_df = dds[,c(2,3)]
rownames(annot_df) = dds[,1]
col = list(Stage=c("G1"="lightblue", "HU45"="blue", "HU90"="darkblue"), Genotype=c("WT"="black", "Rad53K"="#E64B35FF","Mrc1d"="#4DBBD5FF", "Rad9d"="#00A087FF"))
ha <- HeatmapAnnotation(df=annot_df, col=col, which="col")
dict <- read.table("SGD_features.tab", header=F, sep="\t", fill=TRUE, quote="", stringsAsFactor=F)
dict = dict[,c(4, 5)]
colnames(dict) = c("ID", "Name")
write.table(dict, paste0("dict_list.txt"))
for (num in c("count", "rank")) {
    if (num == "count") {
        tdata <- data[,rownames(annot_df)]
        print(dim(tdata))
        print(dim(data))
        rownames(tdata) = rownames(data)
        colnames(tdata) = rownames(annot_df)
        fig_name = "log10(x+1)"
    } else {
        tdata <- data[,rownames(annot_df)]
        tdata <- t(colRanks(tdata))
        print(dim(tdata))
        print(dim(data))
        rownames(tdata) = rownames(data)
        colnames(tdata) = rownames(annot_df)
        fig_name = "Ranking"
    }
    for (f in c("basic", "nod")) {
        a <- read.table(paste0("gene_list_", f, ".txt"), header=F)
        colnames(a) = "Name"
        a = merge(a, dict, by="Name")
        a <- subset(a, !is.na(a[,"ID"]))
        if (dim(a)[1] == 0) next
        #print(a)
        ids = unlist(sapply(a[,"ID"], function(x){if(any(x == rownames(tdata))) return(x); return(NA)}))
        ids = ids[!is.na(ids)]
        temp = tdata[ids,]
        rownames(temp) = unlist(sapply(rownames(temp), function(x){return(a[which(a[,"ID"] == x)[1],"Name"])}))
        print(head(temp))
        print(dim(temp))
        pdf(paste0("yeast_rna_seq_", num, "_", f, ".pdf"))
        draw(Heatmap(temp, name=fig_name, column_title="Samples", row_title="Genes", top_annotation=ha))
        dev.off()
        mat_scaled = t(scale(t(temp)))
        pdf(paste0("yeast_rna_seq_", num, "_", f, "_scaled.pdf"))
        draw(Heatmap(mat_scaled, name=fig_name, column_title="Samples", row_title="Genes", top_annotation=ha))
        dev.off()
    }
}






library(ComplexHeatmap)
library(ggplot2)
library(matrixStats)


data <- read.table('yeast_count_data.tsv', sep="\t", header=T)

dds <- read.table(file.path(ANN_DIR, 'design_matrix.tsv'), sep="\t", header=T)
annot_df = dds[,c(2,3)]
rownames(annot_df) = dds[,1]


rownames(data) = data[,1]
data <- data[,-c(1)]
data <- apply(data, c(1,2), function(x) {log10(x+1)})
col = list(Stage=c("G1"="lightblue", "HU45"="blue", "HU90"="darkblue"), Genotype=c("WT"="black", "Rad53K"="#E64B35FF","Mrc1d"="#4DBBD5FF", "Rad9d"="#00A087FF"))
ha <- HeatmapAnnotation(df=annot_df, col=col, which="col")
tdata <- data[,rownames(annot_df)]
pdf("yeast_rna_seq_count.pdf")
draw(Heatmap(tdata, name="log10(x+1)", column_title="Samples", row_title="Genes", top_annotation=ha))
dev.off()
tdata <- t(colRanks(tdata))
pdf("yeast_rna_seq_rank.pdf")
draw(Heatmap(tdata, name="Ranking", column_title="Samples", row_title="Genes", top_annotation=ha))
dev.off()



a <- t(data[,rownames(annot_df)])
a <- a[rowVars(a) != 0,]
a <- a[,colMeans(a) > 0]
ir.pca <- prcomp(a, center=TRUE, scale=TRUE)
sample_condition = annot_df
scores = as.data.frame(ir.pca$x)

g <- ggplot(data=scores, aes(x=PC1, y=PC2))+geom_point(aes(color=sample_condition$Stage, pch=sample_condition$Genotype), size=4)+scale_color_manual(values=c("lightblue", "blue", "darkblue"))+theme_bw()
pdf("yeast_rna_pca.pdf")
plot(g)
dev.off()

loadings = ir.pca$rotation[,c(1,2)]
loadings = loadings[order(loadings[,1], decreasing=TRUE),]
write.table(loadings, 'yeast_rna_loading.tsv', sep="\t")


ir.pca <- prcomp(t(a), center=TRUE, scale=TRUE)
sample_condition = colnames(a)
scores = as.data.frame(ir.pca$x)
g <- ggplot(data=scores, aes(x=PC1, y=PC2, label=sample_condition))+geom_text(size=2, alpha=0.5)+theme_bw()
pdf("yeast_rna_gene_pca.pdf")
plot(g)
dev.off()
scores = as.data.frame(ir.pca$x)
gene_exp = rowMeans(t(a))
require(viridis)
g <- ggplot(data=scores, aes(x=PC1, y=PC2))+geom_point(aes(color=gene_exp), size=2, alpha=0.5)+scale_color_viridis()+theme_bw()
pdf("yeast_rna_gene_pca_exp.pdf")
plot(g)
dev.off()
loadings = ir.pca$rotation[,c(1,2)]
loadings = loadings[order(loadings[,1], decreasing=TRUE),]
write.table(loadings, 'yeast_rna_gene_loading.tsv', sep="\t")

