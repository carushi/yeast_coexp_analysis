extract.cluster <- function(mat, rcl.list, r.dend) {
    mnames = rownames(mat)
    for (i in 1:length(rcl.list)) {
        clu <- mnames[rcl.list[[as.character(i)]]]
        if (i == 1) {
            out <- cbind(clu, paste("cluster", i, sep=""))
            colnames(out) <- c("GeneID", "Cluster")
        } else {
            clu <- cbind(clu, paste("cluster", i, sep=""))
            out <- rbind(out, clu)
        }
    }
    return(out)
}

filter.mat <- function(mat, colname, rowname) {
    colnames(mat) = colname
    rownames(mat) = rowname
    stopifnot(!any(is.na(rownames(mat))))
    stopifnot(!any(is.na(colnames(mat))))
    stopifnot(all(!duplicated(rownames(mat))))
    mat[is.na(mat)] = 0
    return(mat)
}

clustering.by.complex.heatmap <- function(mat, network, header, bg_gene, tcul_num, ha)
{
    heatmap2 = Heatmap(mat, clustering_distance_rows = "pearson", clustering_distance_columns="pearson", km=tclu_num, right_annotation = ha)
    row_ord = row_order(heatmap2)
    row_dend = row_dend(heatmap2)
    out = extract.cluster(mat, row_ord, row_dend)
    out = out[unlist(sapply(rownames(mat), function(x){return(which(x == out[,1]))})),]
    split = gsub(pattern="^cluster", "", out[,2])
    names(split) = out[,1]
    png(paste0('df_heatmap_', network, '_', header, '.png'))
    draw(heatmap2, row_split=split, column_split=split)
    dev.off()
    out = rbind(out, data.frame(GeneID=unlist(sapply(bg_gene, function(x){if(!any(x == out[,1]))return(x);})), Cluster='insignificant'))
    write.table(out, paste0(network, '_', header, '_cluster_', tclu_num, '.tsv'))
    ha = rowAnnotation(gene=col_list, col = list(gene = col_fun))
    cha = HeatmapAnnotation(gene=col_list, col = list(gene = col_fun))
    heatmap3 = Heatmap(mat2,  row_split=split, column_split=split, cluster_row_slices=FALSE, cluster_column_slices=FALSE, clustering_distance_rows="pearson", clustering_distance_columns="pearson", right_annotation = ha)
    png(paste0('df_heatmap_after_', network, '_', header, '.png'))
    draw(heatmap3)
    dev.off()
}

clustering.by.hierarchical.clustering <- function(mat, network, header, bg_gene)
{
    require(stringr)
    require(dynamicTreeCut)
    m_dist <- dist(mat, diag = FALSE)
    m_hclust <- hclust(m_dist, method= "average")
    for (filter in c('', '_dc')) {
        if (filter == '') {
            groups <- cutree(m_hclust, k=10)
            cluster = 10
        } else {
            groups = cutreeDynamic(m_hclust, distM=as.matrix(m_dist))
            cluster = length(unique(groups))
        }
        stopifnot(length(groups) == dim(mat)[1])
        heatmap2 = Heatmap(mat, cluster_rows=m_hclust, clustering_distance_columns="pearson", right_annotation = ha)
        pdf(paste0('df_heatmap_', network, '_', header, filter, '.pdf'))
        draw(heatmap2)
        dev.off()
        out = data.frame(GeneID=rownames(mat), Cluster=paste0('cluster', groups))
        ha = rowAnnotation(gene=col_list, col = list(gene = col_fun))
        cha = HeatmapAnnotation(gene=col_list, col = list(gene = col_fun))
        split = str_pad(gsub(pattern="^cluster", "", out[,2]), 2, pad="0")
        names(split) = out[,1]
        row_labels = rownames(mat)
        if (dim(mat)[1] > 50) row_labels = rep('', dim(mat)[1])
        col_labels = colnames(mat)
        if (dim(mat)[2] > 50) col_labels = rep('', dim(mat)[2])
        col_cluster = rainbow(length(unique(split)))
        names(col_cluster) = sort(unique(split))
        hac = rowAnnotation(gene=col_list, cluster=split, col = list(gene = col_fun, cluster=col_cluster))
        heatmap2 = Heatmap(mat, cluster_rows=m_hclust, clustering_distance_columns="pearson", right_annotation = hac, row_labels=row_labels, column_labels=col_labels)
        pdf(paste0('df_heatmap_', network, '_', header, '_split', filter, '.pdf'))
        draw(heatmap2)
        dev.off()
        out = rbind(out, data.frame(GeneID=unlist(sapply(bg_gene, function(x){if(!any(x == out[,1]))return(x);})), Cluster='insignificant'))
        write.table(out, paste0(network, '_', header, '_cluster_', tclu_num, filter, '.tsv'))
        heatmap3 = Heatmap(mat2, cluster_rows=m_hclust, cluster_columns=m_hclust, right_annotation = hac, row_labels=row_labels, column_labels=row_labels)
        pdf(paste0('df_heatmap_after_', network, '_', header, filter, '.pdf'))
        draw(heatmap3)
        dev.off()
    }
}



library(rhdf5)
library(ComplexHeatmap)
library(circlize)
library(viridis)

clu_num = 10

# File names of DEG list 
dir = './'
files = list.files(path=dir, recursive = FALSE, pattern = "*.tsv", full.names = TRUE)
# P-value threshold for DEG
thres = 0.05
for (network in c('meta', 'prio')[1]) {
    fname = paste0("yeast_", network, "AggNet.hdf5")
    hdf5 = H5Fopen(fname)
    gene_list = hdf5$col
    coexp_mat = hdf5$agg
    for (deg_file in files) {
        print(deg_file)
        deg <- read.table(deg_file, header=TRUE, sep="\t")
        deg <- deg[!is.na(rownames(deg)),]
        deg <- deg[!is.na(deg[,'padj']),]
        deg <- deg[!grepl(rownames(deg), pattern='^NA*'),]
        header = gsub(pattern = "\\.tsv$", "", basename(deg_file))
        match = unlist(sapply(rownames(deg), function(x) {
            return(any(x == gene_list))
        }))
        deg <- deg[unlist(sapply(rownames(deg), function(x){return(any(x == gene_list))})),]
        bg_gene = rownames(deg)
        deg <- deg[deg[,'padj'] <= thres,]
        print(rownames(deg))
        print(c("->", dim(deg)))
        if (dim(deg)[1] < 10) next
        gene_order = unlist(sapply(rownames(deg), function(x){return(which(x == gene_list)[1])}))
        mat <- filter.mat(coexp_mat[gene_order, ], gene_list, rownames(deg))
        mat2 <- filter.mat(coexp_mat[gene_order, gene_order], rownames(deg), rownames(deg))
        stopifnot(isSymmetric(mat2))
        mat = mat[,colMeans(mat) != unlist(apply(mat, 2, max, na.rm=TRUE))]
        flag = (rowMeans(mat) == unlist(apply(mat, 1, max, na.rm=TRUE)))
        mat = mat[!flag,]
        gene_order = gene_order[!flag]
        tclu_num = min(dim(mat)[1], clu_num)
        col_fun = colorRamp2(c(-2, 0, 2), c("darkblue", "white", "darkred"))
        col_list = unlist(sapply(rownames(mat), function(x){return(deg[which(rownames(deg) == x)[1],'log2FoldChange'])}))
        names(col_list) = rownames(mat)
        ha = rowAnnotation(gene=col_list, col = list(gene = col_fun))
        # clustering.by.complex.heatmap(mat, network, header, bg_gene) built-in clustering in ComplexHeatmap
        clustering.by.hierarchical.clustering(mat, network, header, bg_gene)
    }
    H5Fclose(hdf5)
}
