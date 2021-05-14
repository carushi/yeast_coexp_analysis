source("/home/rkawaguc/ipython/ChIP/rna_seq_script/go_enrichment.R")
# coexpression cluster data
dir = "./"
# DEG data
rna_dir = "./"
for (fname in list.files(rna_dir, pattern="*.tsv", full.names=TRUE)) {
    header = gsub(pattern = "\\.tsv$", "", basename(fname))
    for (filter in c('', '_dc')[2]) {
        for (network in c('prio', 'meta')[2]) {
            data <- read.table(fname, header=T, sep="\t", stringsAsFactor=F)
            coexp_cluster = paste0(dir, network, '_', header, '_cluster_', 10, filter, '.tsv')
            if (file.exists(coexp_cluster)) {
                cluster = read.table(coexp_cluster, stringsAsFactor=F)
            } else {
                print(c('No enough DEGs', coexp_cluster))
                next
            }
            data <- subset(data, !is.na(rownames(data)))
            base_prefix = unlist(strsplit(basename(fname), ".", fixed=TRUE))[1]
            data <- data[unlist(sapply(rownames(data), function(x){return(any(x == cluster[,1]))})),]
            for (cindex in 1:10) {
                pvalue = rep(1, dim(data)[1])
                tcluster = cluster[cluster[,2] == paste0('cluster', cindex),]
                if (dim(tcluster)[1] == 0) 
                    break
                tindex = unlist(sapply(tcluster[,1], function(x){return(which(rownames(data) == x))}))
                pvalue[tindex] = data[tindex, ][,'padj']
                pvalue[is.na(data[,'padj'])] = NA
                names(pvalue) = rownames(data)
                if (length(which(pvalue <= 0.05)) == 0) {
                    exit()
                }
                if (dim(tcluster)[1] <= 3) next
                plot.KEGG.enrichment(data[,2], pvalue, 0.05, rownames(data), paste0("kegg_", network, '_', header, '_', cindex, filter), "/home/rkawaguc/ipython/ChIP/id_list/kegg_id_list.txt")
                plot.GO.enrichment(data[,2], pvalue, 0.05, rownames(data), paste0("go_", network, '_', header, '_', cindex, filter))
            }
        }
    }
}