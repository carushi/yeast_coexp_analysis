plot.KEGG.pathway <- function(map, gene.data, outfile, kegg.code, all.gene.data=NULL) {
    require(pathview)
    require(org.Sc.sgd.db)
    cpd = NULL
    kegg = T
    layer = F
    tail = "png"
    tryCatch({
        pv.out <- pathview(gene.data=gene.data, gene.idtype="ORF", pathway.id=map, species=kegg.code, out.suffix=map, keys.align="y", kegg.native=kegg, match.data=T, same.layer=layer)
        plot.name <- paste0(outfile, "_", map, ".", tail)
        system(paste0("mv ", kegg.code, map, ".", map, ".", tail, " ", plot.name))
    }, error=function(cond) {
        return(NA)
    })
    # print('???')
    if (is.null(all.gene.data)) return()
    tryCatch({
        pv.out <- pathview(gene.data=all.gene.data, gene.idtype="ORF", pathway.id=map, species=kegg.code, out.suffix=map, keys.align="y", kegg.native=kegg, match.data=T, same.layer=layer)
        plot.name <- paste0(outfile, "_", map, "_all.", tail)
        system(paste0("mv ", kegg.code, map, ".", map, ".", tail, " ", plot.name))
    }, error=function(cond) {
        return(NA)
    })
}

orf2common <- function() {
    require(org.Sc.sgd.db)
    orf2gene = org.Sc.sgdENTREZID
    mapped_gene <- mappedkeys(orf2gene)
    xx <- as.list(orf2gene[mapped_gene])
    return(xx)
}
print.KEGG.enrichment <- function(gene.data, all.genes, outfile, kegg.code) {
    require(clusterProfiler)
    dict = orf2common()
    all.genes.ez = unique(unlist(sapply(all.genes, function(x){return(dict[x])})))
    for (change in c("up", "down", "change")) {
        if (change == "up")
            part = gene.data[gene.data > 0]
        else if (change == "down") 
            part = gene.data[gene.data < 0]
        else
            part = gene.data
        if (length(part) == 0)  next
        part.genes.ez = unique(unlist(sapply(names(part), function(x){return(dict[x])})))
        ekegg <- enrichKEGG(unique(names(part)), organism=kegg.code, keyType="kegg", pvalueCutoff=0.05, pAdjustMethod="BH", unique(all.genes), minGSSize=3, maxGSSize=500, qvalueCutoff=0.2, use_internal_data=FALSE)
        print(ekegg)
        fname = paste0(outfile, "_enrich_", change, ".txt")
        write.table(data.frame(ekegg), file=fname)
        if (dim(data.frame(ekegg))[1] == 0) next
        png(paste0(outfile, "_enrich_", change, ".png"), width=700, height=700)
        plot(dotplot(ekegg, showCategory=20))
        dev.off()
        png(paste0(outfile, "_gene_enrich_", change, ".png"), width=500, height=500)
        names(part) = part.genes.ez
        plot(cnetplot(ekegg, foldChange=part))
        dev.off()
    }
}

plot.KEGG.enrichment <- function(fc, pvalue, thres, gene_name, outfile, kegg_id_file) {
    enriched <- gene_name[!is.na(pvalue) & (pvalue < thres)]
    fcs <- fc[!is.na(pvalue) & (pvalue < thres)]
    names(fcs) = enriched
    kegg.code = "sce"
    pathways  = read.table(kegg_id_file, header=F, sep=" ", colClasses=c("character"))[,1]
    df = data.frame(FC=fcs, row.names=enriched)
    adf = data.frame(FC=fc, row.names=gene_name)
    for (i in 1:length(pathways)) {
        plot.KEGG.pathway(pathways[i], df, outfile, kegg.code, adf)
    }
    print.KEGG.enrichment(fcs, gene_name, outfile, kegg.code)
}

plot.GO.enrichment <- function(fc, pvalue, thres, gene_name, outfile) {
    require(clusterProfiler)
    require(org.Sc.sgd.db)
    enriched <- gene_name[!is.na(pvalue) & (pvalue < thres)]
    gene.data <- fc[!is.na(pvalue) & (pvalue < thres)]
    names(gene.data) = enriched
    dict = orf2common()
    all.genes.ez = unique(unlist(sapply(gene_name, function(x){return(dict[x])})))
    for (change in c("up", "down", "change")) {
        if (change == "up")
            part = gene.data[gene.data > 0]
        else if (change == "down") 
            part = gene.data[gene.data < 0]
        else
            part = gene.data
        #print(part[1:10])
        if (length(part) == 0)  next
        part.genes.ez = unique(unlist(sapply(names(part), function(x){return(dict[x])})))
        ego <- enrichGO(gene          = part.genes.ez,
                        universe      = all.genes.ez,
                        OrgDb         = org.Sc.sgd.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)    
        print(ego)
        fname = paste0(outfile, "_enrich_", change, ".txt")
        write.table(data.frame(ego), file=fname)
        if (dim(data.frame(ego))[1] == 0) next
        png(paste0(outfile, "_enrich_", change, ".png"), width=700, height=700)
        plot(dotplot(ego, showCategory=20))
        dev.off()
        png(paste0(outfile, "_gene_enrich_", change, ".png"), width=500, height=500)
        names(part) = part.genes.ez
        plot(cnetplot(ego, foldChange=part))
        dev.off()
    }

}
