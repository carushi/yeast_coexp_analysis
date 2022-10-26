library(ChIPpeakAnno)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(dplyr)
require(motifStack)
require(BSgenome)
require(BSgenome.Scerevisiae.UCSC.sacCer3)
require(org.Sc.sgd.db)
ANN_DIR <<- "../../../ann/"
CHIP_ANN <<- "../../ann/"

txdb <<- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
txdb_genes <<- genes(txdb)
bsgenome <<- getBSgenome("BSgenome.Scerevisiae.UCSC.sacCer3")
origin <<- makeGRangesFromDataFrame(read.table(file.path(CHIP_ANN, "all_origin_sc_roman.tsv"), header=T, sep="\t", quote="\""), keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field=c("chr"))

obtain.freqs.bg <- function() {
    if (file.exists("freqs.rds")) {
        freqs <- readRDS("freqs.rds")
    } else {
        chrms <- names(bsgenome)
        chrms <- chrms[chrms != "chrM"]
        whole_genome <<- NULL
        for (chr in chrms) {
            if (!any(chr == c("chrIV"))) next
            if (is.null(whole_genome)) whole_genome <<- unlist(bsgenome[[chr]])
            else whole_genome <<- c(whole_genome, unlist(bsgenome[[chr]]))
        }
        freqs <- oligoFrequency(whole_genome, MarkovOrder=3)
        saveRDS(freqs, "freqs.rds")
    }
    return(freqs)
}

convert.origin.data <- function(a) {
    data <- a[,c('chr', 'start', 'end', 'name', 'score')]
    gr1 <- toGRanges(data, format="BED", header=FALSE)
    return(gr1)
}

read.idr.output <- function(bed, threshold=NULL) {
    a <- read.table(bed, header=F)
    if (dim(a)[2] == 20) {
        colnames(a) <- c("chr", "ostart", "oend", "name", "score", "strand", "signalValue", "pvalue", "qvalue", "summit", "localIDR", "globalIDR", "rep1_start", "rep1_end", "rep1_signalValue", "rep1_summit", "rep2_start", "rep2_end", "rep2_signalValue", "rep2_summit")
    } else if (dim(a)[2] == 17) {
        colnames(a) <- c("chr", "ostart", "oend", "name", "score", "strand", "signalValue", "pvalue", "qvalue", "localIDR", "globalIDR", "rep1_start", "rep1_end", "rep1_signalValue",  "rep2_start", "rep2_end", "rep2_signalValue")
    } else if (dim(a)[2] == 3) {
        colnames(a) <- c("chr", "ostart", "oend")
        a[,"name"] = 1:dim(a)[1]
        a[,"score"] = 1
    } else {
        colnames(a) <- c("chr", "ostart", "oend", "name", "score")
        a <- a[,1:5]
    }
    a[,"start"] = a[,"ostart"]
    a[,"end"] = a[,"oend"]
    if (!is.null(threshold) && any(colnames(a) == 'globalIDR')) {
        a <- subset(a, a$globalIDR >= threshold)
    }
    a <- a[!duplicated(a[,c("start", "end")]),]
    gr1 <- toGRanges(a, format="BED", header=FALSE) 
    return(gr1)
}

orf2common <- function() {
    orf2gene = org.Sc.sgdENTREZID
    mapped_gene <- mappedkeys(orf2gene)
    xx <- as.list(orf2gene[mapped_gene])
    return(xx)
}

plot.enrich.terms <- function(header, gene_list) {
    dict = orf2common()
    part.genes.ez = unique(unlist(sapply(names(gene_list), function(x){return(dict[x])})))
    print(head(gene_list))
    print(head(part.genes.ez))
    ego <- enrichGO(gene          = part.genes.ez,
                    OrgDb         = org.Sc.sgd.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)    
    fname = paste0(header, "_enrich_go.txt")
    write.table(data.frame(ego), file=fname)
    if (dim(data.frame(ego))[1] > 0) {
        png(paste0(header, "_genrich.png"), width=700, height=700)
        plot(dotplot(ego, showCategory=20))
        dev.off()
    }
    kegg.code = "sce"
    ekegg <- enrichKEGG(names(gene_list), organism=kegg.code, keyType="kegg", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, use_internal_data=FALSE)
    fname = paste0(header, "_enrich_kegg.txt")
    write.table(data.frame(ekegg), file=fname)
    if (dim(data.frame(ekegg))[1] > 0) {
        png(paste0(header, "_kenrich.png"), width=700, height=700)
        plot(dotplot(ekegg, showCategory=20))
        dev.off()
    }
    if (dim(data.frame(ekegg))[1] > 0) {
        png(paste0(header, "_kgene_enrich.png"), width=500, height=500)
        plot(cnetplot(ekegg, foldChange=gene_list))
        dev.off()
    }
    names(gene_list) = part.genes.ez
    if (dim(data.frame(ego))[1] > 0) {
        png(paste0(header, "_ggene_enrich.png"), width=500, height=500)
        plot(cnetplot(ego, foldChange=gene_list))
        dev.off()
    }
}

count.profile.around.tss <- function(peaks, header, sp, stop_plot=FALSE, MERGED=FALSE) {
    print(peaks)
    if (length(peaks) <= 2) 
        return(NULL)
    pdf(paste0(header, "_tss.pdf"))
    p <- binOverFeature(peaks, annotationData=txdb_genes, featureSite=c("FeatureStart"),
            radius=5000, nbins=20, FUN=length, errFun=0,
            ylab="count", main="Distribution of aggregated peak numbers around TSS")
    dev.off()
    p <- data.frame(x=as.numeric(rownames(p)), y=as.numeric(p[,1]), sample=sp, type=rep('TSS', length(rownames(p))))
    pdf(paste0(header, "_tts.pdf"))
    s <- binOverFeature(peaks, annotationData=txdb_genes, featureSite=c("FeatureEnd"),
            radius=5000, nbins=20, FUN=length, errFun=0,
            ylab="count", main="Distribution of aggregated peak numbers around TTS")
    dev.off()
    s <- data.frame(x=as.numeric(rownames(s)), y=as.numeric(s[,1]), sample=sp, type=rep('TTS', length(rownames(s))))
    o <- binOverFeature(peaks, annotationData=origin,
            radius=5000, nbins=20, FUN=length, errFun=0, select='nearest',
            ylab="count", main="Distribution of aggregated peak numbers around origin")
    if (!is.null(o)) {
        o <- data.frame(x=as.numeric(rownames(o)), y=as.numeric(o[,1]), sample=sp, type=rep('origin', length(rownames(o))))
    }
    dev.off()
    result = list(tss=p, tts=s, origin=o)
    write.table(rbind(rbind(p, s), o), file=paste0(header, '_binoverfeature_mat.tsv'), sep="\t")
    if (stop_plot) return(result)
    require(ggplot2)
    require(ggsci)
    aCR<-assignChromosomeRegion(peaks, nucleotideLevel=FALSE, proximal.promoter.cutoff = c(upstream=500, downstream=50),
                                    immediate.downstream.cutoff = c(upstream=50, downstream=500),
                                    TxDb=txdb)
    write.table(aCR, file=paste0(header, '_annot.csv'), sep=",")
    pdf(paste0(header, "_annot.pdf"))
    barplot(aCR$percentage, las=3)+
    dev.off()
    tdf = data.frame(aCR)
    print(head(tdf))
    tdf = cbind(region=factor(rownames(tdf), levels=c('Promoters', 'Exons', 'Introns', 'immediateDownstream', 'Intergenic.Region')), tdf)
    pdf(paste0(header, '_annot_bar.pdf'))
    g <- ggplot(tdf, aes(x=region, y=percentage, fill=region))+geom_bar(stat='identity')+scale_fill_npg()+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    plot(g)
    dev.off()
    pdf(paste0(header, '_annot_pie.pdf'))
    g <- ggplot(tdf, aes(x="", y=percentage, fill=region))+geom_bar(stat='identity')+scale_fill_npg()+theme_bw()
    g <- g+coord_polar("y", start=0)
    plot(g)
    dev.off()
    overlaps.anno <- annotatePeakInBatch(peaks, 
                                        AnnotationData=txdb_genes, 
                                        output="nearestBiDirectionalPromoters",
                                        bindingRegion=c(-2000, 500))
    df <- data.frame(overlaps.anno)
    colnames(df)[colnames(df) == "distance"] = "dist"
    df <- group_by(df, peak) %>% dplyr::slice(which.min(dist))
    if (MERGED)
        df = df[,c("peak", "signalValue", "gene_id")]
    else
        df = df[,c("peak", "score", "gene_id")]
    vec = unlist(df[,2])
    names(vec) = unlist(df[,3])
    write.csv(data.frame(unname(overlaps.anno)), paste0(header, ".csv"))
    return(result)
}


motif.detection <- function(header, peaks, MERGED=FALSE) {
    seq <- getAllPeakSequence(peaks, upstream=20, downstream=20, genome=bsgenome)
    if (MERGED) {
        seq <- seq[order(mcols(peaks)[,"signalValue"], decreasing=TRUE)]
    } else {
        seq <- seq[order(mcols(peaks)[,"score"], decreasing=TRUE)]
    }
    write2FASTA(seq, paste0(header, ".fa"))
    if (length(seq) <= 20) return()
    oligoLength = 7
    freqs <- obtain.freqs.bg()
    tryCatch({
        os <- oligoSummary(seq, oligoLength=oligoLength, MarkovOrder=3, 
                        quickMotif=TRUE, freq=freqs)
        print(os)
        ## plot the results
        zscore <- sort(os$zscore)
        png(paste0(header, "_motif_freq.png"))
        h <- hist(zscore, breaks=100, xlim=c(-50, 50), main="Histogram of Z-score")
        text(zscore[length(zscore)], max(h$counts)/10, 
            labels=names(zscore[length(zscore)]), adj=1)
        dev.off()
        saveRDS(os, "test_os.rds")
        pfms <- mapply(function(.ele, id)
            new("pfm", mat=.ele, name=paste("motif for ", header)), 
            os$motifs, 1:length(os$motifs))
        png(paste0(header, "_motif_logo.png"))
        motifStack(pfms[[1]])    
        dev.off()
    })
    set.seed(0)
    require(BCRANK)
    require(seqLogo)
    BCRANKout <- bcrank(paste0(header, ".fa"), restarts = 25, 
        use.P1 = TRUE, use.P2 = TRUE)
    toptable(BCRANKout)
    topMotif <- toptable(BCRANKout, 1)
    weightMatrix <- pwm(topMotif, normalize = FALSE)
    weightMatrixNormalized <- pwm(topMotif, normalize = TRUE)
    png(paste0(header, "_seqlogo_weight.png"))
    seqLogo(weightMatrixNormalized)
    dev.off()
    if (!grepl('narrow', header, fixed=TRUE) && !grepl('summits', header, fixed=TRUE)) return()
    seq$sequence = unlist(sapply(seq$sequence, function(x){substring(x, nchar(x)/2-20, nchar(x)/2+20)}))
    letter_counts <- consensusMatrix(seq$sequence)
    probs <- prop.table(letter_counts[1:4,], 2)
    png(paste0(header, "_seqlogo.png"))
    seqLogo(probs, ic.scale=FALSE)
    dev.off()
}

plot.tss_bin.data <- function(stack, tss_bin, header) {
    require(ggplot2)
    require(ggsci)
    pdf(paste0(header, "bin.pdf"))
    tss_bin$sample = factor(tss_bin$sample)
    g <- ggplot(tss_bin, aes(x=x, y=y, color=sample))+geom_line()+theme_bw()    
    if (grepl('time', header, fixed=TRUE)) {
        g <- g+scale_color_manual(values=c("lightblue", "blue", "darkblue"))
    } else {
        g <- g+scale_color_manual(values=c("black","#E64B35FF","#4DBBD5FF"))
    }
    plot(g)
    dev.off()
}

plot.basic.annotation.for.merged.peaks <- function(dir, bg_bed) {
    gr1 <- read.idr.output(paste0(dir, bg_bed))
    print(head(gr1))
    result <- count.profile.around.tss(gr1, bg_bed, 0, FALSE)
    write.table(result, file='output_annotation.csv', sep=",", quote=F)
}
plot.basic.annotation.for.origin <- function() {
    a <- read.table(file.path(CHIP_ANN, "origin_with_dna_amount_ann.tsv"), header=T, sep="\t")
    a[,"name"] = 1:dim(a)[1]
    a[,"score"] = 1
    for (target in c('all', 'valid', 'early', 'late')) {
        if (target == 'valid')
            temp <- a[a[,"origin"],]
        else if (target == 'early')
            temp <- a[a[,"origin"] & !a[,"late_origin"],]
        else if (target == 'late')
            temp <- a[a[,"late_origin"],]
        else
            temp <- a[,]
        gr1 <- convert.origin.data(temp)
        header = paste0(target, '_origin_5000')
        print(head(gr1))
        result <- count.profile.around.tss(gr1, header, 0, FALSE)
        write.table(result, file=paste0(header, '_annotation.csv'), sep=",", quote=F)
    }
}

summarize.called.peaks <- function() {
    require(reshape2)
    require(ggsci)
    require(ggplot2)
    dir <- "./"
    for (peak in c('narrow', 'broad')[1]) {
        for (wt in c('', '1', '2', '3')) {
            a <- read.table(paste0("merged_peaks_", peak, "_sorted.bed_annot.csv"), sep=",")
            a <- a[,1, drop=FALSE]
            colnames(a)[1] = 'all'
            print(head(a))
            for (i in 0:8) {
                bed <- paste0(dir, 'peak_pe_', i, '_', peak, 'Peak_annot.csv')
                b <- read.table(bed, sep=" ", header=T, stringsAsFactors=FALSE)
                a <- cbind(a, b[,1])
                colnames(a)[dim(a)[2]] = i
            }
            print(head(a))
            a <- cbind(region=rownames(a), a)
            print(head(a))
            data = melt(a)
            print(head(data))
            png(paste0(peak, '_annotation_ratio.png'))
            g <- ggplot(data, aes(fill=region, y=value, x=variable)) + 
            geom_bar(position="stack", stat="identity")+scale_fill_npg()
            plot(g)
            dev.off()
        }
    }
}
summarize.called.origin <- function() {
    require(reshape2)
    require(ggsci)
    require(ggplot2)
    target = 'all'
    dir = "annotation_script/"
    a <- read.table(paste0(dir, target, "_origin_5000_annot.csv"), sep=",", header=T)
    a <- a[,1, drop=FALSE]
    colnames(a)[1] = 'all'
    print(head(a))
    for (i in c('valid', 'early', 'late')) {
        bed <- paste0(dir, i, '_origin_5000_annot.csv')
        b <- read.table(bed, sep=",", header=T, stringsAsFactors=FALSE)
        a <- cbind(a, b[,1])
        colnames(a)[dim(a)[2]] = i
    }
    a <- cbind(region=rownames(a), a)
    print(head(a))
    data = melt(a)
    print(head(data))
    png(paste0('origin', '_annotation_ratio.png'))
    g <- ggplot(data, aes(fill=region, y=value, x=variable)) + 
    geom_bar(position="stack", stat="identity")+scale_fill_npg()
    plot(g)
    dev.off()
}

for (i in c('', '1', '2', '3')) {
    plot.basic.annotation.for.merged.peaks(paste0('all_peaks_pe_narrow_wt', i, '_merged.bed'))
}
plot.basic.annotation.for.origin()
summarize.called.origin()
summarize.called.peaks()
