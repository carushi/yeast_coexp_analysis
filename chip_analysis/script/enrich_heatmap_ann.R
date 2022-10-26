require(GenomicRanges)
require(GenomicFeatures)
require(EnrichedHeatmap)

LOCAL <<- FALSE
HYPER_CHIP <<- c('wohc', 'ohc', 'pho', '')[1]
# without-hyperchipable, overlapped-with-hyperchipable, pho84, or all TSSs

ANN_DIR <<- '../../../ann/'
CHIP_ANN <<- "../../ann/"
CLUSTER_FILE <<- file.path(CHIP_ANN, 'cluster', 'cluster_file.txt')
ORI_TIME_FILE <<- file.path(CHIP_ANN, 'origin', 'ori_db_time.tsv')
if (HYPER_CHIP == 'wohc') {
    TSS_FILE <<- file.path(CHIP_ANN, 'gene', "unoverlapped_gene.gtf")
} else if (HYPER_CHIP == 'ohc') {
    TSS_FILE <<- file.path(CHIP_ANN, "gene", "overlapped_gene.gtf")
} else if (HYPER_CHIP == 'pho') {
    TSS_FILE <<- file.path(CHIP_ANN, "gene", "pho84_gene.gtf")
} else {
    TSS_FILE <<- file.path(ANN_DIR, "Saccharomyces_cerevisiae.R64-1-1.98_chr_sorted.gtf")
}
PEAK_FILE <<- file.path(CHIP_ANN, "peak", "all_peaks_pe_narrow_wt_merged.bed")


DEFAULT <<- "_"
VALID_STAT <<- TRUE
ORI_STAT <<- FALSE
OPTIONAL <<- TRUE

relative <<- TRUE
global <<- FALSE
statistics <<- FALSE

set.genotype.stage <- function(dataset) {
    if (dataset == 'TF') {
        GENOT <<- c('WT', 'ixr1d', 'swi4d', 'swi6d')
        STAGE <<- c('G1', 'HU45', 'HU90')
    } else if (dataset == 'null') {
        GENOT <<- c('WT', 'sml1d', 'mec1dsml1d', 'rad53dsml1d')
        STAGE <<- c('G1') 
    } else {
        GENOT <<- c('WT', 'Rad53K', 'Mrc1d')
        STAGE <<- c('G1', 'HU45', 'HU90')
    }
}

write.gr.to.bed <- function(gr, output, rangeOutput, upstream, downstream) {
    df <- data.frame(seqnames=seqnames(gr),
        starts=start(gr)-1,
        ends=end(gr),
        names=c(rep(".", length(gr))),
        scores=c(rep(".", length(gr))),
        strands=strand(gr),
        gene_id=names(gr))
    for (i in 1:dim(df)[2]) {
        if (df[i,'strands'] == '+') {
            df[i, 'starts'] = df[i, 'starts']-upstream
            df[i, 'ends'] = df[i, 'ends']+downstream
        } else {
            df[i, 'starts'] = df[i, 'starts']-downstream
            df[i, 'ends'] = df[i, 'ends']+upstream
        }
    }
    write.table(df, file=output, quote=F, sep="\t", row.names=F, col.names=F)
    saveRDS(gr, file=rangeOutput)
}

write.tss.data <- function(target) {
    output_up = paste0(target, '_genomic_range_750_up.bed')
    rangeOutput_up = paste0(target, '_genomic_range_750_up.rds')
    output_down = paste0(target, '_genomic_range_750_down.bed')
    rangeOutput_down = paste0(target, '_genomic_range_750_down.rds')
    fname = TSS_FILE
    dist = 0
    txdb <- makeTxDbFromGFF(fname, format="auto")
    if (grepl('tss', target) ) {
        promoter <- promoters(genes(txdb), upstream = dist, downstream = dist)
        write.gr.to.bed(promoter, output_up, rangeOutput_up, 750, 0)
        write.gr.to.bed(promoter, output_down, rangeOutput_down, 0, 750)
    } else {
        tts <- transcripts(txdb)
        tts <- resize(tts, dist, fix="end")
        tts <- resize(tts, dist+dist+1, fix="start")
        names <- as.vector(unlist(sapply(tts$tx_name, function(x){return(strsplit(x, '_', fixed=TRUE)[[1]][1])})))
        tts$gene_id <- names
        names(tts) <- tts$gene_id
        write.gr.to.bed(tts, output_up, rangeOutput_up, 750, 0)
        write.gr.to.bed(tts, output_down, rangeOutput_down, 0, 750)
    }
}

read.tss.data <- function(target) {
    fname = TSS_FILE
    dist = 0
    txdb <- makeTxDbFromGFF(fname, format="auto")
    if (grepl('tss', target) ) {
        promoter <- promoters(genes(txdb), upstream = dist, downstream = dist)
        if (target == 'tss')
            return(promoter)
        plus <- promoter[strand(promoter) == "+"]
        minus <- promoter[strand(promoter) == "-"]
        if (target == 'tss_plus')
            return(plus)
        else
            return(minus)
    } else {
        tts <- transcripts(txdb)
        tts <- GenomicRanges::resize(tts, dist, fix="end")
        tts <- GenomicRanges::resize(tts, dist+dist+1, fix="start")
        names <- as.vector(unlist(sapply(as.character(tts$tx_name), function(x){return(strsplit(x, '_', fixed=TRUE)[[1]][1])})))
        tts$gene_id <- names
        names(tts) <- names
        return(tts)
    }
}
read.origin.data <- function(target, filt=FALSE, time=TRUE, ori.id=FALSE) {
    fname = file.path(CHIP_ANN, "origin", paste0(target, '_center_wo_hc_time.bed'))
    a <- read.table(fname, header=T, sep="\t", quote="")
    a[,1] = unlist(sapply(a[,1], function(x){if (x == 'chrMito') return('chrM'); return(x)}))
    if (ori.id) {
        a <- a[a[,5] != '',]
    }
    return(makeGRangesFromDataFrame(a, ignore.strand=TRUE, keep.extra.columns=TRUE))
}

plot.cluster.annotated.origin.heatmap <- function(fc_range, target, target_range, tail, ext, btail, col_fun, top_ann, filt_prefix, time_prefix) {
    axis_name = c('-15000', '0', '', '+15000')
    target_range <- target_range[!is.na(target_range$time),]
    time_data <- target_range$time
    row_order = order(time_data)
    mat = list()
    ag <- NULL
    for (genot in GENOT) {
        g <- NULL
        for (stage in STAGE) {
            key = paste0(genot, '_', stage)
            header = paste0(target, '_', key, tail, btail, filt_prefix, time_prefix)
            mat[[key]] = normalizeToMatrix(fc_range, target_range, value_column = key, mean_mode = "w0", w = 25, extend=ext)
            mat[[key]][mat[[key]] < 0] = 0
            mat[[key]][mat[[key]] > 2] = 2
            tg <- EnrichedHeatmap(mat[[key]], name=key, col=col_fun, top_annotation=top_ann, axis_name=axis_name)
            row_order_list = row_order(draw(tg+Heatmap(target_range$time, col = c("white", "orange"), name = "Replication time", show_row_names=FALSE, width=unit(5, 'mm')), row_order=row_order))
            write.table(as.data.frame(target_range)[unlist(row_order_list),], file=paste0('heatmap_', header, '_row_list.tsv'), sep='\t')
            if (is.null(g)) 
                g <- tg
            else
                g <- g+tg
        }
        header = paste0(target, '_', genot, tail, btail, filt_prefix, time_prefix)
        if (is.null(ag)) ag <- g
        else ag <- ag + g
    }
    header = paste0(target, tail, btail, filt_prefix, time_prefix)
    png(paste0('heatmap_', header, '_all.png'), width=1200)
    # row_order_list = row_order(draw(ag))
    draw(ag + Heatmap(target_range$time, col = c("white", "orange"), name = "Replication time", show_row_names=FALSE, width=unit(5, 'mm')), row_order=row_order)
    write.table(as.data.frame(target_range)[unlist(row_order_list),], file=paste0('heatmap_', header, '_all_row_list.tsv'), sep='\t')
    dev.off()
    width = 15
    if (length(GENOT) == 4) width = 20
    pdf(paste0('heatmap_', header, '_all.pdf'), width=width)
    # row_order_list = row_order(draw(ag))
    draw(ag + Heatmap(target_range$time, col = c("white", "orange"), name = "Replication time", show_row_names=FALSE, width=unit(5, 'mm')), row_order=row_order)
    dev.off()
}

extract.cluster <- function(deg, ori_cluster, target_range) {
    if (deg == '') {
        cluster <- ori_cluster[ori_cluster[,2] != 'insignificant',]
        cluster[,2] <- unlist(sapply(cluster[,2], function(x){return(as.integer(sub('^cluster', '', x)))}))
        selected_target = target_range[unlist(sapply(target_range$gene_id, function(x){return(any(x == cluster[,1]))})),]
        cluster = cluster[sapply(selected_target$gene_id, function(x){return(which(x == cluster[,1])[1])}),]
    } else if (deg == '_global') {
        cluster = ori_cluster
        cluster[,2] = ''
        selected_target = target_range[unlist(sapply(target_range$gene_id, function(x){return(any(x == cluster[,1]))})),]
        cluster = cluster[sapply(selected_target$gene_id, function(x){return(which(x == cluster[,1])[1])}),]
    } else { #if (deg == '_all') 
        cluster = ori_cluster
        cluster[cluster[,2] != 'insignificant',2] = '1 significant'
        cluster[cluster[,2] == 'insignificant',2] = '2 insignificant'
        selected_target = target_range[unlist(sapply(target_range$gene_id, function(x){return(any(x == cluster[,1]))})),]
        cluster = cluster[sapply(selected_target$gene_id, function(x){return(which(x == cluster[,1])[1])}),]
    }
    return(list(cluster, selected_target))
}
 extract.each.rep.heatmap <- function(fc_range, target, selected_target, axis_name, tail, btail, filt_prefix, time_prefix, top_ann)
{
    mat = list()
    for (genot in GENOT) {
        print(genot)
        for (stage in STAGE) {
            for (i in 0:1) {
                key = paste0(genot, '_', stage, '_', i)
                header = paste0(target, '_', key, tail, btail, filt_prefix, time_prefix)
                if (any(key == colnames(data.frame(fc_range))))
                    mat[[key]] = normalizeToMatrix(fc_range, selected_target, value_column = key, mean_mode = "w0", w = 25, extend=ext)
            }
        }
    }
    for (genot in GENOT) {
        print(genot)
        for (stage in STAGE) {
            key_1 = paste0(genot, '_', stage, '_', 0)
            key_2 = paste0(genot, '_', stage, '_', 1)
            if (!any(names(mat) == key_2)) next
            x = as.vector(mat[[key_1]])
            y = as.vector(mat[[key_2]])
            x <- x[is.finite(y)]
            y <- y[is.finite(y)]
            y <- y[is.finite(x)]
            x <- x[is.finite(x)]
            png(paste0('correlation_rep_', genot, '_', stage, '.png'))
            plot(x, y)
            dev.off()
            print(paste('cor', genot, stage, cor(x, y)))
        }
    }
    return(list(mat))
}

add.cluster.heatmap <- function(fc_range, target, selected_target, axis_name, tail, btail, filt_prefix, time_prefix, top_ann)
{
    mat = list()
    ag <- NULL
    for (genot in GENOT) {
        g <- NULL
        for (stage in STAGE) {
            key = paste0(genot, '_', stage)
            header = paste0(target, '_', key, tail, btail, filt_prefix, time_prefix)
            mat[[key]] = normalizeToMatrix(fc_range, selected_target, value_column = key, mean_mode = "w0", w = 25, extend=ext)
            trim_mat = mat[[key]]
            trim_mat[trim_mat < 0] = 0
            trim_mat[trim_mat > 2] = 2
            g <- g + EnrichedHeatmap(trim_mat, name=key, col=col_fun, top_annotation=top_ann, axis_name=axis_name)
        }
        header = paste0(target, '_', genot, tail, btail, filt_prefix, time_prefix)
        ag <- ag + g
    }
    return(list(mat, ag))
}

extract.signals.around.tss <- function(mat, target_range, keys, clusters, header, genot='WT') {
    # 1000bp - 40 window
    all_mat <- list()
    MAX = 80
    for (key in keys) {
        tmat <- NULL
        for (i in 0:2) {
            if (!any(names(mat) == paste0(key, '_', i))) next
            if (i == 1) tmat <- mat[[paste0(key, '_', i)]]
            else tmat <- cbind(tmat, mat[[paste0(key, '_', i)]])
        }
        for (label in c('upstream', 'downstream')) {
            if (grepl(label, 'downstream')) {
                dists = c(0, 100, 200, 300, 400)
                pixels = c(0, 4, 8, 12, 16)
                max_pixels = 20 # 500bp
                for (i in 1:length(dists)) {
                    left_index = (MAX/2+1)
                    dist = 500-dists[i]
                    print(c(max_pixels-pixels[i], dist))
                    signals = rowSums(tmat[,(left_index+pixels[i]):(left_index+max_pixels)])/dist
                    result_label = paste0(key, '_', label, '_', dists[i])
                    if (is.null(all_mat[[result_label]])) {
                        all_mat[[result_label]] <- signals
                    } else {
                        all_mat[[result_label]] <- cbind(all_mat[[result_label]], signals)
                    }
                }
            } else {
                dists = c(500, 750, 1000)
                pixels = c(20, 30, 40)
                for (i in 1:length(dists)) {
                    left_index = (MAX/2)
                    signals = rowSums(tmat[,(left_index-pixels[i]+1):left_index])/dists[i]
                    result_label = paste0(key, '_', label, '_', dists[i])
                    if (is.null(all_mat[[result_label]])) {
                        all_mat[[result_label]] <- signals
                    } else {
                        all_mat[[result_label]] <- cbind(all_mat[[result_label]], signals)
                    }
                }
            }
        }
    }
    matrix <- do.call(cbind, all_mat)
    colnames(matrix) <- names(all_mat)
    rownames(matrix) <- rownames(mat[[paste0(keys[1], '_', 0)]])
    m <- merge(matrix, clusters, by.x=0, by.y='GeneID')
    for (udist in c(500, 750, 1000)) {
        for (ddist in c(-1, 0, 100, 200, 300, 400)) {
            x = obtain.upstream.signal(m, udist, paste0(genot, '_G1'), ddist)
            for (time in c('HU45', 'HU90')) {
                y = obtain.upstream.signal(m, udist, paste0(genot, '_', time), ddist)
                if (is.null(x) || is.null(y)) next
                df = data.frame(x=x, y=y, gene=m[,'Row.names'])
                m <- cbind(m, comp.resid(df))
                label = paste0(genot, '_', time, '_', udist)
                if (!is.null(ddist)) label = paste0(label, '_', ddist)
                colnames(m)[dim(m)[2]] <- label
            }
        }
    }
    for (time in c('HU45', 'HU90')) {
        if (!any(colnames(m) == paste0(genot, '_', time, '_upstream_', 500))) next
        png(paste0('correlation_upstream_', header, '_', time, '.png'))
        plot(m[,paste0(genot, '_', 'G1', '_upstream_', 500)], m[,paste0(genot, '_', time, '_upstream_', 500)])
        dev.off()
        print(paste('signal_cor', genot, time, cor(m[,paste0(genot, '_', 'G1', '_upstream_', 500)], m[,paste0(genot, '_', time, '_upstream_', 500)])))
    }
    m <- cbind(m, merge.peak.existence(target_range, m[,'Row.names']))
    m <- cbind(m, merge.origin.existence(target_range, m[,'Row.names']))
    write.table(m, stat.file.name(header), quote=F, sep="\t")
}

stat.file.name <- function(header, mod='') {
    if (VALID_STAT) {
        return(paste0('statistics_residual_', header, '_valid', mod, '.tsv'))
    } else if (ORI_STAT) {
        return(paste0('statistics_residual_', header, '_ori', mod, '.tsv'))
    } else {
        return(paste0('statistics_residual_', header, mod, '.tsv'))
    }
}
read.origin.bed <- function(valid=TRUE) {
    a <- read.table(file.path(CHIP_ANN, "origin", "origin_with_dna_amount_ann.tsv"), header=T, sep='\t', quote="")
    if (valid) {
        if (VALID_STAT)
            a <- subset(a, a[,'origin'] == TRUE)
        else # ori_id annotation
            a <- subset(a, a[,'gene'] != '')
    }
    origin_label <- rep('all', dim(a)[1])
    origin_label[a[,'origin'] == TRUE] = 'early'
    origin_label[a[,'late_origin'] == TRUE] = 'late'
    a <- cbind(a, origin_label=origin_label)
    a <- a[,c('chr', 'start', 'end', 'origin_label', 'gene', 'gene2')]
    bed <- makeGRangesFromDataFrame(a, keep.extra.columns=TRUE, ignore.strand=TRUE)
    return(bed)
}

merge.origin.existence <- function(target_range, genes) {
    selected_range = extract.grange(target_range, genes)
    bed <- read.origin.bed((VALID_STAT | ORI_STAT))
    peak <- nearest(selected_range, bed, select=c("arbitrary"), ignore.strand=TRUE)

    df <- cbind(data.frame(selected_range), data.frame(bed[peak,]))
    return(df)
}

extract.tss.and.tts.nearest.origin <- function(header) {
    for (genot in GENOT) {
        df <- NULL
        if (genot == 'WT')
            tail = ''
        else
            tail = paste0('_', genot)
        m <- read.table(stat.file.name(paste0(header, tail)), header=T, sep="\t", quote="")
        genes <- m[,'Row.names']
        for (target in c('tss', 'tts')) {
            target_range <- read.tss.data(target)
            selected_range <- target_range[unlist(sapply(names(target_range), function(x){return(any(x == genes))})),]
            bed <- read.origin.bed()
            overlap <- nearest(selected_range, bed, select=c("arbitrary"), ignore.strand=TRUE)
            dbed <- data.frame(bed[overlap,])
            dbed <- dbed[,c('start', 'end', 'width', 'origin_label', 'gene', 'gene2')]
            colnames(dbed) <- paste0(c('start_', 'end_', 'width_', '', 'gene_', 'gene2_'), 'origin_', target)
            df <- data.frame(selected_range)
            df <- df[,c('seqnames', 'start', 'end', 'strand', 'gene_id')]
            colnames(df) <- c(paste0(colnames(df)[1:4], '_gene_', target), paste0('gene_id_', target))
            df <- cbind(df, dbed)
            p <- t(apply(df, c(1), function(x){
                gstrand <- x[paste0('strand_gene_', target)]
                if (gstrand == '+')
                    gstart <- as.numeric(x[paste0('start_gene_', target)])
                else
                    gstart <- as.numeric(x[paste0('end_gene_', target)])
                ostart <- as.numeric(x[paste0('start', '_origin_', target)])
                oend <- as.numeric(x[paste0('end', '_origin_', target)])
                ocenter <- (ostart+oend)/2
                if (ostart <= gstart && gstart < oend)
                    dist <- 0
                else
                    dist <- min(abs(gstart-c(ostart, oend)))
                label = "NA"
                if (gstrand == '+') {
                    if (ocenter < gstart)
                        label = "codirectional"
                    else if (gstart < ocenter)
                        label = "head on"
                } else {
                    if (gstart < ocenter)
                        label = "codirectional"
                    else if (ocenter < gstart)
                        label = "head on"
                }
                return(c(dist, label))
            }))
            dist <- as.matrix(p)
            colnames(dist) <- paste0(c('dist_', 'location_'), target)
            df <- cbind(df, dist)
            m <- merge(m, df, by.x='Row.names', by.y=paste0('gene_id_', target), all.x=TRUE)
        }
        write.table(m, stat.file.name(paste0(header, tail), '_mod'))
    }
}

extract.grange <- function(target_range, genes) {
    selected_range <- target_range[unlist(sapply(names(target_range), function(x){return(any(x == genes))})),]
    selected_range <- selected_range[unlist(sapply(genes, function(x) {
        return(which(x == selected_range$gene_id)[1])
    })),]
    stopifnot(all(unlist(sapply(1:length(genes), function(x){return(genes[x] == names(selected_range)[x])}))))
    return(selected_range)
}

merge.peak.existence <- function(target_range, genes) {
    require(genomation)
    selected_range <- extract.grange(target_range, genes)
    bed <- readBed(PEAK_FILE)
    peak <- nearest(selected_range, bed, select=c("arbitrary"), ignore.strand=TRUE)
    df <- cbind(data.frame(selected_range), data.frame(bed[peak,]))
    id <- read.table(file.path(ANN_DIR, 'SGD_orf_id.tab'), header=F, sep="\t", quote="\"")
    colnames(id) = c('gene_id', 'gene_name')
    return(merge(df, id, by='gene_id', all.x=T))
}
comp.resid <- function(df) {
    require(ramify)
    linear <- lm(y ~ x, data=df)
    return(resid(linear))
}
obtain.upstream.signal <- function(mat, udist, label, ddist) {
    if (!any(colnames(mat) == paste0(label, '_upstream_', udist))) return(NULL)
    signal <- mat[,paste0(label, '_upstream_', udist)]
    if (ddist >= 0) {
        signal <- signal-mat[,paste0(label, '_downstream_', ddist)]
    }
    return(signal)
}

clip.heatmap.mat <- function(mat, keys, average=TRUE) {
    require(ramify)
    all_mat <- list(ave_sig=NULL, ave_caseNULL)
    for (key in keys) {
        print(key)
        tmat = mat[[key]]
        tmat = clip(tmat, 0, max(tmat))
        MAX = 80
        dist = 19
        stopifnot(dim(tmat)[2] == MAX)
        if (average) {
            ave_sig = rowSums(tmat[,(MAX/2-dist):MAX/2])/dist-rowSums(tmat[,(MAX/2+1):(MAX/2+1+dist)])/dist
            ave_case = rowSums(tmat[,(MAX/2-dist):MAX/2])/dist
        } else {
            ave_sig = rowSums(tmat[,(MAX/2-dist):MAX/2])-rowSums(tmat[,(MAX/2+1):(MAX/2+1+dist)])
            ave_sig = clip(ave_sig, 0, max(ave_sig))
        }
        if (is.null(all_mat[[1]])) {
            all_mat[[1]] <- ave_sig
            all_mat[[2]] <- ave_case
        } else {
            all_mat[[1]] <- cbind(all_mat[[1]], ave_sig)
            all_mat[[2]] <- cbind(all_mat[[2]], ave_case)
        }
    }
    rownames(all_mat[[1]]) = rownames(mat[[keys[1]]])
    rownames(all_mat[[2]]) = rownames(mat[[keys[1]]])
    return(all_mat)

}

read.fc.data <- function(fprefix, cluster) {
    rna_table = read.table(file.path(CHIP_ANN, "deg", fprefix), header=T, sep="\t", stringsAsFactor=F)
    fc = unlist(sapply(cluster[,1], function(x) {
        index = which(x == rownames(rna_table))
        if (length(index) == 0) return(0)
        else return(rna_table[index[1],2])
    }))
    row_order=order(cluster[,2], -fc, decreasing=F)
    fc[fc < -2] = -2
    fc[fc > 2] = 2
    return(list(fc, row_order))
}

check.tss.cluster.signals <- function(fc_range, target, target_range, tail, ext, btail, col_fun, filt_prefix, MIN, MAX, global=FALSE, add_fc=TRUE)
{
    fc_prefix = ''
    if (add_fc) fc_prefix = '_fc'
    time_prefix = ''
    axis_name = c(paste0('-', ext), '0', paste0('+', ext))
    cluster_deg <- NULL
    if (OPTIONAL) {
        header = paste0(target, tail, btail, filt_prefix, fc_prefix, time_prefix)
        extract.tss.and.tts.nearest.origin(header)
        return()
    }
    for (file in read.table(CLUSTER_FILE, header=F, quote='\"')[,1]) {
        if (!grepl('WT', file) || !grepl('time', file)) next
        ori_cluster = read.table(file.path(CHIP_ANN, 'coexp', file), header=T, sep=" ", stringsAsFactor=F)
        if (is.null(cluster_deg)) {
            cluster_deg <- ori_cluster
            colnames(cluster_deg) <- c('GeneID', paste0(basename(file), '_cluster'))
        } else {
            colnames(ori_cluster) <- c('GeneID', paste0(basename(file), '_cluster'))
            cluster_deg <- merge(cluster_deg, ori_cluster, by='GeneID')
        }
        ctail = sub('\\.tsv$', '', basename(file))
        header = paste0(target, tail, btail, filt_prefix, time_prefix)
        fprefix = gsub("_cluster.*.tsv$", ".tsv", gsub("meta_", "", basename(file)))
        result <- read.fc.data(fprefix, cluster_deg)
        fc <- result[[1]]
        cluster_deg <- cbind(cluster_deg, fc[cluster_deg[,1]])
        colnames(cluster_deg)[dim(cluster_deg)[2]] <- paste0(basename(file), '_fc')
    }
    result <- extract.each.rep.heatmap(fc_range, target, target_range, axis_name, tail, btail, paste0(filt_prefix, fc_prefix), time_prefix, top_ann)
    header = paste0(target, tail, btail, filt_prefix, fc_prefix, time_prefix)
    mat = result[[1]]
    for (genot in GENOT) {
        if (genot == 'WT')
            extract.signals.around.tss(mat, target_range, paste0(genot, '_', STAGE), cluster_deg, header)
        else
            extract.signals.around.tss(mat, target_range, paste0(genot, '_', STAGE), cluster_deg, paste0(header, '_', genot), genot=genot)
    }
}

plot.cluster.annotated.heatmap.global <- function(fc_range, target, target_range, tail, ext, btail, col_fun, filt_prefix, fc_prefix, time_prefix, MIN, MAX, global, add_fc, cluster_row, ori_cluster, file, axis_name) {
    row_order = NULL
    deg <- c('_global')
    if (cluster_row) {
        results = extract.cluster(deg, ori_cluster, target_range)
        cluster <- results[[1]]
        selected_target <- results[[2]]
        max_cluster = length(unique(cluster[,2]))
    } else {
        results = extract.cluster(deg, ori_cluster, target_range)
        selected_target <- results[[2]]
        max_cluster = 1
    }
    if (tail == '_prop')
        top_ann = HeatmapAnnotation(enriched=anno_enriched(ylim=c(MIN, MAX), gp=gpar(col=rainbow(max_cluster))))
    else
        top_ann = HeatmapAnnotation(enriched = anno_enriched(gp=gpar(col=rainbow(max_cluster))))
    if (cluster_row && add_fc) {
        fprefix = gsub("_cluster.*.tsv$", ".tsv", gsub("meta_", "", basename(file)))
        rna_table = read.table(file.path(CHIP_ANN, "deg", fprefix), header=T, sep="\t", stringsAsFactor=F)
        fc = unlist(sapply(cluster[,1], function(x) {
            index = which(x == rownames(rna_table))
            if (length(index) == 0) return(0)
            else return(rna_table[index[1],2])
        }))
        row_order=order(cluster[,2], -fc, decreasing=F)
        fc[fc < -2] = -2
        fc[fc > 2] = 2
    }
    if (cluster_row){
        cluster_ann = Heatmap(cluster[,2], col = structure(rainbow(max_cluster), names = sort(unique(cluster[,2]))), name = "cluster",
                              show_row_names = FALSE, width = unit(3, "mm"), row_split = cluster[,2], row_order=row_order)
    } else {
        cluster_ann = NULL
    }
    if (add_fc) {
        cluster_ann <- cluster_ann+Heatmap(fc,  col = c("blue", "white", "red"), name = "log2FC",
                        show_row_names = FALSE, width = unit(3, "mm"))
    }
    result <- add.cluster.heatmap(fc_range, target, selected_target, axis_name, tail, btail, paste0(filt_prefix, fc_prefix), time_prefix, top_ann)
    ag <- result[[2]]
    width = 12
    if (length(GENOT) == 4) width =15
    if (cluster_row) ctail = sub('\\.tsv$', '', basename(file))
    else ctail = ''
    header = paste0(target, tail, btail, filt_prefix, fc_prefix, time_prefix)
    png(paste0('heatmap_', header, ctail, deg, '_all.png'), width=1200)
    draw(cluster_ann+ag, row_order=row_order)
    dev.off()
    pdf(paste0('heatmap_', header, ctail, deg, '_all.pdf'), width=width)
    draw(cluster_ann+ag, row_order=row_order)
    dev.off()
    return()
}

plot.cluster.annotated.heatmap.deg <- function(fc_range, target, target_range, tail, ext, btail, ctail, col_fun, filt_prefix, fc_prefix, time_prefix, MIN, MAX, global, add_fc, cluster_row, ori_cluster, file, axis_name) {
    degs = c('_all', '')
    for (deg in degs) {
        results = extract.cluster(deg, ori_cluster, target_range)
        cluster <- results[[1]]
        selected_target <- results[[2]]
        max_cluster = length(unique(cluster[,2]))
        if (tail == '_prop')
            top_ann = HeatmapAnnotation(enriched=anno_enriched(ylim=c(MIN, MAX), gp=gpar(col=rainbow(max_cluster))))
        else
            top_ann = HeatmapAnnotation(enriched = anno_enriched(gp=gpar(col=rainbow(max_cluster))))
        row_order = NULL
        if (add_fc) {
            fprefix = gsub("_cluster.*.tsv$", ".tsv", gsub("meta_", "", basename(file)))
            rna_table = read.table(file.path(CHIP_ANN, "deg", fprefix), header=T, sep="\t", stringsAsFactor=F)
            fc = unlist(sapply(cluster[,1], function(x) {
                index = which(x == rownames(rna_table))
                if (length(index) == 0) return(0)
                else return(rna_table[index[1],2])
            }))
            row_order=order(cluster[,2], -fc, decreasing=F)
            fc[fc < -2] = -2
            fc[fc > 2] = 2
        }
        if (cluster_row){
            cluster_ann = Heatmap(cluster[,2], col = structure(rainbow(max_cluster), names = sort(unique(cluster[,2]))), name = "cluster",
                show_row_names = FALSE, width = unit(3, "mm"), row_split = cluster[,2], row_order=row_order)
        } else {
            cluster_ann = NULL
        }
        if (add_fc) {
            cluster_ann <- cluster_ann+Heatmap(fc,  col = c("blue", "white", "red"), name = "log2FC",
                            show_row_names = FALSE, width = unit(3, "mm"),  show_row_dend = FALSE)
        }
        result <- add.cluster.heatmap(fc_range, target, selected_target, axis_name, tail, btail, paste0(filt_prefix, fc_prefix), time_prefix, top_ann)
        ag <- result[[2]]
        width = 12
        if (length(GENOT) == 4)
            width =15
        header = paste0(target, tail, btail, filt_prefix, fc_prefix, time_prefix)
        png(paste0('heatmap_', header, ctail, deg, '_all.png'), width=1200)
        draw(cluster_ann+ag, row_order=row_order)
        dev.off()
        pdf(paste0('heatmap_', header, ctail, deg, '_all.pdf'), width=width)
        # row_order_list = row_order(draw(ag))
        draw(cluster_ann+ag, row_order=row_order)
        dev.off()
        write.table(cluster[row_order,], file=paste0('heatmap_', header, ctail, deg, '_all_row_list.tsv'), sep='\t')
    }
}


plot.cluster.annotated.heatmap <- function(fc_rnage, target, target_range, tail, ext, btail, col_fun, filt_prefix, MIN, MAX, global=FALSE, add_fc=TRUE, cluster_row=TRUE)
{
    fc_prefix = ''
    if (add_fc) fc_prefix = '_fc'
    time_prefix = ''
    axis_name = c(paste0('-', ext), '0', paste0('+', ext))
    for (file in read.table(CLUSTER_FILE, header=F, quote='\"')[,1]) {
        ori_cluster = read.table(file.path(CHIP_ANN, 'coexp', file), header=T, sep=" ", stringsAsFactor=F)
        ctail = sub('\\.tsv$', '', basename(file))
        header = paste0(target, tail, btail, filt_prefix, time_prefix)
        if (global) {
            plot.cluster.annotated.heatmap.global(fc_range, target, target_range, tail, ext, btail, col_fun, filt_prefix, fc_prefix, time_prefix, MIN, MAX, global, add_fc, cluster_row, ori_cluster, file, axis_name)
            if (!cluster_row) break
        } else {
            plot.cluster.annotated.heatmap.deg(fc_range, target, target_range, tail, ext, btail, ctail, col_fun, filt_prefix, fc_prefix, time_prefix, MIN, MAX, global, add_fc, cluster_row, ori_cluster, file, axis_name)
        }
        break
    }
}

if (HYPER_CHIP == 'wohc') {
    btail_list <- c('_wohc')
} else if (HYPER_CHIP == 'ohc') {
    btail_list <- c('_ohc')
} else if (HYPER_CHIP == 'pho') {
    btail_list <- c('_pho')
} else {
    btail_list <- c(DEFAULT)
}

for (btail in btail_list) {
    # fold change
    rep_tail <- ''
    if (statistics) {
        rep_tail <- '_rep'
    }
    if (btail == '_wohc' || btail == '_ohc' || btail == '_pho') {
        fc <- read.table(paste0('chip_case_cont_mean', DEFAULT, rep_tail, '.tsv'), header=T, sep='\t', stringsAsFactor=F)
    } else {
        fc <- read.table(paste0('chip_case_cont_mean', btail, rep_tail, '.tsv'), header=T, sep='\t', stringsAsFactor=F)
    }
    fc_range <- makeGRangesFromDataFrame(fc, keep.extra.columns=TRUE, ignore.strand=TRUE)
    filt = TRUE
    filt_prefix = '_wo_hc'
    time = TRUE
    time_prefix = '_time'
    MIN = 0
    MAX = 1.2
    args = commandArgs(trailingOnly=TRUE)
    dataset = ''
    if (length(args) >= 1) dataset = args[1]
    set.genotype.stage(dataset)
    if (length(args) >= 2) {
        MIN = args[2]
        MAX = args[3]
    }
    if (any(args == 'global')) global <<- TRUE
    if (any(args == 'statistics')) statistics <<- TRUE


    # for (target in c('tss', 'tts', 'all', 'valid', 'early', 'late')) {
    for (target in c('tss', 'valid')) {
        if (relative) {
            require(circlize)
            tail = '_prop'
            if (any(target == c('valid', 'early', 'late', 'all'))) {
                col_fun = colorRamp2(c(-0, 0.5, 1.5), c("blue", "white", "yellow"))
                top_ann = HeatmapAnnotation(enriched=anno_enriched(ylim=c(MIN, MAX)))
            } else {
                col_fun = colorRamp2(c(-0, 0.5, 2.0), c("blue", "white", "yellow"))
                top_ann = HeatmapAnnotation(enriched=anno_enriched(ylim=c(MIN, MAX)))
            }
        } else {
            tail = ''
            col_fun = c("blue", "white", "red")
            top_ann = HeatmapAnnotation(enriched = anno_enriched())
        }
        if (grepl('tss', target) || grepl('tts', target)) {
            target_range = read.tss.data(target)
            ext = 1000
            if (statistics) {
                check.tss.cluster.signals(fc_range, target, target_range, tail, ext, btail, col_fun, filt_prefix, MIN, MAX, global)
                break
            } else {
                plot.cluster.annotated.heatmap(fc_range, target, target_range, tail, ext, btail, col_fun, filt_prefix, MIN, MAX, global, add_fc=TRUE, cluster_row=TRUE)
            }
        } else {
            if (statistics) 
                break
            target_range = read.origin.data(target, filt, time=time, ori.id=ORI_STAT)
            ext = 15000
            plot.cluster.annotated.origin.heatmap(fc_range, target, target_range, tail, ext, btail, col_fun, top_ann, filt_prefix, time_prefix)
        }
    }
}
