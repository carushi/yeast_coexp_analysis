require(pROC)
require(matrixStats)
threshold = 0.05

ScaleResidualMatrix <- function(res.data, time=45, max_value=NA) {
    if (is.na(max_value)) {
        max_value = max(unlist(res.data[, c('WT_G1_upstream_500', paste0('WT_HU', 45, '_upstream_500'), paste0('WT_HU', 90, '_upstream_500'))]), na.rm=TRUE)
    }
    for (col in c('WT_G1_upstream_500', 'WT_HU45_upstream_500', 'WT_HU90_upstream_500')) {
        res.data[,col] = res.data[,col]/max_value
    }
    return(res.data)
}

ComputeAUROC <- function(dbr, gene.set) {
    flag = unlist(sapply(gene.set, function(x) {
        if (any(x == dbr)) return(1)
        return(0)
    }))
    rocScore = roc(flag, 1:length(gene.set))
    print(auc(rocScore))
}

ComputeResidAUROC <- function(top.diff.gene, resid, res.filter) {
    if (!is.null(res.filter)) {
        res.filter[is.na(res.filter)] = 0
        resid[res.filter] = NA
    }
    resid = resid[!is.na(top.diff.gene)]
    top.diff.gene = top.diff.gene[!is.na(top.diff.gene)]
    rocScore = roc(top.diff.gene, resid)
    print(auc(rocScore))
    return(rocScore)
}

TopResidual <- function(res, hu, g1, res.filter, header) {
    require(ggplot2)
    thres = abs(res[!res.filter][order(abs(res[!res.filter]), na.last=TRUE, decreasing=TRUE)[1000]])
    df = data.frame(HU=hu, G1=g1)
    flag=((!res.filter) & (abs(res) >= thres))
    stopifnot(sum(flag) == 1000)
    df = cbind(df, flag=flag)
    g <- ggplot(df, aes(x=G1, y=HU, color=flag))+geom_point()+theme_bw()+ylim(NA, 1)+xlim(NA, 1)
    pdf(paste0('plot_resid_', header, '.pdf'))
    plot(g)
    dev.off()
    return(df[,3])
}

SetSignificanceByBoth <- function(df) {
    vec = unlist(sapply(1:dim(df)[1], function(x) {
        if (is.na(df[x,'significance_residual']) || is.na(df[x,'significance_diffbind'])) {
            return('Non')
        }
        if (df[x,'significance_residual'] == 'Significant') {
            if (df[x,'significance_diffbind'] == 'Significant') {
                return('Both')
            } else {
                return('Residual only')
            }
        } else {
            if (df[x,'significance_diffbind'] == 'Significant') {
                return('Diff only')
            } else {
                return('Non')
            }
        }
    }))
    return(vec)
}

SetSignificanceByOne <- function(vec, threshold=NA) {
    if (is.na(threshold)) {
        significance = unlist(sapply(vec, function(x) {
            if (is.na(x)) return('Non')
            if (x == 1) return('Significant')
            else return('Non')
        }))
    } else {
        significance = unlist(sapply(vec, function(x) {
            if (is.na(x)) return('Non')
            else if (x <= threshold) return('Significant')
            else return('Non')
        }))
    }
    return(significance)
}

SignificanceOverlap <- function(fc, res, pvalue, header, top.res.flag, res.filter, threshold) {
    require(ggplot2)
    df = cbind(fc, pvalue)
    res_df = cbind(res, top.res.flag)
    res_df = cbind(res_df, res.filter)
    rownames(res_df) = names(res)
    df = merge(df, res_df, by=0, all=TRUE)
    colnames(df) = c('gene_id', 'DiffBind', 'pvalue', 'Residual', 'significance_residual', 'res_filtered')
    df[,'significance_residual'] = SetSignificanceByOne(df[,'significance_residual'], NA)
    df = cbind(df, significance_diffbind=SetSignificanceByOne(df[,'pvalue'], threshold))
    df = cbind(df, significance=SetSignificanceByBoth(df))
    write.table(df, paste0('matrix_resid_diffbind', header, '.tsv'), sep="\t")
}

PlotCorrelation <- function(df, header) {
    require(ggplot2)
    cont.table = table(df[,'significance'])
    print(fisher.test(matrix(cont.table[c(1,2, 4, 3)], nrow=2), alternative='greater'))
    df = df[!is.na(df[,'Residual']),]
    df = df[order(df[,'significance_residual'], na.last=TRUE),]
    pdf(paste0('comp_fc_pvalue_', header, '.pdf'))
    g <- ggplot(df, aes(x=DiffBind, y=pvalue, color=significance_residual))+geom_point(alpha=0.3)+theme_bw()+scale_y_log10()
    plot(g)
    dev.off() 
    pdf(paste0('comp_resid_pvalue_', header, '.pdf'))
    g <- ggplot(df, aes(x=Residual, y=pvalue, color=significance_residual))+geom_point(alpha=0.3)+theme_bw()+scale_y_log10()
    plot(g)
    dev.off()
    print(c('correlation', cor(df[,'Residual'], df[,'DiffBind'], method='spearman', use="complete.obs")))
    for (col in c('significance_residual', 'significance_diffbind', 'significance')) {
        g <- ggplot(df, aes(x=Residual, y=DiffBind))+geom_point(aes_string(fill=col, color=col))+theme_bw()
        pdf(paste0('comp_resid_diffbind_', header, '_', col, '.pdf'))
        plot(g)
        dev.off()
    }
}

FilterResidual <- function(res.data, time) {
    return(res.data[,paste0('WT_HU', 45, '_upstream_500')] < -0.075 | res.data[,paste0('WT_G1_upstream_500')] < -0.075 | res.data[,paste0('WT_HU', 90, '_upstream_500')] < -0.075)
}

checkTopResData <- function(res.data, top.res.flag) {
    top_db <- read.table("../../ann/gene/top_db_lu.csv", header=T, sep=",", stringsAsFactors=F)    
    print(min(abs(res.data[top.res.flag, 'WT_HU45_500_.1'])))
    index = unlist(sapply(top_db[,'gene_id'], function(x) {
        return(which(x == res.data[,'gene_id']))
    }))
    print(length(index))
    print(dim(top_db))
    print(min(abs(res.data[index,'WT_HU45_500_.1'])))
    print(min(res.data[index,'WT_G1_upstream_500']))
    print(min(res.data[index,'WT_HU45_upstream_500']))
    print(min(abs(res.data[top.res.flag, 'WT_HU45_500_.1'])))
    result = unlist(sapply(top_db[,'gene_id'], function(x) {
        return(any(x == res.data[top.res.flag,'gene_id']))
    }))
    print(c(min(res.data[top.res.flag,'WT_HU45_500_.1']), max(res.data[top.res.flag, 'WT_HU45_500_.1']), min(abs(res.data[top.res.flag, 'WT_HU45_500_.1']))))
    stopifnot(all(result))
}

extractMaxThres <- function() {
    thres = list(CP=0, TF=0)
    for (data_label in c('TF', 'CP')) {
        if (data_label == 'CP') {
            res.data <- read.table("statistics_residual_tss_prop__wo_hc_fc_mod.tsv", header=T, sep=" ", quote="\"", stringsAsFactors=F)
        } else {
            res.data <- read.table("statistics_residual_tss_prop__wo_hc_fc_mod.tsv", header=T, sep=" ", quote="\"", stringsAsFactors=F)
        }
        new_thres = max(res.data[,paste0('WT_', c('G1', 'HU45', 'HU90'), '_upstream_500')], na.rm=TRUE)
        new_thres = max(new_thres, max(res.data[,paste0('WT_', c('G1', 'HU45', 'HU90'), '_downstream_0')], na.rm=TRUE))
        thres[[data_label]] = new_thres
    }
    return(thres)
}

plotDiffResidual <- function(fc, pvalue, res.data, top.res.flag, res.filter, top.diff.flag, time, data_label) {
    df = cbind(cbind(fc, pvalue), top.diff.flag)
    res =  cbind(y=res.data[,paste0('WT_HU', time, '_upstream_500')], x=res.data[,paste0('WT_G1_upstream_500')])
    rownames(res) = res.data[,'gene_id']
    saveRDS(df, 'test_df.rds')
    res_df = cbind(res, top.res.flag)
    res_df = cbind(res_df, res.filter)
    rownames(res_df) = rownames(res)
    saveRDS(res_df, 'test_res.rds')
    df = merge(df, res_df, by=0, all=TRUE)
    print(head(df))
    df[,'top.diff.flag'] = factor(df[,'top.diff.flag'])
    g <- ggplot(df)+geom_point(aes(x=x, y=y, color=top.diff.flag), alpha=0.3)+theme_classic()+scale_color_manual(values=c('gray', 'red'))
    pdf(paste0('plot_rep_diffbind_', time, '_', data_label, '.pdf'), useDingbats=FALSE)
    plot(g)
    dev.off()
}

max_value = extractMaxThres()



for (filt in c('', '_subt')[1]) {
    mat <- read.table(paste0('final_result_diffbind', filt, '.tsv'), sep="\t", quote="\"", header=T, stringsAsFactors=FALSE)
    for (time in c(45, 90)) {
        max_value = extractMaxThres()
        print(c('max threshold', max_value))
        for (data_label in c('TF', 'CP')) {
            header = paste0(filt, '_', data_label, '_', time)
            if (data_label == 'CP') {
                res.data <- read.table("statistics_residual_tss_prop__wo_hc_fc_mod.tsv", header=T, sep=" ", quote="\"", stringsAsFactors=F)
            } else {
                res.data <- read.table("statistics_residual_tss_prop__wo_hc_fc_mod.tsv", header=T, sep=" ", quote="\"", stringsAsFactors=F)
            }
            res.data <- ScaleResidualMatrix(res.data, time, max_value[[data_label]])
            pvalue = mat[,paste0('WT_HU', time, '_FDR')]
            fc = mat[,paste0('WT_HU', time, '_Fold')]
            names(fc) = mat[,'gene_id']
            resid = res.data[,paste0('WT_HU', time, '_500_.1')]
            names(resid) = res.data[,'gene_id']
            res.filter = FilterResidual(res.data, time)
            top.diff.flag = (pvalue <= threshold)
            top.res.flag = TopResidual(resid, res.data[,paste0('WT_HU', time, '_upstream_500')], res.data[,paste0('WT_G1_upstream_500')], res.filter, paste0('WT_HU', data_label))
            plotDiffResidual(fc, pvalue, res.data, top.res.flag, res.filter, top.diff.flag, time, data_label)
            SignificanceOverlap(fc, resid, pvalue, header, top.res.flag, res.filter, threshold)
        }
    }
}

for (filt in c('', '_subt')[1]) {
    add = 1
    cols = c('black', 'red', 'gray', 'orange')
    time = 45
    for (data_label in c('CP', 'TF')) {
        header = paste0(filt, '_', data_label, '_', time)
        df = read.table(paste0('matrix_resid_diffbind', header, '.tsv'), sep="\t", header=T, stringsAsFactors=FALSE)
        PlotCorrelation(df, header)
    }
    pdf(paste0('roc_all.pdf'))
    for (data_label in c('CP', 'TF')) {
        header = paste0(filt, '_', data_label, '_', time)
        print(header)
        df = read.table(paste0('matrix_resid_diffbind', header, '.tsv'), sep="\t", header=T, stringsAsFactors=FALSE)
        roc = ComputeResidAUROC(df[,'significance_diffbind'] == 'Significant', abs(df[,'Residual']), df[,'res_filtered'])
        plot(roc, col = cols[add], add = (add != 1))
        add = add+1
        roc = ComputeResidAUROC(df[,'significance_residual'] == 'Significant', -log10(df[,'pvalue']), NA)
        plot(roc, col = cols[add], add = (add != 1))
        add = add+1
    }
    dev.off()
}

