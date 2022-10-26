
dir = "./"
oname = paste0(dir, "genot_difference_")
data <- NULL
for (g in c("Rad53K", "Mrc1d", "Rad9d")) {
    ofile = paste0(oname, g, ".tsv")
    if (!file.exists(ofile)) next
    a <- read.table(ofile, header=T, sep="\t")
    a <- a[,c(2, 6)]
    colnames(a) = c(paste0("log2FC_", g), paste0("padj_", g))
    if (is.null(data)) {
        data <- a
    } else {
        stopifnot(rownames(data) == rownames(a))
        data <- cbind(data, a)
    }
}
if (!is.null(data))
    write.table(data, "de_result_mc_genotype.tsv", sep="\t", quote=F)

oname = paste0(dir, "genotype_difference_")
data <- NULL
for (time in c("G1", "HU45", "HU90")) {
    for (g in c("Rad53K", "Mrc1d", "Rad9d")) {
        ofile = paste0(oname, time, "_", g, ".tsv")
        if (!file.exists(ofile)) next
        a <- read.table(ofile, header=T, sep="\t")
        a <- a[,c(2, 6)]
        colnames(a) = c(paste0("log2FC_", time, "_", g), paste0("padj_", time, "_", g))
        if (is.null(data)) {
            data <- a
        } else {
            stopifnot(rownames(data) == rownames(a))
            data <- cbind(data, a)
        }
    }
}
if (!is.null(data))
    write.table(data, "de_result_sc_genotype.tsv", sep="\t", quote=F)

oname = paste0(dir, "time_difference_")
data <- NULL
for (time in c("HU45", "HU90")) {
    for (g in c("WT", "Rad53K", "Mrc1d", "Rad9d")) {
        ofile = paste0(oname, time, "_", g, ".tsv")
        if (!file.exists(ofile)) next
        a <- read.table(ofile, header=T, sep="\t")
        a <- a[,c(2, 6)]
        colnames(a) = c(paste0("log2FC_", time, "_", g), paste0("padj_", time, "_", g))
        if (is.null(data)) {
            data <- a
        } else {
            stopifnot(rownames(data) == rownames(a))
            data <- cbind(data, a)
        }
    }
}
if (!is.null(data))
    write.table(data, "de_result_sc_time.tsv", sep="\t", quote=F)

oname = paste0(dir, "genotype_difference_")
data <- NULL
for (time in c("G1", "HU45", "HU90")) {
    genotypes = c("WT", "Rad53K", "Mrc1d", "Rad9d")
    for (gi in 1:3) {
        for (gh in (gi+1):4) {
            ofile = paste0(oname, time, "_", genotypes[gh], "_", genotypes[gi], ".tsv")
            if (!file.exists(ofile)) next
            a <- read.table(ofile, header=T, sep="\t")
            a <- a[,c(2, 6)]
            colnames(a) = c(paste0("log2FC_", time, "_", genotypes[gh], "_", genotypes[gi]), paste0("padj_", time, "_", genotypes[gh], "_", genotypes[gi]))
            if (is.null(data)) {
                data <- a
            } else {
                stopifnot(rownames(data) == rownames(a))
                data <- cbind(data, a)
            }
        }
    }
}
if (!is.null(data))
    write.table(data, "de_result_sc_each_genot.tsv", sep="\t", quote=F)
