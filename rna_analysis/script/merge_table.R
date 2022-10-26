a <- list.files(path=".", pattern="*Gene.out.tab")
a <- a[order(a)]
all <- NULL
for (f in a) {
    header = unlist(strsplit(f, '-'))[1]
    data <- read.table(f, sep="\t", row.names=1)
    temp = cbind(rownames(data), data[,3,drop=F]) # stranded
    temp = cbind(rownames(data), data[,1,drop=F]) # unstranded
    colnames(temp) = c("ID", header)
    if (is.null(all)) {
        all <- temp
    } else {
        all <- merge(all, temp, all=T, by="ID")
    }
}
write.table(all, file='yeast_count_data_with_metadata.tsv', sep="\t")
for (l in c("N_ambiguous", "N_multimapping", "N_noFeature", "N_unmapped")) {
    all <- subset(all, all[,1] != l)
}
write.table(all, file='yeast_count_data.tsv', sep="\t")

