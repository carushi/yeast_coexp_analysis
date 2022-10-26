library(DESeq2)
library(ggplot2)


data <- read.table('yeast_count_data.tsv', sep="\t", header=T)
annot_df$Unique = as.character(rep(1:12,2))
rownames(data) = data[,1]
data <- data[,-c(1)]
data <- data[,rownames(annot_df)]
fact_ann <- annot_df
fact_ann[,1] <- factor(fact_ann[,1], levels=c("G1", "HU45", "HU90"))
fact_ann[,2] <- factor(fact_ann[,2], levels=c("WT", "Rad53K", "Mrc1d", "Rad9d"))


d2s <- DESeqDataSetFromMatrix(countData=data, colData=fact_ann, design=~ Stage + Genotype)
d2s <- estimateSizeFactors(d2s)
d2s <- estimateDispersions(d2s)
d2s <- nbinomWaldTest(d2s)
for (x in 2:4) {
    res <- results(d2s, contrast=c("Genotype", as.character(annot_df$Genotype[x]), "WT"), independentFiltering=F, cooksCutoff=Inf)
    write.table(res, paste0("genot_difference_", annot_df$Genotype[x], ".tsv"), sep="\t")
}

for (x in c("HU45", "HU90")) {
    res <- results(d2s, contrast=c("Stage", x, "G1"),independentFiltering=F, cooksCutoff=Inf)
    write.table(res, paste0("stage_difference_", x, ".tsv"), sep="\t")
}

d2s <- DESeqDataSetFromMatrix(countData=data, colData=fact_ann, design=~ Unique)
d2s <- estimateSizeFactors(d2s)
d2s <- estimateDispersions(d2s)
d2s <- nbinomWaldTest(d2s)
png(paste0("dispersion_rna_seq_unique.png"))
plotDispEsts(d2s)
dev.off()

for (x in 1:4) {
    res <- results(d2s, contrast=c("Unique", x+4, x), independentFiltering=F, cooksCutoff=Inf)
    write.table(res, paste0("time_difference_HU45_", annot_df$Genotype[x], ".tsv"), sep="\t")
    res <- results(d2s, contrast=c("Unique", x+8, x), independentFiltering=F, cooksCutoff=Inf)
    write.table(res, paste0("time_difference_HU90_", annot_df$Genotype[x], ".tsv"), sep="\t")
}


for (i in 1:3) {
    start = seq(1, 12, 4)[i]
    time = c("G1", "HU45", "HU90")[i]
    for (j in 1:3) {
        res <- results(d2s, contrast=c("Unique", start+j, start), independentFiltering=F, cooksCutoff=Inf)
        ofname = paste0("genotype_difference_", time, "_", annot_df$Genotype[1+j], ".tsv")
        write.table(res, ofname, sep="\t")
        print(annot_df[start,])
        print(annot_df[start+j,])
        print(c(i, ofname))
    }
}

print(annot_df)
for (i in 1:3) {
    start = seq(1, 12, 4)[i]
    time = c("G1", "HU45", "HU90")[i]
    for (j in 1:3) {
        res <- results(d2s, contrast=c("Unique", start+j, start), independentFiltering=F, cooksCutoff=Inf)
        ofname = paste0("genotype_difference_", time, "_", annot_df$Genotype[1+j], ".tsv")
        write.table(res, ofname, sep="\t")
        print(annot_df[start,])
        print(annot_df[start+j,])
        print(c(i, ofname))
    }
}

for (i in 1:3) {
    start = seq(1, 12, 4)[i]
    time = c("G1", "HU45", "HU90")[i]
    for (j in 0:2) {
        for (h in (j+1):3) {
            case <- data[,annot_df$Unique == start+h]
            cont <- data[,annot_df$Unique == start+j]
            case <- rowMeans(case)
            cont <- rowMeans(cont)
            cont <- log10(cont+1)
            case <- log10(case+1)
            res <- results(d2s, contrast=c("Unique", start+h, start+j), independentFiltering=F, cooksCutoff=Inf)
            ofname = paste0("genotype_difference_", time, "_", annot_df$Genotype[1+h], "_", annot_df$Genotype[1+j], ".tsv")
            write.table(res, ofname, sep="\t")
            a <- data.frame(cbind(case=case, cont=cont[names(case)]))
            ofname = paste0("scatter_genotype_difference_", time, "_", annot_df$Genotype[1+h], "_", annot_df$Genotype[1+j], ".png")
            png(ofname)
            g <- ggplot(a, aes(x=cont, y=case))+geom_point()
            g <- g+xlab("Cont log10(mean(x)+1)")+ylab("Case log10(mean(x)+1)")
            plot(g)
            dev.off()
        }
    }
}

