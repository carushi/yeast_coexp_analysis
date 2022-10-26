require(GGally)

ANN_DIR="../../../ann/"

tss <- read.table('transcript_promoter_sorted.bed', sep="\t", header=F, quote='\"', stringsAsFactors=F)
colnames(tss) <- c('chr', 'start', 'end', '.', 'strand', 'gene_id', 'gene_name')
tss <- tss[,-c(4)]

dir = '../data_cp/' # diffbind results for CP  dataset
wt_rep1 = read.table(file.path(dir, 'dba_result_1_2.tsv'), header=T, sep="\t", stringsAsFactor=F)
raw_rep1 = read.table(file.path(dir, 'chip_count_integrated.tsv'), header=T, sep="\t", stringsAsFactor=F)
signal_gene = wt_rep1[which(wt_rep1[,'FDR'] < 0.05),]
dds <- read.table(file.path(ANN_DIR, 'design_matrix_cp.tsv'), sep="\t", header=T)
annot_df = dds[,c(2,3)]
rownames(annot_df) = paste0('X', dds[,1])
annot_df = annot_df[annot_df[,'Genotype'] == 'WT',]
mat = raw_rep1[,rownames(annot_df)]
colnames(mat) = paste0(annot_df[,'Genotype'], '_', annot_df[,'Stage'], '_', as.integer((1:dim(mat)[2]-1)/3+1))
pdf(paste0('cor_', 'tf', '.pdf'), width=15, height=18)
ggpairs(mat, titel="correlation")
dev.off()
rmat = mat[,1:3]+mat[,4:6]
pdf(paste0('cor_', 'tf', '_rem.pdf'), width=15, height=18)
ggpairs(rmat, titel="correlation")
dev.off()
colnames(wt_rep1)[1] = 'chr'
mat <- mat[,-c(4:5)]
# colnames(mat)[4:dim(mat)[2]] = paste0(prefix, '_', colnames(mat)[4:dim(mat)[2]])
m <- merge(tss, wt_rep1, by=c('chr', 'start', 'end'), all=FALSE)
m[,'FDR'] <= 0.05
        


dir = '../data_tf/' # diffbind results for TF dataset
wt_rep1 = read.table(file.path(dir, 'dba_result_1_2.tsv'), header=T, sep="\t", stringsAsFactor=F)
raw_rep1 = read.table(file.path(dir, 'chip_count_integrated.tsv'), header=T, sep="\t", stringsAsFactor=F)

dds <- read.table(file.path(ANN_DIR, 'design_matrix_tf.tsv'), sep="\t", header=T)
annot_df = dds[,c(2,3)]
rownames(annot_df) = paste0('X', dds[,1])
annot_df = annot_df[annot_df[,'Genotype'] == 'WT',]
mat = raw_rep1[,rownames(annot_df)]
colnames(mat) = paste0(annot_df[,'Genotype'], '_', annot_df[,'Stage'], '_',  as.integer((1:dim(mat)[2]-1)/3+1))
pdf(paste0('cor_', 'cp', '.pdf'), width=15, height=18)
ggpairs(mat, titel="correlation")
dev.off()
rmat = mat[,1:3]+mat[,4:6]
pdf(paste0('cor_', 'cp', '_rem.pdf'), width=15, height=18)
ggpairs(rmat, titel="correlation")
dev.off()

