require(ggplot2)
require(cowplot)
require(ggsci)
require(scales)
library(gglorenz)
library(DescTools)
library(ineq)
library(lawstat)

files <- c("Rad_rep1_pe.txt", "Rad_rep2_pe.txt", "Swi6_rep1_pe.txt")
# total_counts <- c(4196755+250020, 4202122+242632, 473120+191002)*2 # input reads, STAR output
total_counts <- c(9395176, 9603578, 2003838) # obtained from samtools output
all_data <- NULL
for (i in 1:3) {
    all_data <- rbind(all_data, cbind(rep(i, dim(mat)[1]), mat[,4]/(total_counts[i])))
}
colnames(all_data) <- c('Dataset', 'Coverage')
all_data <- as.data.frame(all_data)
all_data[,1] <- factor(all_data[,1], labels=c('Rad_rep1', 'Rad_rep2', 'Swi6'))
g <- ggplot(all_data, aes(x=Coverage, group=Dataset, color=Dataset))+geom_density()+theme_minimal_grid(12)+scale_color_npg()+xlim(0, 0.001)
pdf(paste0('cov_rep_density_norm_single.pdf'))
plot(g)
dev.off()
g <- ggplot(all_data, aes(x=Coverage, group=Dataset, color=Dataset))+geom_density()+theme_minimal_grid(12)+scale_color_npg()+xlim(0, 0.00025)
pdf(paste0('cov_rep_density_zoom_norm_single.pdf'))
plot(g)
dev.off()


files <- c("Rad_rep1_pe_uniq.txt", "Rad_rep2_pe_uniq.txt", "Swi6_rep1_pe_uniq.txt")
genot <- c('Rad_rep1', 'Rad_rep2', 'Swi6')
all_data <- NULL
for (i in 1:3) {
    mat <- read.table(files[i], header=F, sep="\t", stringsAsFactor=F)
    print(head(mat))
    print(Gini(mat[,4]))
    print(max(mat[,4]))
    print(mean(mat[,4]))
    pdf(paste0('lc_', i, '.pdf'))
    plot(Lc(mat[,4]))
    dev.off()
    all_data <- rbind(all_data, cbind(rep(genot[i], dim(mat)[1]), mat[,4]))
}


colnames(all_data) <- c('genot', 'coverage')
all_data <- as.data.frame(all_data)
all_data[,1] <- factor(all_data[,1], levels=c('Swi6', 'Rad_rep1', 'Rad_rep2'))
all_data[,2] <- as.numeric(as.character(all_data[,2]))
g <- ggplot(all_data, aes(x=coverage, color=genot))+geom_histogram(alpha=0.5)+scale_x_continuous(trans = log10_trans())
png('test_hist.png')
plot(g)
dev.off()

g <- ggplot(all_data, aes(x=coverage, color=genot, linetype=genot)) +
    stat_lorenz(desc=TRUE, lwd=1.6)+coord_fixed()+geom_abline(linetype = "dashed")+
    hrbrthemes::scale_x_percent() +
    hrbrthemes::scale_y_percent() +
    theme(legend.title = element_blank())+theme_cowplot(25)+
    scale_color_manual(values=c(do.call(rgb, as.list(c(248,118,108)/255)), do.call(rgb, as.list(c(81,182,76)/255)), do.call(rgb, as.list(c(40,70,255)/255))))+
    scale_linetype_manual(values=c(1,1,4))
pdf('lorenz_curve_asc.pdf')
plot(g)
dev.off()
g <- ggplot(all_data, aes(x=coverage, color=genot, linetype=genot)) +
    stat_lorenz(desc=FALSE, lwd=1.6)+coord_fixed()+geom_abline(linetype = "dashed")+
    hrbrthemes::scale_x_percent() +
    hrbrthemes::scale_y_percent() +
    theme(legend.title = element_blank())+theme_cowplot(25)+
    scale_color_manual(values=c(do.call(rgb, as.list(c(248,118,108)/255)), do.call(rgb, as.list(c(81,182,76)/255)), do.call(rgb, as.list(c(40,70,255)/255))))+
    scale_linetype_manual(values=c(1,1,4))
pdf('lorenz_curve.pdf')
plot(g)
dev.off()

for (i in 1:3) {
    tdata <- all_data[all_data[,1] == genot[i],]
    print(Gini(tdata[,2]))
    print(max(tdata[,2]))
    print(mean(tdata[,2]))
    g <- ggplot(tdata, aes(x=coverage, color=genot)) +
        stat_lorenz(desc=TRUE, alpha=0.35)+coord_fixed()+geom_abline(linetype = "dashed")+
        hrbrthemes::scale_x_percent() +
        hrbrthemes::scale_y_percent() +
        theme(legend.title = element_blank())
    pdf(paste0('lorenz_curve_', i, '.pdf'))
    plot(g)
    dev.off()
}