# Pipeline used for yeast coexp analysis

## RNA-seq analysis
* rna_analysis/
  * read processing
    * alignment.sh
  * summarize read count data
    * merge_table.R
    * plot_raw_count_data.R
  * DE-seq2 analysis
    * plot_de_result.R
    * merge_de_result.R
    * sort_and_join.sh
    * (go_enrichment.R)

## ChIP-seq analysis
* chip_analysis/
  * read processing
    * alignment.sh
  * peak calling
    * peak_calling.sh
    * extract_wt_peaks.sh
    * merge_all_peaks_dist.sh
  * coverage comparison
    * normalize.py 
    * enrich_heatmap_ann.R
  * promoter coverage comparison
    * draw_lorenz_curve.R
    * (lorenz_curve.R)
  * diffbind analysis
    * raw_corr_diffbind.R
    * diffbind_comparison_tf.R

## Co-expression cluster analysis
* coexp_analysis/
  * coexp_deg.R
    * extract the coexp clusters by dynamicTreeCut
    * input: coexp matrix from CoCoCoNet and output of DESeq2
  * coexp_enrichment.R
    * apply an enrichment analysis for each coexp cluster
    * input: output of coexp_deg.R and output of DESeq2

## Reference
* [Prevalent and Dynamic Binding of the Cell Cycle Checkpoint Kinase Rad53 to Gene Promoters.](https://www.biorxiv.org/content/10.1101/2021.05.23.445333v3) 
  * Yi-Jun Sheu, Risa Karakida Kawaguchi, Jesse Gillis, Bruce Stillman.
  * bioRxiv 2021.05.23.445333; doi: https://doi.org/10.1101/2021.05.23.445333