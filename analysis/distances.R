wd = "~/cage/analysis/"
setwd(wd)

library("data.table") # %like%
library("rtracklayer")
library("GenomicRanges")
library("GenomicFeatures")

library("dplyr")
library("readr")

grp1_file = "txrevise/scripts/processed/Homo_sapiens.GRCh38.96_regular/txrevise_regular.grp_1.upstream.gff3"
grp2_file = "txrevise/scripts/processed/Homo_sapiens.GRCh38.96_regular/txrevise_regular.grp_2.upstream.gff3"
promoters_file = "../qtlmap_prep/FANTOM5_promoter_annotations.tsv"

promoter_annots = read_tsv(promoters_file, col_types="ccciiicii")

upstream1 = GenomicFeatures::makeTxDbFromGFF(grp1_file)
upstream2 = GenomicFeatures::makeTxDbFromGFF(grp2_file)
exons_list1 = GenomicFeatures::exonsBy(upstream1, by = "tx", use.names = TRUE)
exons_list2 = GenomicFeatures::exonsBy(upstream2, by = "tx", use.names = TRUE)

genes_col = c()
promoters_col = c()
distances_col = c() # distances from each promoter to its closest transcript, -1 if no transcript

#for every gene

start_time = Sys.time()
for (gene in unique(promoter_annots$gene_name)) {

  promoters = promoter_annots %>%
    filter(gene_name == gene)

  annot1 = exons_list1[names(exons_list1) %like% paste0(gene, ".grp_1.upstream")]
  annot2 = exons_list2[names(exons_list2) %like% paste0(gene, ".grp_2.upstream")]
  exons = c(annot1, annot2)

  #for every promoter belonging to the gene
  for (i in 1 : dim(promoters)[1]) {

    peak_start = promoters$peak_start[i]
    peak_end = promoters$peak_end[i]
    strand = promoters$strand[i]
    chr = promoters$chr[i]
    tss_id = promoters$tss_id[i]

    genes_col = c(genes_col, gene)
    promoters_col = c(promoters_col, tss_id)

    promoter = GRanges(seqnames=chr, IRanges(start=peak_start, end=peak_end))

    # if no exons, distance is -1
    if (length(exons) < 1) {
      distances_col = c(distances_col, -1)
      next
    }

    #for every exon in transcripts

    #if promoter overlaps with an exon, distance is 0
    #otherwise, for every transcript belonging to the gene
    #record the smallest distance from either promoter end to any exon's strand-dependent start
    overlaps_with_exon = FALSE
    promoter_distances = c()
    for (j in 1 : length(exons@unlistData@ranges)) {

      exon_start = exons@unlistData@ranges@start[j]
      exon_end = exon_start + exons@unlistData@ranges@width[j]

      if ((peak_start <= exon_start & peak_end >= exon_end) | # peak surrounds transcript
          (peak_start >= exon_start & peak_start <= exon_end) | # peak start overlaps with transcript
          (peak_end >= exon_start & peak_end <= exon_end)) # peak end overlaps with transcript
      {
        overlaps_with_exon = TRUE
        break
      }

      exon_real_start = exon_start
      if (strand == "+" & peak_start < exon_start) {
        promoter_distances = c(promoter_distances, exon_start - peak_end)
      }
      if (strand == "-" & peak_end > exon_end) {
        promoter_distances = c(promoter_distances, peak_start - exon_end)
      }

    }
    if (overlaps_with_exon) {
      distances_col = c(distances_col, 0)
    }
    else if (length(promoter_distances) == 0) {
      distances_col = c(distances_col, -1)
    }
    else {
      distances_col = c(distances_col, min(promoter_distances))
    }
  }
}

Sys.time() - start_time

gathered = tibble(gene=genes_col, promoter=promoters_col, distance=distances_col)

length(gathered$distance) # 93 554 #compare to n of transcripts
sum(gathered$distance == -1) / length(gathered$distance) # 10.4% # no transcript at all or no transcript upstream
sum(gathered$distance == 0) / length(gathered$distance) # 63.5%
sum(gathered$distance > 0) / length(gathered$distance) # 26.1%
to_plot = gathered[gathered$distance > 0,] # 24 393

sum(gathered$distance > 1000) / length(gathered$distance) # 5.0%
sum(gathered$distance > 1000) / sum(gathered$distance >= 0) # 5.6%

library("ggplot2")
ggplot(to_plot) + geom_histogram(aes(x=distance)) + scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000))
