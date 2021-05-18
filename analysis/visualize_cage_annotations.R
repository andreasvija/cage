library("data.table") # %like%
library("rtracklayer")
library("GenomicRanges")
library("GenomicFeatures")

library("dplyr")
library("readr")
library("optparse")

library("wiggleplotr")

setwd("cage/analysis")

promoter_annots = read_tsv("../qtlmap_prep/FANTOM5_promoter_annotations.tsv", col_types="ccciiicii")
new_transcripts = readRDS("txrevise/scripts/processed/Homo_sapiens.GRCh38.96_CAGE-20/new_transcripts_20.rds")

grp1 = "txrevise/scripts/processed/Homo_sapiens.GRCh38.96_regular/txrevise_regular.grp_1.upstream.gff3"
grp2 = "txrevise/scripts/processed/Homo_sapiens.GRCh38.96_regular/txrevise_regular.grp_2.upstream.gff3"
exons_list1 = GenomicFeatures::exonsBy(GenomicFeatures::makeTxDbFromGFF(grp1), by = "tx", use.names = TRUE)
exons_list2 = GenomicFeatures::exonsBy(GenomicFeatures::makeTxDbFromGFF(grp2), by = "tx", use.names = TRUE)

all_genes = unique(promoter_annots$gene_name)
ann_over_cage_genes = unique(read_tsv("ann_over_cage.tsv")$gene)

#for every gene
for (gene in all_genes) {

  promoters = promoter_annots %>%
    filter(gene_name == gene)
  chromosome = promoter_annots$chr[1]

  annot1 = exons_list1[names(exons_list1) %like% paste0(gene, ".grp_1.upstream")]
  annot2 = exons_list2[names(exons_list2) %like% paste0(gene, ".grp_2.upstream")]
  exons = c(annot1, annot2)

  if (length(exons) == 0) {next}

  gene_new_transcripts = new_transcripts[names(new_transcripts) %like% paste0(gene, ".new.upstream")]
  if (length(gene_new_transcripts) == 0) {next}


  spec = promoters

  spec_wide = spec %>%
    mutate(peak_start = peak_start - 150, peak_end = peak_end + 150)

  rangeslists = c()
  for (i in 1:dim(spec)[1]) {
    row = spec[i,]
    rangeslists = c(rangeslists,
                    GRanges(seqnames=chromosome,
                            ranges=IRanges(as.numeric(spec[i,4]), as.numeric(spec[i,5])),
                            strand=as.character(spec[i,7])))
  }
  names(rangeslists) = spec$tss_id

  rangeslists_wide = c()
  for (i in 1:dim(spec_wide)[1]) {
    row = spec_wide[i,]
    rangeslists_wide = c(rangeslists_wide,
                         GRanges(seqnames=chromosome,
                                 ranges=IRanges(as.numeric(spec_wide[i,4]), as.numeric(spec_wide[i,5])),
                                 strand=as.character(spec_wide[i,7])))
  }
  names(rangeslists_wide) = spec_wide$tss_id


  filter_start = min(spec$peak_start) - 100000
  filter_end = max(spec$peak_end) + 100000
  region_filter = GRanges(seqnames=chromosome, IRanges(start=filter_start, end=filter_end))

  filtered_exons = lapply(exons, pintersect, region_filter, drop.nohit.ranges=TRUE)
  filtered_exons = filtered_exons[lapply(filtered_exons, length) > 0]
  #filtered_exons = pintersect(GRanges(seqnames=chromosome, transcript), region_filter, drop.nohit.ranges=TRUE)
  if (length(filtered_exons) == 0) {next}

  light_blue = c(filtered_exons, rangeslists_wide, as.list(gene_new_transcripts))
  dark_blue = c(filtered_exons, rangeslists)

  print(plotTranscripts(exons=light_blue, cdss=dark_blue))
}
