library("dplyr")
library("readr")
library("tidyr")

setwd("~/cage/qtlmap_prep")
set.seed(123)



# https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_annotation/hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt.gz (26.06.2020)
promoterAnnots = read_tsv("hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt", col_types="ccnccccccc")
colnames(promoterAnnots) = c("complex", "transcript", "distance", "gene_id", "hgnc_mgi_id",
                             "uniprot_id", "gene_name", "gene_symbol", "gene_synonyms", "gene_source")
# 210 250 annotations



# remove useless columns and rows

sum(is.na(promoterAnnots$gene_symbol) & !is.na(promoterAnnots$gene_id)) # 0
sum(is.na(promoterAnnots$gene_symbol) & !is.na(promoterAnnots$hgnc_mgi_id)) # 0

sum(is.na(promoterAnnots$gene_symbol) & !is.na(promoterAnnots$uniprot_id)) # 439
sum(is.na(promoterAnnots$uniprot_id) & !is.na(promoterAnnots$gene_symbol)) # 41 765

sum(is.na(promoterAnnots$gene_symbol) & !is.na(promoterAnnots$gene_name)) # 0
sum(is.na(promoterAnnots$gene_symbol) & !is.na(promoterAnnots$gene_synonyms)) # 0
sum(is.na(promoterAnnots$gene_symbol) & !is.na(promoterAnnots$gene_source)) # 0

# gene_symbol is obviously the identificator with the biggest coverage
promoterAnnots = promoterAnnots[!is.na(promoterAnnots$gene_symbol),]
# 99 953 annotations, over 22 591 genes meaning 110 297 promoters had no gene name
promoterAnnots = promoterAnnots %>%
  select(complex, gene_symbol) %>%
  rename(gene_name=gene_symbol)



# https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz (26.06.2020)
promoterCounts = read_tsv("hg38_fair+new_CAGE_peaks_phase1and2.bed", col_types="cnncncnnc")
colnames(promoterCounts) = c("chr", "peak_start", "peak_end", "complex",
                             "peak_score", "strand", "tss_start", "tss_end", "rgb")

promoterAnnots = promoterAnnots %>%
  inner_join(promoterCounts, by="complex") %>% # 99 880, meaning 73 annots had no counts in bed file
  separate(complex, c("old_tss_id", "tss_id"), sep = ";") %>%
  select(tss_id, gene_name, chr, peak_start, peak_end, peak_score, strand, tss_start, tss_end)

fun_chr <- function(s) {
  return(strsplit(s, "chr")[[1]][2])
}
promoterAnnots$chr = sapply(promoterAnnots$chr, FUN=fun_chr)



apply(is.na(promoterAnnots), 2, sum)

# remove genes not on chromosomes 1-22
promoterAnnots = promoterAnnots[promoterAnnots$chr %in% as.character(1:22),] # 96 562 over 21 696

# split gene names from spaces as different genes
promoterAnnots$gene_name = sapply(promoterAnnots$gene_name, strsplit, split=" ")
promoterAnnots = unnest(promoterAnnots, cols=c(gene_name)) # 97 320



# give genes ID-s instead of names

# mapping version 1 - file

geneMapping = read_tsv("transcript_usage_Ensembl_96_phenotype_metadata.tsv", col_types="ccccciicccci") %>%
  select(gene_id, gene_name) %>% # 207 749
  unique() # 58 434

# remove gene names that have multiple ID-s from gene mapping
uniqueGeneNames = geneMapping %>%
  group_by(gene_name) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  filter(n == 1)
uniqueGeneNames = uniqueGeneNames$gene_name

geneMapping = geneMapping[geneMapping$gene_name %in% uniqueGeneNames,] # 56 835


# mapping version 2 - gprofiler2 (not used - as good as from file, but more of a black box)
'
library("gprofiler2")

geneMapping = unique(promoterAnnots$gene_name) %>%
  gconvert(mthreshold=1, filter_na=FALSE) %>%
  mutate(gene_name=input, gene_id=target) %>%
  select(gene_name, gene_id) # 21 564

geneMapping$gene_id[geneMapping$gene_id == "None"] = NA
geneMapping = geneMapping[!is.na(geneMapping$gene_id),] # 20 481
'

# make gene column ENSEMBL ID from source, discard row if impossible

joined = inner_join(promoterAnnots, geneMapping, by=c("gene_name")) # 94 928 / 94 946
fantom_genes = unique(promoterAnnots$gene_name) # 21 535
overlap_genes = unique(joined$gene_name[!is.na(joined$gene_id)]) # 20 456 / 20 454

1 - length(overlap_genes)/length(fantom_genes) # 5% gene names not mapped (or not uniquely)
missing = fantom_genes[!(fantom_genes %in% overlap_genes)] # most LOCXXXXX

promoterAnnots = joined

# remove promotors mapping to multiple (non-na) gene id-s

dim(table(promoterAnnots$tss_id)[table(promoterAnnots$tss_id) != 1]) / length(unique(promoterAnnots$tss_id)) # <0.7% removed
not_dupes = promoterAnnots %>%
  select(tss_id) %>%
  group_by(tss_id) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  filter(n==1) %>%
  select(tss_id)
promoterAnnots = inner_join(promoterAnnots, not_dupes, by="tss_id") # 93 663 over 20 201 / 93 651

# remove promotors whose genes have promotors on multiple chromosomes

combinations = unique(promoterAnnots[c("gene_id", "chr")])
multi_chr_genes = names(table(combinations$gene_id)[table(combinations$gene_id) > 1])
promoterAnnots = promoterAnnots[!(promoterAnnots$gene_id %in% multi_chr_genes),] # 93 554 over 20 193



# save new annotations
promoterAnnots = promoterAnnots %>%
  mutate(gene_name = gene_id) %>%
  select(-gene_id)
write.table(promoterAnnots, file="FANTOM5_promoter_annotations.tsv", sep="\t", row.names=FALSE, quote=FALSE)
