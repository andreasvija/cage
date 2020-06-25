library("dplyr")
library("readr")
library("tidyr")
library("gprofiler2")

setwd("~/cage/qtlmap_prep")
set.seed(123)



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
# 99 953 annotations, meaning 110 297 promoters had no gene name
promoterAnnots = promoterAnnots %>% select(complex, gene_symbol)



# better format
# tss_id, gene_name, chr(without "chr"), peak_start, peak_end, strand
# hg19::chr1:564571..564600,+;hg_1.1
# XXXX::chr[chr]:[peak_start]..[peak_end],[strand];[tss_id] and gene_symbol is gene_name

fun_tss_id <- function(s) {
  return(strsplit(s, ";")[[1]][2])
}
promoterAnnots$tss_id = sapply(promoterAnnots$complex, FUN=fun_tss_id)
promoterAnnots$gene_name = promoterAnnots$gene_symbol
fun_chr <- function(s) {
  return(strsplit(strsplit(s, "::chr")[[1]][2], ":")[[1]][1])
}
promoterAnnots$chr = sapply(promoterAnnots$complex, FUN=fun_chr)

fun_peak_start <- function(s) {
  return(strsplit(strsplit(s, "\\.\\.")[[1]][1], "chr[0-9MXY]+:")[[1]][2])
}
promoterAnnots$peak_start = sapply(promoterAnnots$complex, FUN=fun_peak_start)
fun_peak_end <- function(s) {
  return(strsplit(strsplit(s, "\\.\\.")[[1]][2], ",")[[1]][1])
}
promoterAnnots$peak_end = sapply(promoterAnnots$complex, FUN=fun_peak_end)
sum(promoterAnnots$peak_start > promoterAnnots$peak_end) # 0

fun_strand <- function(s) {
  return(strsplit(strsplit(s, ",")[[1]][2], ";")[[1]][1])
}
promoterAnnots$strand = sapply(promoterAnnots$complex, FUN=fun_strand)

promoterAnnots = promoterAnnots %>% select(-complex, -gene_symbol)
apply(is.na(promoterAnnots), 2, sum)

# remove genes not on chromosomes 1-22
promoterAnnots = promoterAnnots[!(promoterAnnots$chr %in% c("X", "Y", "M")),] # 96 631

# split gene names from spaces as different genes
promoterAnnots$gene_name = sapply(promoterAnnots$gene_name, strsplit, split=" ")
promoterAnnots = unnest(promoterAnnots, cols=c(gene_name)) # 97 422



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
geneMapping = unique(promoterAnnots$gene_name) %>%
  gconvert(mthreshold=1, filter_na=FALSE) %>%
  mutate(gene_name=input, gene_id=target) %>%
  select(gene_name, gene_id) # 21 564

geneMapping$gene_id[geneMapping$gene_id == "None"] = NA
geneMapping = geneMapping[!is.na(geneMapping$gene_id),] # 20 481
'

# make gene column ENSEMBL ID from source, discard row if impossible

joined = inner_join(promoterAnnots, geneMapping, by=c("gene_name")) # 95 023 / 95 042
fantom_genes = unique(promoterAnnots$gene_name) # 21 564
overlap_genes = unique(joined$gene_name[!is.na(joined$gene_id)]) # 20 481

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
promoterAnnots = inner_join(promoterAnnots, not_dupes, by="tss_id") # 93 718 / 93 706



# save new annotations
promoterAnnots = promoterAnnots %>%
  mutate(gene_name = gene_id) %>%
  select(-gene_id)
write.table(promoterAnnots, file="FANTOM5_promoter_annotations.tsv", sep="\t", row.names=FALSE, quote=FALSE)
