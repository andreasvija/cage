wd = "~/cage/analysis/"
setwd(wd)

library("GenomicRanges")
library("rtracklayer")
library("wiggleplotr")
library("cowplot")
library("data.table")
library("GenomicFeatures")

library("dplyr")
library("readr")
library("stringr")
library("ggplot2")

set.seed(123)


sample_mapping = read_tsv("../qtlmap_prep/sampleMetadata.tsv") %>%
  select(sample_id, genotype_id)
sample_mapping_ = read_tsv("GEUVADIS_EUR.tsv") %>%
  select(sample_id, genotype_id)
overlap = intersect(sample_mapping$genotype_id, sample_mapping_$genotype_id) # n = 78


genotypes = read_tsv("variantinfo.vcf")[c(3,10:163)] %>%
  reshape2::melt(id.vars=c("ID"))
names(genotypes) = c("variant", "genotype_id", "alleles")
genotypes = genotypes %>%
  merge(sample_mapping) %>%
  select(-genotype_id) %>%
  mutate(track_id="CAGE")

genotypes_ = read_tsv("variantinfo_geuvadis.vcf")[c(3,10:454)] %>%
  reshape2::melt(id.vars=c("ID"))
names(genotypes_) = c("variant", "genotype_id", "alleles")
genotypes_ = genotypes_ %>%
  filter(genotype_id %in% sample_mapping_$genotype_id) %>%
  merge(sample_mapping_) %>%
  select(-genotype_id) %>%
  mutate(track_id="RNA-seq")

genotypes = rbind(genotypes, genotypes_)
mapping = list("0", "1", "1", "2")
names(mapping) = c("0|0", "0|1", "1|0", "1|1")
genotypes$alleles = mapping[genotypes$alleles]
genotypes$alleles = unlist(genotypes$alleles)


summ = read_tsv("../qtlmap_prep/bwa.counts.summary", col_types = paste0(c("c", rep("n", 154)), collapse=""))
renamer <- function(name) {
  name = str_replace(name, "../align/results/", "")
  name = str_replace(name, ".sorted.bam", "")
  return(name)
}

colnames(summ) = unname(sapply(colnames(summ), renamer))
summ = t(summ[1,2:length(summ)])
summ = data_frame(sample_id=rownames(summ), N=summ)

summ_ = read_tsv("txrev-25_merged_gene_counts.txt")
summ_ = colSums(summ_[,-1])
summ_ = data_frame(sample_id=labels(summ_), N=summ_) %>%
  filter(sample_id %in% sample_mapping_$sample_id) # let's only look at european individuals

summ = rbind(summ, summ_)

sample_data = data_frame(sample_id=summ$sample_id, scaling_factor=summ$N/1000000) %>%
  mutate(bigWig = paste0("~/cage/align/results/", sample_id, ".bw"))

upstream1 = GenomicFeatures::makeTxDbFromGFF(
  "txrevise/scripts/processed/Homo_sapiens.GRCh38.96_regular/txrevise_regular.grp_1.upstream.gff3")
exons_list1 = GenomicFeatures::exonsBy(upstream1, by = "tx", use.names = TRUE)


observed_gene = "ENSG00000151694"
observed_variant = "chr2_9555777_A_G"

promoterMetadata = read_tsv("../qtlmap_prep/FANTOM5_promoter_annotations.tsv", col_types="ccciiicii") %>%
  filter(gene_name == observed_gene)

new_transcripts = readRDS("txrevise/scripts/processed/Homo_sapiens.GRCh38.96_CAGE-20/new_transcripts_20.rds")
gene_new_transcripts = new_transcripts[names(new_transcripts) %like% paste0(observed_gene, ".new.upstream")]


spec = promoterMetadata %>%
  filter(gene_name==observed_gene)

spec_wide = spec %>%
  mutate(peak_start = peak_start - 150, peak_end = peak_end + 150)


rangeslists = c()
for (i in 1:dim(spec)[1]) {
  row = spec[i,]
  rangeslists = c(rangeslists,
                  GRanges(seqnames=2,
                          ranges=IRanges(as.numeric(spec[i,4]), as.numeric(spec[i,5])),
                          strand=as.character(spec[i,7])))
}
names(rangeslists) = spec$tss_id


annot1 = exons_list1[names(exons_list1) %like% paste0(observed_gene, ".grp_1.upstream")]

filter_start = min(spec$peak_start) - 1000
filter_end = max(spec$peak_end) + 1000
region_filter = GRanges(seqnames=2, IRanges(start=filter_start, end=filter_end))

filtered_exons = lapply(annot1, pintersect, region_filter, drop.nohit.ranges=TRUE)
filtered_exons = filtered_exons[lapply(filtered_exons, length) > 0]

filtered_new_transcripts = lapply(gene_new_transcripts, pintersect, region_filter, drop.nohit.ranges=TRUE)
filtered_new_transcripts = filtered_new_transcripts[lapply(filtered_new_transcripts, length) > 0]


light_blue = c(filtered_exons, rangeslists, as.list(filtered_new_transcripts))
dark_blue = c(filtered_exons, rangeslists)


genotypes_sub = genotypes %>%
  filter(variant == observed_variant) %>%
  select(sample_id, alleles, track_id)

sample_data_this = sample_data %>%
  merge(genotypes_sub) %>%
  mutate(colour_group=alleles) %>%
  select(-alleles)

sample_data_this_cage = sample_data_this %>%
  filter(track_id=="CAGE")
sample_data_this_geuvadis = sample_data_this %>%
  filter(track_id=="RNA-seq")

fill_palette = c("#225ea8", "#884072", "#ef233c")

a = plotCoverage(exons=light_blue, cdss=dark_blue, track_data=sample_data_this_cage,
                 coverage_type="line", fill_palette=fill_palette,
                 plot_fraction=1, rescale_introns=FALSE, heights=c(0.3, 0.7))

b = plotCoverage(exons=light_blue, cdss=dark_blue, track_data=sample_data_this_geuvadis,
                 coverage_type="line", fill_palette=fill_palette,
                 plot_fraction=1, rescale_introns=FALSE, heights=c(0.3, 0.7))

ggsave(paste0("plots_one/cage.pdf"), a)
ggsave(paste0("plots_one/rna-seq.pdf"), b)
