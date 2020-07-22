# TODO: <simpleError: object 'sample_id' not found> : no individual -> variant mapping found?
wd = "~/cage/analysis/"
setwd(wd)

library("GenomicRanges")
library("rtracklayer")
library("wiggleplotr")
library("cowplot")
library("data.table")
library("GenomicFeatures")

library("dplyr")
library("ggplot2")
library("readr")
library("stringr")
library("reshape2")
library("tidyr")
library("optparse")

set.seed(123)


option_list <- list(
  make_option(c("--kind"), type="character", default="cage_over_tx",
              help="What kind of visualization to make - cage_over_tx, tx_over_cage or ann_over_cage", metavar = "path")
)
opt <- parse_args(OptionParser(option_list=option_list))


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

summ_ = read_tsv("../pipeline/txrevise/results/featureCounts/merged_gene_counts.txt")
summ_ = colSums(summ_[,-1])
summ_ = data_frame(sample_id=labels(summ_), N=summ_) %>%
  filter(sample_id %in% sample_mapping_$sample_id) # let's only look at european individuals

summ = rbind(summ, summ_)

sample_data = data_frame(sample_id=summ$sample_id, scaling_factor=summ$N/1000000) %>%
  mutate(bigWig = paste0("~/cage/align/results/", sample_id, ".bw"))


upstream1 = GenomicFeatures::makeTxDbFromGFF("txrevise/scripts/processed/Homo_sapiens.GRCh38.96/merged/txrevise.grp_1.upstream.gff3")
upstream2 = GenomicFeatures::makeTxDbFromGFF("txrevise/scripts/processed/Homo_sapiens.GRCh38.96/merged/txrevise.grp_2.upstream.gff3")
exons_list1 = GenomicFeatures::exonsBy(upstream1, by = "tx", use.names = TRUE)
exons_list2 = GenomicFeatures::exonsBy(upstream2, by = "tx", use.names = TRUE)


to_visualize = read_tsv(paste0(opt$kind, ".tsv"))
promoterMetadata = read_tsv("../qtlmap_prep/FANTOM5_promoter_annotations.tsv", col_types="ccciiicii") %>%
  filter(gene_name %in% to_visualize$gene)


start_time = Sys.time()
for (n in c(1: (dim(to_visualize)[1]) )) { # c(1: (dim(to_visualize)[1]) )
temp = tryCatch({

  observed_gene = to_visualize$gene[n]
  observed_variant = to_visualize$top_variant[n]


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

  rangeslists_wide = c()
  for (i in 1:dim(spec_wide)[1]) {
    row = spec_wide[i,]
    rangeslists_wide = c(rangeslists_wide,
                         GRanges(seqnames=2,
                                 ranges=IRanges(as.numeric(spec_wide[i,4]), as.numeric(spec_wide[i,5])),
                                 strand=as.character(spec_wide[i,7])))
  }
  names(rangeslists_wide) = spec_wide$tss_id


  annot1 = exons_list1[names(exons_list1) %like% paste0(observed_gene, ".grp_1.upstream")]
  annot2 = exons_list2[names(exons_list2) %like% paste0(observed_gene, ".grp_2.upstream")]


  filter_start = min(spec$peak_start) - 1000
  filter_end = max(spec$peak_end) + 1000
  region_filter = GRanges(seqnames=2, IRanges(start=filter_start, end=filter_end))

  filtered_exons_g1 = lapply(annot1, pintersect, region_filter, drop.nohit.ranges=TRUE)
  filtered_exons_g1 = filtered_exons_g1[lapply(filtered_exons_g1, length) > 0]
  filtered_exons_g2 = lapply(annot2, pintersect, region_filter, drop.nohit.ranges=TRUE)
  filtered_exons_g2 = filtered_exons_g2[lapply(filtered_exons_g2, length) > 0]

  # coding regions must be subset of exons i think
  # g2 should be dark blue, g1 light blue
  wide_prom_and_g1 = c(filtered_exons_g2, filtered_exons_g1, rangeslists_wide)
  prom_and_g2 = c(filtered_exons_g2, rangeslists)


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

  a = plotCoverage(exons=wide_prom_and_g1, cdss=prom_and_g2, track_data=sample_data_this_cage,
                   coverage_type="line", fill_palette=fill_palette, transcript_label=FALSE,
                   plot_fraction=0.2, return_subplots_list=TRUE, rescale_introns=TRUE,
                   new_intron_length=100, heights=c(0.5, 0.5))

  a = plotCoverage(exons=wide_prom_and_g1, cdss=prom_and_g2, track_data=sample_data_this_geuvadis,
                   coverage_type="line", fill_palette=fill_palette, transcript_label=FALSE,
                   plot_fraction=0.2, return_subplots_list=TRUE, rescale_introns=TRUE,
                   new_intron_length=100, heights=c(0.5, 0.5))

  combo = plot_grid(a[[1]], b[[1]], b[[2]], ncol=1, nrow=3, align="v")

  ggsave(paste0("plots/", opt$kind, "/", n, "_", observed_gene, ".pdf"), combo)

}, error = function(e){print(e)})
}
Sys.time() - start_time
