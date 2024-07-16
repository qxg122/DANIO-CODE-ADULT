library(readr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(pbapply)
library(ggplot2)
library(ggsci)
library(stringr)
library(GenomicFeatures)
library(ChIPseeker)
library(clusterProfiler)
library(LOLA)
library(UpSetR)
library(networkD3)
library(rtracklayer)
library(genomation)
library(RColorBrewer)
library(BSgenome.Drerio.UCSC.danRer11)


# Figure 1
# DANIO-CODE datasets are obtained from https://github.com/DANIO-CODE/DANIO-CODE_Data_analysis/blob/master/Figures/Figure1/DCC-annotations-2021-05.csv
# This adult zebrafish track hub can be accessed through the URL https://data.cyverse.org/dav-anon/iplant/home/guanguan/Adult_zebrafish_hub/hub.txt in the UCSC Genome Browser 


# Figure 2
# Read the ChromHMM results
setwd("./ChromHMM_file/")
tissue_names <- c("blood", "brain", "colon", "heart", "intestine", "kidney", "liver", "muscle", "skin", "spleen", "testis")
dir_paths <- paste0("./ChromHMM_file/", tissue_names, '/ChromHMM_out')
emission_data <- list()
for (dir_path in dir_paths) {
  file_path <- file.path(dir_path, "emissions_10.txt")
  tissue_name <- basename(dirname(dir_path))
  emission_df <- read.table(file_path, header = TRUE)
  emission_df$tissue <- tissue_name
  emission_data[[dir_path]] <- emission_df
}

# Annotate the states according to the emission of different histone modification marks
State <- c('1_TssA', '2_TssFlank1', '3_TssFlank2', '4_EnhA', '5_EnhFlank', '6_Openchrom', '7_Bivalent', '8_cHC', '9_fHC','10_Quies')
State_indice <- list(
  blood = c(5, 4, 1, 2, 6, 10, 7, 9, 8, 3),
  brain = c(8, 7, 10, 9, 5, 4, 6, 2, 1, 3),
  colon = c(8, 7, 9, 10, 6, 4, 5, 2, 1, 3),
  heart = c(7, 8, 9, 10, 6, 4, 5, 1, 2, 3),
  intestine = c(6, 4, 10, 5, 2, 1, 3, 9, 7, 8),
  kidney = c(8, 7, 9, 10, 5, 4, 6, 2, 1, 3),
  liver = c(8, 7, 10, 9, 5, 4, 6, 2, 1, 3),
  muscle = c(8, 7, 9, 10, 6, 2, 1, 3, 4, 5),
  skin = c(6, 4, 1, 3, 2, 5, 10, 9, 7, 8),
  spleen = c(7, 9, 8, 10, 5, 4, 6, 2, 1, 3),
  testis = c(8, 7, 9, 10, 4, 6, 2, 1, 3, 5)
)
emission_data <- lapply(seq_along(emission_data), function(x) {
  emission_data[[x]]$State <- State[State_indice[[x]]]
  emission_data[[x]]
})
names(emission_data) <- tissue_names
# for (i in seq_along(emission_data)) {
#   emission_data[[i]]$State <- State[State_indice[[i]]]
# }

# Plot the emission of marks with annotation for all adult tissues
pivoted_data <- list()
for (i in seq_along(emission_data)) {
  emi <- emission_data[[i]]
  tissue_name <- unique(emi$tissue)[1]
  emi_pivoted <- emi %>%
    pivot_longer(!c(State,tissue,State), names_to = "Mark", values_to = "Emission")
  pivoted_data[[i]] <- emi_pivoted
}
all_pivoted_data <- do.call(rbind, pivoted_data)
all_pivoted_data$State <- factor(all_pivoted_data$State,
                                 levels = c('1_TssA', '2_TssFlank1','3_TssFlank2', '4_EnhA','5_EnhFlank', '6_Openchrom','7_Bivalent', '8_cHC','9_fHC', '10_Quies'),
                                 ordered = TRUE)
ggplot(all_pivoted_data, aes(State, Mark, fill = Emission)) + geom_tile() +
  theme_classic() + coord_flip() + facet_grid(~ tissue) +
  theme(axis.text.x = element_text(hjust = 1, angle = 90, face = 'bold'),
        axis.text.y = element_text(face = 'bold'),
        axis.title = element_text(face = 'bold'),
        axis.title.x = element_text(face = 'bold'),
        axis.title.y = element_text(face = 'bold'),
        aspect.ratio = 2,
        strip.text.x = element_text(face = 'bold', size = 14)
  ) +
  scale_fill_distiller(palette = "Blues", direction = 1)

# Assign the studied states from ChromHMM to open chromatin regions
atacPeakFiles <- list.files("./Adult_Zebrafish_CREs_Annotation", pattern = ".+.atac.rep1_vs_rep2.idr0.05.narrowPeak", full.names = TRUE)
names(atacPeakFiles) <- tissue_names
importUniqueNarrowPeak <- function(file){
  tmp <- import(file)
  tmp <- tmp[order(tmp$pValue, decreasing = TRUE)]
  unique(tmp)
}
allPeaks <- lapply(atacPeakFiles, importUniqueNarrowPeak)
allPeaks <- as(allPeaks, Class = "GRangesList")
chromHMMFiles <- list.files(path = paste0("./ChromHMM_file/", tissue_names, "/ChromHMM_out"), 
                            pattern = ".+_10_dense\\.bed", full.names = T)
names(chromHMMFiles) <- tissue_names
segs <- lapply(chromHMMFiles, rtracklayer::import)
chromHMMSegmentsList <- lapply(segs, function(x) split(x, x$name))
fWrap <- function(peaks, segments){
  function(x, y){
    x[peaks %over% segments[[y]] ] <- y
    x
  }
}
segsStage <- lapply(names(chromHMMSegmentsList), function(z) {
  purrr::reduce(names(chromHMMSegmentsList[[z]]),
                fWrap(allPeaks[[z]], chromHMMSegmentsList[[z]]),
                .init = vector(mode = "character",
                               length = length(allPeaks[[z]])))
})
names(segsStage) <- tissue_names
for (i in names(segsStage)) {
  allPeaks[[i]]$name <- segsStage[[i]]
}
allPeaks_annotation <- lapply(names(allPeaks), function(x) {
  allPeaks[[x]][allPeaks[[x]]$name==allPeaks[[x]]$name]$name <- State[State_indice[[x]][as.numeric(allPeaks[[x]]$name)]]
  allPeaks[[x]]
})
names(allPeaks_annotation) <- tissue_names

# Prepare the CRE annotation track hub 
paper_colors <- c("#1F78B4", "#33A02C", "#B2DF8A", "#E31A1C", "#FB9A99", 
                  "#FFA500", "#654321","#7B68EE", "#800080", "#A1A2A3")
names(paper_colors) <- c('1_TssA', '2_TssFlank1','3_TssFlank2', '4_EnhA', '5_EnhFlank', 
                         '6_Openchrom','7_Bivalent', '8_cHC','9_fHC', '10_Quies')
tissueExport <- lapply(allPeaks_annotation, function(x){
  tmp <- x
  mcols(tmp) <- DataFrame(name = tmp$name,
                          score = tmp$score,
                          itemRgb = paper_colors[tmp$name],
                          thick = ranges(tmp),
                          row.names = seqnames(tmp))
  tmp$score <- ifelse(tmp$score > 1000, 1000, tmp$score)
  tmp
})
lapply(names(tissueExport), function(name){
  export(tissueExport[[name]], 
         str_c("./Adult_Zebrafish_CREs_Annotation/", name, "_PADREs_formal_idr_github_10_chrstart.bed"))
})
system("cd ./Adult_Zebrafish_CREs_Annotation/; ls *_PADREs_formal_idr_github_10_chrstart.bed | while read id; do sort -k1,1 -k2n $id > $id.sorted.bed; done")

# Plot the number of different CREs across all adult tissues
segments <- list()
for(tissue in tissue_names) {
  segmentFiles <- str_c("./Adult_Zebrafish_CREs_Annotation/", tissue, "_PADREs_formal_idr_github_10_chrstart.bed.sorted.bed")
  segments_Files <- rtracklayer::import(segmentFiles)
  segments[[tissue]] <- segments_Files
}
segmentsNumber <- lapply(segments, function(x) {
  x %>% as.data.frame %>% group_by(name) %>% summarize(n = n())
})
segmentsNumberDf <- purrr::reduce(segmentsNumber, full_join, by = "name")
plotSegmentNumbers <- function(segments, 
                               title = "Number of segments per tissues",
                               include.y.axis = T,
                               show_legend = T){
  segmentsNumber <- lapply(segments, function(x) {
    x %>% as.data.frame %>% group_by(name) %>% summarize(n = n())
  })
  segmentsNumberDf <- purrr::reduce(segmentsNumber, full_join, by = "name")
  names(segmentsNumberDf) <- c("name", names(segmentsNumber))
  segmentsNumberDf %>%  pivot_longer(cols = -1, 
                                     names_to = "tissues", 
                                     values_to = "number") %>%
    mutate(name = factor(name, levels = rev(stringr::str_sort(unique(name), 
                                                              numeric = T))),
           tissues = factor(tissues, 
                            levels = 
                              rev(c("brain", "blood", "colon", 
                                    "heart", "intestine", "kidney", 
                                    "liver", "muscle", "skin", 
                                    "spleen", "testis")))) %>%
    ggplot(aes(name, number, fill = tissues)) + 
    geom_col(position = "dodge", show.legend = show_legend) + 
    geom_text(aes(label = number), position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 1.5) +  
    scale_fill_brewer(palette = "Set3", direction = -1) + 
    theme_bw() + {if(include.y.axis){
      theme()}else{
        theme(axis.text.y = 
                element_blank(),
              axis.title.y = element_blank())
      }} +
    coord_flip() +
    guides(fill = guide_legend(reverse = TRUE)) +
    ggtitle(title)
}
plotSegmentNumbers(segments, show_legend = T) + ylim(0, 20000)

# Validate CREs annotation by the distribution of annotated promoter on Ensembl reference genome
txDb <- makeTxDbFromGFF('https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/genes/danRer11.ensGene.gtf.gz')
brain_promoters <- segments$brain[segments$brain$name == '1_TssA']
peakAnno <- annotatePeak(brain_promoters, tssRegion=c(-3000, 3000),TxDb=txDb)
plotAnnoPie(peakAnno)
promoterAnno <- annotatePeak(brain_promoters, tssRegion=c(-500, 500),TxDb=txDb)
annotatedPeaksGR <- as.GRanges(promoterAnno)
annotatedPeaksGR[annotatedPeaksGR$annotation == 'Promoter'] %>% length()
annotatedPeaksGR$geneId <- sub("\\..*$", "", annotatedPeaksGR$geneId)
genelist <- bitr(annotatedPeaksGR$geneId, fromType="ENSEMBL",
                 toType="SYMBOL", OrgDb='org.Dr.eg.db')
tmp <- as.list(txDb)
tmp$genes$gene_id %>% unique() %>% length()
geneid <- sub("\\..*$", "", tmp$genes$gene_id)
genesymbol <- bitr(geneid, fromType="ENSEMBL",
                   toType="SYMBOL", OrgDb='org.Dr.eg.db')
intersect(genesymbol$SYMBOL, genelist$SYMBOL)
annotatedPeaksGR$geneSymbol <- ifelse(annotatedPeaksGR$geneId %in% genelist$ENSEMBL, genelist$SYMBOL[match(annotatedPeaksGR$geneId, genelist$ENSEMBL)], NA)
annotatedPeaksGR$geneSymbol %>% unique() %>% length()
annotatedPeaksGR[!is.na(annotatedPeaksGR$geneSymbol) & annotatedPeaksGR$annotation == 'Promoter']


# Figure 3
# LOLA: The enrichment of embryonic elements in adult elements
regionDB = loadRegionDB("./Adult_Zebrafish_CREs_Annotation/LOLA_github/")
regionSetA = readBed('./Adult_Zebrafish_CREs_Annotation/Hpf12_PADREs.bb.bed')
regionSetB = readBed('./Adult_Zebrafish_CREs_Annotation/Prim5_PADREs.bb.bed')
regionSetC = readBed('./Adult_Zebrafish_CREs_Annotation/LongPec_PADREs.bb.bed')
userSets = GRangesList(regionSetA, regionSetB, regionSetC)
Universe = unlist(userSets)
locResultsRestricted = runLOLA(userSets, Universe, regionDB, cores=1)
plotTopLOLAEnrichments(locResultsRestricted)
locResultsRestricted$cols <- with(locResultsRestricted, factor(ifelse(pValueLog < 10,"<10", 
                                                                      ifelse(pValueLog >= 10 & pValueLog <= 100, "10-100", ">100")),
                                                               levels = c("<10", "10-100", ">100")))
tissue_order <- c("blood_elements", "brain_elements", "colon_elements", "heart_elements", "intestine_elements", "kidney_elements", "liver_elements", "muscle_elements", "skin_elements", "spleen_elements", "testis_elements")
locResultsRestricted$description <- factor(locResultsRestricted$description, levels = desired_order)
locResultsRestricted$userSet <- ifelse(locResultsRestricted$userSet == 1, "Hpf12",
                                       ifelse(locResultsRestricted$userSet == 2, "Prim5",
                                              ifelse(locResultsRestricted$userSet == 3, "LongPec", locResultsRestricted$userSet)))
stage_reorder <- c("Hpf12", "Prim5", "LongPec")
locResultsRestricted$userSet <- factor(locResultsRestricted$userSet, levels = stage_reorder)
ggplot(locResultsRestricted, aes(x = support/size, y = description)) +
  geom_point(aes(size = support, color = as.factor(userSet))) +
  scale_colour_manual("Embryonic stages", values = c("gray", "dark blue", "red")) +
  xlim(0, 1.02) +
  labs(size = "Number of shared elements", shape = "Embryonic stages", x = 'EnrichmentRatio', y = 'Adult tissue', title = 'Enrichment of embryonic CREs in adult CREs') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = 'bold'),
        axis.title.y = element_text(size = 12, face = 'bold'),
        title = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(size = 10, face = 'bold'),
        axis.title.x = element_text(size = 12, face = 'bold'),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text =  element_text(size = 10)
  )

# UpsetR: The intersections of reused embryonic elements across adult tissues.
peakFiles <- list.files("./Adult_Zebrafish_CREs_Annotation", pattern = '_PADREs_formal_idr_github_10_chrstart.bed.sorted.bed', full.names = T)
peakGranges <- lapply(peakFiles, rtracklayer::import)
names(peakGranges) <- tissue_names
embryo_PADREs_file <- list.files('./Adult_Zebrafish_CREs_Annotation', pattern = '_PADREs.bb.bed', full.names = TRUE)
names(embryo_PADREs_file) <- c('Hpf12', 'Prim5', 'LongPec')
embryo_PADREs <- lapply(embryo_PADREs_file, rtracklayer::import)
upsethpf12 <- lapply(peakGranges, function(x) which(embryo_PADREs$Hpf12 %over% x))
upsetprim5 <- lapply(peakGranges, function(x) which(embryo_PADREs$Prim5 %over% x))
upsetlongpec <- lapply(peakGranges, function(x) which(embryo_PADREs$LongPec %over% x))
upset(fromList(upsethpf12), nsets = 11, nintersects = 11, order.by = "freq", 
      mainbar.y.max = 11500, set_size.scale_max = 35000, show.numbers = "yes",
      keep.order = T, sets = rev(c("brain", "blood", "colon", "heart", "intestine", "kidney", "liver", "muscle", "skin", "spleen", "testis")))
upset(fromList(upsetprim5), nsets = 11, nintersects = 11, order.by = "freq", 
      mainbar.y.max = 11500, set_size.scale_max = 35000, show.numbers = "yes",
      keep.order = T, sets = rev(c("brain", "blood", "colon", "heart", "intestine", "kidney", "liver", "muscle", "skin", "spleen", "testis")))
upset(fromList(upsetlongpec), nsets = 11, nintersects = 11, order.by = "freq", 
      mainbar.y.max = 11500, set_size.scale_max = 35000, show.numbers = "yes",  set_size.show = FALSE,
      keep.order = T, sets = rev(c("brain", "blood", "colon", "heart", "intestine", "kidney", "liver", "muscle", "skin", "spleen", "testis")))

# visualize the functional transition of CREs from embyros to adults using sankey plot.
# Read the CREs data of embryonic stages (Obtained from DANIO-CODE track hub)
segments_1 <- list()
stage_names <- c("Hpf12", "Prim5", "LongPec")
library(stringr)
for(stage in stage_names) {
  segmentFiles <- str_c("./Adult_Zebrafish_CREs_Annotation/", stage, "_PADREs.bb.bed")
  segments_Files <- rtracklayer::import(segmentFiles)
  segments_1[[stage]] <- segments_Files
}
sankey_df <- data.frame(Longpec_Element = character(), Brain_Element = character(), Count = numeric(), stringsAsFactors = FALSE)
for (element in unique(segments_1$LongPec$name)) {
  longpec_overlap <- subsetByOverlaps(segments_1$LongPec, segments$brain)
  longpec_overlap_element <- longpec_overlap[longpec_overlap$name == element]
  brain_overlap_element <- subsetByOverlaps(segments$brain, longpec_overlap_element)
  element_counts <- table(brain_overlap_element$name)
  selected_elements <- c("1_TssA", "2_TssFlank1", "3_TssFlank2", "4_EnhA", "5_EnhFlank", "6_Openchrom", "7_Bivalent", "8_cHC", "9_fHC", "10_Quiescent")
  selected_counts <- element_counts[selected_elements]
  df <- data.frame(Longpec_Element = element, Brain_Element = names(selected_counts), Value = as.numeric(selected_counts), stringsAsFactors = FALSE)
  sankey_df <- rbind(sankey_df, df)
}
df <- data.frame(Longpec_Element = element, Brain_Element = names(selected_counts), Value = as.numeric(selected_counts), stringsAsFactors = FALSE)
sankey_df <- rbind(sankey_df, df)
nodes <- data.frame(name = c(as.character(segments_1$LongPec$name), 
                             as.character(segments$brain$name)) %>% unique())
edges <- sankey_df
edges1 <- edges[!is.na(edges$Value), ]
rownames(edges1) <- 1:nrow(edges1)
edges1$IDsource <- match(edges1$Longpec_Element, nodes$name)-1
edges1$IDtarget <- match(edges1$Brain_Element, nodes$name)-1
sankey <- sankeyNetwork(Links = edges1,
              Nodes = nodes,
              Source = "IDsource",
              Target = "IDtarget",
              Value = "Value",
              NodeID = "name",
              LinkGroup = 'Longpec_Element',
              sinksRight = F,
              nodeWidth = 10,
              fontSize = 18,
              nodePadding = 10)
htmlwidgets::onRender(sankey, jsCode = '
  function(el, x) {
    d3.select(el).selectAll(".node text")
      .attr("x", function(d) { return d.x === 0 ? -6 : 15; })
      .attr("text-anchor", function(d) { return d.x === 0 ? "end" : "start"; });
  }
')


# Figure 4
# Plot the epigenomic feature of CREs by heatmap (Promoter, Enhancer, Open chromatin)
Enh_Granges <- lapply(seq_along(peakGranges), function(x) {
  tmp <- peakGranges[[x]][peakGranges[[x]]$name == '4_EnhA']
  tmp <- tmp[!is.na(tmp)]
  tmp
})
names(Enh_Granges) <- tissue_names
Enh_Granges_type <- lapply(seq_along(Enh_Granges), function(x) {
  Enh_Granges[[x]]$type <- 'Enhancer'
  Enh_Granges[[x]]
})
names(Enh_Granges_type) <- tissue_names
Promo_Granges <- lapply(seq_along(peakGranges), function(x) {
  tmp <- peakGranges[[x]][peakGranges[[x]]$name == '1_TssA']
  tmp
})
names(Promo_Granges) <- tissue_names
Promo_Granges_type <- lapply(seq_along(Promo_Granges), function(x) {
  Promo_Granges[[x]]$type <- 'Promoter'
  Promo_Granges[[x]]
})
names(Promo_Granges_type) <- tissue_names
Openchromatin_Granges <- lapply(seq_along(peakGranges), function(x) {
  tmp <- peakGranges[[x]][peakGranges[[x]]$name == '6_Openchrom']
  tmp
})
names(Openchromatin_Granges) <- tissue_names
Openchromatin_Granges_type <- lapply(seq_along(Openchromatin_Granges), function(x) {
  Openchromatin_Granges[[x]]$type <- 'OpenChro'
  Openchromatin_Granges[[x]]
})
names(Openchromatin_Granges_type) <- tissue_names
chromHMM_Enh_Promo_Open <- lapply(names(Promo_Granges_type), function(x) {
  tmp <- c(Promo_Granges_type[[x]], Enh_Granges_type[[x]], Openchromatin_Granges_type[[x]])
  tmp
}) 
names(chromHMM_Enh_Promo_Open) <- tissue_names
fishPlotInputs <- list(
  blood = list(
    peakRanges = peakGranges$blood,
    assayFiles =  c(ATAC = "./atac/*/call-macs2_signal_track_pooled/execution/rep.pooled.pval.signal.bigwig",
                    H3K27ac = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig",
                    H3K4me3 = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig")
  ),
  brain = list(
    peakRanges = peakGranges$brain,
    assayFiles = c(ATAC = "./atac/*/call-macs2_signal_track_pooled/execution/rep.pooled.pval.signal.bigwig",
                   H3K27ac = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_SRR9662084_1.fastq.nodup.pval.signal.bigwig",
                   H3K4me3 = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_SRR9662084_1.fastq.nodup.pval.signal.bigwig",
                   H3K27ac_SE = "./Adult_Zebrafish_CREs_Annotation/GSE75734_adult_brain_drerio_H3K27ac_normalized_lift.wig.bigwig",
                   CpG_Methylation = "./Adult_Zebrafish_CREs_Annotation/Bog_brain_trimmed_SE_danRer11_deduplicated_CpG.bigwig")
  ),
  colon = list(
    peakRanges = peakGranges$colon,
    assayFiles =  c(ATAC = "./atac/*/call-macs2_signal_track_pooled/execution/rep.pooled.pval.signal.bigwig",
                    H3K27ac = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig",
                    H3K4me3 = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig")
  ),
  heart = list(
    peakRanges = peakGranges$heart,
    assayFiles = c(ATAC = "./atac/*/call-macs2_signal_track_pooled/execution/rep.pooled.pval.signal.bigwig",
                   H3K27ac = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_SRR9662087_1.fastq.nodup.pval.signal.bigwig",
                   H3K4me3 = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_SRR9662087_1.fastq.nodup.pval.signal.bigwig",
                   H3K27ac_SE = "./Adult_Zebrafish_CREs_Annotation/GSE75734_adult_1year_heart_drerio_H3K27ac_normalized_lift.wig.bigwig")
  ),
  intestine = list(
    peakRanges = peakGranges$intestine,
    assayFiles = c(ATAC = "./atac/*/call-macs2_signal_track_pooled/execution/rep.pooled.pval.signal.bigwig",
                   H3K27ac = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig",
                   H3K4me3 = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig",
                   H3K27ac_SE = "./Adult_Zebrafish_CREs_Annotation/GSE75734_adult_1year_intestine_drerio_H3K27ac_normalized_lift.wig.bigwig")
  ),
  kidney = list(
    peakRanges = peakGranges$kidney,
    assayFiles =  c(ATAC = "./atac/*/call-macs2_signal_track_pooled/execution/rep.pooled.pval.signal.bigwig",
                    H3K27ac = "./chip/*/call-macs2_signal_track/shard-0/execution/SRR9662045_1.fastq.nodup_x_ctl.pooled.pval.signal.bigwig",
                    H3K4me3 = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig")
  ),
  liver = list(
    peakRanges = peakGranges$liver,
    assayFiles =  c(ATAC = "./atac/*/call-macs2_signal_track_pooled/execution/rep.pooled.pval.signal.bigwig",
                    H3K27ac = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig",
                    H3K4me3 = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig")
  ),
  muscle = list(
    peakRanges = peakGranges$muscle,
    assayFiles =  c(ATAC = "./atac/*/call-macs2_signal_track_pooled/execution/rep.pooled.pval.signal.bigwig",
                    H3K27ac = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig",
                    H3K4me3 = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig",
                    CpG_Methylation = "./Adult_Zebrafish_CREs_Annotation/Pot_muscle_danRer11_deduplicated_CpG.bigwig")
  ),
  skin = list(
    peakRanges = peakGranges$skin,
    assayFiles =  c(ATAC = "./atac/*/call-macs2_signal_track_pooled/execution/rep.pooled.pval.signal.bigwig",
                    H3K27ac = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig",
                    H3K4me3 = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig")
  ),
  spleen = list(
    peakRanges = peakGranges$spleen,
    assayFiles =  c(ATAC = "./atac/*/call-macs2_signal_track_pooled/execution/rep.pooled.pval.signal.bigwig",
                    H3K27ac = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig",
                    H3K4me3 = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig")
  ),
  testis = list(
    peakRanges = peakGranges$testis,
    assayFiles =  c(ATAC = "./atac/*/call-macs2_signal_track_pooled/execution/rep.pooled.pval.signal.bigwig",
                    H3K27ac = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_SRR9662100_1.fastq.nodup.pval.signal.bigwig",
                    H3K4me3 = "./chip/*/call-macs2_signal_track_pooled/execution/rep.pooled_x_SRR9662100_1.fastq.nodup.pval.signal.bigwig")
  )
)
emmisionWindows <- lapply(names(fishPlotInputs), function(x) {
  ovlp <- findOverlaps(fishPlotInputs[[x]]$peakRanges, chromHMM_Enh_Promo_Open[[x]])
  tmp <- fishPlotInputs[[x]]$peakRanges[queryHits(ovlp)]
  tmp$type <- chromHMM_Enh_Promo_Open[[x]]$type[subjectHits(ovlp)]
  tmp
})
names(emmisionWindows) <- tissue_names
Emissions <- lapply(names(fishPlotInputs), function(x) {
  targets <- lapply(fishPlotInputs[[x]]$assayFiles, import)
  list(heatmaps = ScoreMatrixList(targets,
                                  emmisionWindows[[x]] + 2000,
                                  # strand.aware = T,
                                  bin.num = 400,
                                  weight.col = "score"),
       in_data = emmisionWindows)
})
names(Emissions) <- tissue_names
lapply("brain", function(x) {
  tmp <- intersectScoreMatrixList(Emissions[[x]]$heatmaps)
  pdf(str_c(x, "_enh_promo_open_formal_10_heatmap_github_ranked.pdf"))
  grp <- emmisionWindows[[x]]$type[as.numeric(rownames(tmp$ATAC))]
  multiHeatMatrix(tmp, xcoords = c(-2000, 2000),
                  group = factor(grp),
                  order = TRUE,
                  group.col = c('blue', 'yellow', 'red'),
                  winsorize = c(0, 95), col = list(brewer.pal(6, "Blues"),
                                                   brewer.pal(6, "Reds"),
                                                   brewer.pal(6, "Greens"),
                                                   brewer.pal(6, "Oranges"),
                                                   brewer.pal(6, "Purples")))
  dev.off()
})
# Plot the aggregated signal of epigenomic feature of CREs (Promoter, Enhancer, Open chromatin)
customWinsorize <- function(x, probs = c(0, .99)){
  vals <- quantile(x, probs = probs)
  tmp <- x
  tmp[x < vals[1]] <- vals[1]
  tmp[x > vals[2]] <- vals[2]
  tmp
}
metaplot <- lapply(names(Emissions$brain$heatmaps), function(x) {
      Emissions$brain$heatmaps[[x]]@.Data %>% customWinsorize(probs = c(0, .95)) %>% as.data.frame() %>% 
        mutate(elements = emmisionWindows$brain$type[as.numeric(rownames(Emissions$brain$heatmaps[[x]]))]) %>%
        group_by(elements) %>% summarize(across(everything(), mean)) %>%
        ungroup %T>% { colnames(.)[2:ncol(.)] <- 1:(ncol(.) - 1) } %>%
        pivot_longer(-elements, names_to = "position", values_to = "score") %>%
        mutate(position = as.numeric(position),
               mark = x)
}) %>% { do.call(rbind, .) }
MetaPlt <- metaplot %>% 
  mutate(mark = factor(mark, 
                       levels = c("ATAC", "H3K27ac", "H3K4me3", "CpG_Methylation")),
         elements = factor(elements, 
                           levels = c("Enhancer", "Promoter", "OpenChro"))) %>%
  ggplot(aes(position, score, group = elements, color = elements)) + 
  geom_line(linewidth=0.5) +  facet_wrap(vars(mark), scales = "free_y", nrow = 1) + 
  theme_bw() + scale_color_manual(values = c("red", "blue", "yellow"), 
                                  labels = c("Enhancer",
                                             "Promoter", 
                                             "OpenChro")) +
  scale_x_continuous(labels = c(-1000, -500, 0, 500, 1000)) + 
  labs(color = "Elements") + theme(axis.ticks.length.y = unit(-0.1, "cm"),
                                   axis.text.y = element_text(hjust = 0,
                                                              margin = 
                                                                margin(l = 10,
                                                                       r = -20)),
                                   axis.text.x = element_blank(),
                                   axis.title.x = element_blank())
MetaPlt +
  theme(aspect.ratio = 0.8)
dev.off()

# Calculate the PhastCons conservation score of adult CREs
phastCons <- rtracklayer::import("https://research.nhgri.nih.gov/manuscripts/Burgess/zebrafish/downloads/NHGRI-1/danRer11/danRer11Tracks/ZF_GC_CC_GF.danRer11.bw")
seqlevels(phastCons) <- seqlevels(BSgenome.Drerio.UCSC.danRer11)
tmpCons <- coverage(phastCons, weight = "score")
segments <- lapply(segments, function(x){
  seqlevels(x) <- seqlevels(BSgenome.Drerio.UCSC.danRer11)
  x
})
conservation <- pblapply(segments, function(x){
  binnedAverage(x, tmpCons, varname = "phastCons")
})
# boxplot
conservation %>% names() %>% lapply(function(x){
  tmp <- conservation[[x]]
  data.frame(tissue = x, name = tmp$name, phastCons = tmp$phastCons)
}) %>% {do.call(rbind, .)} %>%
  mutate(name = factor(name, levels = stringr::str_sort(unique(name), 
                                                        numeric = T)),
         tissue = factor(tissue, 
                         levels = 
                           c("blood", "brain", "colon", "heart", 
                             "intestine", "kidney", "liver", "muscle", 
                             "skin", "spleen", "testis"))) %>%
  ggplot(aes(name, phastCons, fill = tissue)) + geom_boxplot(outlier.shape = NA) + theme_bw() +
  scale_fill_brewer(palette = "Set3") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Average cyprinides PhastCons conservation")
# heatmap
lapply(names(conservation), function(x) {
  tmp <- conservation[[x]]
  data.frame(tissue = x, name = tmp$name, phastCons = tmp$phastCons)
}) %>% {do.call(rbind, .)} %>%
  mutate(name = factor(name, levels = stringr::str_sort(unique(name), 
                                                        numeric = T)),
         tissue = factor(tissue, 
                         levels = 
                           c("blood", "brain", "colon", "heart", 
                             "intestine", "kidney", "liver", "muscle", 
                             "skin", "spleen", "testis"))) %>% group_by(tissue,name) %>% summarize(median_phastCons = median(phastCons)) %>% 
  ggplot(aes(x = tissue, y = name, fill = median_phastCons)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "skyblue") +
  theme_minimal() + theme(axis.text.x = element_text(hjust = 1, face = 'bold', size = 10),
                          axis.text.y = element_text(face = 'bold'),
                          axis.title = element_text(face = 'bold'),
                          axis.title.x = element_text(face = 'bold'),
                          axis.title.y = element_text(face = 'bold')) +
  labs(title = "Average cyprinides PhastCons conservation",
       x = "Tissue",
       y = "Name",
       fill = "phastCons")
