library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(stringr)

# Figure 1
# DANIO-CODE datasets are obtained from https://github.com/DANIO-CODE/DANIO-CODE_Data_analysis/blob/master/Figures/Figure1/DCC-annotations-2021-05.csv
# This adult zebrafish track hub can be accessed through the URL https://data.cyverse.org/dav-anon/iplant/home/guanguan/Adult_zebrafish_hub/hub.txt in the UCSC Genome Browser 

# Figure 2
setwd("./ChromHMM_file/")
tissue_names <- c("blood", "brain", "colon", "heart", "intestine", "kidney", "liver", "muscle", "skin", "spleen", "testis")
dir_paths <- paste0("./ChromHMM_file/", tissue_names, '/ChromHMM_out')
emission_data <- list()
for (dir_path in dir_paths) {
  file_path <- file.path(dir_path, "emissions_10.txt")
  tissue_name <- basename(dirname(dir_path))
  emission_df <- read.table(file_path, header = TRUE)
  emission_df$tissue <- tissue_names
  emission_data[[dir_path]] <- emission_df
}

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
emission_data

# for (i in seq_along(emission_data)) {
#   emission_data[[i]]$State <- State[State_indice[[i]]]
# }

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

# Emission for adult brain only
brain_emission <- emission_data$brain
brain_emission
brain_emission_pivoted <- brain_emission %>%
  pivot_longer(!c(State,tissue,State), names_to = "Mark", values_to = "Emission")
brain_emission_pivoted$State <- factor(brain_emission_pivoted$State,
                                       levels = rev(c('1_TssA', '2_TssFlank1','3_TssFlank2', '4_EnhA','5_EnhFlank', '6_Openchrom','7_Bivalent', '8_cHC','9_fHC', '10_Quies')),
                                       ordered = TRUE)
ggplot(brain_emission_pivoted, aes(State, Mark, fill = Emission)) + geom_tile() +
  theme_classic() + coord_flip() + facet_grid(~ tissue) +
  theme(axis.text.x = element_text(hjust = 1, angle = 90, face = 'bold', size = 14),
        axis.text.y = element_text(face = 'bold', size = 14),
        axis.title = element_text(face = 'bold', size = 14),
        axis.title.x = element_text(face = 'bold', size = 14),
        axis.title.y = element_text(face = 'bold', size = 14),
        strip.text.x = element_text(face = 'bold', size = 14)
  ) +
  scale_fill_distiller(palette = "Blues", direction = 1)
dev.off()

# CRE annotation with ChromHMM results
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
segsStage
names(segsStage) <- tissue_names
for (i in names(segsStage)) {
  allPeaks[[i]]$name <- segsStage[[i]]
}

allPeaks_annotation <- lapply(names(allPeaks), function(x) {
  allPeaks[[x]][allPeaks[[x]]$name==allPeaks[[x]]$name]$name <- State[State_indice[[x]][as.numeric(allPeaks[[x]]$name)]]
  allPeaks[[x]]
})
names(allPeaks_annotation) <- tissue_names
lapply(allPeaks_annotation, function(x) {
  unique(x$name)
})

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

segments <- list()
for(tissue in tissue_names) {
  segmentFiles <- str_c("/rds/projects/m/muellerf-cage/Qianhong/Adult_Zebrafish_CREs_Annotation/", tissue, "_PADREs_formal_idr_github_10_chrstart.bed.sorted.bed")
  segments_Files <- rtracklayer::import(segmentFiles)
  segments[[tissue]] <- segments_Files
}
