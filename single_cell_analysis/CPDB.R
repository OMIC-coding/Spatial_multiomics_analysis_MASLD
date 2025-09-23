devtools::install_github("zktuong/ktplots", dependencies = TRUE)
library(schard)

scRNA_masld = schard::h5ad2seurat(filename='/data/project/AI4Omic/MASLD/results/scRNA/CellphoneDB/adata_major.h5ad')
#scRNA_masld$cell_type_lvl1 <- gsub(" ", '.', scRNA_masld$cell_type_lvl1)
library(ktplots)
data(cpdb_output2) # legacy reasons
pvals <- read.delim("/data/project/AI4Omic/MASLD/results/scRNA/CellphoneDB/statistical_analysis_pvalues_MASH.txt", check.names = FALSE)
pvals <- pvals %>% filter(!classification %in% c('Adhesion by Collagen/Integrin'))
means <- read.delim("/data/project/AI4Omic/MASLD/results/scRNA/CellphoneDB/statistical_analysis_means_MASH.txt", check.names = FALSE)
decon = read.csv("/data/project/AI4Omic/MASLD/results/scRNA/CellphoneDB/statistical_analysis_deconvoluted_MASH.txt", sep="\t")
interaction_annotation = data.frame('interacting_pair'=pvals$interacting_pair, 'classification'=pvals$classification)
pdf('/data/project/AI4Omic/MASLD/results/scRNA/CellphoneDB/circos-Central Vein EC-HSC.pdf')
cell_type1 <- 'LAM'
cell_type2 <- 'Treg'



devtools::install_github("Sarah145/CCPlotR")

library(CCPlotR)
library(dplyr)
data("toy_data")

interaction <- read.delim("/data/project/AI4Omic/MASLD/results/scRNA/CellphoneDB/statistical_analysis_interaction_scores_MASH.txt", check.names = FALSE)
CCplot_data <- interaction %>% select('gene_a', 'gene_b', 'Central_Vein_EC|HSC', 'HSC|Central_Vein_EC') %>%
  rename('ligand' = 'gene_a', 'receptor' = 'gene_b') %>%
  filter(receptor != '', ligand != '') %>% 
  tidyr::pivot_longer(cols = c('Central_Vein_EC|HSC', 'HSC|Central_Vein_EC'), names_to = 'interaction', values_to = 'score') %>%
  tidyr::separate(interaction, into = c('source', 'target'), sep='\\|') %>%
  filter(score > 0)

genes <- unique(c(CCplot_data$ligand, CCplot_data$receptor))
library(Seurat)
AE <- AverageExpression(scRNA_masld, features = genes, return.seurat = FALSE, group.by = 'cell_type_lvl2_noblank')
CCplot_exp <- AE %>% as.data.frame() %>% tibble::rownames_to_column('gene') %>% tidyr::pivot_longer(cols = -gene, names_to = 'cell_type', values_to = 'mean_exp')
CCplot_exp$cell_type <- CCplot_exp$cell_type %>% gsub('RNA.', '', .) %>% gsub('\\.', '_', .)
CCplot_exp %>% filter(cell_type %in% c('Central Vein EC', 'HSC'))
pdf('/data/project/AI4Omic/MASLD/results/scRNA/CellphoneDB/circos-RSPO3-LGR6.pdf', width = 7, height = 6)
cc_circos(CCplot_data, exp_df = CCplot_exp, option = 'C', n_top_ints = 15, cell_cols = c(`HSC` = "#9467bd", `Central_Vein_EC` = "#E2E8A7"), palette = 'PuRd')
dev.off()



plot_cpdb2(
  scdata = scRNA_masld %>% as.SingleCellExperiment(),
  cell_type1 = cell_type1,
  cell_type2 = cell_type2,
  celltype_key = "cell_type_lvl2_noblank", # column name where the cell ids are located in the metadata
  means = means,
  pvals = pvals,
  keep_significant_only = TRUE,
  frac = 0.2,
  plot_score_as_thickness = TRUE,
  deconvoluted = decon, # new options from here on specific to plot_cpdb2
  desiredInteractions = list(
    c(cell_type1, cell_type2)
  ),
  interaction_grouping = interaction_annotation,
  node_group_colors = c(
    cell_type1 = "#E2E8A7",
    cell_type2 = "#9467bd"
  ),
  return_df = FALSE,
  remove_self = TRUE,
)
dev.off()

pvals[pvals$interacting_pair %in% df[[1]]$pair, ]$classification %>% table()

df[[1]]$
edge_group_colors = c(
  "Activating" = "#e15759",
  "Chemotaxis" = "#59a14f",
  "Inhibitory" = "#4e79a7",
  "Intracellular trafficking" = "#9c755f",
  "DC_development" = "#B07aa1",
  "Unknown" = "#e7e7e7"
),


pdf('/data/project/AI4Omic/MASLD/results/scRNA/CellphoneDB/circos-HGF-MET.pdf')
plot_cpdb4(
  scdata = scRNA_masld %>% as.SingleCellExperiment(),
  interaction = "RSPO3-LGR6",
  cell_type1 = "HSC",
  cell_type2 = "HSC",
  celltype_key = "cell_type_lvl2_noblank",
  means = means,
  pvals = pvals,
  deconvoluted = decon,
  keep_significant_only = TRUE
)
