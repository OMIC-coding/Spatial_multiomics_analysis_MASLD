# CellChat
install.packages('BiocManager')
install.packages('devtools')
BiocManager::install('SummarizedExperiment')
BiocManager::install('SingleCellExperiment')
devtools::install_github('JiekaiLab/dior')
devtools::install_github("jinworks/CellChat", force = TRUE)
devtools::install_github('immunogenomics/presto')
install.packages('NMF')

library(reticulate)
use_condaenv('GeneTrajectory', required = TRUE)
use_python(python = '/home/bailab/miniconda3/envs/GeneTrajectory/bin/python', required = TRUE)


devtools::install_github("cellgeni/schard")

# Load the relevant packages
library(dplyr)
library(CellChat)
library(Seurat)
library(SummarizedExperiment)

### Set the ligand-receptor database. Here we will use the "Secreted signalling" database for cell-cell communication (let's look at ECM-receptor in the future )
CellChatDB <- CellChatDB.human # The othe roption is CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Other options include ECM-Receptor and Cell-Cell Contact

#adata_seurat = schard::h5ad2seurat(filename='/data/project/AI4Omic/MASLD/results/scRNA/CellphoneDB/adata_major.h5ad')
adata_seurat = schard::h5ad2seurat(filename='/data/project/AI4Omic/MASLD/results/scRNA/CellphoneDB/adata_major_29samples.h5ad')
adata_seurat$samples <- adata_seurat$sample

# Data preparation
data.input <- GetAssayData(adata_seurat, assay = 'RNA', layer = 'data')
meta.data <- data.frame(adata_seurat@meta.data)

#### First by merged cell types ####
#burkhardt.identity <- data.frame(group = burkhardt.meta.data$Type, row.names = row.names(burkhardt.meta.data))

ctrl.use <- rownames(meta.data)[meta.data$disease_status == "CTRL"] # extract the cell names from Control data
masl.use <- rownames(meta.data)[meta.data$disease_status == "MASL"] # extract the cell names from MASL data
mash.use <- rownames(meta.data)[meta.data$disease_status == "MASH"] # extract the cell names from MASH data

ctrl.data.input <- data.input[, ctrl.use]
#masld.data.input <- data.input[, masld.use]
masl.data.input <- data.input[, masl.use]
mash.data.input <- data.input[, mash.use]

ctrl.meta <- meta.data[meta.data$disease_status == "CTRL", ]
masl.meta <- meta.data[meta.data$disease_status == 'MASL', ]
mash.meta <- meta.data[meta.data$disease_status == "MASH", ]

# Third Revision from cell_type_lvl2 to cell_type_lvl1
ctrl.identity <- data.frame(group = ctrl.meta$cell_type_lvl2, row.names = row.names(ctrl.meta), samples = ctrl.meta$sample)
masl.identity <- data.frame(group = masl.meta$cell_type_lvl2, row.names = row.names(masl.meta), samples = masl.meta$sample)
mash.identity <- data.frame(group = mash.meta$cell_type_lvl2, row.names = row.names(mash.meta), samples = mash.meta$sample)


# Create the cellchat objects
ctrl.cc <- createCellChat(object = ctrl.data.input, do.sparse = T, meta = ctrl.identity, group.by = "group")
levels(ctrl.cc@idents) # show factor levels of the cell labels

masl.cc <- createCellChat(object = masl.data.input, do.sparse = T, meta = masl.identity, group.by = "group")
levels(masl.cc@idents) # show factor levels of the cell labels

mash.cc <- createCellChat(object = mash.data.input, do.sparse = T, meta = mash.identity, group.by = "group")
levels(mash.cc@idents) # show factor levels of the cell labels

ctrlGroupSize <- as.numeric(table(ctrl.cc@idents)) # Get the number of cells in each group
maslGroupSize <- as.numeric(table(masl.cc@idents)) # Get the number of cells in each group
mashGroupSize <- as.numeric(table(mash.cc@idents)) # Get the number of cells in each group

# Set the databases
ctrl.cc@DB <- CellChatDB.use 
masl.cc@DB <- CellChatDB.use 
mash.cc@DB <- CellChatDB.use


## We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
### Control
ctrl.cc <- subsetData(ctrl.cc) # We subset the expression data of signalling genes to save on computational cost
ctrl.cc <- identifyOverExpressedGenes(ctrl.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
ctrl.cc <- identifyOverExpressedInteractions(ctrl.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
ctrl.cc <- smoothData(ctrl.cc, adj=PPI.human) # Other option includes PPI.human We're told that we may have to comment these out.

ctrl.cc <- computeCommunProb(ctrl.cc, raw.use = FALSE, population.size = FALSE)
ctrl.cc <- computeCommunProbPathway(ctrl.cc) # Calculate the probabilities at the signalling level
ctrl.cc <- aggregateNet(ctrl.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

### masl
masl.cc <- subsetData(masl.cc) # We subset the expression data of signalling genes to save on computational cost
masl.cc <- identifyOverExpressedGenes(masl.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
masl.cc <- identifyOverExpressedInteractions(masl.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
masl.cc <- smoothData(masl.cc, adj=PPI.human) # Other option includes PPI.human We're told that we may have to comment these out.

masl.cc <- computeCommunProb(masl.cc, raw.use = FALSE, population.size = FALSE)
masl.cc <- computeCommunProbPathway(masl.cc) # Calculate the probabilities at the signalling level
masl.cc <- aggregateNet(masl.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

### mash
mash.cc <- subsetData(mash.cc) # We subset the expression data of signalling genes to save on computational cost
mash.cc <- identifyOverExpressedGenes(mash.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
mash.cc <- identifyOverExpressedInteractions(mash.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
mash.cc <- smoothData(mash.cc, adj=PPI.human) # Other option includes PPI.human We're told that we may have to comment these out.

mash.cc <- computeCommunProb(mash.cc, raw.use = FALSE, population.size = FALSE)
mash.cc <- computeCommunProbPathway(mash.cc) # Calculate the probabilities at the signalling level
mash.cc <- aggregateNet(mash.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities


ctrl.net <- subsetCommunication(ctrl.cc)
masl.net <- subsetCommunication(masl.cc)
mash.net <- subsetCommunication(mash.cc)

### save results
write.csv(ctrl.net, file = "/data/project/AI4Omic/MASLD/results/scRNA/flowsig/cellchat_communications_ctrl.csv")
write.csv(masl.net, file = "/data/project/AI4Omic/MASLD/results/scRNA/flowsig/cellchat_communications_masl.csv")
write.csv(mash.net, file = "/data/project/AI4Omic/MASLD/results/scRNA/flowsig/cellchat_communications_mash.csv")