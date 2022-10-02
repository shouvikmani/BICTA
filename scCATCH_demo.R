library(scCATCH)
load(paste0(system.file(package = "scCATCH"), "/extdata/mouse_kidney_203.rda"))

# demo_geneinfo
demo_geneinfo()

# revise gene symbols
mouse_kidney_203 <- rev_gene(data = mouse_kidney_203, data_type = "data", species = "Mouse", geneinfo = geneinfo)

## ----createscCATCH, echo=TRUE-------------------------------------------------
obj <- createscCATCH(data = mouse_kidney_203, cluster = mouse_kidney_203_cluster)

## ----findmarkergene, echo=TRUE------------------------------------------------
# demo_geneinfo
demo_marker()

# find highly expressed genes
obj <- findmarkergene(object = obj, species = "Mouse", marker = cellmatch, tissue = "Kidney")

## ----findcelltype, echo=TRUE--------------------------------------------------
obj <- findcelltype(object = obj)

# Results is stored in obj
obj@celltype

## ----findcelltype_new, echo=TRUE----------------------------------------------
# The most strict condition to identify marker genes
obj <- findmarkergene(object = obj, species = "Mouse", marker = cellmatch,tissue = "Kidney", use_method = "1")

# The most loose condition to identify marker genes
obj <- findmarkergene(object = obj, species = "Mouse", marker = cellmatch, tissue = "Kidney", use_method = "2")

# Other conditions to identify marker genes
obj <- findmarkergene(object = obj,species = "Mouse", marker = cellmatch, tissue = "Kidney", use_method = "1", comp_cluster = 1)

## ----findcelltype1, echo=TRUE-------------------------------------------------
# Example
cellmatch_new <- cellmatch[cellmatch$species == "Mouse" & cellmatch$tissue %in% c("Kidney", "Liver", "Lung", "Brain"), ]
obj <- findmarkergene(object = obj, if_use_custom_marker = TRUE, marker = cellmatch_new)
obj <- findcelltype(obj)

# Example
cellmatch_new <- cellmatch[cellmatch$species == "Mouse" & cellmatch$cancer %in% c("Lung Cancer", "Lymph node", "Renal Cell Carcinoma", "Prostate Cancer"), ]
obj <- findmarkergene(object = obj, if_use_custom_marker = TRUE, marker = cellmatch_new)
obj <- findcelltype(obj)

# Example
cellmatch_new <- cellmatch[cellmatch$species == "Mouse", ]
cellmatch_new <- cellmatch[cellmatch$cancer %in% c("Lung Cancer", "Lymph node", "Renal Cell Carcinoma", "Prostate Cancer") | cellmatch$tissue %in% c("Kidney", "Liver", "Lung", "Brain"), ]
obj <- findmarkergene(object = obj, if_use_custom_marker = TRUE, marker = cellmatch_new)
obj <- findcelltype(obj)

## ----findcelltype2, echo=TRUE-------------------------------------------------
# Example

# cellmatch_new is provided by users
# cellmatch_new <- rbind(cellmatch, cellmatch_new)

# Then use the new cellmatch
# a. define the species, tissue, and cancer
obj <- findmarkergene(object = obj, species = "Mouse", marker = cellmatch_new, tissue = "Kidney")
obj <- findcelltype(obj)

# b. directly use custom cellmatch
obj <- findmarkergene(object = obj, if_use_custom_marker = TRUE, marker = cellmatch_new)
obj <- findcelltype(obj)

## ----findcelltype3, echo=TRUE-------------------------------------------------
# Please refer to demo_marker to build a marker data.frame (new_cellmatch) for another species, e.g., rat
# Then use the new marker
obj <- findmarkergene(object = obj, species = "Rat", if_use_custom_marker = TRUE, marker = cellmatch_new, tissue = "Kidney")
obj <- findcelltype(obj)
