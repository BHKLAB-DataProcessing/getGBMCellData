library(stringr)
library(openxlsx)

input_dir <- "/pfs/downloadGBMData/"
out_dir <- "/pfs/out/"
# input_dir <- "~/Documents/pfs/downloadGBMData/"
# out_dir <- "~/Documents/pfs/getGBMCellData/" 

cell <- read.xlsx(paste(input_dir, "mmc2.xlsx", sep=""),rowNames = TRUE)#Cell_line names
cell$Patient_id <- rownames(cell)

saveRDS(cell, paste0(out_dir, "cell.rds"))