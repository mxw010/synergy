rm(list=ls(all=T))
# Load neccessrary packages
require("RColorBrewer")
require("dplyr")
require("ggplot2")
require("reshape")
require("scales")
require("data.table")
require("stringr")
library(tidyr)
#library(tidyext)


# Get working directory of current file
# current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(current_working_dir)
setwd("/Users/mxw010/Dropbox/Stanton Lab Synergy/2D sensitivities_Emma Chory and Meng Wang synergy analysis")
source('Meng/get_compiled_data.R')
experiments = list.files("data/")
#experiments = experiments[-grep("data", experiments)]

dose_cutoff = 5

# experiments = experiments[5:length(experiments)]
exp = experiments[1]

compiled_plates <- NULL
average_dose <- NULL
individual_doses <- NULL

for (exp in experiments) {
  # load data for each experiment
  results = fread(paste0("data/", exp, "/responses.csv"), stringsAsFactors = FALSE, fill = TRUE)
  metadata = fread(paste0("data/", exp, "/metadata.csv"), stringsAsFactors = FALSE, fill = TRUE)
  
  #metadata %>% filter(BlockId == 6)
  
  rows = as.data.frame(str_split_fixed(metadata$RowConcs, pattern = ",", n = 10), stringsAsFactors = FALSE)
  rows = mutate_all(rows, function(x) as.numeric(as.character(x)))
  colnames(rows) = 1:10
  rows$Drug = metadata$RowName
  rows$BlockId = metadata$BlockId
  rows = rows %>% pivot_longer(as.character(1:10), names_to = "Row")
  
  colnames(rows) = c("Name1", "BlockId", "Row", "Drug1")
  rows$Row = as.numeric(rows$Row)
  
  cols = as.data.frame(str_split_fixed(metadata$ColConcs, pattern = ",", n = 10), stringsAsFactors = FALSE)
  cols = mutate_all(cols, function(x) as.numeric(as.character(x)))
  colnames(cols) = 1:10
  cols$Drug = metadata$ColName
  cols$BlockId = metadata$BlockId
  cols = cols %>% pivot_longer(as.character(1:10), names_to = "Row")
  colnames(cols) = c("Name2", "BlockId", "Col", "Drug2")
  cols$Col = as.numeric(cols$Col)
  cols = cols %>% filter(!is.na(Drug2))
  
  cols = data.frame(cols)
  rows = data.frame(rows)
  results = data.frame(results)
  
  results = merge(results, rows, all.x = TRUE)
  results = merge(results, cols, all.x = TRUE)
  
  drugcombos = unique(results %>% select(BlockId, Name1, Name2)) %>% arrange(BlockId)
  rownames(drugcombos) = 1:nrow(drugcombos)
  
  # organize based on combo with less drugs
  count1 = length(unique(drugcombos$Name1))
  count2 = length(unique(drugcombos$Name2))
  if (count2 < count1) {
    blockids = drugcombos$BlockId
    drugcombos = drugcombos %>% select(-BlockId)
    colnames(drugcombos) = rev(colnames(drugcombos))
    drugcombos$BlockId = blockids
    colnames(results) = gsub(pattern = "1", "XX", colnames(results))
    colnames(results) = gsub(pattern = "2", "1", colnames(results))
    colnames(results) = gsub(pattern = "XX", "2", colnames(results))
  }
  for (id in unique(drugcombos$BlockId)) {
    y <- get_compiled_data(results, id, exp)
    if (!is.null(y)) {
      compiled_plates <- rbind(compiled_plates, y$compiled_plates)
      average_dose <- rbind(average_dose, y$average_dose)
      individual_doses <- rbind(individual_doses, y$individual_doses)
    }
  }
}

write.csv(compiled_plates,'Meng/compiled_plates.csv',quote=F,row.names=F)
write.csv(average_dose,'Meng/average_dose.csv',quote=F,row.names=F)
write.csv(individual_doses,'Meng/individual_doses.csv',quote=F,row.names=F)
