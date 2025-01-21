##############################################################
#            Name: One-step metminer
#            Prj: MetMiner pipeline 
#            Assignment: export data
#            Author: Xiao Wang (shawnwang2016@126.com)
#            Date: Dec 31, 2024
##############################################################

# getopts -----------------------------------------------------------------

suppressMessages(library(getopt))
suppressMessages(library(crayon))
msg_yes = green$bold$italic
msg_no = red$bold$italic
msg_run = blue$bold$italic$underline
msg_warning = yellow$bold$italic

command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'wd', 'w', 1, 'character', 'working directory '
),byrow = T, ncol = 5)

args = getopt(command)

##> defaults

if (is.null(args$wd)){
  message(msg_no("-w error\n need set working directory"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}


wd = args$wd |> as.character()
setwd(wd)
# load package ------------------------------------------------------------

options(stringsAsFactors = F)
options(warn = -1)

suppressMessages(library(tidymass))
suppressMessages(library(tidyverse))
suppressMessages(library(tidymassshiny))
suppressMessages(library(plantmdb))
suppressMessages(library(progressr))
suppressMessages(library(patchwork))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
handlers(global = TRUE)
handlers("progress")


# export  -----------------------------------------------------------------

##> merge object
load("massdataset/07.object_pos_anno.rda")
load("massdataset/07.object_neg_anno.rda")
dir.create("Data_export/POS",showWarnings = F,recursive = T)
dir.create("Data_export/NEG",showWarnings = F,recursive = T)
export_ms2_data(object = object_pos_anno,file_type = 'mgf',path = "Data_export/POS/")
export_ms2_data(object = object_neg_anno,file_type = 'mgf',path = "Data_export/NEG/")

object_merge <- 
  merge_mass_dataset_fix(
    x = object_neg_anno,
    y = object_pos_anno,
    sample_direction = 'inner',
    variable_direction = 'full'
  ) %>%
  activate_mass_dataset("sample_info") %>% 
  filter(class != "QC")

expmat <- object_merge %>% extract_expression_data() %>% 
  rownames_to_column("variable_id")

variable_info.pos <- object_pos_anno %>% 
  extract_variable_info() %>% 
  dplyr::select(variable_id,mz,rt)

variable_info.neg <- object_neg_anno %>% 
  extract_variable_info() %>% 
  dplyr::select(variable_id,mz,rt)

variable_info = rbind(variable_info.pos,variable_info.neg)
save(object_merge,file = "massdataset/08.object_merge.rda")
writexl::write_xlsx(expmat,"Data_export/expmat.xlsx")
writexl::write_xlsx(variable_info,"Data_export/variable_info.xlsx")