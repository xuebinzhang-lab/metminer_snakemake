##############################################################
#            Name: One-step metminer
#            Prj: MetMiner pipeline 
#            Assignment: create massdataset
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
  'wd', 'w', 1, 'character', 'working directory ',
  'format', 'f', 1, 'character', 'input data format: rda or csv, default: rda',
  'metadata', 'm', 1, 'character', 'An .csv file of sample information'
),byrow = T, ncol = 5)

args = getopt(command)

##> defaults

if (is.null(args$wd)){
  message(msg_no("-w error\n need set working directory"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}


if (is.null(args$metadata)){
  message(msg_warning("-w error\need upload sample information in .csv format"))
  q(status = 1)
}


if (is.null(args$format)){
  args$format = 'rda'
}

wd = args$wd |> as.character()
metadata = args$metadata |> as.character()
format = args$format |> as.character()

setwd(wd)

# load package ------------------------------------------------------------

options(stringsAsFactors = F)
options(warn = -1)

suppressMessages(library(tidymass))
suppressMessages(library(tidyverse))
suppressMessages(library(progressr))

handlers(global = TRUE)
handlers("progress")

# start object -------------------------------------------------------

sample_info <- read.csv(metadata)

##> start with massdataset
if(format == 'rda') {
  if(!file.exists(paste0(wd,"/massdataset/object_neg_raw.rda")) | !file.exists(paste0(wd,"/massdataset/object_pos_raw.rda"))) {
    message(msg_no("\nError: no .rda file detected, please make sure your input is .rda format"))
    q(status = 1)
  }
  load(paste0(wd,"/massdataset/object_pos_raw.rda"))
  object_pos_raw <-
    object_pos_raw %>% 
    activate_mass_dataset('sample_info') %>% 
    dplyr::select(sample_id) %>% 
    dplyr::left_join(sample_info)
  load(paste0(wd,"/massdataset/object_neg_raw.rda"))
  object_neg_raw <-
    object_neg_raw %>% 
    activate_mass_dataset('sample_info') %>% 
    dplyr::select(sample_id) %>% 
    dplyr::left_join(sample_info)
} 

if(format == "csv") {
  if(!file.exists(paste0(wd,"/peakpicking_tbl.csv"))){
    message(msg_no("\nError: no peakpicking_tbl.csv file detected, please make sure your input file is renamed as peakpicking_tbl.csv"))
    q(status = 1)
  }
  temp_tbl = read.csv(paste0(wd,"/peakpicking_tbl.csv"))
  variable_info = temp_tbl %>% 
    dplyr::select(variable_id,mz,rt,ion) %>%
    dplyr::mutate(ion = dplyr::case_when(
      str_detect(ion,regex("\\+|pos",ignore_case = T)) ~ "pos",
      str_detect(ion,regex("\\-|neg",ignore_case = T)) ~ "neg"))
  expression_data = temp_tbl %>% 
    dplyr::select(-mz,-rt,-ion) %>%
    tibble::column_to_rownames("variable_id") %>%
    dplyr::mutate_if(is.character,as.numeric)
  ##> pos
  variable_info_pos = 
    variable_info %>%
    dplyr::filter(ion == 'pos')
  expression_data_pos <- 
    expression_data %>% 
    tibble::rownames_to_column("variable_id") %>%
    dplyr::inner_join(variable_info_pos,by = 'variable_id') %>%
    tibble::column_to_rownames("variable_id") %>%
    dplyr::select(sample_info %>% dplyr::pull(sample_id))
  object_pos_raw <- 
    create_mass_dataset(
      expression_data = expression_data_pos,
      variable_info = variable_info_pos,
      sample_info = sample_info
    )
  ##> neg
  variable_info_neg = 
    variable_info %>%
    dplyr::filter(ion == 'neg')
  expression_data_neg <- 
    expression_data %>% 
    tibble::rownames_to_column("variable_id") %>%
    dplyr::inner_join(variable_info_neg,by = 'variable_id') %>%
    tibble::column_to_rownames("variable_id") %>%
    dplyr::select(sample_info %>% dplyr::pull(sample_id))
  object_neg_raw <- 
    create_mass_dataset(
      expression_data = expression_data_neg,
      variable_info = variable_info_neg,
      sample_info = sample_info
    )
}


# save mass_dataset -----------------------------------------------------

mass_dataset_pos_file = paste0(wd,"/massdataset/01.object_pos_raw.rda")
mass_dataset_neg_file = paste0(wd,"/massdataset/01.object_neg_raw.rda")

save(object_neg_raw,file = mass_dataset_neg_file)
save(object_pos_raw,file = mass_dataset_pos_file)

