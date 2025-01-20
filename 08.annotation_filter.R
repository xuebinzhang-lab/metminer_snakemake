##############################################################
#            Name: One-step metminer
#            Prj: MetMiner pipeline 
#            Assignment: deduplication
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
  'parameters', 's', 1, 'character','An excel file with default parameters',
  'threads', 't', 2, 'integer', 'threads for peak picking steps'
),byrow = T, ncol = 5)

args = getopt(command)

##> defaults

if (is.null(args$wd)){
  message(msg_no("-w error\n need set working directory"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}

if (is.null(args$parameters)){
  message(msg_no("-w error\n need provide parameters tabels"))
  q(status = 1)
}


if (is.null(args$threads)){
  args$threads = 1
}

if (is.null(args$customized_db)){
  args$customized_db = NULL
}
wd = args$wd |> as.character()
parameters = args$parameters |>  as.character()
threads = args$threads |> as.character()
customized_db = args$customized_db |>as.character()
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
handlers(global = TRUE)
handlers("progress")


# run_annotation ----------------------------------------------------------

##> parameters

message(msg_run("\nStep2. Annotation filtering"))
temp_para = readxl::read_xlsx(parameters) 
para_anno_filter_tbl = temp_para |> 
  dplyr::filter(process == "annotation filtering") |> 
  dplyr::select(para,Value)

para_anno_filter = para_anno_filter_tbl |> pull(Value) |> setNames(para_anno_filter_tbl$para) |> as.list()

Adduct_pos = stringr::str_split(
  string = para_anno_filter$Adduct_pos, pattern = "\\|",n = Inf,simplify = FALSE
) |> unlist()

Adduct_neg = stringr::str_split(
  string = para_anno_filter$Adduct_neg, pattern = "\\|",n = Inf,simplify = FALSE
) |> unlist()

multi_anno = para_anno_filter$multi_anno %>% as.character()
redundancy = para_anno_filter$redundancy %>% as.character()
column = para_anno_filter$column %>% as.character()

load("massdataset/07.object_neg_anno.rda")
load("massdataset/07.object_pos_anno.rda")

annotation_tbl_pos <- 
  object_pos_anno %>% extract_annotation_table()

