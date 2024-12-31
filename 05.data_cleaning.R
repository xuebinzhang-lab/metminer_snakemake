##############################################################
#            Name: One-step metminer
#            Prj: MetMiner pipeline 
#            Assignment: data cleaning
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
  'outlier', 'o', 1, 'character', 'method for outlier remove: tidymass, parameters_table or ignore',
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
  args$threads = NA
}

if (is.null(args$outlier)){
  args$outlier = 'ignore'
}

wd = args$wd |> as.character()
parameters = args$parameters |>  as.character()
threads = args$threads |> as.character()
outlier_method = args$outlier |> as.character()

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

handlers(global = TRUE)
handlers("progress")

# import mass_dataset -------------------------------------------------------
data_clean_out_file = paste0(wd,"/data_cleaning/")
if(!file.exists(data_clean_out_file)) {dir.create(data_clean_out_file,showWarnings = F,recursive = T)}
mass_dataset_pos_file = paste0(wd,"/massdataset/01.object_pos_raw.rda")
mass_dataset_neg_file = paste0(wd,"/massdataset/01.object_neg_raw.rda")

load(mass_dataset_pos_file)
load(mass_dataset_neg_file)

temp_para = readxl::read_xlsx(parameters) 

# noise remove ------------------------------------------------------------
para_mv_tbl = temp_para |> 
  dplyr::filter(process == "remove noise") |> 
  dplyr::select(para,Value)
para_mv = para_mv_tbl |> pull(Value) |> setNames(para_mv_tbl$para) |> as.list()
##> positive 

object_pos_mv.list <- tidymassshiny::find_noise_multiple(
  object = object_pos_raw,
  tag = para_mv$tag |> as.character(),
  qc_na_freq = para_mv$qc_na_freq |> as.numeric(),
  S_na_freq = para_mv$S_na_freq |> as.numeric()
)
object_pos_mv = object_pos_mv.list$object_mv
p1.pos <- object_pos_raw |> massqc::show_sample_missing_values(
  color_by = para_mv$color_by |> as.character(),
  order_by = para_mv$order_by|> as.character(),
  percentage = para_mv$percentage |> as.logical(),
  show_x_text = para_mv$show_x_text |> as.logical(),
  show_x_ticks = para_mv$show_x_ticks |> as.logical(),
  desc = para_mv$desc |> as.logical()
)

p2.pos <- object_pos_mv |> massqc::show_sample_missing_values(
  color_by = para_mv$color_by |> as.character(),
  order_by = para_mv$order_by|> as.character(),
  percentage = para_mv$percentage |> as.logical(),
  show_x_text = para_mv$show_x_text |> as.logical(),
  show_x_ticks = para_mv$show_x_ticks |> as.logical(),
  desc = para_mv$desc |> as.logical()
)

p.pos.mv = (p1.pos + ggtitle('raw data')) / (p2.pos + ggtitle('remove noise'))

ggsave(filename = paste0(data_clean_out_file,"01.mv_variable_pos.pdf"),plot = p.pos.mv,width = 12,height = 10)
ggsave(filename = paste0(data_clean_out_file,"01.mv_variable_pos.png"),plot = p.pos.mv,width = 12,height = 10,dpi = 300)

##> negative 
object_neg_mv.list <- tidymassshiny::find_noise_multiple(
  object = object_neg_raw,
  tag = para_mv$tag |> as.character(),
  qc_na_freq = para_mv$qc_na_freq |> as.numeric(),
  S_na_freq = para_mv$S_na_freq |> as.numeric()
)

object_neg_mv = object_neg_mv.list$object_mv
p1.neg <- object_neg_raw |> massqc::show_sample_missing_values(
  color_by = para_mv$color_by |> as.character(),
  order_by = para_mv$order_by|> as.character(),
  percentage = para_mv$percentage |> as.logical(),
  show_x_text = para_mv$show_x_text |> as.logical(),
  show_x_ticks = para_mv$show_x_ticks |> as.logical(),
  desc = para_mv$desc |> as.logical()
)

p2.neg <- object_neg_mv |> massqc::show_sample_missing_values(
  color_by = para_mv$color_by |> as.character(),
  order_by = para_mv$order_by|> as.character(),
  percentage = para_mv$percentage |> as.logical(),
  show_x_text = para_mv$show_x_text |> as.logical(),
  show_x_ticks = para_mv$show_x_ticks |> as.logical(),
  desc = para_mv$desc |> as.logical()
)

p.neg.mv = (p1.neg + ggtitle('raw data')) / (p2.neg + ggtitle('remove noise'))

ggsave(filename = paste0(data_clean_out_file,"01.mv_variable_neg.pdf"),plot = p.neg.mv,width = 12,height = 10)
ggsave(filename = paste0(data_clean_out_file,"01.mv_variable_neg.png"),plot = p.neg.mv,width = 12,height = 10,dpi = 300)

