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
message(msg_run("\nStep1. remove noise..."))
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

noise_tbl_pos = object_pos_mv.list$noisy_tbl
noise_tbl_neg = object_neg_mv.list$noisy_tbl
writexl::write_xlsx(x = list(pos = noise_tbl_pos,
                             neg = noise_tbl_neg),
                    path = paste0(data_clean_out_file,"01.noise_information.xlsx"))
save(object_pos_mv,file = paste0(wd,"/massdataset/02.object_pos_mv.rda"))
save(object_neg_mv,file = paste0(wd,"/massdataset/02.object_neg_mv.rda"))
# remove outlier ----------------------------------------------------------
message(msg_run("Step2. detect outlier..."))

##> parameters
para_outlier_tbl = temp_para |> 
  dplyr::filter(process == "remove outlier") |> 
  dplyr::select(para,Value)
para_outlier = para_outlier_tbl |> pull(Value) |> setNames(para_outlier_tbl$para) |> as.list()

##> detect outlier by masscleaner
outlier_samples.neg <-
  object_neg_mv %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples.pos <-
  object_pos_mv %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_pos = extract_outlier_table(outlier_samples.pos)
outlier_neg = extract_outlier_table(outlier_samples.neg)

if(para_outlier$method %>% as.character() == "tidymass") {
  outlier_sid.pos =
    outlier_pos %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(
      !sample_id,values_to = 'judge',names_to = 'condition'
    ) %>%
    dplyr::filter(str_detect(condition,para_outlier$by_witch)) %>%
    group_by(sample_id) %>%
    summarise(is_outlier = all(judge == TRUE)) %>%
    filter(is_outlier) %>%
    ungroup() %>%
    pull(sample_id)
  if(length(outlier_sid.pos) == 0 ) {
    outlier_sid.pos <- NA
  }
  outlier_sid.neg =
    outlier_neg %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(
      !sample_id,values_to = 'judge',names_to = 'condition'
    ) %>%
    dplyr::filter(str_detect(condition,para_outlier$by_witch)) %>%
    group_by(sample_id) %>%
    summarise(is_outlier = all(judge == TRUE)) %>%
    filter(is_outlier) %>%
    ungroup() %>%
    pull(sample_id)
  if(length(outlier_sid.neg) == 0 ) {
    outlier_sid.neg <- NA
  }
  out_para_record = paste0("method: tidymass\nby: ", para_outlier$by_witch)
} else if(para_outlier$method == 'ignore') {
  outlier_sid.pos <- NA
  outlier_sid.neg <- NA
} else if(para_outlier$method == 'customized') {
  if(any(is.na(para_outlier$pos_sid %>% as.character()))) {
    outlier_sid.pos <- NA
  } else {
    outlier_sid.pos <- para_outlier$pos_sid %>% as.character() %>% str_split(string = .,pattern = "\\|",n = Inf) %>% unlist()
  }
  if(any(is.na(para_outlier$neg_sid %>% as.character()))) {
    outlier_sid.neg <- NA
  } else {
    outlier_sid.neg <- para_outlier$neg_sid %>% as.character() %>% str_split(string = .,pattern = "\\|",n = Inf)  %>% unlist()
  }
  out_para_record = paste0("method: customized")
}

if(any(is.na(outlier_sid.pos))) {
  message(msg_yes("no outlier detected in positive model"))
} 

if(any(is.na(outlier_sid.neg))) {
  message(msg_yes("no outlier detected in negative model"))
} 

if(all(!is.na(outlier_sid.pos)) & all(outlier_sid.pos %in% rownames(outlier_neg))) {
  object_pos_outlier =
    object_pos_mv %>% activate_mass_dataset(what = 'sample_info') %>%
    filter(!sample_id %in% outlier_sid.pos)
  message(msg_yes(paste0(paste0(outlier_sid.pos,collapse = ", ")," were labeled as outlier and removed from your positive dataset! according: ",para_outlier$method)))
  pca_before.pos = massqc::massqc_pca(
    object = object_pos_mv %>% +1 %>% log(2) %>% scale(),
    color_by = para_outlier$color_by %>% as.character()
  )
  pca_after.pos = massqc::massqc_pca(
    object = object_pos_outlier %>% +1 %>% log(2) %>% scale(),
    color_by = para_outlier$color_by %>% as.character()
  )
  pca_plt_out_pos = (pca_before.pos + ggtitle(label = "Before")) + (pca_after.pos+ ggtitle(label = "After")) + 
    patchwork::plot_layout(guides = 'collect') + patchwork::plot_annotation(tag_levels = "a")
  ggsave(filename = paste0(data_clean_out_file,"/02.PCA_before_after_remove_outlier_pos.pdf"),plot = pca_plt_out_pos,width = 12,height = 5.8)
  ggsave(filename = paste0(data_clean_out_file,"/02.PCA_before_after_remove_outlier_pos.png"),plot = pca_plt_out_pos,width = 12,height = 5.8)
  res_outlier_record_pos = c("Outlier in positive model:",
                         outlier_sid.pos,
                         out_para_record)
  writeLines(res_outlier_record_pos,paste0(data_clean_out_file,"/02.outlier_record_pos.txt"))
} else {
  message(msg_no("error: The sample id seems not included in this object!"))
}

if(all(!is.na(outlier_sid.neg)) & all(outlier_sid.neg %in% rownames(outlier_neg))) {
  object_neg_outlier =
    object_neg_mv %>% activate_mass_dataset(what = 'sample_info') %>%
    filter(!sample_id %in% outlier_sid.neg)
  message(msg_yes(paste0(paste0(outlier_sid.neg,collapse = ", ")," were labeled as outlier and removed from your positive dataset! according: ",para_outlier$method)))
  pca_before.neg = massqc::massqc_pca(
    object = object_neg_mv %>% +1 %>% log(2) %>% scale(),
    color_by = para_outlier$color_by %>% as.character()
  )
  pca_after.neg = massqc::massqc_pca(
    object = object_neg_outlier %>% +1 %>% log(2) %>% scale(),
    color_by = para_outlier$color_by %>% as.character()
  )
  pca_plt_out_neg = (pca_before.neg + ggtitle(label = "Before")) + (pca_after.neg+ ggtitle(label = "After")) + 
    patchwork::plot_layout(guides = 'collect') + patchwork::plot_annotation(tag_levels = "a")
  ggsave(filename = paste0(data_clean_out_file,"/02.PCA_before_after_remove_outlier_neg.pdf"),plot = pca_plt_out_neg,width = 12,height = 5.8)
  ggsave(filename = paste0(data_clean_out_file,"/02.PCA_before_after_remove_outlier_neg.png"),plot = pca_plt_out_neg,width = 12,height = 5.8)
  res_outlier_record_neg = c("Outlier in negative model:",
                         outlier_sid.neg,
                         out_para_record)
  writeLines(res_outlier_record_neg,paste0(data_clean_out_file,"/02.outlier_record_neg.txt"))
} else {
  message(msg_no("error: The sample id seems not included in this object!"))
}

if(!exists("object_pos_outlier")) {
  object_pos_outlier <- object_pos_mv
}

if(!exists("object_neg_outlier")) {
  object_neg_outlier <- object_neg_mv
}

save(object_pos_outlier,file = paste0(wd,"/massdataset/03.object_pos_outlier.rda"))
save(object_neg_outlier,file = paste0(wd,"/massdataset/03.object_neg_outlier.rda"))


# missing value imputation ------------------------------------------------------------------










