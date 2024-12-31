##############################################################
#            Name: One-step metminer
#            Prj: MetMiner pipeline 
#            Assignment: peak picking
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
  'peak_picking', 'o',2, 'character', 'An excel file which paramounter_part2 generated',
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

if (is.null(args$peak_picking)){
  message(msg_warning("-w error\no peak picking optimize parameters upload, use defalut."))
  args$peak_picking = "default"
}

if (is.null(args$threads)){
  args$threads = NA
}

wd = args$wd |> as.character()
parameters = args$parameters |>  as.character()
peak_picking = args$peak_picking |> as.character()
threads = args$threads |> as.character()
setwd(wd)
##> file check
folder_check.neg = paste0(wd,"/MS1/NEG/")
folder_check.pos = paste0(wd,"/MS1/POS/")

if(!file.exists(folder_check.neg) | !file.exists(folder_check.pos)) {
  message(msg_no('Please check your working directory, make sure ms1 files were correctly set'))
  q(status = 1)
}

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



# parameters -------------------------------------------------------------


temp_para = readxl::read_xlsx(parameters) %>% 
  dplyr::filter(process == 'peak picking') %>% 
  dplyr::select(para,Value)
if(peak_picking != 'default') {
  pp_para.pos <- readxl::read_xlsx(peak_picking) %>% 
    dplyr::select(para,Value) %>% setNames(c("para","Value_opt")) %>% 
    mutate(Value_opt = as.character(Value_opt))
  pp_para.neg <- readxl::read_xlsx(peak_picking,sheet = 2) %>% 
    dplyr::select(para,Value) %>% setNames(c("para","Value_opt")) %>% 
    mutate(Value_opt = as.character(Value_opt))
  temp_para.pos <-
    temp_para %>% left_join(pp_para.pos,by = "para") %>% 
    dplyr::mutate(Value = dplyr::case_when(
      !is.na(Value_opt) ~ Value_opt,
      TRUE ~ Value
    )) %>% 
    dplyr::select(para,Value)
  temp_para.neg <-
    temp_para %>% left_join(pp_para.neg,by = "para") %>% 
    dplyr::mutate(Value = dplyr::case_when(
      !is.na(Value_opt) ~ Value_opt,
      TRUE ~ Value
    )) %>% 
    dplyr::select(para,Value)
} else {
  temp_para.pos <- temp_para;
  temp_para.neg <- temp_para
}

para.pos = temp_para.pos %>% pull(Value) %>% setNames(temp_para$para) %>% as.list()
para.neg = temp_para.neg %>% pull(Value) %>% setNames(temp_para$para) %>% as.list()

if(!is.na(threads)) {
  para.neg$threads <- as.character(threads);
  para.pos$threads <- as.character(threads);
}

# run peak picking progress -----------------------------------------------

message(msg_run('Start peak picking progress in positive model'))

process_data(
  path = paste0(wd,"/MS1/POS/"),
  polarity = "positive",
  ppm = para.pos$ppm %>% as.numeric(),
  peakwidth = c(para.pos$p_min %>% as.numeric(),para.pos$p_max %>% as.numeric()),
  snthresh = para.pos$snthresh %>% as.numeric(),
  prefilter = c(para.pos$pre_left %>% as.numeric(),para.pos$pre_right %>% as.numeric()),
  fitgauss = para.pos$fitgauss %>% as.logical(),
  integrate = para.pos$integrate %>% as.numeric(),
  mzdiff = para.pos$mzdiff %>% as.numeric(),
  noise = para.pos$noise %>% as.numeric(),
  threads = para.pos$threads %>% as.numeric(),
  binSize = para.pos$binSize %>% as.numeric(),
  bw = para.pos$bw %>% as.numeric(),
  output_tic = para.pos$out_put_peak %>% as.logical(),
  output_bpc = para.pos$out_put_peak %>% as.logical(),
  output_rt_correction_plot = para.pos$out_put_peak %>% as.logical(),
  min_fraction = para.pos$min_fraction %>% as.numeric(),
  fill_peaks = para.pos$fill_peaks %>% as.logical()
)

message(msg_run('Start peak picking progress in positive model'))

process_data(
  path = paste0(wd,"/MS1/NEG/"),
  polarity = "negative",
  ppm = para.neg$ppm %>% as.numeric(),
  peakwidth = c(para.neg$p_min %>% as.numeric(),para.neg$p_max %>% as.numeric()),
  snthresh = para.neg$snthresh %>% as.numeric(),
  prefilter = c(para.neg$pre_left %>% as.numeric(),para.neg$pre_right %>% as.numeric()),
  fitgauss = para.neg$fitgauss %>% as.logical(),
  integrate = para.neg$integrate %>% as.numeric(),
  mzdiff = para.neg$mzdiff %>% as.numeric(),
  noise = para.neg$noise %>% as.numeric(),
  threads = para.neg$threads %>% as.numeric(),
  binSize = para.neg$binSize %>% as.numeric(),
  bw = para.neg$bw %>% as.numeric(),
  output_tic = para.neg$out_put_peak %>% as.logical(),
  output_bpc = para.neg$out_put_peak %>% as.logical(),
  output_rt_correction_plot = para.neg$out_put_peak %>% as.logical(),
  min_fraction = para.neg$min_fraction %>% as.numeric(),
  fill_peaks = para.neg$fill_peaks %>% as.logical()
)

##> save mass_dataset

mass_dataset_dir = paste0(wd,"/massdataset/")

if(!file.exists(mass_dataset_dir)) {dir.create(mass_dataset_dir,showWarnings = F,recursive = T)}
load(paste0(wd,"/MS1/POS/Result/object"))
object_pos_raw <- object
save(object_pos_raw,file = paste0(mass_dataset_dir,"/object_pos_raw.rda"))
load(paste0(wd,"/MS1/NEG/Result/object"))
object_neg_raw <- object
save(object_neg_raw,file = paste0(mass_dataset_dir,"/object_neg_raw.rda")) 

message(msg_yes(paste0("\nPeak picking step finish!\nCheck your result at:",mass_dataset_dir)))
