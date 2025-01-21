##############################################################
#            Name: One-step metminer
#            Prj: MetMiner pipeline 
#            Assignment: remove wrong samples
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
  'outlier', 'o', 1, 'character','sample id delimed by |'
),byrow = T, ncol = 5)

args = getopt(command)

##> defaults

if (is.null(args$wd)){
  message(msg_no("-w error\n need set working directory"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}

if (is.null(args$outlier)){
  message(msg_no("-w error\n need provide outlier"))
  q(status = 1)
}


wd = args$wd |> as.character()
outlier = args$outlier |>  as.character()

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


# remove wrong samples with large NA rate----------------------------------------------------
load("massdataset/01.object_neg_raw.rda")
load("massdataset/01.object_pos_raw.rda")
outlier_sid = outlier %>% str_split(string = .,pattern = "\\|",n = Inf)  %>% unlist()
object_pos_raw =
  object_pos_raw %>% activate_mass_dataset(what = 'sample_info') %>%
  filter(!sample_id %in% outlier_sid)
object_neg_raw =
  object_neg_raw %>% activate_mass_dataset(what = 'sample_info') %>%
  filter(!sample_id %in% outlier_sid)

save(object_neg_raw,file = "massdataset/01.object_neg_raw.rda")
save(object_pos_raw,file = "massdataset/01.object_pos_raw.rda")