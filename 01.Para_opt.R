##############################################################
#            Name: One-step metminer
#            Prj: MetMiner pipeline 
#            Assignment: parameter optimize
#            Author: Xiao Wang (shawnwang2016@126.com)
#            Date: Dec 30, 2024
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
  'massSDrange', 'm', 1, 'integer','massSDrange, default = 2',
  'smooth', 's', 1, 'integer', 'smooth, default = 0',
  'cutoff', 'c', 1, 'double', 'cutoff, default = 0.95',
  'filenum','n',1, 'integer', 'default 3',
  'thread', 't', 1, 'integer', 'default 5'
),byrow = T, ncol = 5)

args = getopt(command)

##> defaults

if (is.null(args$wd)){
  message(msg_no("-w error\n need set working directory"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}

if (is.null(args$massSDrange)){
  args$massSDrange = 2
}

if (is.null(args$smooth)){
  args$smooth = 0
}

if (is.null(args$cutoff)){
  args$cutoff = 0.95
}

if (is.null(args$thread)){
  args$thread = 5
}

if (is.null(args$filenum)){
  args$filenum = 3
}

wd = args$wd |> as.character()
massSDrange = args$massSDrange |> as.numeric()
smooth = args$smooth |> as.numeric()
cutoff = args$cutoff |> as.numeric()
thread = args$thread |> as.numeric()
filenum = args$filenum |> as.numeric()
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

# part1. calculate ppm-cutoff --------------------------------------------------------
message(msg_run("Step1. find best ppm cutoff"))

##> record runing time
with_progress({
  p <- progressor(along = 1:2)
  ##> positive model
  p("Processing Positive Model (step1)") 
  res_part1_pos = paramounter_part1(
    directory = paste0(wd,"/MS1/POS/QC/"),
    smooth = smooth,
    cutoff = cutoff,
    thread = thread,
    filenum = filenum
  )
  
  ##> positive model
  p("Processing Negative Model (step1)") 
  res_part1_neg = paramounter_part1(
    directory = paste0(wd,"/MS1/NEG/QC/"),
    smooth = smooth,
    cutoff = cutoff,
    thread = thread,
    filenum = filenum
  )
})

p_pos = res_part1_pos$plot + ggtitle("Positive")
p_neg = res_part1_neg$plot + ggtitle("Negative")
p = p_pos+p_pos

output_dir = paste0(wd,"/Results/parameter_optimize/")

if(!file.exists(output_dir)) {
  dir.create(path = output_dir,showWarnings = F,recursive = T)
}
ggsave(filename = paste0(output_dir,"ppmcut.pdf"),plot = p,width = 12,height = 5.8)
save.image(paste0(output_dir,"01.para_opt.rda"))
# part2. calculate new parameters --------------------------------------------------------

message(msg_run("Step2. find best parameters"))

##> record runing time
start_time <- Sys.time() 
with_progress({
  p <- progressor(along = 1:2)
  ##> positive model
  p("Processing Positive Model (step1)") 
  res_part2_pos = paramounter_part2(
    directory = paste0(wd,"/MS1/POS/QC/"),
    massSDrange = massSDrange,
    smooth = smooth,
    cutoff = cutoff,
    thread = thread,
    filenum = filenum,
    ppmCut = res_part1_pos$ppmCut %>% as.numeric()
  )
  
  ##> positive model
  p("Processing Negative Model (step1)") 
  res_part2_neg = paramounter_part2(
    directory = paste0(wd,"/MS1/NEG/QC/"),
    smooth = smooth,
    cutoff = cutoff,
    thread = thread,
    filenum = filenum,
    massSDrange = massSDrange,
    ppmCut = res_part1_neg$ppmCut %>% as.numeric()
  )
})


writexl::write_xlsx(list(
  pos = res_part2_pos,
  neg = res_part2_neg
),path = paste0(output_dir,"parameters.xlsx"))
save.image(paste0(output_dir,"01.para_opt.rda"))
message(msg_yes("\nAll done!"))