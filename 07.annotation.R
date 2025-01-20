##############################################################
#            Name: One-step metminer
#            Prj: MetMiner pipeline 
#            Assignment: feature annotation
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
  'threads', 't', 2, 'integer', 'threads for peak picking steps',
  'customized_db','d',2,'character','customized database folder path (metid required in .rda format)'
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

message(msg_run("\nStep2. Annotation"))
temp_para = readxl::read_xlsx(parameters) 
para_anno_tbl = temp_para |> 
  dplyr::filter(process == "annotation") |> 
  dplyr::select(para,Value)

para_anno = para_anno_tbl |> pull(Value) |> setNames(para_anno_tbl$para) |> as.list()


ms1.match.ppm = para_anno$ms1.match.ppm |> as.numeric()
ms2.match.ppm = para_anno$ms2.match.ppm |> as.numeric()
mz.ppm.thr = para_anno$mz.ppm.thr |> as.numeric()
ms2.match.tol = para_anno$ms2.match.tol |> as.numeric()
fraction.weight = para_anno$fraction.weight |> as.numeric()
dp.forward.weight = para_anno$dp.forward.weight |> as.numeric()
dp.reverse.weight = para_anno$dp.reverse.weight |> as.numeric()
remove_fragment_intensity_cutoff	 = para_anno$remove_fragment_intensity_cutoff	 |> as.numeric()
rt.match.tol = para_anno$rt.match.tol |> as.numeric()
ce = para_anno$ce |> as.character()
column = para_anno$column |> as.character()
rt.match.weight = para_anno$rt.match.weight |> as.numeric()
ms2.match.weight = para_anno$ms2.match.weight |> as.numeric()
candidate.num = para_anno$candidate.num |> as.numeric()
threads = para_anno$threads |> as.numeric()
total.score.tol = para_anno$total.score.tol |>as.numeric()
ms1.match.weight = para_anno$ms1.match.weight |>as.numeric()
db_name = stringr::str_split(
  string = para_anno$database, pattern = "\\|",n = Inf,simplify = FALSE
) |> unlist()

##> run anno 
buildin_db <- list(
  MoNA = plantmdb::mona_database0.0.4,
  Massbank = plantmdb::massbank_database0.0.4,
  ReSpect = plantmdb::respect_database0.0.1,
  PlaSMA = plantmdb::plasma_database0.0.1,
  KEGG = plantmdb::kegg_plant_database0.0.1,
  KNApSAcK = plantmdb::knapsack_agri_database0.0.1,
  Ath_Cyc = plantmdb::ath_plantcyc.database0.0.1,
  Zma_Cyc = plantmdb::zma_plantcyc_database0.0.1,
  MetaboBASE = plantmdb::metabobase_database0.0.1,
  Orbitrap = plantmdb::orbitrap_database0.0.3
)

# 根据 db_name 筛选内置数据库
if (length(db_name) == 0) {
  buildin_db <- NULL
} else {
  temp_anno_idx <- match(db_name, names(buildin_db))
  buildin_db <- buildin_db[temp_anno_idx]
}

## Customized ms database
if (!is.null(customized_db)) {
  # 捕获自定义数据库目录中的 .rda 文件
  temp_file_name <- tryCatch({
    dir(customized_db, pattern = ".*\\.rda$")
  }, error = function(e) {
    warning("Error reading the specified customized_db directory. Returning NULL.")
    return(NULL)
  })
  
  # 如果没有找到 .rda 文件或发生错误，直接使用内置数据库
  if (is.null(temp_file_name) || length(temp_file_name) == 0) {
    db <- buildin_db
  } else {
    cuz_db <- list() # 初始化自定义数据库列表
    
    for (i in seq_along(temp_file_name)) {
      # 加载 .rda 文件到新环境以避免变量冲突
      env <- new.env()
      load(file = file.path(customized_db, temp_file_name[[i]]), envir = env)
      cuz_db[[i]] <- mget(ls(env), envir = env) # 提取所有加载的对象
    }
    
    # 使用文件名（去除 .rda 后缀）作为数据库名称
    cuz_name <- stringr::str_remove(temp_file_name, "\\.rda$")
    names(cuz_db) <- cuz_name
    
    # 合并自定义数据库和内置数据库
    db <- if (is.null(buildin_db)) cuz_db else c(buildin_db, cuz_db)
  }
} else {
  # 如果没有自定义数据库路径，仅使用内置数据库
  db <- buildin_db
}

dir.create("Feature_annotation/Database",showWarnings = F,recursive = T)
save(db,file = "Feature_annotation/Database/auto_saved.dblist")

# run annotation ----------------------------------------------------------

load("massdataset/06.object_pos_ms2.rda")
load("massdataset/06.object_neg_ms2.rda")
tags = names(db)

# 初始化变量
object_pos_anno <- NULL
object_neg_anno <- NULL

# 使用 map
map(.x = 1:length(tags), .f = function(.x) {
  # 判断当前循环是第几次，并分配初始值
  if (.x == 1) {
    object_pos_temp <- object_pos_ms2
    object_neg_temp <- object_neg_ms2
  } else {
    object_pos_temp <- object_pos_anno
    object_neg_temp <- object_neg_anno
  }
  
  # 执行正极性代谢物注释
  object_pos_anno <<- annotate_metabolites_mass_dataset(
    object = object_pos_temp,
    polarity = "positive",
    database = db[[.x]],
    ms1.match.ppm = ms1.match.ppm,
    ms2.match.ppm = ms2.match.ppm,
    rt.match.tol = rt.match.tol,
    candidate.num = candidate.num,
    column = column,
    threads = threads,
    mz.ppm.thr = mz.ppm.thr,
    ms2.match.tol = ms2.match.tol,
    fraction.weight = fraction.weight,
    dp.forward.weight = dp.forward.weight,
    dp.reverse.weight = dp.reverse.weight,
    remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff,
    ce = ce,
    ms1.match.weight = ms1.match.weight,
    rt.match.weight = rt.match.weight,
    ms2.match.weight = ms2.match.weight,
    total.score.tol = total.score.tol
  )
  
  # 执行负极性代谢物注释
  object_neg_anno <<- annotate_metabolites_mass_dataset(
    object = object_neg_temp,
    polarity = "negative",
    database = db[[.x]],
    ms1.match.ppm = ms1.match.ppm,
    ms2.match.ppm = ms2.match.ppm,
    rt.match.tol = rt.match.tol,
    candidate.num = candidate.num,
    column = column,
    threads = threads,
    mz.ppm.thr = mz.ppm.thr,
    ms2.match.tol = ms2.match.tol,
    fraction.weight = fraction.weight,
    dp.forward.weight = dp.forward.weight,
    dp.reverse.weight = dp.reverse.weight,
    remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff,
    ce = ce,
    ms1.match.weight = ms1.match.weight,
    rt.match.weight = rt.match.weight,
    ms2.match.weight = ms2.match.weight,
    total.score.tol = total.score.tol
  )
})


save(object_neg_anno,file = "massdataset/07.object_neg_anno.rda")
save(object_pos_anno,file = "massdataset/07.object_pos_anno.rda")