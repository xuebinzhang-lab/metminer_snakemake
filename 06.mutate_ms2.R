##############################################################
#            Name: One-step metminer
#            Prj: MetMiner pipeline 
#            Assignment: Add ms2
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
  'parameters', 's', 1, 'character','An excel file with default parameters'
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
suppressMessages(library(furrr))

handlers(global = TRUE)
handlers("progress")


# functions ---------------------------------------------------------------

progressr::handlers("progress")
mda_mutate_ms2 = function (object, 
                           column = c("rp", "hilic"), 
                           polarity = c("positive", "negative"), 
                           ms1.ms2.match.mz.tol = 15, 
                           ms1.ms2.match.rt.tol = 30,
                           path = ".",
                           n_session = 30) 
{
  mda_read_ms2 <- function(x){
    oplan = future::plan("multisession", workers = n_session)
    on.exit(future::plan(oplan), add = TRUE)
    p <- progressr::progressor(steps = length(x))
    round_n = length(x)
    ms2_data <- furrr::future_map(.x = 1:length(x), function(.x) {
      p(sprintf("%s",paste0('mgf file :',ms2_data_name[.x],"(",.x,"/",round_n,")")))
      temp_ms2_type <- stringr::str_split(string = x[.x], 
                                          pattern = "\\.")[[1]]
      temp_ms2_type <- temp_ms2_type[length(temp_ms2_type)]
      if (temp_ms2_type == "mzXML" | temp_ms2_type == "mzxml") {
        data <- masstools::read_mzxml(file = x[.x])
        data <- convert_ms2_mzxml2mgf(data)
      }
      if (temp_ms2_type == "mzML" | temp_ms2_type == "mzml") {
        data <- masstools::read_mzxml(file = x[.x])
        data <- convert_ms2_mzxml2mgf(data)
      }
      if (temp_ms2_type == "mgf") {
        data <- suppressMessages(masstools::read_mgf(file = x[.x]))
      }
      return(data)
    })
    return(ms2_data)
  }
  
  mda_ms2_reform <- function(x1) {
    
    oplan = future::plan(
      list(
        future::tweak(
          future::multiprocess, 
          workers = 8), 
        future::tweak(
          future::multicore,
          workers = 10)
      )
    )
    
    on.exit(future::plan(oplan), add = TRUE)
    p <- progressr::progressor(steps = length(x1))
    round_n = length(x1)
    out = furrr::future_map(.x = 1:length(x1), .f = function(.x) {
      p(sprintf("%s",paste0('ms2 name :',ms2_data_name[.x],"(",.x,"/",round_n,")")))
      ms2_data_tmp = ms2_data[[.x]]
      tmp_out <- furrr::future_map( 1:length(ms2_data_tmp), .f = function(.y) {
        ms2_detail_tmp = ms2_data_tmp[[.y]]
        info <- ms2_detail_tmp$info
        info <- data.frame(name = paste("mz", info[1],
                                        "rt", info[2], sep = ""), mz = info[1], rt = info[2],
                           file = ms2_data_name[.x], stringsAsFactors = FALSE,
                           check.names = FALSE)
        rownames(info) <- NULL
        ms2_detail_tmp$info <- info
        return(ms2_detail_tmp)
      })
      return(tmp_out)
    })
    return(out)
  }
  
  check_object_class(object = object, class = "mass_dataset")
  column = match.arg(column)
  polarity = match.arg(polarity)
  variable_info = object@variable_info
  ms2_list = list.files(path = path, pattern = "mgf|mzXML|mzxml|mzML|mzml", 
                        all.files = TRUE, full.names = TRUE, recursive = TRUE)
  if (length(ms2_list) == 0) {
    stop("No MS2 in ", path)
  }
  ms2_data_name = basename(ms2_list)
  
  
  
  message(crayon::green("Reading mgf files..."))
  
  progressr::with_progress({
    ms2_data <- mda_read_ms2(x = ms2_list)
  })
  
  names(ms2_data) <- ms2_data_name
  
  message(crayon::green("Reform ms2 data"))
  
  progressr::with_progress({
    ms2_data <- mda_ms2_reform(x1 = ms2_data)
  })

  ms2_data <- do.call(what = c, args = ms2_data)
  ms1.info <- lapply(ms2_data, function(x) {
    x[[1]]
  }) %>% dplyr::bind_rows()
  rownames(ms1.info) <- NULL
  duplicated.name <- unique(ms1.info$name[duplicated(ms1.info$name)])
  if (length(duplicated.name) > 0) {
    lapply(duplicated.name, function(x) {
      ms1.info$name[which(ms1.info$name == x)] <- paste(x, 
                                                        c(seq_len(sum(ms1.info$name == x))), sep = "_")
    })
  }
  ms2.info <- lapply(ms2_data, function(x) {
    x[[2]]
  })
  names(ms2.info) <- ms1.info$name
  match.result <- masstools::mz_rt_match(data1 = variable_info[, 
                                                               c(2, 3)], data2 = ms1.info[, c(2, 3)], mz.tol = ms1.ms2.match.mz.tol, 
                                         rt.tol = ms1.ms2.match.rt.tol, rt.error.type = "abs")
  if (is.null(match.result)) {
    message(crayon::red("No variable are matched with MS2 spectra."))
    return(object)
  }
  if (nrow(match.result) == 0) {
    message(crayon::red("No variable are matched with MS2 spectra."))
    return(object)
  }
  message(crayon::green(length(unique(match.result[, 1])), 
                        "out of", nrow(variable_info), "variable have MS2 spectra."))
  message(crayon::green("Selecting the most intense MS2 spectrum for each peak..."))
  temp.idx <- unique(match.result[, 1])
  match.result <- lapply(temp.idx, function(idx) {
    idx2 <- match.result[which(match.result[, 1] == idx), 
                         2]
    if (length(idx2) == 1) {
      return(c(idx, idx2))
    }
    else {
      temp.ms2.info <- ms2.info[idx2]
      return(c(idx, idx2[which.max(unlist(lapply(temp.ms2.info, 
                                                 function(y) {
                                                   y <- y[order(y[, 2], decreasing = TRUE), , 
                                                          drop = FALSE]
                                                   if (nrow(y) > 5) y <- y[seq_len(5), ]
                                                   sum(y[, 2])
                                                 })))]))
    }
  })
  match.result <- do.call(rbind, match.result) %>% as.data.frame()
  colnames(match.result) <- c("Index1", "Index2")
  match.result <- data.frame(match.result, variable_info$variable_id[match.result$Index1], 
                             ms1.info$name[match.result$Index2], stringsAsFactors = FALSE, 
                             check.names = FALSE)
  colnames(match.result) <- c("Index1.ms1.data", "Index.ms2.spectra", 
                              "MS1.peak.name", "MS2.spectra.name")
  ms1.info <- ms1.info[unique(match.result[, 2]), , drop = FALSE]
  ms2.info <- ms2.info[unique(match.result[, 2])]
  match.result$Index.ms2.spectra <- match(match.result$MS2.spectra.name, 
                                          ms1.info$name)
  ms2_data = new(Class = "ms2_data", column = column, polarity = polarity, 
                 variable_id = match.result$MS1.peak.name, ms2_spectrum_id = match.result$MS2.spectra.name, 
                 ms2_mz = ms1.info$mz[match.result$Index.ms2.spectra], 
                 ms2_rt = ms1.info$rt[match.result$Index.ms2.spectra], 
                 ms2_file = ms1.info$file[match.result$Index.ms2.spectra], 
                 ms2_spectra = ms2.info[match.result$Index.ms2.spectra], 
                 mz_tol = ms1.ms2.match.mz.tol, rt_tol = ms1.ms2.match.rt.tol)
  if (length(object@ms2_data) == 0) {
    name = paste(sort(ms2_data_name), collapse = ";")
    ms2_data = list(name = ms2_data)
    names(ms2_data) = name
    object@ms2_data = ms2_data
  }
  else {
    name = paste(sort(ms2_data_name), collapse = ";")
    if (any(names(object@ms2_data) == name)) {
      object@ms2_data[[match(name, names(object@ms2_data))]] = ms2_data
    }
    else {
      object@ms2_data = c(object@ms2_data, ms2_data)
      names(object@ms2_data)[length(names(object@ms2_data))] = name
    }
  }
  process_info = object@process_info
  parameter <- new(Class = "tidymass_parameter", pacakge_name = "massdataset", 
                   function_name = "mutate_ms2()", parameter = list(column = column, 
                                                                    polarity = polarity, ms1.ms2.match.mz.tol = ms1.ms2.match.mz.tol, 
                                                                    ms1.ms2.match.rt.tol = ms1.ms2.match.rt.tol, path = path), 
                   time = Sys.time())
  if (all(names(process_info) != "mutate_ms2")) {
    process_info$mutate_ms2 = parameter
  }
  else {
    process_info$mutate_ms2 = c(process_info$mutate_ms2, 
                                parameter)
  }
  object@process_info = process_info
  return(object)
}

# parameters -----------------------------------------------------------------



message(msg_run("\nStep1. Add MS2"))
temp_para = readxl::read_xlsx(parameters) 
para_mv_tbl = temp_para |> 
  dplyr::filter(process == "mutate_ms2") |> 
  dplyr::select(para,Value)
para_mv = para_mv_tbl |> pull(Value) |> setNames(para_mv_tbl$para) |> as.list()

column = para_mv$column |> as.character()
ms1.ms2.match.mz.tol = para_mv$ms1.ms2.match.mz.tol |> as.numeric()
ms1.ms2.match.rt.tol = para_mv$ms1.ms2.match.rt.tol |> as.numeric()
n_session = para_mv$n_session |> as.numeric()


# load object -------------------------------------------------------------
load("massdataset/05.object_pos_norm.rda")
load("massdataset/05.object_neg_norm.rda")

object_pos_ms2 <-
  object_pos_norm %>%
  mda_mutate_ms2(
    object = .,
    column = column,
    polarity = "positive",
    ms1.ms2.match.mz.tol = ms1.ms2.match.mz.tol,
    ms1.ms2.match.rt.tol = ms1.ms2.match.rt.tol,
    path = "MS2/POS/",
    n_session = n_session
  )
save(object_pos_ms2,file = "massdataset/06.object_pos_ms2.rda")

object_neg_ms2 <-
  object_neg_norm %>%
  mda_mutate_ms2(
    object = .,
    column = column,
    polarity = "negative",
    ms1.ms2.match.mz.tol = ms1.ms2.match.mz.tol,
    ms1.ms2.match.rt.tol = ms1.ms2.match.rt.tol,
    path = "MS2/NEG/",
    n_session = n_session
  )
save(object_neg_ms2,file = "massdataset/06.object_neg_ms2.rda")