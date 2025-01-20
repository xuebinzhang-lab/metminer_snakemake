##############################################################
#            Name: One-step metminer
#            Prj: MetMiner pipeline 
#            Assignment: raw data quality control
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
  'metadata', 'm', 1, 'character', 'An .csv file of sample information',
  'col_by','c',1,'character', 'pca color by'
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

if (is.null(args$col_by)){
  args$col_by = 'group'
}

wd = args$wd |> as.character()
metadata = args$metadata |> as.character()
threads = args$threads |> as.character()
outlier_method = args$outlier |> as.character()
pca_col_by = args$col_by |> as.character()

setwd(wd)
##> file check
folder_check.neg = paste0(wd,"/MS1/NEG/Result/object")
folder_check.pos = paste0(wd,"/MS1/POS/Result/object")

if(!file.exists(folder_check.neg) | !file.exists(folder_check.pos)) {
  message(msg_no('Please check your working directory, make sure peak picking step has finished'))
  q(status = 1)
}


folder_check_qc = paste0(wd,"/Raw_data_Quanlity_Control/")
if(!file.exists(folder_check_qc)) {dir.create(folder_check_qc,showWarnings = F,recursive = T)}

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


# functions ---------------------------------------------------------------

mass_qc_modified = 
  function (object, path = ".", type = c("html", "pdf", "all"),pca_col_by = 'group',pca_scale = T) 
  {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    type = match.arg(type)
    if (!is(object = object, class2 = "mass_dataset")) {
      stop("obejct should be mass_dataset class.\n")
    }
    options(warn = -1)
    if (length(grep("Report", dir(path))) > 0) {
      output_path = file.path(path, paste("Report", length(grep("Report", 
                                                                dir(path))) + 1, sep = "_"))
    }
    else {
      output_path = file.path(path, "Report")
    }
    message("Get report template.")
    rmarkdown::draft(file = output_path, template = "massqc", 
                     package = "massqc", create_dir = TRUE, edit = FALSE)
    message("Parameters.")
    parameters <- massdataset::extract_process_info(object) %>% 
      lapply(function(x) {
        if (length(x) == 1) {
          massdataset::parse_tidymass_parameter(object = x)
        }
        else {
          x %>% lapply(function(y) {
            massdataset::parse_tidymass_parameter(object = y)
          }) %>% dplyr::bind_rows()
        }
      }) %>% dplyr::bind_rows() %>% dplyr::arrange(time)
    save(parameters, file = file.path(output_path, "parameters.rda"))
    message("Sample infromation.")
    save(object, file = file.path(output_path, "object.rda"))
    if (nrow(object) >= 1000) {
      hex = TRUE
    }
    else {
      hex = FALSE
    }
    plot <- massdataset::show_mz_rt_plot(object = log(object), 
                                         hex = hex)
    ggplot2::ggsave(filename = file.path(output_path, "mz_rt_plot.png"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    ggplot2::ggsave(filename = file.path(output_path, "mz_rt_plot.pdf"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    message("Missing values in dataset.")
    plot <- massdataset::show_missing_values(object = object, 
                                             show_row_names = ifelse(nrow(object) < 20, TRUE, FALSE), 
                                             show_column_names = ifelse(ncol(object) < 20, TRUE, FALSE), 
                                             percentage = TRUE, return_as_ggplot = TRUE)
    ggplot2::ggsave(filename = file.path(output_path, "mv_plot.png"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    ggplot2::ggsave(filename = file.path(output_path, "mv_plot.pdf"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    message("Missing values in all variables.")
    if (nrow(object) > 10000) {
      plot <- massdataset::show_variable_missing_values(object = object[seq_len(10000), 
      ], order_by = "rt", show_x_text = ifelse(nrow(object) < 
                                                 20, TRUE, FALSE), show_x_ticks = ifelse(nrow(object) < 
                                                                                           20, TRUE, FALSE), percentage = TRUE) + labs(title = "Only first 10000 variables.")
    }
    else {
      plot <- massdataset::show_variable_missing_values(object = object, 
                                                        order_by = "rt", show_x_text = ifelse(nrow(object) < 
                                                                                                20, TRUE, FALSE), show_x_ticks = ifelse(nrow(object) < 
                                                                                                                                          20, TRUE, FALSE), percentage = TRUE)
    }
    ggplot2::ggsave(filename = file.path(output_path, "variable_mv_plot.png"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    ggplot2::ggsave(filename = file.path(output_path, "variable_mv_plot.pdf"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    message("Missing values in all samples.")
    plot <- massdataset::show_sample_missing_values(object = object, 
                                                    color_by = "class", order_by = "injection.order", percentage = TRUE)
    ggplot2::ggsave(filename = file.path(output_path, "sample_mv_plot.png"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    ggplot2::ggsave(filename = file.path(output_path, "sample_mv_plot.pdf"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    message("RSD for variables.")
    if (sum(object@sample_info$class == "QC") >= 3) {
      plot1 <- massqc_rsd_plot(object = object %>% massdataset::activate_mass_dataset(what = "sample_info") %>% 
                                 dplyr::filter(class == "QC"), order_by = "rt", show_x_text = ifelse(nrow(object) < 
                                                                                                       20, TRUE, FALSE), show_x_ticks = ifelse(nrow(object) < 
                                                                                                                                                 20, TRUE, FALSE)) + labs(title = "All QC samples")
      plot2 = massqc_cumulative_rsd_plot(object = object %>% 
                                           massdataset::activate_mass_dataset(what = "sample_info") %>% 
                                           dplyr::filter(class == "QC"), rsd_cutoff = 30) + 
        labs(title = "All QC samples")
    }
    else {
      plot1 <- massqc_rsd_plot(object = object %>% massdataset::activate_mass_dataset(what = "sample_info") %>% 
                                 dplyr::filter(class == "Subject"), order_by = "rt", 
                               show_x_text = ifelse(nrow(object) < 20, TRUE, FALSE), 
                               show_x_ticks = ifelse(nrow(object) < 20, TRUE, FALSE)) + 
        labs(title = "All Subject samples")
      plot2 <- massqc_cumulative_rsd_plot(object = object %>% 
                                            massdataset::activate_mass_dataset(what = "sample_info") %>% 
                                            dplyr::filter(class == "Subject"), rsd_cutoff = 30) + 
        labs(title = "All Subject samples")
    }
    ggplot2::ggsave(filename = file.path(output_path, "variable_rsd.pdf"), 
                    plot = plot1, width = 8, height = 6, dpi = 600)
    ggplot2::ggsave(filename = file.path(output_path, "variable_cumulative_rsd.pdf"), 
                    plot = plot2, width = 8, height = 6, dpi = 600)
    plot <- plot1 + plot2 + patchwork::plot_layout(nrow = 1)
    ggplot2::ggsave(filename = file.path(output_path, "rsd.png"), 
                    plot = plot, width = 12, height = 6, dpi = 600)
    message("Intensity for all the variables in all samples.")
    if (sum(object@sample_info$class == "QC") >= 3) {
      plot <- massqc_sample_boxplot(object = log(object) %>% 
                                      massdataset::activate_mass_dataset(what = "sample_info") %>% 
                                      dplyr::filter(class == "QC"), color_by = ifelse(any(colnames(object@sample_info) == 
                                                                                            "batch"), "batch", "class"))
    }
    else {
      plot <- massqc_sample_boxplot(object = log(object)[, 
                                                         seq_len(30)], color_by = ifelse(any(colnames(object@sample_info) == 
                                                                                               "batch"), "batch", "class")) + labs(title = ifelse(ncol(object) > 
                                                                                                                                                    30, "First 30 samples", ""))
    }
    if (ncol(object) > 100) {
      plot <- plot + ggplot2::theme(axis.text.x = element_blank(), 
                                    axis.ticks.x = element_blank())
    }
    else {
      plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                                        hjust = 1, vjust = 1))
    }
    ggplot2::ggsave(filename = file.path(output_path, "sample_box.pdf"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    ggplot2::ggsave(filename = file.path(output_path, "sample_box.png"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    message("Sample correlation.")
    if (sum(object@sample_info$class == "QC") >= 3) {
      plot <- massqc_sample_correlation(object = object %>% 
                                          massdataset::activate_mass_dataset(what = "sample_info") %>% 
                                          dplyr::filter(class == "QC"), method = "square")
    }
    else {
      plot <- massqc_sample_correlation(object = object[seq_len(100)], 
                                        method = "square") + labs(title = ifelse(ncol(object) > 
                                                                                   100, "First 100 samples", ""))
    }
    ggplot2::ggsave(filename = file.path(output_path, "sample_correlation.png"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    ggplot2::ggsave(filename = file.path(output_path, "sample_correlation.pdf"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    message("PCA score plot.")
    sample_info <- extract_sample_info(object)
    if (any(colnames(sample_info) == "batch")) {
      sample_info$batch <- as.character(sample_info$batch)
      methods::slot(object, "sample_info") <- sample_info
    }
    temp_obj <- object %>% +1 %>% log(2)
    if(isTRUE(pca_scale)) {
      temp_obj = object %>% +1 %>% log(2) %>% scale()
    }
    plot <- temp_obj %>% massqc_pca(color_by = pca_col_by)
    message("PCA done.")
    ggplot2::ggsave(filename = file.path(output_path, "pca.pdf"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    ggplot2::ggsave(filename = file.path(output_path, "pca.png"), 
                    plot = plot, width = 8, height = 6, dpi = 600)
    message("Render report.")
    if (type == "html" | type == "all") {
      rmarkdown::render(file.path(output_path, "massqc.template.Rmd"), 
                        rmarkdown::html_document())
      file.rename(from = file.path(output_path, "massqc.template.html"), 
                  to = file.path(output_path, "massqc_report.html"))
    }
    if (type == "pdf" | type == "all") {
      rmarkdown::render(file.path(output_path, "massqc.template.Rmd"), 
                        rmarkdown::pdf_document())
      file.rename(from = file.path(output_path, "massqc.template.pdf"), 
                  to = file.path(output_path, "massqc_report.pdf"))
    }
    message("Remove some files.")
    file = dir(output_path)
    remove_file = grep("png|Rmd|parameters|rda", file, value = TRUE)
    unlink(x = file.path(output_path, remove_file), recursive = TRUE, 
           force = TRUE)
  }

# run_qc ------------------------------------------------------------------

##> positive model

load(folder_check.pos)

sample_info <- read.csv(metadata)
object <- object %>% activate_mass_dataset('sample_info') %>% 
  dplyr::select(sample_id) %>% 
  dplyr::left_join(sample_info)
dir.create(paste0(folder_check_qc,"/POS/"),showWarnings = F,recursive = T)
mass_qc_modified(
  object = object,path = paste0(folder_check_qc,"/POS/"),type = 'html',pca_col_by = pca_col_by,
)

load(folder_check.neg)

sample_info <- read.csv(metadata)
object <- object %>% activate_mass_dataset('sample_info') %>% 
  dplyr::select(sample_id) %>% 
  dplyr::left_join(sample_info)
dir.create(paste0(folder_check_qc,"/NEG/"),showWarnings = F,recursive = T)
mass_qc_modified(
  object = object,path = paste0(folder_check_qc,"/NEG/"),type = 'html',pca_col_by = pca_col_by,
)