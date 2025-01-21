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
suppressMessages(library(stringr))
handlers(global = TRUE)
handlers("progress")


# functions ---------------------------------------------------------------

re_form_reg = function(string){
  string %>% str_replace_all("\\(","\\\\(") %>% str_replace_all("\\)","\\\\)") %>% str_replace_all("\\-","\\\\-")%>% str_replace_all("\\+","\\\\+")
}
# run_annotation ----------------------------------------------------------

##> parameters

message(msg_run("\nStep2. Annotation filtering"))
temp_para = readxl::read_xlsx(parameters) 
para_anno_filter_tbl = temp_para |> 
  dplyr::filter(process == "annotation filtering") |> 
  dplyr::select(para,Value)

para_anno_filter = para_anno_filter_tbl |> pull(Value) |> setNames(para_anno_filter_tbl$para) |> as.list()

Adduct_pos = para_anno_filter$Adduct_pos

Adduct_neg = para_anno_filter$Adduct_neg

multi_anno = para_anno_filter$multi_anno %>% as.character()
redundancy = para_anno_filter$redundancy %>% as.character()
unclass = para_anno_filter$unclass %>% as.logical()


load("massdataset/07.object_neg_anno.rda")
load("massdataset/07.object_pos_anno.rda")

annotation_tbl_pos <- 
  object_pos_anno %>% extract_annotation_table()
annotation_tbl_neg <- 
  object_neg_anno %>% extract_annotation_table()

# run filtering -----------------------------------------------------------
##> step1 Level 3 filtering by adduct

L3_annotation_neg <- annotation_tbl_neg %>% 
  dplyr::filter(Level == 3) %>% 
  dplyr::filter(stringr::str_detect(Adduct,re_form_reg(Adduct_neg)))
  
L3_annotation_pos <- annotation_tbl_pos %>% 
  dplyr::filter(Level == 3) %>% 
  dplyr::filter(stringr::str_detect(Adduct,re_form_reg(Adduct_pos)))

L12_annotation_neg =  annotation_tbl_neg %>% 
  dplyr::filter(Level != 3)

L12_annotation_pos <- annotation_tbl_pos %>% 
  dplyr::filter(Level != 3)
temp_step1_anno <- rbind(L3_annotation_neg,L3_annotation_pos,L12_annotation_neg,L12_annotation_pos) %>% tibble() %>% 
  mutate(Adduct_level = case_when(
    stringr::str_detect(Adduct,re_form_reg("(M+H)+|(M-H)-")) ~ 1,
    TRUE ~ 2
  ))

##> step2 filter multiple duplications
temp_step2_anno <-
  temp_step1_anno %>% 
  dplyr::group_by(variable_id) %>% 
  dplyr::slice_min(order_by = Level) %>% 
  dplyr::slice_max(order_by = Total.score) %>% 
  dplyr::slice_min(order_by = Adduct_level)

if(multi_anno == "keep highest total score") {
  temp_step2_anno = temp_step2_anno
} else if(multi_anno == "keep the first one") {
  temp_step2_anno = 
    temp_step2_anno %>% dplyr::slice_head(n = 1)
} else {
  temp_step2_anno = temp_step1_anno
}

##> step2 remove redundancy

if(redundancy == "keep highest total score") {
  temp_step3_anno = temp_step2_anno %>% 
    ungroup() %>% 
    group_by(Compound.name) %>% 
    dplyr::slice_min(order_by = Level) %>% 
    dplyr::slice_max(order_by = Total.score) %>% 
    dplyr::slice_min(order_by = Adduct_level)
} else if(redundancy == "keep the first one") {
  temp_step3_anno = 
    temp_step2_anno %>% 
    ungroup() %>% 
    group_by(Compound.name) %>% 
    dplyr::slice_min(order_by = Level) %>% 
    dplyr::slice_max(order_by = Total.score) %>% 
    dplyr::slice_min(order_by = Adduct_level) %>% 
    dplyr::slice_head(n = 1)
} else {
  temp_step3_anno = temp_step2_anno
}


output = list(
  Annotation_Full = rbind(annotation_tbl_pos,annotation_tbl_neg),
  Annotation_clean = temp_step3_anno
)


# KEGG and classification -------------------------------------------------

temp_step4_classification <- 
  temp_step3_anno %>% 
  ungroup() %>% 
  dplyr::select(variable_id,Compound.name,Lab.ID) %>% 
  dplyr::left_join(plantmdb::class.database,by = "Lab.ID")

temp_unclass <- temp_step4_classification %>% 
  dplyr::filter(is.na(superclass))
temp_class <- temp_step4_classification %>% 
  dplyr::filter(!is.na(superclass))
if(isTRUE(unclass)) {
  query1 = temp_unclass %>% dplyr::pull(Compound.name) %>% unique()
  n2cid = MDAtoolkits::mda_get_cid_fast(query = query1,core_num = 10)
  query2 = n2cid %>% dplyr::pull(InChIKey) %>% unique()
  Inchi2class = MDAtoolkits::cfb_crawler(query = query2,delay_max = 0.6,ssl.verifypeer = FALSE)
  n2kegg = MDAtoolkits::mda_name2kegg(query = query1,core_num = 10)
  temp_unclass_final <- 
  temp_unclass %>% 
    dplyr::select(variable_id,Compound.name,Lab.ID) %>% 
    dplyr::left_join(n2cid %>% dplyr::rename("Compound.name" = "query") )%>% 
    dplyr::left_join(Inchi2class) %>% 
    dplyr::left_join(n2kegg) %>% 
    tibble() %>% 
    dplyr::select(colnames(temp_unclass)) %>% 
    mutate(KEGG.ID = case_when(
      KEGG.ID == "" ~ NA,
      TRUE ~ KEGG.ID
    ))
  temp_step4_classification <- rbind(temp_class,temp_unclass_final)
    
}

output = list(
  Annotation_Full = rbind(annotation_tbl_pos,annotation_tbl_neg),
  Annotation_clean = temp_step3_anno,
  Annotation_class = temp_step4_classification
)

writexl::write_xlsx(output,"Feature_annotation/Compound_annotation.xlsx")


