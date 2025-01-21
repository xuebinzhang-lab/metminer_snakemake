                  ##############################################################
#            Name: One-step metminer
#            Prj: MetMiner pipeline 
#            Assignment: DAM analysis
#            Author: Xiao Wang (shawnwang2016@126.com)
#            Date: Dec 31, 2024
##############################################################
options(stringsAsFactors = F)
options(warn = -1)
suppressMessages(if (!require('getopt')) BiocManager::install('getopt'))
suppressMessages(if (!require('crayon')) BiocManager::install('crayon'))
# functions ---------------------------------------------------------------

msg_yes = green$bold$italic
msg_no = red$bold$italic
msg_run = blue$bold$italic$underline
msg_warning = yellow$bold$italic


# Setting R script parameters ---------------------------------------------

command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'wd', 'w', 1, 'character', 'working directory',
  'compare_group', 'c', 1, 'character','input group file: \n                     1. first column is sample names correspounding to peak matrix, \n                     2. second column is sample group of samples',
  'kegg_database','k', 1,'character', 'kegg database',
  'parameters', 's', 1, 'character','An excel file with default parameters'
),byrow = T, ncol = 5)
args = getopt(command)


## help information
if (!is.null(args$help)) {
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status=1)
}

## default value
if (is.null(args$wd)){
  message(msg_no("Need working directory! please read help carefully!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}
if (is.null(args$compare_group)){
  message(msg_no("Need group file! please read help carefully!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}


suppressMessages(library(MDAtoolkits))
suppressMessages(library(tidyverse))
suppressMessages(library(ropls))
suppressMessages(library(PCAtools))
suppressMessages(library(devtools))
suppressMessages(library(circlize))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ggprism))
suppressMessages(library(clusterProfiler))
suppressMessages(library(splitstackshape))
suppressMessages(library(rstatix))
suppressMessages(library(bruceR))
suppressMessages(library(tidymass))
# nameing parameters ------------------------------------------------------
wd = args$wd |> as.character()
compare_group = args$compare_group |>  as.character()
annotation_file = args$annotation_file |> as.character()
kegg_database = args$kegg_database |> as.character()
parameters = args$parameters |> as.character()
load(kegg_database)
setwd(wd)

load("massdataset/08.object_merge.rda")
sample_info <- object_merge %>% extract_sample_info() %>% 
  dplyr::select(sample_id,group) %>% distinct()
compare_groups = import(compare_group)

temp_para = readxl::read_xlsx(parameters) 
para_dam_tbl = temp_para |> 
  dplyr::filter(process == "DAM") |> 
  dplyr::select(para,Value)

para_dam = para_dam_tbl |> pull(Value) |> setNames(para_dam_tbl$para) |> as.list()

VIP_cut <- para_dam$VIP %>% as.numeric()
pval_cut <- para_dam$pvalue %>% as.numeric()
qval_cut <- para_dam$qvalue %>% as.numeric()
log2fc_cut <- para_dam$log2fc %>% as.numeric()
test.method <- para_dam$test.method %>% as.character()
pls.method <- para_dam$pls.method %>% as.character()
# run ---------------------------------------------------------------------

select = dplyr::select
filter = dplyr::filter
left_join = dplyr::left_join
unique = base::unique
## get peak area file ==============
message(msg_run("Step1. File format check ... "))

peak_raw <- object_merge %>% extract_expression_data() %>% rownames_to_column('CompoundID')
## check whether chr in peak file
peak_check = peak_raw %>% pivot_longer(
  !CompoundID,names_to = "sample",values_to = "peak"
)
check_tmp1 = class(peak_check$peak)
if (check_tmp1 != "numeric") {
  suppressWarnings({
    check_tmp2 = peak_check %>%
      dplyr::mutate(peak = as.numeric(peak)) %>%
      dplyr::filter(is.na(peak))
  })
  
  error_compound = unique(check_tmp2$CompoundID)
  message(msg_no("Peak file check: Failed! The following cells has non-numeric charaters!\n",paste0(capture.output(as.data.frame(check_tmp2)), collapse = "\n")))
  q(status = 1)
}

## check sample number
if(ncol(peak_raw)< 4) {
  message(msg_no("Peak file check: Failed! \nSample number smaller than 4, we need at least two samples, each sample have at least two biological repeats!"))
  q(status = 1)
}
## NA value check
NA_check =  nrow(peak_raw[is.na(peak_raw[,-1]),-1])
if (NA_check == 0) {
  message(msg_yes("Peak file check: Pass!"))
} else {
  message(msg_no("Peak file check: Failed! \nNA value was not allowed in peak area matrix, please check your input file!"))
  q(status = 1)
}

##  get group file =========
group_df = sample_info %>% setNames(c("Sample","Group"))
## check sample number
if(ncol(group_df) != 2) {
  message(msg_no("Group file check: Failed! \nGroup file ONLY allowed two columns, \n   1. the 1st column is sampleIDs correspounding to the normalization peak file.\n   2. the 2nd column is group labels of the samples."))
  q(status = 1)
}
## NA value check
NA_check =  nrow(group_df[is.na(group_df),-1])
if (is.null(NA_check)) {
  message(msg_yes("Group file check: Pass!" ))
} else {
  message(msg_no("Group file check: Failed! \nBlank was not allowed in input file, maybe at the bottom of the input file have blank rows, remove them under softwares like textmate or vscode not recommand in excel, please check your input file!"))
  q(status = 1)
}

if(any(isFALSE(compare_groups %>% pull(left) %in% group_df$Group)) | any(isFALSE(compare_groups %>% pull(right)  %in% group_df$Group))) {
  unlink(store_path, recursive = TRUE)
  message(msg_no("Error: Wrong group labels, please check the prameter\n<-a ...> and <-b ...>"))
  q(status = 1)
}

##  get annotation file =========
annotation <- import("Feature_annotation/Compound_annotation.xlsx",sheet = 2) %>% 
  dplyr::select(variable_id,Compound.name)
classification <- import("Feature_annotation/Compound_annotation.xlsx",sheet = 3)
annot_df <- left_join(annotation,classification) %>% 
  dplyr::rename("CompoundID" = "variable_id")
## check sample number
if(ncol(annot_df) < 2) {
  message(msg_no("Compound annotation file should be more than two columns, with the 1st column as Compound ID. Please check your file"))
  q(status = 1)
} else {
  message(msg_yes("Annotation file check: Pass!" ))
}


# PCA analysis ------------------------------------------------------------
library(furrr)
future::plan("multicore", workers = 40)
furrr::future_map(.x = 1:nrow(compare_groups),.f = function(.x) {
  left = compare_groups[.x,1] %>% as.character()
  right = compare_groups[.x,2] %>% as.character()
  message(msg_run("Step2. Start PCA analysis ..."))
  store_path = paste0("DAM/",left,"_vs_",right,"/")
  dir.create(store_path,showWarnings = F,recursive = T)
  left_index = group_df %>%
    filter(Group == left) %>%
    select(Sample) %>%
    .$Sample
  l_num = length(left_index)
  right_index = group_df %>%
    filter(Group == right) %>%
    select(Sample) %>%
    .$Sample
  r_num = length(right_index)
  pca_mat = peak_raw %>%
    column_to_rownames("CompoundID") %>%
    select(all_of(left_index),all_of(right_index))
  
  colnames(pca_mat) = c(paste0(left,"_rep",str_pad(string = c(1:l_num),width = 2,side = "left",pad = "0")),
                        paste0(right,"_rep",str_pad(string = c(1:r_num),width = 2,side = "left",pad = "0"))) %>% gsub(" ","",.)
  
  pca_group = data.frame(
    row.names = colnames(pca_mat),
    group = rep(c(left,right),times = c(l_num,r_num))
  )
  
  message(msg_run("running pca."))
  pca_result = PCAtools::pca(
    mat = pca_mat,
    metadata = pca_group,
    scale = T,removeVar = 0.1
  )
  if(sum(l_num,r_num) < 20) {p_size = 3} else {p_size = 1}
  
  message(msg_run("running biplot"))
  
  pdf(file = paste0(store_path,"/01.PCA1.pdf"),height = 8,width = 8)
  suppressMessages(
    biplot(
      pcaobj = pca_result,
      colby = "group",
      hline = 0,vline = 0,pointSize = p_size,
      title = "biplot",
      subtitle = paste0(left, " vs " ,right),encircle = T,legendPosition = "top"
    )
  )
  
  invisible(dev.off())
  pdf(file = paste0(store_path,"/01.PCA_clean.pdf"),height = 8,width = 8)
  biplot(
    pcaobj = pca_result,
    colby = "group",
    hline = 0,vline = 0,lab = "",legendPosition = "top",pointSize = p_size
  )
  invisible(dev.off())
  
  png(file = paste0(store_path,"/01.PCA1.png"),height = 2000,width = 2000,res = 300)
  suppressMessages(
    biplot(
      pcaobj = pca_result,
      colby = "group",
      hline = 0,vline = 0,pointSize = p_size,
      title = "biplot",
      subtitle = paste0(left, " vs " ,right),encircle = T,legendPosition = "top"
    )
  )
  
  invisible(dev.off())
  png(file = paste0(store_path,"/01.PCA_clean.png"),height = 2000,width = 2000,res = 300)
  biplot(
    pcaobj = pca_result,
    colby = "group",
    hline = 0,vline = 0,lab = "",legendPosition = "top",pointSize = p_size
  )
  invisible(dev.off())
  
  if(is.null(pca_result)) {q(status = 1)} else {
    message(msg_yes(paste0("Congratulation PCA analysis done! \nplease change directory to \n",getwd(),gsub("./","/",store_path),"\nchecking the result.")))
  }
  
  # DAM analysis ------------------------------------------------------------
  message(msg_run("Step3. Start DAM analysis ... "))
  
  message(msg_run(paste("Run DAM analysis, \nsignificant test method:",test.method,"\npls method:",pls.method)))
  DAM_all = MDAtoolkits::DAM_analysis(
    x = peak_raw,
    left_index = left_index,
    right_index = right_index,
    left = left,
    right = right,
    method = test.method,
    method2 = pls.method
  )
  
  message(msg_run(paste("Export",pls.method,"result!")))
  ## output opls-da figure
  if(!is.null(DAM_all$opls_out)) {
    pdf(paste0(store_path,"/01.opls_da_result.pdf"),width = 10,height = 10)
    plot(DAM_all$opls_out,parPaletteVc = c("purple", "gold"))
    invisible(dev.off())
    
    png(paste0(store_path,"/01.opls_da_result.png"),width = 2000,height = 2000,res = 300)
    plot(DAM_all$opls_out,parPaletteVc = c("purple", "gold"))
    invisible(dev.off())
  }
  
  
  DAM_tbl_all = DAM_all$DAM_tbl
  DAM_tbl_all.anno = DAM_tbl_all %>%
    left_join(.,annot_df,by = "CompoundID")
  write.csv(x = DAM_tbl_all,file = paste0(store_path,"/02.All_compound_DAM_result.csv"),row.names = F)
  write.csv(x = DAM_tbl_all.anno,file = paste0(store_path,"/02.All_compound_DAM_result_annotation.csv"),row.names = F)
  
  message(msg_run("Upset plot..."))
  ## upset plot
  
  FvsH_upset = DAM_tbl_all %>%
    mutate(
      log2fc = case_when(
        abs(log2fc) >= 1 ~ 1,
        TRUE ~ 0
      ),
      pvalue = case_when(
        pvalue <= 0.05 ~ 1,
        TRUE ~ 0
      ),
      FDR = case_when(
        FDR <= 0.05 ~ 1,
        TRUE ~ 0
      ),
      VIP = case_when(
        VIP >= 1 ~ 1,
        TRUE ~ 0
      )
    ) %>%
    column_to_rownames("CompoundID") %>%
    as.matrix()
  
  m = make_comb_mat(FvsH_upset)
  ss = set_size(m)
  cs = comb_size(m)
  
  ht = UpSet(
    m,
    set_order = order(ss),
    comb_order = order(comb_degree(m), -cs),
    top_annotation = HeatmapAnnotation(
      "DM Intersections" = anno_barplot(
        cs,
        ylim = c(0, max(cs)*1.1),
        border = FALSE,
        gp = gpar(fill = "salmon",alpha = 0.6),
        height = unit(4, "cm")
      ),
      annotation_name_side = "left",
      annotation_name_rot = 90
    ),
    left_annotation = rowAnnotation(
      "DM number" = anno_barplot(-ss,
                                 baseline = 0,
                                 border = FALSE,
                                 gp = gpar(fill = "cyan",alpha = 0.6),
                                 width = unit(4, "cm")
      ),
      set_name = anno_text(set_name(m),
                           location = 0.5,
                           just = "center",
                           width = max_text_width(set_name(m)) + unit(4, "mm"))
    ),
    right_annotation = NULL,
    show_row_names = FALSE
  )
  
  ht = draw(ht)
  pdf(paste0(store_path,"/03.DAM_upset.pdf"),width = 6,height = 4)
  ht
  invisible(dev.off())
  
  png(paste0(store_path,"/03.DAM_upset.png"),width = 2500,height = 2000,res = 300)
  ht
  invisible(dev.off())
  
  
  message(msg_run("generate output tables ..."))
  ## output DAM
  
  DAM_compounds = DAM_tbl_all %>%
    filter(
      abs(log2fc) >= log2fc_cut
    ) %>%
    filter(
      pvalue <= pval_cut
    ) %>%
    filter(
      FDR <= qval_cut
    ) %>%
    filter(
      VIP >= VIP_cut
    )
  DAM_compounds.anno = DAM_compounds %>% left_join(.,annot_df,by = "CompoundID")
  if(nrow(DAM_compounds) == 0) {
    message(msg_warning("Warning: There is no DAM compounds under the selected cutoff, \nmaybe you should set a more flexible threshold such as ") %+% magenta$bold$italic$underline("< -q 1 -l 0 >"))
    q(status = 1)
  } else {
    write.csv(x = DAM_compounds,file = paste0(store_path,"/03.DAM_compounds.csv"),row.names = F)
    write.csv(x = DAM_compounds.anno,file = paste0(store_path,"/03.DAM_compounds_annot.csv"),row.names = F)
  }
  
  
  message(msg_run("volcan plot ..."))
  # volcano plot ------------------------------------------------------------
  
  
  
  pdf(paste0(store_path,"/04.volcano_plot.pdf"),width = 8,height = 5)
  DAM_volcano(x = DAM_tbl_all,title = paste0(left," vs ",right),
              pval_cut = pval_cut,log2fc_cut = log2fc_cut,
              qval_cut = qval_cut,VIP_cut = VIP_cut)
  invisible(dev.off())
  
  png(paste0(store_path,"/04.volcano_plot.png"),width = 3000,height = 2200,res = 300)
  DAM_volcano(x = DAM_tbl_all,title = paste0(left," vs ",right),
              pval_cut = pval_cut,log2fc_cut = log2fc_cut,
              qval_cut = qval_cut,VIP_cut = VIP_cut)
  invisible(dev.off())
  
  # heatmap -----------------------------------------------------------------
  
  message(msg_run("heatmap ..."))
  
  DAM_compounds %>%
    select("CompoundID") %>%
    inner_join(.,pca_mat %>% rownames_to_column("CompoundID"),by = "CompoundID") %>%
    column_to_rownames("CompoundID")-> heat_mat
  
  heat_mat=
    heat_mat %>% as.data.frame() %>%
    rownames_to_column("CompoundID") %>%
    arrange(CompoundID) %>%
    column_to_rownames("CompoundID") %>%
    as.matrix()
  ht_DAM_main = Heatmap(
    matrix = t(scale(t(heat_mat))),
    show_row_names = F,
    cluster_columns = F,
    name = "z-scored peak area",
    border_gp = gpar(color = "black",size = 1)
  )
  
  ht_DAM_log2area = Heatmap(
    matrix = log2(heat_mat+1),
    name = "log2(peak area+1)",
    show_row_names = F,cluster_columns = F,
    col = colorRamp2(c(min(log2(heat_mat+1)), max(log2(heat_mat+1))), c("yellow", "#FF1493")),
    border_gp = gpar(color = "black",size = 1)
  )
  
  ht_DAM_log2fc =
    DAM_compounds %>%
    select(CompoundID,log2fc) %>%
    arrange(CompoundID) %>%
    column_to_rownames("CompoundID") %>%
    Heatmap(
      show_row_names = F,cluster_rows = F,
      #row_order = eko_ht_order,
      name = "log2 (fold channge)",
      col = colorRamp2(c(-4,-2, 0, 2,4) ,c("blue","green", "white", "purple","red")),
      border_gp = gpar(color = "black",size = 1)
    )
  
  ht_DAM_pval =
    DAM_compounds %>%
    select(CompoundID,pvalue) %>%
    arrange(CompoundID) %>%
    column_to_rownames("CompoundID") %>%
    Heatmap(
      show_row_names = F,cluster_rows = F,
      name = "pvalue",
      #row_order = eko_ht_order,
      col = colorRamp2(c(0,0.05), c("salmon", "white")),
      border_gp = gpar(color = "black",size = 1),
    )
  ht_DAM_VIP =
    DAM_compounds %>%
    select(CompoundID,VIP) %>%
    arrange(CompoundID) %>%
    column_to_rownames("CompoundID") %>%
    Heatmap(
      show_row_names = F,cluster_rows = F,
      name = "VIP",
      #row_order = eko_ht_order,
      col = colorRamp2(c(0,2), c("white", "cyan")),
      border_gp = gpar(color = "black",size = 1),
    )
  ht_DAM_FDR =
    DAM_compounds %>%
    select(CompoundID,FDR) %>%
    arrange(CompoundID) %>%
    column_to_rownames("CompoundID") %>%
    Heatmap(
      show_row_names = F,cluster_rows = F,
      name = "FDR",
      col = colorRamp2(c(0,0.05), c("#CC6699", "white")),
      border_gp = gpar(color = "black",size = 1),
    )
  ht_DAM = ht_DAM_main+ht_DAM_log2area+ht_DAM_log2fc+ht_DAM_pval+ht_DAM_FDR+ht_DAM_VIP
  pdf(paste0(store_path,"/05.DAM_heatmap.pdf"),width = 12,height = 15)
  draw(ht_DAM)
  invisible(dev.off())
  
  png(paste0(store_path,"/05.DAM_heatmap.png"),width = 3000,height = 3300,res = 300)
  draw(ht_DAM)
  invisible(dev.off())
  
  message(msg_yes(paste0("Congratulation DAM analysis done! \nplease change directory to \n",getwd(),gsub("./","/",store_path),"\nchecking the result.")))
  
  ## running enrichment analysis... ====
  
  message(msg_run("Step4. Start enrichment analysis ..."))
  DAM_compounds.anno =
    DAM_compounds.anno %>%
    dplyr::rename("Name" = "Compound.name",
                  "mf" = "Formula",
                  "KEGG" = "KEGG.ID")
  
  DAM_clean = DAM_compounds.anno %>%
    select(CompoundID,Name,mf,KEGG)
  
  ## for KEGG enrichment
  kegg_t2c  = 
  kegg_t2c %>% left_join(annot_df %>% dplyr::select(CompoundID,KEGG.ID) %>% 
                           setNames(c("CompoundID","KEGGID"))) %>% 
    filter(!is.na(CompoundID)) %>% 
    dplyr::select(PathwayID,CompoundID)
  kegg_t2n <- 
    kegg_t2n %>% inner_join(kegg_t2c %>% select(PathwayID))
  clist = DAM_clean %>%
    dplyr::select(CompoundID,KEGG) %>% 
    filter(is.na(KEGG)) %>% 
    pull(CompoundID)
  
  
  message(msg_run("Running KEGG enrichment ..."))
  
  eko = enricher(
    gene = clist,pvalueCutoff = 0.05,qvalueCutoff = 1,pAdjustMethod = "BH",TERM2GENE = kegg_t2c,TERM2NAME = kegg_t2n
  )
  eko_tbl = as.data.frame(eko)
  if(nrow(eko_tbl) == 0) {
    message(msg_no("There are no significant enriched KEGG pathway according to detected DAMs!"))
  } else {
    eko_bar = barplot(eko,showCategory = 20)
    eko_bar =
      eko_bar + ggtitle(paste0("KEGG enrichment of ", left,"_vs_",right))
    eko_bubble = dotplot(eko,showCategory = 20)
    eko_bubble =
      eko_bubble + ggtitle(paste0("KEGG enrichment of ", left,"_vs_",right))
    MostLenDes = nchar(x = eko_tbl$Description)
    if (max(MostLenDes) >= 150) {
      fig_l = 20
    } else {
      fig_l = 12
    }
    ggsave(eko_bar,filename = paste0(store_path,"/06.KEGG_barplot.pdf"),height = 10,width = fig_l)
    ggsave(eko_bar,filename = paste0(store_path,"/06.KEGG_barplot.png"),height = 10,width = fig_l)
    ggsave(eko_bubble,filename = paste0(store_path,"/06.KEGG_bubbleplot.pdf"),height = 10,width = fig_l)
    ggsave(eko_bubble,filename = paste0(store_path,"/06.KEGG_bubbleplot.png"),height = 10,width = fig_l)
    ## output enrich table
    write.csv(eko_tbl,file = paste0(store_path,"/07.KEGG_enrichment.csv"),row.names = F)
    eko_tbl %>%
      select(Description,geneID) %>%
      cSplit(.,"geneID","/","long") %>%
      rename("KEGG" = "geneID") %>%
      inner_join(.,DAM_clean,by = "KEGG") -> eko_p2c
    eko_p2c %>%
      left_join(.,DAM_compounds,"CompoundID") ->eko_out_tbl
    write.csv(eko_out_tbl,file = paste0(store_path,"/07.Pathway2Compound.csv"),row.names = F)
    ## main matrix
    eko_heat_mat1 = heat_mat %>% as.data.frame() %>%
      rownames_to_column("CompoundID") %>%
      inner_join(.,data.frame(
        CompoundID = unique(eko_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    ## log2fc
    eko_heat_mat2 = DAM_tbl_all %>%
      select(CompoundID,log2fc) %>%
      inner_join(.,data.frame(
        CompoundID = unique(eko_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    ## pvalue
    eko_heat_mat3 = DAM_tbl_all %>%
      select(CompoundID,pvalue) %>%
      inner_join(.,data.frame(
        CompoundID = unique(eko_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    eko_heat_mat4 = DAM_tbl_all %>%
      select(CompoundID,VIP) %>%
      inner_join(.,data.frame(
        CompoundID = unique(eko_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    eko_heat_mat5 = DAM_tbl_all %>%
      select(CompoundID,FDR) %>%
      inner_join(.,data.frame(
        CompoundID = unique(eko_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    ## pvalue tag
    eko_heat_tag = eko_heat_mat3 %>%
      as.data.frame() %>%
      mutate(
        pvalue = case_when(
          pvalue <= 0.0001 ~ "****",
          pvalue > 0.0001 & pvalue <= 0.001 ~ "***",
          pvalue > 0.001 & pvalue <= 0.01 ~ "**",
          pvalue > 0.01 & pvalue <= 0.05 ~ "*",
          TRUE ~ ""
        )
      )
    eko_ht_logmat = Heatmap(
      matrix = log2(eko_heat_mat1),
      name = "log2(peak area)",
      show_row_names = F,cluster_columns = F,
      col = colorRamp2(c(min(log2(eko_heat_mat1)), max(log2(eko_heat_mat1))), c("yellow", "#FF1493")),
      border_gp = gpar(color = "black",size = 1)
    )
    eko_ht_main = Heatmap(
      matrix = t(scale(t(eko_heat_mat1))),
      name = "zscored peak area",
      show_row_names = F,cluster_columns = F,
      border_gp = gpar(color = "black",size = 1)
    )
    #tmp_1 = draw(eko_ht_main)
    #eko_ht_order = row_order(tmp_1)
    eko_ht_log2 = Heatmap(
      matrix = eko_heat_mat2,
      show_row_names = F,cluster_rows = F,
      #row_order = eko_ht_order,
      name = "log2 (fold channge)",
      col = colorRamp2(c(-4,-2, 0, 2,4) ,c("blue","green", "white", "purple","red")),
      border_gp = gpar(color = "black",size = 1)
    )
    
    eko_ht_pval = Heatmap(
      matrix = eko_heat_mat3,
      show_row_names = F,cluster_rows = F,
      name = "pvalue",
      #row_order = eko_ht_order,
      col = colorRamp2(c(0.05,0), c("white", "salmon")),
      border_gp = gpar(color = "black",size = 1),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(eko_heat_tag[i, j], x, y, gp = gpar(fontsize = 10))
      }
    )
    eko_ht_FDR = Heatmap(
      matrix = eko_heat_mat5,
      show_row_names = F,cluster_rows = F,
      name = "FDR",
      #row_order = eko_ht_order,
      col = colorRamp2(c(0.05,0), c("white", "#CC6699")),
      border_gp = gpar(color = "black",size = 1)
    )
    
    eko_ht_VIP = Heatmap(
      matrix = eko_heat_mat4,
      show_row_names = F,cluster_rows = F,
      name = "VIP",
      #row_order = eko_ht_order,
      col = colorRamp2(c(2,0), c("cyan", "white")),
      border_gp = gpar(color = "black",size = 1)
    )
    
    eko_ht = eko_ht_main+eko_ht_logmat+eko_ht_log2+eko_ht_pval+eko_ht_FDR+eko_ht_VIP
    pdf(paste0(store_path,"/08.KEGG_Enrichment_DAM_heatmap.pdf"),width = 10,height = 18)
    draw(eko_ht)
    invisible(dev.off())
    png(paste0(store_path,"/08.KEGG_Enrichment_DAM_heatmap.png"),width = 2200,height = 2600,res = 300 )
    draw(eko_ht)
    invisible(dev.off())
  }
  
  message(msg_yes("KEGG annotation finish!"))
  
  message(msg_run("Start KEGG pathway analysis"))
  invisible(dev.off())
  ## pathway analysis
  eko = enricher(
    gene = clist,pvalueCutoff = 1,qvalueCutoff = 1,pAdjustMethod = "BH",TERM2GENE = kegg_t2c,TERM2NAME = kegg_t2n
  )
  eko_all_tbl = as.data.frame(eko)
  if(nrow(eko_all_tbl) == 0) {
    message(msg_no("No compounds can be matched to KEGG database!"))
  } else  {
    write.csv(eko_all_tbl,file = paste0(store_path,"/09.KEGG_Pathway_analysis_result.csv"),row.names = F)
    pathway_blt = barplot(eko,showCategory = 20)+ggtitle(paste0("Pathway analysis of ",left," vs ",right))
    ggsave(plot = pathway_blt,filename = paste0(store_path,"/09.KEGG_Pathway_bar.pdf"),width = 15,height = 17)
    eko_all_tbl %>%
      select(Description,geneID) %>%
      cSplit(.,"geneID","/","long") %>%
      rename("KEGG" = "geneID") %>%
      inner_join(.,DAM_clean,by = "KEGG") -> eko_p2c
    eko_p2c %>%
      left_join(.,DAM_compounds,"CompoundID") ->eko_out_tbl
    write.csv(eko_out_tbl,file = paste0(store_path,"/09.Pathway2Compound_all.csv"),row.names = F)
    ## main matrix
    eko_heat_mat1 = heat_mat %>% as.data.frame() %>%
      rownames_to_column("CompoundID") %>%
      inner_join(.,data.frame(
        CompoundID = unique(eko_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    ## log2fc
    eko_heat_mat2 = DAM_tbl_all %>%
      select(CompoundID,log2fc) %>%
      inner_join(.,data.frame(
        CompoundID = unique(eko_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    ## pvalue
    eko_heat_mat3 = DAM_tbl_all %>%
      select(CompoundID,pvalue) %>%
      inner_join(.,data.frame(
        CompoundID = unique(eko_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    eko_heat_mat4 = DAM_tbl_all %>%
      select(CompoundID,VIP) %>%
      inner_join(.,data.frame(
        CompoundID = unique(eko_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    eko_heat_mat5 = DAM_tbl_all %>%
      select(CompoundID,FDR) %>%
      inner_join(.,data.frame(
        CompoundID = unique(eko_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    ## pvalue tag
    eko_heat_tag = eko_heat_mat3 %>%
      as.data.frame() %>%
      mutate(
        pvalue = case_when(
          pvalue <= 0.0001 ~ "****",
          pvalue > 0.0001 & pvalue <= 0.001 ~ "***",
          pvalue > 0.001 & pvalue <= 0.01 ~ "**",
          pvalue > 0.01 & pvalue <= 0.05 ~ "*",
          TRUE ~ ""
        )
      )
    eko_ht_logmat = Heatmap(
      matrix = log2(eko_heat_mat1),
      name = "log2(peak area)",
      show_row_names = F,cluster_columns = F,
      col = colorRamp2(c(min(log2(eko_heat_mat1)), max(log2(eko_heat_mat1))), c("yellow", "#FF1493")),
      border_gp = gpar(color = "black",size = 1)
    )
    eko_ht_main = Heatmap(
      matrix = t(scale(t(eko_heat_mat1))),
      name = "zscored peak area",
      show_row_names = F,cluster_columns = F,
      border_gp = gpar(color = "black",size = 1)
    )
    #tmp_1 = draw(eko_ht_main)
    #eko_ht_order = row_order(tmp_1)
    eko_ht_log2 = Heatmap(
      matrix = eko_heat_mat2,
      show_row_names = F,cluster_rows = F,
      #row_order = eko_ht_order,
      name = "log2 (fold channge)",
      col = colorRamp2(c(-4,-2, 0, 2,4) ,c("blue","green", "white", "purple","red")),
      border_gp = gpar(color = "black",size = 1)
    )
    
    eko_ht_pval = Heatmap(
      matrix = eko_heat_mat3,
      show_row_names = F,cluster_rows = F,
      name = "pvalue",
      #row_order = eko_ht_order,
      col = colorRamp2(c(0.05,0), c("white", "salmon")),
      border_gp = gpar(color = "black",size = 1),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(eko_heat_tag[i, j], x, y, gp = gpar(fontsize = 10))
      }
    )
    eko_ht_FDR = Heatmap(
      matrix = eko_heat_mat5,
      show_row_names = F,cluster_rows = F,
      name = "FDR",
      #row_order = eko_ht_order,
      col = colorRamp2(c(0.05,0), c("white", "#CC6699")),
      border_gp = gpar(color = "black",size = 1)
    )
    
    eko_ht_VIP = Heatmap(
      matrix = eko_heat_mat4,
      show_row_names = F,cluster_rows = F,
      name = "VIP",
      #row_order = eko_ht_order,
      col = colorRamp2(c(2,0), c("cyan", "white")),
      border_gp = gpar(color = "black",size = 1)
    )
    
    eko_ht = eko_ht_main+eko_ht_logmat+eko_ht_log2+eko_ht_pval+eko_ht_FDR+eko_ht_VIP
    pdf(paste0(store_path,"/09.KEGG_Pathway_DAM_heatmap.pdf"),width = 10,height = 18)
    draw(eko_ht)
    invisible(dev.off())
    png(paste0(store_path,"/09.KEGG_pathway_DAM_heatmap.png"),width = 2200,height = 2600,res = 300 )
    draw(eko_ht)
    invisible(dev.off())
  }
  
  message(msg_yes("Pathway analysis finish."))
  
  # new ---------------------------------------------------------------------
  message(msg_run("Running ClassyFire enrichment ..."))
  
  DAM_clean = DAM_compounds.anno %>%
    select(CompoundID,Name,mf,KEGG)
  
  ## for classifire enrichment
  
  clist = DAM_clean %>% pull(CompoundID)
  
  ## classyfire database
  met2NAME = 
    annot_df %>% 
    dplyr::select(CompoundID,superclass,class,subclass,parent_levels) %>% 
    pivot_longer(!CompoundID,names_to = 'category',values_to = 'NAME') %>% 
    separate_longer_delim(cols = NAME,delim = "|") %>% 
    drop_na() %>% 
    filter(NAME != "NA") %>% 
    filter(NAME != "") %>% 
    mutate(NAME = str_trim(NAME))
  TERM2NAME = met2NAME %>% 
    dplyr::select(category,NAME) %>% 
    distinct() %>% 
    group_by(category) %>% 
    mutate(n = n()) %>% 
    mutate(Level = case_when(
      category == "superclass" ~ 'A',
      category == "class" ~ 'B',
      category == "subclass" ~ 'C',
      category == "parent_levels" ~ 'D'
    )) %>% 
    ungroup() %>% 
    dplyr::select(-category) %>% 
    arrange(Level) %>% 
    group_by(Level) %>% 
    mutate(TERM = paste0(Level,str_pad(string = 1:n,width = 5,side = "left",pad = "0"))) %>% 
    ungroup %>% 
    dplyr::select(TERM,NAME) %>% distinct()
  TERM2GENE = TERM2NAME %>% left_join(met2NAME,by = "NAME") %>% 
    dplyr::select(TERM,CompoundID) %>% distinct()
  
  eclass = enricher(
    gene = clist,pvalueCutoff = 0.05,qvalueCutoff = 1,pAdjustMethod = "BH",TERM2GENE = TERM2GENE,TERM2NAME = TERM2NAME
  )
  eclass_tbl = as.data.frame(eclass)
  if(nrow(eclass_tbl) == 0) {
    message(msg_no("There are no significant enriched ClassyFire pathway according to detected DAMs!"))
  } else {
    eclass_bar = barplot(eclass,showCategory = 20)
    eclass_bar =
      eclass_bar + ggtitle(paste0("ClassyFire enrichment of ", left,"_vs_",right))
    eclass_bubble = dotplot(eclass,showCategory = 20)
    eclass_bubble =
      eclass_bubble + ggtitle(paste0("ClassyFire enrichment of ", left,"_vs_",right))
    MostLenDes = nchar(x = eclass_tbl$Description)
    if (max(MostLenDes) >= 150) {
      fig_l = 20
    } else {
      fig_l = 12
    }
    ggsave(eclass_bar,filename = paste0(store_path,"/10.ClassyFire_barplot.pdf"),height = 10,width = fig_l)
    ggsave(eclass_bar,filename = paste0(store_path,"/10.ClassyFire_barplot.png"),height = 10,width = fig_l)
    ggsave(eclass_bubble,filename = paste0(store_path,"/10.ClassyFire_bubbleplot.pdf"),height = 10,width = fig_l)
    ggsave(eclass_bubble,filename = paste0(store_path,"/10.ClassyFire_bubbleplot.png"),height = 10,width = fig_l)
    ## output enrich table
    write.csv(eclass_tbl,file = paste0(store_path,"/11.ClassyFire_enrichment.csv"),row.names = F)
    eclass_tbl %>%
      select(Description,geneID) %>%
      cSplit(.,"geneID","/","long") %>%
      rename("CompoundID" = "geneID") %>%
      inner_join(.,DAM_clean,by = "CompoundID") -> eclass_p2c
    eclass_p2c %>%
      left_join(.,DAM_compounds,"CompoundID") ->eclass_out_tbl
    write.csv(eclass_out_tbl,file = paste0(store_path,"/11.Pathway2Compound.csv"),row.names = F)
    ## main matrix
    eclass_heat_mat1 = heat_mat %>% as.data.frame() %>%
      rownames_to_column("CompoundID") %>%
      inner_join(.,data.frame(
        CompoundID = unique(eclass_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    ## log2fc
    eclass_heat_mat2 = DAM_tbl_all %>%
      select(CompoundID,log2fc) %>%
      inner_join(.,data.frame(
        CompoundID = unique(eclass_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    ## pvalue
    eclass_heat_mat3 = DAM_tbl_all %>%
      select(CompoundID,pvalue) %>%
      inner_join(.,data.frame(
        CompoundID = unique(eclass_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    eclass_heat_mat4 = DAM_tbl_all %>%
      select(CompoundID,VIP) %>%
      inner_join(.,data.frame(
        CompoundID = unique(eclass_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    eclass_heat_mat5 = DAM_tbl_all %>%
      select(CompoundID,FDR) %>%
      inner_join(.,data.frame(
        CompoundID = unique(eclass_p2c$CompoundID)
      ),"CompoundID") %>%
      arrange(CompoundID) %>%
      column_to_rownames("CompoundID") %>% as.matrix()
    ## pvalue tag
    eclass_heat_tag = eclass_heat_mat3 %>%
      as.data.frame() %>%
      mutate(
        pvalue = case_when(
          pvalue <= 0.0001 ~ "****",
          pvalue > 0.0001 & pvalue <= 0.001 ~ "***",
          pvalue > 0.001 & pvalue <= 0.01 ~ "**",
          pvalue > 0.01 & pvalue <= 0.05 ~ "*",
          TRUE ~ ""
        )
      )
    eclass_ht_logmat = Heatmap(
      matrix = log2(eclass_heat_mat1),
      name = "log2(peak area)",
      show_row_names = F,cluster_columns = F,
      col = colorRamp2(c(min(log2(eclass_heat_mat1)), max(log2(eclass_heat_mat1))), c("yellow", "#FF1493")),
      border_gp = gpar(color = "black",size = 1)
    )
    eclass_ht_main = Heatmap(
      matrix = t(scale(t(eclass_heat_mat1))),
      name = "zscored peak area",
      show_row_names = F,cluster_columns = F,
      border_gp = gpar(color = "black",size = 1)
    )
    #tmp_1 = draw(eclass_ht_main)
    #eclass_ht_order = row_order(tmp_1)
    eclass_ht_log2 = Heatmap(
      matrix = eclass_heat_mat2,
      show_row_names = F,cluster_rows = F,
      #row_order = eclass_ht_order,
      name = "log2 (fold channge)",
      col = colorRamp2(c(-4,-2, 0, 2,4) ,c("blue","green", "white", "purple","red")),
      border_gp = gpar(color = "black",size = 1)
    )
    
    eclass_ht_pval = Heatmap(
      matrix = eclass_heat_mat3,
      show_row_names = F,cluster_rows = F,
      name = "pvalue",
      #row_order = eclass_ht_order,
      col = colorRamp2(c(0.05,0), c("white", "salmon")),
      border_gp = gpar(color = "black",size = 1),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(eclass_heat_tag[i, j], x, y, gp = gpar(fontsize = 10))
      }
    )
    eclass_ht_FDR = Heatmap(
      matrix = eclass_heat_mat5,
      show_row_names = F,cluster_rows = F,
      name = "FDR",
      #row_order = eclass_ht_order,
      col = colorRamp2(c(0.05,0), c("white", "#CC6699")),
      border_gp = gpar(color = "black",size = 1)
    )
    
    eclass_ht_VIP = Heatmap(
      matrix = eclass_heat_mat4,
      show_row_names = F,cluster_rows = F,
      name = "VIP",
      #row_order = eclass_ht_order,
      col = colorRamp2(c(2,0), c("cyan", "white")),
      border_gp = gpar(color = "black",size = 1)
    )
    
    eclass_ht = eclass_ht_main+eclass_ht_logmat+eclass_ht_log2+eclass_ht_pval+eclass_ht_FDR+eclass_ht_VIP
    pdf(paste0(store_path,"/08.ClassyFire_Enrichment_DAM_heatmap.pdf"),width = 10,height = 18)
    draw(eclass_ht)
    invisible(dev.off())
    png(paste0(store_path,"/08.ClassyFire_Enrichment_DAM_heatmap.png"),width = 2200,height = 2600,res = 300 )
    draw(eclass_ht)
    invisible(dev.off())
  }
  
  message(msg_yes("ClassyFire annotation finish!"))
  
  message(msg_yes(paste0("Congratulation ALL analysis done! \nplease change directory to \n",getwd(),gsub("./","/",store_path),"\nchecking the result.")))
  
  
})


