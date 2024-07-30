library(CrossTalkeR)

# the method always consider the first path as control: the multiple control case will be handle soon
paths <- c(
  'Pre' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/pre_cellphonedb_results/filtered_corrected.csv',
  'During' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/during_cellphonedb_results/filtered_corrected.csv',
  'JustFinish' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/justfinish_cellphonedb_results/filtered_corrected.csv',
  'After' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/after_cellphonedb_results/filtered_corrected.csv'
)


# Generating the report and the object
data <- CrossTalkeR::generate_report(
  lrpaths=paths,
  threshold=0,
  out_path='/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/', 
  out_file='All_time.html', 
  output_fmt = "html_document", 
 # report = FALSE
)


# compare -----------------------------------------------------------------------------
# pre vs during
paths <- c(
  'Pre' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/pre_cellphonedb_results/filtered_corrected.csv',
  'During' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/during_cellphonedb_results/filtered_corrected.csv'
)

# Selected gene list     
genes <- c('Type II IFNR|R', "CTLA4|L", "SIRPA|R", "CD274|L", "CD28|L", "TGFB1|L", "TIGIT|L", "CXCL10|L", "CXCL9|L", "CXCL13|L", "PVR|R", "CD244|R", "BTLA|R")

# Generating the report and the object
data <- CrossTalkeR::generate_report(
  genes = genes,
  lrpaths=paths,
  threshold=0,
  out_path='/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/', 
  out_file='pre_vs_during.html', 
  output_fmt = "html_document", 
  # report = FALSE
)

png(paste(output_dir,"During.Pre.network.png",sep = "/"), width = 3800, height =3800, res = 500)　
plot_cci(
  graph = data@graphs$During_x_Pre, coords = data@coords, colors = data@colors, plt_name = "", pg = data@rankings$During_x_Pre$Pagerank,  
  vnames = F, leg = T
)
dev.off() #irisplot_2.png　出力
#ggsave(plot = network_plot, filename = paste(output_dir,"During.Pre.network.png",sep = "/"),height = 5, width = 5, dpi = 500)

# during vs justfinish
paths <- c(
  'During' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/during_cellphonedb_results/filtered_corrected.csv', 
  'JustFinish' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/justfinish_cellphonedb_results/filtered_corrected.csv'
)

# Selected gene list     
genes <- c('Type II IFNR|R', "CTLA4|L", "SIRPA|R", "CD274|L", "CD28|L", "TGFB1|L", "TIGIT|L", "CXCL10|L", "CXCL9|L", "CXCL13|L", "PVR|R", "CD244|R", "BTLA|R")

# Generating the report and the object
data <- CrossTalkeR::generate_report(
  genes = genes,
  lrpaths=paths,
  threshold=0,
  out_path='/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/', 
  out_file='during_vs_justfinish.html', 
  output_fmt = "html_document", 
  # report = FALSE
)

png(paste(output_dir,"Justfinish.Pre.network.png",sep = "/"), width = 3800, height =3800, res = 500)　
plot_cci(
  graph = data@graphs$JustFinish_x_During, coords = data@coords, colors = data@colors, plt_name = "", pg = data@rankings$JustFinish_x_During$Pagerank,  
  vnames = T, leg = T
)
dev.off() #irisplot_2.png　出力


# during vs justfinish
paths <- c(
  'JustFinish' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/justfinish_cellphonedb_results/filtered_corrected.csv', 
  'After' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/after_cellphonedb_results/filtered_corrected.csv'
)

# Selected gene list     
genes <- c('Type II IFNR|R', "CTLA4|L", "SIRPA|R", "CD274|L", "CD28|L", "TGFB1|L", "TIGIT|L", "CXCL10|L", "CXCL9|L")

# Generating the report and the object
data <- CrossTalkeR::generate_report(
  genes = genes,
  lrpaths=paths,
  threshold=0,
  out_path='/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/', 
  out_file='justfinish_vs_after.html', 
  #output_fmt = "html_document", 
  # report = FALSE
)

# pre vs justfinish
paths <- c(
  'Pre' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/pre_cellphonedb_results/filtered_corrected.csv', 
  'JustFinish' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/justfinish_cellphonedb_results/filtered_corrected.csv'
)

# Selected gene list     
genes <- c('Type II IFNR|R', "CTLA4|L", "SIRPA|R", "CD274|L", "CD28|L", "TGFB1|L", "TIGIT|L", "CXCL10|L", "CXCL9|L", "CXCL13|L", "PVR|R", "CD244|R", "BTLA|R", "HVEM|L")

# Generating the repo rt and the object
data <- CrossTalkeR::generate_report(
  genes = genes,
  lrpaths=paths,
  threshold=0,
  out_path='/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/', 
  #out_file='pre_vs_justfinish.html', 
  out_file = 'pre_vs_justfinish.pdf',
  output_fmt = "pdf_document",
  #output_fmt = "html_document", 
  report = TRUE
)

png(paste(output_dir,"Justfinish.Pre.network.png",sep = "/"), width = 3800, height =3800, res = 500)　
options(repr.plot.width=20, repr.plot.height=20)
plot_cci(
  graph = data@graphs$JustFinish_x_Pre, coords = data@coords[V(data@graphs$JustFinish_x_Pre), ], colors = data@colors[V(data@graphs$JustFinish_x_Pre)], 
  plt_name = "", pg = data@rankings$JustFinish_x_Pre$Pagerank[V(data@graphs$JustFinish_x_Pre)],  
  vnames = T, leg = T, ignore_alpha = FALSE,
  log = FALSE
)
dev.off() #irisplot_2.png　出力


# pre vs after
paths <- c(
  'Pre' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/pre_cellphonedb_results/filtered_corrected.csv', 
  'After' = '/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/after_cellphonedb_results/filtered_corrected.csv'
)

# Selected gene list     
genes <- c('Type II IFNR|R', "CTLA4|L", "SIRPA|R", "CD274|L", "CD28|L", "TGFB1|L", "TIGIT|L", "CXCL10|L", "CXCL9|L")

# Generating the report and the object
data <- CrossTalkeR::generate_report(
  genes = genes,
  lrpaths=paths,
  threshold=0,
  out_path='/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/', 
  out_file='pre_vs_after.html', 
  #output_fmt = "html_document", 
  # report = FALSE
)

