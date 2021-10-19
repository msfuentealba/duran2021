library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
library(GeneOverlap)
library(ggvenn)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(ggpubr)
library(gridExtra)

data <- read_tsv("./raw_data/linear_model_results.tsv")
xbp_wt <- data %>% filter(control=="wt"&mutant=="xbp") %>% group_by(uniprot) %>% summarise(score = mean(score))
fad_wt <- data %>% filter(control=="wt"&mutant=="fad") %>% group_by(uniprot) %>% summarise(score = mean(score))
fx_fad <- data %>% filter(control=="fad"&mutant=="fx") %>% group_by(uniprot) %>% summarise(score = mean(score))

gl_fad_wt <- fad_wt$score
names(gl_fad_wt) <- fad_wt$uniprot
gl_fad_wt <- gl_fad_wt %>% sort %>% rev

gl_fx_fad <- fx_fad$score
names(gl_fx_fad) <- fx_fad$uniprot
gl_fx_fad <- gl_fx_fad %>% sort %>% rev

set.seed(10)
g_fad_wt <- gseGO(
  gl_fad_wt,
  ont = "BP",
  OrgDb = org.Mm.eg.db,
  keyType = "UNIPROT",
  minGSSize = 20,
  maxGSSize = 200,
  pvalueCutoff = 1,
  pAdjustMethod = "BH")@result

set.seed(10)
g_fx_fad <- gseGO(
  gl_fx_fad,
  ont = "BP",
  OrgDb = org.Mm.eg.db,
  keyType = "UNIPROT",
  minGSSize = 20,
  maxGSSize = 200,
  pvalueCutoff = 1,
  pAdjustMethod = "BH")@result

fad_wt_up <- g_fad_wt$Description[g_fad_wt$NES>0&g_fad_wt$p.adjust<0.05]
fad_wt_dn <- g_fad_wt$Description[g_fad_wt$NES<0&g_fad_wt$p.adjust<0.05]
fx_fad_up <- g_fx_fad$Description[g_fx_fad$NES>0&g_fx_fad$p.adjust<0.05]
fx_fad_dn <- g_fx_fad$Description[g_fx_fad$NES<0&g_fx_fad$p.adjust<0.05]

testGeneOverlap(newGeneOverlap(fad_wt_up,fx_fad_dn,genome.size = length(unique(c(g_fad_wt$Description,g_fx_fad$Description)))))
testGeneOverlap(newGeneOverlap(fad_wt_dn,fx_fad_up,genome.size = length(unique(c(g_fad_wt$Description,g_fx_fad$Description)))))

c <- brewer.pal(11, "RdBu")

p1 <- ggvenn(list(`5xFAD /\nWT\n(up)` = fad_wt_up, `TgXBP1s-5xFAD /\n5xFAD\n(down)` = fx_fad_dn),
       show_percentage = FALSE,
       fill_color = c("white","white"),
       stroke_color =c(rep(c[2],100),rep(c[10],100)),
       stroke_size = 2,
       text_size = 8,
       set_name_size = 8) 

p2 <- ggvenn(list(`5xFAD /\nWT\n(down)` = fad_wt_dn, `TgXBP1s-5xFAD /\n5xFAD\n(up)` = fx_fad_up),
             show_percentage = FALSE,
             fill_color = c("white","white"),
             stroke_color =c(rep(c[10],100),rep(c[2],100)),
             stroke_size = 2,
             text_size = 8,
             set_name_size = 8) 

pdf(file = paste0("./output/panel_b.pdf"), height = 9, width = 10) 
ggarrange(p1,p2)
dev.off()

uni2sym <- data %>% dplyr::select(uniprot,symbol) %>% unique
up_dn_nes <- g_fad_wt[g_fad_wt$Description%in%intersect(fad_wt_up,fx_fad_dn),c(2,5,11)] %>% left_join(g_fx_fad[g_fx_fad$Description%in%intersect(fad_wt_up,fx_fad_dn),c(2,5,11)], by = c("Description"="Description"))
up_dn_nes <- up_dn_nes %>% group_by(Description) %>% summarise(mean = mean(abs(NES.x),abs(NES.y)), inter = list(intersect(unlist(strsplit(core_enrichment.x,split = "\\/")),unlist(strsplit(core_enrichment.y,split = "\\/"))))) %>% arrange(desc(mean))
up_dn_nes <- up_dn_nes[1:10,] %>% unnest(inter) %>% set_names(c("go","nes","uniprot")) %>% left_join(uni2sym, by = c("uniprot"="uniprot")) %>% na.omit
up_dn_nes$on <- 1
up_dn_nes <- reshape2::acast(up_dn_nes, symbol~go, value.var = "on")
up_dn_nes[is.na(up_dn_nes)] <- 0
up_dn_nes <- up_dn_nes[rowSums(up_dn_nes)>5,]

pdf(file = paste0("./output/panel_b_table1.pdf"), height = 7, width = 4) 
grid.table(as.data.frame(matrix(sort(colnames(up_dn_nes)), ncol = 1)), theme = ttheme_minimal(), rows = NULL, cols = NULL)
dev.off()

dn_up_nes <- g_fad_wt[g_fad_wt$Description%in%intersect(fad_wt_dn,fx_fad_up),c(2,5,11)] %>% left_join(g_fx_fad[g_fx_fad$Description%in%intersect(fad_wt_dn,fx_fad_up),c(2,5,11)], by = c("Description"="Description"))
dn_up_nes <- dn_up_nes %>% group_by(Description) %>% summarise(mean = mean(abs(NES.x),abs(NES.y)), inter = list(intersect(unlist(strsplit(core_enrichment.x,split = "\\/")),unlist(strsplit(core_enrichment.y,split = "\\/"))))) %>% arrange(desc(mean))
dn_up_nes <- dn_up_nes[1:10,] %>% unnest(inter) %>% set_names(c("go","nes","uniprot")) %>% left_join(uni2sym, by = c("uniprot"="uniprot")) %>% na.omit
dn_up_nes$on <- 1
dn_up_nes <- reshape2::acast(dn_up_nes, symbol~go, value.var = "on")
dn_up_nes[is.na(dn_up_nes)] <- 0
dn_up_nes <- dn_up_nes[rowSums(dn_up_nes)>5,]

pdf(file = paste0("./output/panel_b_table2.pdf"), height = 7, width = 4) 
grid.table(as.data.frame(matrix(sort(colnames(dn_up_nes)), ncol = 1)), theme = ttheme_minimal(), rows = NULL, cols = NULL)
dev.off()

