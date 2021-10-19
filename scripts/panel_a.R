library(tidyverse)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

#data <- read_rds("./output/linear_model.rds") %>% unnest(c("results"))
#write_tsv(data, file = "./raw_data/linear_model_results.tsv")

data <- read_tsv("./raw_data/linear_model_results.tsv")
xbp_wt <- data %>% filter(control=="wt"&mutant=="xbp") %>% group_by(symbol) %>% summarise(score = mean(score))
xbp_wt$group <- "xbp_wt"
fad_wt <- data %>% filter(control=="wt"&mutant=="fad") %>% group_by(symbol) %>% summarise(score = mean(score))
fad_wt$group <- "fad_wt"
fx_fad <- data %>% filter(control=="fad"&mutant=="fx") %>% group_by(symbol) %>% summarise(score = mean(score))
fx_fad$group <- "fx_fad"
all <- rbind(xbp_wt,fad_wt,fx_fad)
mat <- reshape2::dcast(all,symbol~group,value.var = "score")
rownames(mat) <- mat$symbol
mat$symbol <- NULL
mat <- mat %>% select(xbp_wt,fad_wt,fx_fad) %>% arrange(fad_wt) %>%na.omit

cor(mat,use = "pairwise.complete.obs",method = "s")
cor.test(mat[,1],mat[,2],method = "s",use="pairwise.complete.obs")
cor.test(mat[,1],mat[,3],method = "s",use="pairwise.complete.obs")
cor.test(mat[,2],mat[,3],method = "s",use="pairwise.complete.obs")

col_fun <- colorRamp2(seq(-5,5,1),rev(brewer.pal(11, "RdBu")))
colnames(mat) <- c("TgXBP1s /\nWT","5xFAD /\nWT","TgXBP1s-5xFAD /\n5xFAD")
pdf(file = paste0("./output/panel_a.pdf"), height = 5, width = 2) 
Heatmap(as.matrix(mat),
        border = TRUE,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        heatmap_legend_param = list(direction = "vertical", title = "Significance\nscore", title_position = "topcenter"),
        col = col_fun)
dev.off()
