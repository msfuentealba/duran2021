library(tidyverse)

grid <- read_rds("./output/linear_model.rds")
comb <- grid$results[[1]][,c(1,7)] %>% full_join(grid$results[[2]][,c(1,7)], by = c("uniprot"="uniprot")) %>% full_join(grid$results[[3]][,c(1,7)], by = c("uniprot"="uniprot"))
colnames(comb) <- c("uniprot","WT vs TgXBP1","WT vs 5xFAD","5xFAD vs TgXBP1/5xFAD")
c(grid$results[[1]]$symbol,grid$results[[2]]$symbol,grid$results[[3]]$symbol) %>% unique %>% length #confirm length

r <- cor(comb[,2:4], use = "pairwise.complete.obs", method = "spearman")
r
cor.test(comb$`WT vs TgXBP1`,comb$`WT vs 5xFAD`, method = "spearman")
cor.test(comb$`WT vs TgXBP1`,comb$`5xFAD vs TgXBP1/5xFAD`, method = "spearman")
cor.test(comb$`WT vs 5xFAD`,comb$`5xFAD vs TgXBP1/5xFAD`, method = "spearman")$p.value


#previous analysis validation
#xbp_wt <- readxl::read_xlsx("./info/1.xlsx", col_names = TRUE, skip = 3, sheet = 2)[,c(1,48)] %>% set_names("uniprot","wt_xbp")
#fad_wt <- readxl::read_xlsx("./info/1.xlsx", col_names = TRUE, skip = 3, sheet = 4)[,c(1,74)] %>% set_names("uniprot","wt_fad")
#fx_fad <- readxl::read_xlsx("./info/1.xlsx", col_names = TRUE, skip = 3, sheet = 6)[,c(1,76)] %>% set_names("uniprot","fad_fx")

#comb2 <- xbp_wt %>% left_join(fad_wt) %>% left_join(fx_fad)
#r2 <- cor(comb2[,2:4], use = "pairwise.complete.obs", method = "spearman")
#r2

#heatmap
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

col_fun <- colorRamp2(seq(-1,1,2/10),rev(brewer.pal(11, "RdBu")))

ano_col <- columnAnnotation(Control = c("WT","WT","5xFAD"), 
                            Mutant = c("TgXBP1","5xFAD","TgXBP1/5xFAD"), 
                            col = list(Control = c("WT" = "#fefefe", "5xFAD" = "#00007a"),
                                       Mutant = c("TgXBP1" = "#4d4c4b", "5xFAD" = "#00007a",  "TgXBP1/5xFAD" = "#3a7d22")), 
                            gp = gpar(col = "black", lwd = 1),
                            gap = unit(1, "mm"))

ano_row <- rowAnnotation(Control = c("WT","WT","5xFAD"), 
                         Mutant = c("TgXBP1","5xFAD","TgXBP1/5xFAD"), 
                         col = list(Control = c("WT" = "#fefefe", "5xFAD" = "#00007a"),
                                    Mutant = c("TgXBP1" = "#4d4c4b", "5xFAD" = "#00007a",  "TgXBP1/5xFAD" = "#3a7d22")), 
                         gp = gpar(col = "black", lwd = 1),
                         gap = unit(1, "mm"))

pdf(file = paste0("./output/figures/correlations.pdf"), height = 3, width = 5) 
Heatmap(r,
        col = col_fun,
        rect_gp = gpar(col = "white", lwd = 1),
        #rect_gp = gpar(type = "none")
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(round(r[i,j],2), x = x, y = y, gp=gpar(fontsize=10, fontface = "bold" ))
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA, lty = 1, lwd = 1))
        },
        heatmap_legend_param = list(direction = "vertical", title = "Spearman's\nrho", title_position = "topcenter"),
        top_annotation = ano_col,
        left_annotation = ano_row,
        column_dend_side = "bottom",
        row_dend_side = "left",
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_row_names=FALSE,
        show_column_names=FALSE,
        border = TRUE)
dev.off()        

col_fun_or <- colorRamp2(seq(-5,5,1),rev(brewer.pal(11, "RdBu")))
pdf(file = paste0("./output/figures/heatmap.pdf"), height = 7, width = 5) 
Heatmap(as.matrix(comb[,2:4] %>% na.omit),
        top_annotation = ano_col,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        heatmap_legend_param = list(direction = "vertical", title = "Score", title_position = "topcenter"),
        col = col_fun_or)
dev.off()
