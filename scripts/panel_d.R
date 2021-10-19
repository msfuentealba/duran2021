library(tidyverse)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

data <- read_tsv("./raw_data/human_datasets.tsv")
mat <- reshape2::acast(data, symbol~group, value.var = "logfc", fun.aggregate = mean)
mat[is.nan(mat)] <- NA
dim(mat)
r <- cor(mat, use = "pairwise.complete.obs", method = "spearman")

#calculate p-values
r2 <- r
r2[upper.tri(r2)] <- NA
diag(r2) <- NA
grid <- reshape2::melt(r2) %>% na.omit
grid$value2 <- apply(grid, 1, function(x) cor.test(mat[,x['Var1']],mat[,x['Var2']], use = "complete.obs", method = "spearman")$estimate)
grid$pvalue <- apply(grid, 1, function(x) cor.test(mat[,x['Var1']],mat[,x['Var2']], use = "complete.obs", method = "spearman")$p.value)
grid$fdr <- p.adjust(grid$pvalue, method = "BH")
grid$n <- apply(grid, 1, function(x) sum(rowSums(!is.na(mat[,c(x['Var1'],x['Var2'])]))==2))

colnames(r)
colnames(r) <- c("Human - Age 62\n(Xu et al. 2016)",
                 "Human - Age 84\n(Xu et al. 2016)",
                 "Human - Age 95\n(Xu et al. 2016)",
                 "Human - Alzheimer's disease\n(Xu et al. 2019)",
                 "Mouse - 5xFAD model - 10 months\n(Kim et al. 2019)",
                 "Mouse - 5xFAD model - 5 months\n(Kim et al. 2019)",
                 "Mouse - 5xFAD model\n(This study)",
                 "Mouse - Xbp1 model\n(This study)")
rownames(r) <- colnames(r)
r_out <- reshape2::melt(r) %>% set_names(c("group1","group2","r"))
r_out$group1 <- gsub("\n"," ",r_out$group1)
r_out$group2 <- gsub("\n"," ",r_out$group2)

col_fun <- colorRamp2(seq(-1,1,2/10),rev(brewer.pal(11, "RdBu")))
pdf(file = paste0("./output/panel_c.pdf"), height = 7, width = 8) 
Heatmap(r,
        col = col_fun,
        rect_gp = gpar(col = "white", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(round(r[i,j],2), x = x, y = y, gp=gpar(fontsize=10, fontface = "bold" ))
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA, lty = 1, lwd = 1))
        },
        heatmap_legend_param = list(direction = "vertical", title = "Spearman's\nrho", title_position = "topcenter"),
        clustering_method_rows = "average",
        clustering_method_columns = "average",
        column_dend_side = "top",
        row_dend_side = "left",
        column_names_rot = 45,
        border = TRUE)
dev.off()

write.csv(r, file = "./output/raw_panel_d.csv", row.names = TRUE)


