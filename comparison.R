library(tidyverse)
grid <- read_rds("./output/linear_model.rds")

xbp <- grid$results[[1]]
xbp$group <- "thisxbp"
xbp <- xbp[,c(2:3,8)]

fad <- grid$results[[2]]
fad$group <- "thisfad"
fad <- fad[,c(2:3,8)]

kim5 <- readxl::read_xlsx("./data/datasets.xlsx", sheet = 1, skip = 3)[,c(5,10)] %>% set_names("symbol","logfc")
kim5$symbol <- gsub(";.*", "",kim5$symbol)
kim5$logfc <- kim5$logfc*-1
kim5 <- kim5 %>% na.omit
kim5$group <- "kim5"

kim10 <- readxl::read_xlsx("./data/datasets.xlsx", sheet = 2, skip = 3)[,c(5,10)] %>% set_names("symbol","logfc")
kim10$symbol <- gsub(";.*", "",kim10$symbol)
kim10$logfc <- kim10$logfc*-1
kim10 <- kim10 %>% na.omit
kim10$group <- "kim10"

xu <- readxl::read_xlsx("./data/datasets.xlsx", sheet = 3, skip = 1)[,c(1,34)] %>% set_names("human","logfc")
xu <- xu %>% na.omit
mouse2human <- read_tsv("./data/mouse2human.txt") %>% set_names("mouse","human") %>% na.omit
xu <- xu %>% left_join(mouse2human) %>% na.omit
xu <- xu[,c(3,2)] %>% set_names("symbol","logfc")
xu$group <- "human_alz"

age <- readxl::read_xlsx("./data/datasets.xlsx", sheet = 4)[,c(1,6:8)] %>% set_names("uniprot","b","c","d")
age$uniprot <- gsub("-.*", "",age$uniprot)
uniprot2genes <- read_tsv("./data/uniprot2genes.txt") %>% set_names("uniprot","human") %>% na.omit
age <- age %>% left_join(uniprot2genes) %>% left_join(mouse2human)
age_62 <- age[,c(6,2)] %>% set_names("symbol","logfc") %>% na.omit
age_62$group <- "age_62"
age_84 <- age[,c(6,3)] %>% set_names("symbol","logfc") %>% na.omit
age_84$group <- "age_84"
age_95 <- age[,c(6,4)] %>% set_names("symbol","logfc") %>% na.omit
age_95$group <- "age_95"

all <- rbind(xbp,fad,kim5,kim10,xu,age_62,age_84,age_95)
mat <- reshape2::acast(all, symbol~group, value.var = "logfc", fun.aggregate = mean)
mat[is.nan(mat)] <- NA
dim(mat)
r <- cor(mat, use = "pairwise.complete.obs", method = "spearman")

r2 <- r
r2[upper.tri(r2)] <- NA
diag(r2) <- NA
grid <- reshape2::melt(r2) %>% na.omit
grid$value2 <- apply(grid, 1, function(x) cor.test(mat[,x['Var1']],mat[,x['Var2']], use = "complete.obs", method = "spearman")$estimate)
grid$pvalue <- apply(grid, 1, function(x) cor.test(mat[,x['Var1']],mat[,x['Var2']], use = "complete.obs", method = "spearman")$p.value)
grid$fdr <- p.adjust(grid$pvalue, method = "BH")
grid$n <- apply(grid, 1, function(x) sum(rowSums(!is.na(mat[,c(x['Var1'],x['Var2'])]))==2))
#grid$label <- paste0(round(grid$value2,2),"\n(",format(grid$fdr, scientific = TRUE, digits = 3),")")
# r2 <- reshape2::acast(grid, Var1~Var2, value.var = "value2")
# r2 <- r2[rownames(r),colnames(r)]
# labelmat <- reshape2::acast(grid, Var1~Var2, value.var = "label")
# labelmat <- labelmat[rownames(r),colnames(r)]

colnames(r)
dim(r)
colnames(r)
colnames(r) <- c("Human - Age 62\n(Xu et al. 2016)",
                 "Human - Age 84\n(Xu et al. 2016)",
                 "Human - Age 95\n(Xu et al. 2016)",
                 "Human - Alzheimer's disease\n(Xu et al. 2019)",
                 #"Human - Alzheimer's disease (Stepler et al. 2020)",
                 "Mouse - 5xFAD model - 10 months\n(Kim et al. 2019)",
                 "Mouse - 5xFAD model - 5 months\n(Kim et al. 2019)",
                 "Mouse - 5xFAD model\n(This study)",
                 "Mouse - Xbp1 model\n(This study)")
rownames(r) <- colnames(r)
library(circlize)
library(RColorBrewer)

col_fun <- colorRamp2(seq(-1,1,2/10),rev(brewer.pal(11, "RdBu")))
library(ComplexHeatmap)
pdf(file = paste0("./output/figures/comparison.pdf"), height = 7, width = 8) 
Heatmap(r,
        col = col_fun,
        rect_gp = gpar(col = "white", lwd = 1),
        #rect_gp = gpar(type = "none")
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(round(r[i,j],2), x = x, y = y, gp=gpar(fontsize=10, fontface = "bold" ))
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA, lty = 1, lwd = 1))
        },
        heatmap_legend_param = list(direction = "vertical", title = "Spearman's\nrho", title_position = "topcenter"),
        clustering_method_rows = "average",
        clustering_method_columns = "average",
        #top_annotation = ano_col,
        #left_annotation = ano_row,
        column_dend_side = "top",
        row_dend_side = "left",
        column_names_rot = 45,
        border = TRUE)
dev.off()
# cog <- readxl::read_xlsx("./data/cogn.xlsx", sheet = 1, skip=3, col_names = TRUE)[,c(1,3,5)] %>% set_names("human","logfc","pvalue")
# cog$score <- (-log10(cog$pvalue)*sign(cog$logfc))
# cog$human <- sapply(cog$human, function(x) strsplit(x, split = "\\|")[[1]][1]) %>% as.character
# cog <- cog %>% left_join(mouse2human) %>% na.omit
# cog <- cog[,c(5,2:4)] %>% set_names("symbol","logfc","pvalue","score")
# cog$group <- "human_cognitive1"
# 
# cog2 <- readxl::read_xlsx("./data/cogn.xlsx", sheet = 2, skip=3, col_names = TRUE)[,c(1,3,5)] %>% set_names("human","logfc","pvalue")
# cog2$score <- (-log10(cog2$pvalue)*sign(cog2$logfc))
# cog2$human <- sapply(cog2$human, function(x) strsplit(x, split = "\\|")[[1]][1]) %>% as.character
# cog2 <- cog2 %>% left_join(mouse2human) %>% na.omit
# cog2 <- cog2[,c(5,2:4)] %>% set_names("symbol","logfc","pvalue","score")
# cog2$group <- "human_cognitive2"
