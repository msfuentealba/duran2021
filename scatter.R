library(tidyverse)
grid <- read_rds("./output/linear_model.rds")

wt_xbp <- grid$results[[1]]
wt_fad <- grid$results[[2]]
fad_fx <- grid$results[[3]]

#wt_fad vs fad_fx opposite
inter <- intersect(wt_fad$uniprot,fad_fx$uniprot)
wt_fad <- wt_fad %>% filter(uniprot%in%inter) %>% dplyr::select(uniprot,symbol,score)
fad_fx <- fad_fx %>% filter(uniprot%in%inter) %>% dplyr::select(uniprot,symbol,score)

comb <- wt_fad %>% left_join(fad_fx, by = c("uniprot"="uniprot")) %>% unique 
comb$dir <- ifelse(comb$score.x<0&comb$score.y>0,"dn_up",
            ifelse(comb$score.x>0&comb$score.y<0,"up_dn",NA))
sum(!is.na(comb$dir))/nrow(comb) # 76% going in opposite directions

#source("./backup/rankprodbounds.R")

get_genes <- function(d){
  #rank wt_fad
  subset <- comb %>% filter(dir==d)
  subset <- subset %>% arrange(desc(abs(score.x)))
  subset$rank.x <- 1:nrow(subset)
  
  #rank fad_fx
  subset <- subset %>% arrange(desc(abs(score.y)))
  subset$rank.y <- 1:nrow(subset)
  
  n <- 1000000
  sim_prod <- replicate(n,prod(c(sample(1:nrow(subset),1),sample(1:nrow(subset),1))))
  subset$prod <- sapply(1:nrow(subset), function(x) prod(c(subset$rank.x[x],subset$rank.y[x])))
  subset$pval <- sapply(subset$prod, function(x) sum(x>=sim_prod)/n)
  subset <- subset %>% arrange(pval)
  return(subset)
}

dn_up <- get_genes("dn_up")
head(dn_up)
nrow(dn_up)
dn_up$score.z <- (-log10(dn_up$pval)*1)
hist(dn_up$score.z)
up_dn <- get_genes("up_dn")
nrow(up_dn)
up_dn$score.z <- (-log10(up_dn$pval)*-1)
hist(up_dn$score.z)
intersect(up_dn$uniprot,dn_up$uniprot)

scater <- comb %>% dplyr::select(symbol.x,score.x,score.y) %>% set_names("symbol","WT vs 5xFAD", "5xFAD vs TgXBP1/5xFAD")
cor.test(scater$`WT vs 5xFAD`,scater$`5xFAD vs TgXBP1/5xFAD`, method = "spearman") # double check

#scater$label <- ifelse(scater$symbol%in%c(dn_up$symbol.x[1:10],up_dn$symbol[1:10]),scater$symbol,"")
scater$dir <- ifelse(scater$symbol%in%c(dn_up$symbol.x,up_dn$symbol.x),"reverse","mimic")

#pairs_cols <- brewer.pal(10, "Paired")
pairs_cols <- c("grey80","grey10")

library(ggrepel)
pdf(file = paste0("./output/figures/scatter.pdf"), height = 6, width = 6) 
ggplot(scater, aes(x = `WT vs 5xFAD`, y = `5xFAD vs TgXBP1/5xFAD`, color = dir)) +
  geom_point() +
  scale_color_manual(values = pairs_cols[1:2])+
  theme_bw()+
  #geom_label_repel(aes(label = label), size = 3)+
  theme(legend.position = "none",
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14))
dev.off()

gsea <- rbind(dn_up[,c(1,11)],up_dn[,c(1,11)])
length(gsea$uniprot)/length(unique(gsea$uniprot))
genelist <- gsea$score.z
names(genelist) <- gsea$uniprot
genelist <- sort(genelist) %>% rev
head(genelist)

library(clusterProfiler)
#kegg <- gseKEGG(genelist,
#        organism = "mmu",
#        keyType = "uniprot",
#        minGSSize = 10,
#        pvalueCutoff = 1,
#        maxGSSize = 50)@result

library(org.Mm.eg.db)
go <- gseGO(genelist,
            ont = "BP",
            OrgDb = org.Mm.eg.db,
            keyType = "UNIPROT",
            minGSSize = 10,
            maxGSSize = 200,
            by = 'fgsea')@result
colnames(go)

go <- go %>% filter(p.adjust<0.05)
go_up <- go %>% filter(NES>0)
go_up <- go_up[1:10,c(1:2,5:7,11)]
rownames(go_up) <- NULL
#xlsx::write.xlsx(go_up[,1:5], file = "./output/go_up.xlsx", row.names = FALSE)
#write_csv2(go_up[,1:5], file = "./output/go_up.csv")
library(gridExtra)
library(grid)
go_up$NES <- round(go_up$NES,2)
go_up$pvalue <- formatC(go_up$pvalue, format = "e", digits = 2)
go_up$p.adjust <- formatC(go_up$p.adjust, format = "e", digits = 2)
colnames(go_up) <- c("GO ID","Name","NES","p-value","adjusted\np-value","core_enrichment")
pdf(file = paste0("./output/figures/go_up.pdf"), height = 4, width = 10) 
grid.table(go_up[,1:5], rows = NULL)
dev.off()

go_up$core_enrichment <- sapply(go_up$core_enrichment, function(x) strsplit(x, split = "\\/"))
top_go_up <- go_up[,c(2,6)] %>% unnest(cols = c(core_enrichment))
top_go_up <- top_go_up %>% left_join(comb[,1:2], by = c("core_enrichment"="uniprot"))
top_go_up$value <- 1
table(top_go_up$symbol.x)
top_go_up <- reshape2::acast(top_go_up, symbol.x~Name, value.var = "value")
top_go_up[is.na(top_go_up)] <- 0
top_genes_up <- rownames(top_go_up)[rowSums(top_go_up)==9]

xbp_wt <- readxl::read_xlsx("./info/1_edited.xlsx", col_names = TRUE, sheet = 1)
fx_fad <- readxl::read_xlsx("./info/1_edited.xlsx", col_names = TRUE, sheet = 3)

raw <- xbp_wt[,c(2,4:13)] %>% left_join(fx_fad[,c(2,4:13)])
raw <- raw %>% filter(symbol%in%top_genes_up)
raw <- reshape2::melt(raw) %>% na.omit
raw$group <- gsub("1|2|3|4|5","",raw$variable)
raw$group <- factor(raw$group, levels = c("wt","xbp","fad","fx"))
pdf(file = paste0("./output/figures/top_genes_up.pdf"), height = 5, width = 7) 
ggplot(raw, aes(x=group, y=as.numeric(value), fill = group)) + 
  #geom_boxplot(aes(fill=group)) + 
  stat_summary(fun.y = mean, na.rm = TRUE, 
               geom = "point", shape = "diamond",
               size = 3, alpha = 0.5,
               position = position_dodge(width = .7)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", width = .2, alpha = 0.5,
               position = position_dodge(width = .7)) +
  geom_point(position = position_jitterdodge(jitter.width = .2, dodge.width = .7), alpha = 1, shape = 21) +
  scale_fill_manual(values = c("#fefefe","#4d4c4b","#00007a","#3a7d22"), labels = c("WT", "TgXBP1", "5xFAD","TgXBP1/5xFAD")) + 
  facet_wrap(. ~ symbol, scales = "free", nrow = 4, ncol = 4) +
  labs(fill = "Group")+
  ylab("Normalized protein abundance")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

library(circlize)
library(ComplexHeatmap)
col_fun <- colorRamp2(c(0,1),c("white","grey50"))

pdf(file = paste0("./output/figures/top_go_up.pdf"), height = 15, width = 4) 
Heatmap(top_go_up,
        col = col_fun,
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA, lty = 3, lwd = 1))
          grid.circle(x = x, y = y, r = top_go_up[i, j]*2 * min(unit.c(width, height)),gp = gpar(fill = col_fun(top_go_up[i, j]), col = NA))
        },
        border = TRUE)
dev.off()

go_dn <- go %>% filter(NES<0)
go_dn <- go_dn[1:10,c(1:2,5:7,11)]
rownames(go_dn) <- NULL
go_dn$NES <- round(go_dn$NES,2)
go_dn$pvalue <- formatC(go_dn$pvalue, format = "e", digits = 2)
go_dn$p.adjust <- formatC(go_dn$p.adjust, format = "e", digits = 2)
colnames(go_dn) <- c("GO ID","Name","NES","p-value","adjusted\np-value","core_enrichment")
pdf(file = paste0("./output/figures/go_dn.pdf"), height = 4, width = 10) 
grid.table(go_dn[,1:5], rows = NULL)
dev.off()

#xlsx::write.xlsx2(go_dn, file = "./output/go_dn.xlsx", row.names = FALSE)
go_dn$core_enrichment <- sapply(go_dn$core_enrichment, function(x) strsplit(x, split = "\\/"))
top_go_dn <- go_dn[,c(2,6)] %>% unnest(cols = c(core_enrichment))
top_go_dn <- top_go_dn %>% left_join(comb[,1:2], by = c("core_enrichment"="uniprot"))
top_go_dn$value <- 1
table(top_go_dn$symbol.x)
top_go_dn <- reshape2::acast(top_go_dn, symbol.x~Name, value.var = "value")
top_go_dn[is.na(top_go_dn)] <- 0
top_genes_dn <- rownames(top_go_dn)[rowSums(top_go_dn)>=8]
raw <- xbp_wt[,c(2,4:13)] %>% left_join(fx_fad[,c(2,4:13)])
raw <- raw %>% filter(symbol%in%top_genes_dn)
raw <- reshape2::melt(raw) %>% na.omit
raw$group <- gsub("1|2|3|4|5","",raw$variable)
raw$group <- factor(raw$group, levels = c("wt","xbp","fad","fx"))
pdf(file = paste0("./output/figures/top_genes_dn.pdf"), height = 5, width = 7) 
ggplot(raw, aes(x=group, y=as.numeric(value), fill = group)) + 
  #geom_boxplot(aes(fill=group)) + 
  stat_summary(fun.y = mean, na.rm = TRUE, 
               geom = "point", shape = "diamond",
               size = 3, alpha = 0.5,
               position = position_dodge(width = .7)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", width = .2, alpha = 0.5,
               position = position_dodge(width = .7)) +
  geom_point(position = position_jitterdodge(jitter.width = .2, dodge.width = .7), alpha = 1, shape = 21) +
  scale_fill_manual(values = c("#fefefe","#4d4c4b","#00007a","#3a7d22"), labels = c("WT", "TgXBP1", "5xFAD","TgXBP1/5xFAD")) + 
  facet_wrap(. ~ symbol, scales = "free", nrow = 4, ncol = 4) +
  labs(fill = "Group")+
  ylab("Normalized protein abundance")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()


col_fun <- colorRamp2(c(0,1),c("white","grey50"))

pdf(file = paste0("./output/figures/top_go_dn.pdf"), height = 15, width = 4) 
Heatmap(top_go_dn,
        col = col_fun,
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA, lty = 3, lwd = 1))
          grid.circle(x = x, y = y, r = top_go_dn[i, j]*2 * min(unit.c(width, height)),gp = gpar(fill = col_fun(top_go_dn[i, j]), col = NA))
        },
        border = TRUE)
dev.off()



