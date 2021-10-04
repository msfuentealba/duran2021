library(tidyverse)
library(limma)
library(qvalue)

run = rbind(tibble(run = "r1", sample_id = c("wt1","xbp1","fad1","fx1","wt2")),
            tibble(run = "r2", sample_id = c("xbp2","fad2","fx2","wt3","xbp3")),
            tibble(run = "r3", sample_id = c("fad3","fx3","wt4","xbp4","fad4")),
            tibble(run = "r4", sample_id = c("fx4","wt5","xbp5","fad5","fx5")))

xbp_wt <- readxl::read_xlsx("./info/1_edited.xlsx", col_names = TRUE, sheet = 1)
fad_wt <- readxl::read_xlsx("./info/1_edited.xlsx", col_names = TRUE, sheet = 2)
fx_fad <- readxl::read_xlsx("./info/1_edited.xlsx", col_names = TRUE, sheet = 3)

dpa <- function(control, mutant){
  print(paste0(control," - ",mutant))
  if(control=="wt"&mutant=="fad"){
    data <- fad_wt  
  } else if (control=="wt"&mutant=="xbp"){
    data <- xbp_wt  
  } else if (control=="fad"&mutant=="fx"){
    data <- fx_fad
  }
  #data <- data[rowSums(!is.na(data[,4:8]))>=3&rowSums(!is.na(data[,9:13]))>=3,]
  #data <- data[!grepl("contaminant",data$uniprot),]
  uni2sym <- data[,1:2]
  data <- data[,c(1,4:13)] %>% data.frame
  rownames(data) <- data$uniprot
  data$uniprot <- NULL
  annot <- tibble(sample_id = colnames(data), group = c(rep("control",5),rep("mutant",5)))
  annot <- annot %>% left_join(run)
  annot$group <- factor(annot$group)
  annot$run <- factor(annot$run)
  design <- model.matrix(~0+group+run,annot)
  colnames(design) <- c("control","mutant","run2","run3","run4")
  fit <- lmFit(data, design)
  contrast_matrix <- makeContrasts(contrasts = "mutant-control",levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  results <- topTable(fit2,  coef = "mutant-control", number = Inf, adjust.method = "BH")[,c(1,4,5)]
  results$uniprot <- row.names(results)
  rownames(results) <- NULL
  results <- results %>% left_join(uni2sym)
  results$qvalue <- qvalue(results$P.Value)$qvalues
  results <- results[,c(4,5,1:3,6)] %>% arrange(P.Value)
  results$score <- -log10(results$P.Value)*sign(results$logFC)
  results <- results %>% arrange(qvalue)
  colnames(results) <- c("uniprot","symbol","logfc","pvalue","fdr","qvalue","score")
  return(results)
}

grid <- tibble(control = c("wt","wt","fad"), mutant = c("xbp","fad","fx"))
grid <- grid %>% mutate(results = map2(control,mutant,dpa))
write_rds(grid, file = "./output/linear_model.rds")
