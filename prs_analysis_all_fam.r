#!/usr/bin/Rscript

# arguments: list of fams text file, profile, best, splitvalidate.results.txt, output directory
args = commandArgs(trailingOnly=TRUE)

## load necessary packages
.libPaths("/cluster/home/dkim/y/envs/park_prs/lib/R/library")
library(ggplot2)
library(dplyr)

.libPaths("/cluster/home/aanderson/R_4.1")
library(ggpubr)

# reading all arguments
fam_list <- readLines(args[1])

profile <- read.delim(args[2], header=TRUE, sep="")

best <- read.delim(args[3], header=TRUE, sep="")

results <- read.delim(args[4], header=TRUE, sep="")
colnames(results) <- c("FID", "IID", "Phenotype", "order", "best.pgs", "split")

setwd(args[5])

# analyze & plot prs scores for each software w/ all desired fam files (holds phenotype)
for (fam_path in fam_list) {
  out_stem <- sub('\\.fam$', '', basename(fam_path))

  fam <- read.delim(fam_path, header=FALSE)
  colnames(fam) <- c("FID", "IID", "Father", "Mother", "Sex", "Phenotype")
  
  # get all prs scores into one data frame
  all <- select(profile, FID, IID, SCORE)
  all <- all %>% rename_at('SCORE', ~'PLINK')
  all <- cbind(all, "PRSice-2"=best[,4])
  all <- cbind(all, "lassosum"=results[,5])
  all <- cbind(all, "Phenotype"=fam[,6])
  
  # write scores to a file
  write.table(all, file=paste(out_stem, '_all_prs.tsv', sep=""), row.names=FALSE, quote=FALSE, sep='\t')
  
  # statistical analysis
  plink_ttest <- t.test(PLINK ~ Phenotype, data = all)
  prsice_ttest <- t.test(`PRSice-2` ~ Phenotype, data = all)
  lasso_ttest <- t.test(lassosum ~ Phenotype, data = all)
  
  ## density plots ##
  all <- all %>% rename_at('PRSice-2', ~'PRSice')
  all$Phenotype[all$Phenotype == "0"] <- "Control"
  all$Phenotype[all$Phenotype == "1"] <- "Case"
  
  all_filtered <- all %>% filter(Phenotype %in% c("Control", "Case"))
  
  # PLINK
  plink_dp <- ggdensity(all_filtered, x="PLINK", add="median", 
                        color="Phenotype", fill="Phenotype", palette=c("coral", "steelblue"), 
                        alpha=0.25, xlab="PRS", ylab="Density") + 
    theme(legend.position = "top", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16), 
        axis.title=element_text(size=20),
        axis.text=element_text(size=14)) + 
    labs(title = "PLINK PRS Density by Phenotype", fill = "Phenotype") + 
    annotate("text", x = quantile(all$PLINK, 0.75), y = max(density(all$PLINK)$y) + 1, 
             label = paste("p-value: ", signif(plink_ttest$p.value, digits=3)), size=6, hjust=0)
  
  # PRSice-2
  prsice_dp <- ggdensity(all_filtered, x="PRSice", add="median", 
                        color="Phenotype", fill="Phenotype", palette=c("coral", "steelblue"), 
                        alpha=0.25, xlab="PRS", ylab="Density") + 
    theme(legend.position = "top", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16), 
        axis.title=element_text(size=20),
        axis.text=element_text(size=14)) + 
    labs(title = "PRSice-2 PRS Density by Phenotype", fill = "Phenotype") + 
    annotate("text", x = quantile(all$PRSice, 0.75), y = max(density(all$PRSice)$y) + 1, 
             label = paste("p-value: ", signif(prsice_ttest$p.value, digits=3)), size=6, hjust=0)

  # lassosum
  lassosum_dp <- ggdensity(all_filtered, x="lassosum", add="median", 
                         color="Phenotype", fill="Phenotype", palette=c("coral", "steelblue"), 
                         alpha=0.25, xlab="PRS", ylab="Density") +
    theme(legend.position = "top", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16), 
        axis.title=element_text(size=20),
        axis.text=element_text(size=14)) + 
    labs(title = "lassosum PRS Density by Phenotype", fill = "Phenotype") + 
    annotate("text", x = quantile(all$lassosum, 0.75), y = max(density(all$lassosum)$y) + 1, 
             label = paste("p-value: ", signif(lasso_ttest$p.value, digits=3)), size=6, hjust=0)
  
  # save plots in png format
  ggsave(paste(out_stem, "_plink_dp.png", sep=""), plot = plink_dp)
  
  ggsave(paste(out_stem, "_prsice2_dp.png", sep=""), plot = prsice_dp)
  
  ggsave(paste(out_stem, "_lassosum_dp.png", sep=""), plot = lassosum_dp)
}
