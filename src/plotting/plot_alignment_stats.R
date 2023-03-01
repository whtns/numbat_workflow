
library(glue)
library(dplyr)
library(tidyverse)

test <- read.table("output/multiqc_data_1/multiqc_hisat2.txt", header = T)



align_cats <- c("paired_aligned_one", "paired_aligned_discord_one", "unpaired_aligned_one", "paired_aligned_multi", "unpaired_aligned_multi", "paired_aligned_none")

#### Rearrange the data frame for ggplot 

test2 <- dplyr::select(test, c("Sample", align_cats)) %>% 
	dplyr::mutate(Sample = gsub(".*_S", "S", Sample)) %>%     
	dplyr::mutate(Sample = gsub("_L.*", "", Sample)) %>%        
#	dplyr::mutate(Sample = gsub(".hisat2", "", Sample)) %>% 
	dplyr::group_by(Sample) %>% 
	tidyr::gather("alignment_category", "reads", align_cats) %>% 
	dplyr::group_by(Sample,alignment_category) %>% dplyr::summarise(reads=sum(reads)) %>% 
	identity()

sample_order <- dplyr::filter(test2, alignment_category == "paired_aligned_one") %>% 
	dplyr::arrange(reads) %>% 
	dplyr::pull(Sample)

test2$Sample <- factor(test2$Sample, sample_order)

# 
# data$cell<- ordered(data$cell, levels = unique(data$cell))

plot_read_counts <- function(df, n_samples){
	ggplot(df, aes(x = Sample, y = reads, fill = alignment_category)) + 
		geom_bar(stat = "identity") + 
		ggtitle(glue(n_samples, " samples in dataset")) +
		# scale_y_log10() +
		theme(axis.text.x = element_text(size=rel(0.15),angle = 90, vjust = 1, hjust=1)) +
		NULL	
}


pdf("results/alignment_stats.pdf")
plot_read_counts(test2, "all")


test3 <- dplyr::filter(test2, Sample %in% sample_order[1:300])

plot_read_counts(test3, "bottom 300")
dev.off()


