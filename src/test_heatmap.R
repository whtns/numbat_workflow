source("packages.R")
source("functions.R")


test_heatmap <- make_numbat_heatmaps_test("output/numbat_sridhar_mini/SRR13884242/done.txt", p_min = 0.5, line_width = 0.1)

joint_post <- read_tsv("output/numbat_sridhar_mini/SRR13884242/joint_post_2.tsv") %>%
  mutate(
  cnv_state = ifelse(cnv_state == 'neu', NA, cnv_state))


test_heatmap +
  scale_colour_gradientn(colours = c("purple", "orange")) +
  geom_segment(
  aes(x = seg_start, xend = seg_end, y = cell_index, yend = cell_index, color = avg_entropy, alpha = p_cnv),
  size = 0.1) +
  # labs(colour = "Purple-Orange") +
  NULL


pdf("~/tmp/test.pdf")
test_heatmap +
  aes(color = p_cnv) +
  scale_color_gradient2() +
  # scale_color_manual(
  #   values = c('amp' = 'red', 'del' = 'darkblue', 'bamp' = 'pink', 'loh' = 'darkgreen', 'bdel' = 'blue'),
  #   labels = c('amp' = 'AMP', 'del' = 'DEL', 'bamp' = 'BAMP', 'loh' = 'CNLoH', 'bdel' = 'BDEL'),
  #   limits = force,
  #   na.translate = F
  # ) +
  # scale_alpha_continuous(range = c(0,1), limits = c(0.5, 1), oob = scales::squish) +
  NULL

test_heatmap
dev.off()

browseURL("~/tmp/test.pdf")


ggplot(pd[pd$score1 != 0,], aes(x=x, y=species)) +
  geom_tile(aes(fill  =score1)) +
  scale_fill_gradient2("Score 1", limits = c(0, 4),
                       low = "#762A83", mid = "white", high = "#1B7837") +

  new_scale("fill") +

  geom_tile(aes(fill = score2), data = subset(pd, score2 != 0)) +
  scale_fill_gradient2("Score 2", limits = c(0, 3),
                       low = "#1B7837", mid = "white", high = "#762A83") +

  geom_text(data=pd, aes(label = letters, color = factor(change))) +
  scale_color_manual("Change", values = c("black", "#F2A11F"),
                     labels = c("None", "Some")) +
  coord_fixed(ratio = 1.5, xlim=c(0.5,16.5), ylim=c(0.5, 3.5))


make_numbat_heatmaps_test <- function(done_file, p_min = 0.9, line_width = 0.1){
  # browser()

  sample_id <- path_file(path_dir(done_file))

  numbat_dir = fs::path_split(done_file)[[1]][[2]]

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(glue("output/{numbat_dir}/{sample_id}_numbat.rds"))

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  myannot <- mynb$clone_post[,c("cell", "GT_opt")]
  # myannot <- mynb$clone_post[,c("cell", "GT_opt", "clone_opt")]

  ## numbat ------------------------------
  numbat_heatmap <- plot_phylo_heatmap_new(mynb$gtree, mynb$joint_post, mynb$segs_consensus, mynb$clone_post,
    show_phylo = FALSE,
    annot_bar_width = 1,
  )

  return(numbat_heatmap)
}


