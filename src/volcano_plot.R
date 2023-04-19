library(EnhancedVolcano)

volcano_plot_clones <- function(clone_diffex, sample_id){

  myres <- clone_diffex %>%
    split(.$clone_comparison) %>%
    map(tibble::column_to_rownames, "symbol") %>%
    identity()

  make_volcano <- function(myres, clone_comparison){
    EnhancedVolcano(myres,
                    lab = rownames(myres),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    FCcutoff = 0.5) +
      labs(title = sample_id, subtitle = clone_comparison)
  }

  test0 <- imap(myres, make_volcano)

  return(test0)
}

clone_comparison_volcanos <- imap(total_diffex_clones, volcano_plot_clones)

pdf("~/tmp/clone_comparisons.pdf")
clone_comparison_volcanos
dev.off()

volcano_plot_clone_clusters <- function(clone_diffex, sample_id){
  # browser()
  myres <- clone_diffex %>%
    group_by(clone_comparison, cluster)

  myres <-
    myres %>%
    group_split() %>%
    map(tibble::column_to_rownames, "symbol") %>%
    identity()

  make_volcano <- function(myres){
    # browser()
    mytitle = sample_id
    mysubtitle = glue("{unique(myres$clone_comparison)}_{unique(myres$cluster)}")

    EnhancedVolcano(myres,
                    lab = rownames(myres),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    FCcutoff = 0.5) +
      labs(title = mytitle, subtitle = mysubtitle)
  }

  test0 <- map(myres, make_volcano)

  return(test0)
}

clone_cluster_comparison_volcanos <- imap(cluster_diffex_clones, volcano_plot_clone_clusters)

pdf("~/tmp/clone_cluster_comparisons.pdf")
clone_cluster_comparison_volcanos
dev.off()


