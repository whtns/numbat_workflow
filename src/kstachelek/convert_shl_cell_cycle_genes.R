cell_cycle_genes_file <- "/home/sunlee/Cell_cycle_genes_updated_06.20.2018.csv"
out_file = gsub(".csv", "_w_ensembl_geneids.csv", cell_cycle_genes_file)

cc_genes <- read.csv(cell_cycle_genes_file, stringsAsFactors = FALSE)

library(cataract)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

head(cc_genes)

cc_genes$human.gene <- trimws(cc_genes$human.gene)

map_features_to_symbols <- function(symbols) {
  geneids <- ensembldb::genes(edb, filter = SymbolFilter(symbols), columns = c("gene_id"))
  trxids <- ensembldb::transcripts(edb, filter = SymbolFilter(symbols), columns = c("tx_id"))
  geneids <- as.data.frame(geneids)[c("gene_id", "symbol")]
  trxids <- as.data.frame(trxids)[c("tx_id", "symbol")]
  fullids <- merge(geneids, trxids, by = "symbol")
  # trxids <- mapIds(edb, keys = symbols, column = "TXID", keytype = "GENENAME",
                   # multiVals = "list")
  return(fullids)
}

feature_map <- map_features_to_symbols(cc_genes$human.gene)

anno_cc_genes <- dplyr::left_join(cc_genes, feature_map, by = c("human.gene" = "symbol"))

# anno_cc_genes <- anno_cc_genes[!grepl("LRG", anno_cc_genes$tx_id),]

write.csv(anno_cc_genes, out_file)


