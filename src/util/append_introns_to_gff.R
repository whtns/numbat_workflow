#!/usr/bin/Rscript


# import required librareis -----------------------------------------------

library(GenomicRanges)
library(rtracklayer)
library(parallel)


gff <- import.gff2("Mus_musculus.GRCm38.71.DEXSeq.gff", asRangedData=F)
#Add a level for "intronic_part"
elementMetadata(gff)$type <- factor(elementMetadata(gff)$type, levels=c(levels(elementMetadata(gff)$type), "intronic_part"))
#Fix the exonic_part_number to be a properly formatted character
USE <- which(!is.na(elementMetadata(gff)$exonic_part_number))
exonic_parts <- sprintf("%03i", elementMetadata(gff)$exonic_part_number[USE])
elementMetadata(gff)$exonic_part_number <- as.character(elementMetadata(gff)$exonic_part_number)
elementMetadata(gff)$exonic_part_number[USE] <- exonic_parts

#Split by gene_id
grl <- split(gff, elementMetadata(gff)$gene_id)
#Add the introns as "intronic_parts
add_introns <- function(gr) {
  exons <- gr[which(elementMetadata(gr)$type=="exonic_part"),]
  if(length(exons) > 1) {
    seqname <- seqnames(exons)[-1]
    starts <- end(exons)+1
    starts <- starts[-length(starts)]
    ends <- start(exons)-1
    ends <- ends[-1]
    bounds <- IRanges(start=starts, end=ends)
    strand <- strand(exons)[-1]
    introns <- GRanges(seqnames=seqname, ranges=bounds, strand=strand)
    intron_ids <- sprintf("%03i", c(1:length(introns)))
    #Remove 0-width introns
    DISCARD <- which(width(introns) <= 0)
    if(length(DISCARD) > 0) {
      introns <- introns[-DISCARD]
      intron_ids <- intron_ids[-DISCARD] #Set intron numbers so they follow their respective exonic parts
    }
    if(length(introns) > 0) {
      #create the meta-data
      df <- as.data.frame(elementMetadata(exons))
      nrows <- length(introns)
      metadf <- df[1:nrows,] #does this need to deal with gene_id and transcripts differently?
      metadf <- transform(metadf, gene_id=as.character(gene_id), transcripts=as.character(transcripts))
      metadf$transcripts <- as.character(c(rep(NA, nrows)))
      metadf$type <- factor(c(rep("intronic_part", nrows)), levels=levels(metadf$type))
      metadf$exonic_part_number <- intron_ids
      elementMetadata(introns) <- metadf
      #Merge the GRanges
      gr <- append(gr, introns)
      gr <- gr[order(start(gr), elementMetadata(gr)$type),] #resort
    }
  }
  return(gr)
}
with_introns <- endoapply(grl, add_introns) 
#reorder things
chroms <- sapply(with_introns, function(x) as.factor(seqnames(x))[1])
starts <- sapply(with_introns, function(x) start(x)[1])
o <- order(chroms, starts)
with_introns2 <- with_introns[o]
##Merge into a GRange
#with_introns2 <- unlist(with_introns2, use.names=F, recursive=T)
#Create GFF formatted output
asGFF2 <- function(x) {
  df <- as.data.frame(x)
  aggregates <- which(df$type == "aggregate_gene")
  meta <- character(nrow(df))
  meta[aggregates] <- sprintf("gene_id \"%s\"", df$gene_id[aggregates])
  #This gives introns a transcript "NA" field, which may not be ideal
  meta[-aggregates] <- sprintf("transcripts \"%s\"; exonic_part_number \"%s\"; gene_id \"%s\"", df$transcripts[-aggregates], df$exonic_part_number[-aggregates], df$gene_id[-aggregates])
  paste(df$seqnames, "dexseq_prepare_annotation.py", df$type, df$start, df$end, ".", df$strand, ".", meta, sep="\t")
}
outputGFF <- unlist(lapply(with_introns2, asGFF2))
write.table(outputGFF, file="Mus_musculus.GRCm38.71.DEXSeq.introns.gff", row.names=F, col.names=F, quote=F)