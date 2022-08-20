#' Methylation Plot
#'
#' Represents the methylation levels and the detected DMRs present in a defined region of the genome for a set of samples.
#' @author Raul Sanz, \email{raulsanzr@gmail.com}
#' @param genome Available reference genomes: hg19, hg38 and mm39.
#' @param chr Chromosome.
#' @param start First nucleotide position (included).
#' @param end Last nucleotide position (included).
#' @param sites GRanges containing the CpG sites and its methylation values associated.
#' @param regions GRanges or data.frame containing the coordinates of the detected DMRs.
#' @param enhancers (optional) data.frame containing the coordinates of enhancers.
#' @param group (optional) Variable to group the samples. 
#' @examples 
#' genome <- "hg19"
#' chr <- "chr11"
#' start <- 27015473
#' end <- 27015991
#' group <- c("control", "cond_A", "cond_A", "cond_B", "control", "cond_A", ...)
#' @export
methplot <- function(genome, chr, start, end, sites, regions, enhancers, group){
  # genomic coordinates
  gtrack <- Gviz::GenomeAxisTrack()
  
  # chromosome representation
  itrack <- Gviz::IdeogramTrack(genome=genome, chromosome=chr)
  
  # building the gene model based on the selected reference genome.
  if(genome=="hg19"){ # human reference genome hg19 release.
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    gene_model <- Gviz::GeneRegionTrack(txdb, genome=genome, chromosome=chr, showId=TRUE, geneSymbol=TRUE, name="UCSC")
    symbols <- unlist(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, Gviz::gene(gene_model), "SYMBOL", "ENTREZID", multiVals="first"))
    Gviz::symbol(gene_model) <- symbols[Gviz::gene(gene_model)]
  } else if(genome=="hg38"){ # human reference genome GRCh38 release.
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    gene_model <- Gviz::GeneRegionTrack(txdb, genome=genome, chromosome=chr, showId=TRUE, geneSymbol=TRUE, name="UCSC")
    symbols <- unlist(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, Gviz::gene(gene_model), "SYMBOL", "ENTREZID", multiVals="first"))
    Gviz::symbol(gene_model) <- symbols[Gviz::gene(gene_model)]
  } else if(genome=="mm39"){
    txdb <- TxDb.Mmusculus.UCSC.mm39.refGene::TxDb.Mmusculus.UCSC.mm39.refGene
    # gene model
    gene_model <- Gviz::GeneRegionTrack(txdb, genome=genome, chromosome=chr, showId=TRUE, geneSymbol=TRUE, name="UCSC")
    symbols <- unlist(AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, Gviz::gene(gene_model), "SYMBOL", "ENTREZID", multiVals="first"))
    Gviz::symbol(gene_model) <- symbols[Gviz::gene(gene_model)]
  } else{
    stop(paste0(genome, " is not a valid reference genome. Check the compatible genomes."))
  }
  
  # detected DMRs
  if(class(regions) == "data.frame"){
    DMR_track <- Gviz::AnnotationTrack(start=c(regions$start), end=c(regions$end), chromosome=chr, name="DMRs", col="green", fill="green") 
  } else if(class(regions) == "GRanges"){
    DMR_track <- Gviz::AnnotationTrack(start=c(regions@ranges@start), end=c(regions@ranges@end), chromosome=chr, name="DMRs", col="green", fill="green") 
  }
  
  # heatmap of the methylation values at every CpG site
  heatmap <- Gviz::DataTrack(sites, name=" ",chromosome = chr, type="heatmap", showSampleNames=T, cex.sampleNames=0.7, 
                       gradient=c(colorRampPalette(c("blue", "white", "red"))(n = 500)), separator=2)
  
  # average methylation level per group
  if(missing(group)){ # if no group specified
    methylation <- Gviz::DataTrack(sites, name="Methylation", chromosome=chr, type="a", groups=colnames(sites@elementMetadata))
  } else{
    methylation <- Gviz::DataTrack(sites, name="Methylation", chromosome=chr, type="a", groups=group)
  }
  
  if(missing(enhancers)){ # no enhancers data
    Gviz::plotTracks(list(itrack, gtrack, gene_model, DMR_track, heatmap, methylation), from=start, to=end, 
                     extend.left=0.1, extend.right=0.1, sizes=c(2,2,5,2,10,5))
  } else{ # including enhancers track
    enh <- AnnotationTrack(start=c(enhancers$start), end=c(enhancers$end), chromosome=chr, name="enh", fill="darkgreen", col="darkgreen")
    Gviz::plotTracks(list(itrack, gtrack, gene_model, enh, DMR_track, heatmap, methylation), from=start, to=end, 
                     extend.left=0.1, extend.right=0.1, sizes=c(2,2,5,2,2,10,5))
  }
}