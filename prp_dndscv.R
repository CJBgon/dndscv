#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))


######### input arguments #############
parser <- ArgumentParser()

parser$add_argument(
    "-c", "--cds", type="string", default=NA, 
    help='Path to the CDS file.'
)
parser$add_argument(
    "-f", "--fasta", type="string", default=NA,
    help= 'Path to the reference genome associated with the variants.'
)
parser$add_argument(
    "-o", "--out", type="string", default=NA,
    help= 'Path to the output directory.'
)
parser$add_argument(
    "-e", "--exclude_contig", type="string", default='MT',
    help= 'contigs to exclude frmo the buildref (note, Mitochondrial genomes are unsuited for global dNdS calculations.'
)

args <- parser$parse_args()

## main ##
# expected CDS format:
# [gene.id, gene.name, cds.id, chr, chr.coding.start, chr.coding.end. cds.start, cds.end, length, strand]
library(dndscv)
path_cds_table = args$cds
path_genome_fasta = args$fasta
buildref(
    cdsfile=path_cds_table,
    genomefile=path_genome_fasta,
    outfile = paste0(args$out,"/output_refcds.rda"),
    excludechrs=args$exclude_contig)
