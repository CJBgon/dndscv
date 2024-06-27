#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparser"))


######### input arguments #############
parser <- ArgumentParser()

parser$add_argument(
    "-f", "--input", type="character",
    help ='Path to the variant aggregate.'
)
parser$add_argument(
    "-o", "--out", type="character", default='/out/',
    help = 'Path to the output directory.'
)
parser$add_argument(
    "-r", "--ref", type="character",
    help = '"hg38, hg19: or path to the reference cds file.',
    default = "hg19"
)
parser$add_argument(
    "-s", "--subm", type="character",
    help = 'Substitution model (precomputed models are available in the data directory) (submod_192r_3w.rda,12r and 2r are available too.)',
    default ="192r_3w"
)
p <- add_argument(
    p,
    "--known_cancer_genes", type="character", nargs='+',
    help = 'List of a-priori known cancer genes (to be excluded from the indel background model) - references the cancer gene census v81.  ("cgc81" = cancergenes_cgc81.rda)',
    default = "cgc81"
)
parser$add_argument(
    "-c", "--covariates", type="character",
    help = 'Covariates (a matrix of covariates -columns- for each gene -rows-) [default: reference covariates] [cv=NULL runs dndscv without covariates] (hg19_hg38_epigenome_pcawg.rda)',
    default ="hg19"
)
parser$add_argument(
    "-m", "--maxmuts", type="integer",
    help = 'If n<Inf, arbitrarily the first n mutations by chr position will be kept',
    default = 3
)
parser$add_argument(
    "-t", "--maxtmb", type="integer",
    help = 'to filter out hypermutated samples.',
    default = 3000
)
parser$add_argument(
    "-i", "--indels", action="store_true", default=FALSE,
    help = 'Use unique indel sites instead of the total number of indels (it tends to be more robust)',
)
parser$add_argument(
    "-p", "--mindels", type="integer",
    help = 'Minimum number of indels required to run the indel recurrence module.',
    default = 5
)
parser$add_argument(
    "-a", "--maxcov", type="integer",
    help = 'Maximum number of covariates that will be considered (additional columns in the matrix of covariates will be excluded)',
    default = 20
)
parser$add_argument(
    "-w", "--constrain", action="store_false", default=TRUE
    help = 'This constrains wnon==wspl (this typically leads to higher power to detect selection)',
)
parser$add_argument(
    "-d", "--output_type", type="integer", default=3
    help = 'Output: 1 = Global dN/dS values; 2 = Global dN/dS and dNdSloc; 3 = Global dN/dS, dNdSloc and dNdScv',
)
parser$add_argument(
    "-g", "--gennum", type="integer", default=1
    help = 'NCBI genetic code number (default = 1; standard genetic code). To see the list of genetic codes supported use: ? seqinr::translate. Note that the same genetic code must be used in the dndscv and buildref functions.',
)
parser$add_argument(
    "-n", "--runname", type="character",
    help = 'name of the run - given to output file',
    default='dndscv'
)
# parser$add_argument(
#     "-v", "--genelist", type="character", nargs='+',
#     help = 'limit the run to this vector of genes'
# )

args <- parser$parse_args()

#### Wrangle input #####
library("dndscv")

# data("dataset_simbreast", package="dndscv") # Example dataset
# dndsout = dndscv(mutations, refdb="example_output_refcds.rda", cv=NULL)
if (args$covariates == 'covariates_hg19_hg38_epigenome_pcawg.rda') {
    load('/usr/src/app/data/covariates_hg19_hg38_epigenome_pcawg.rda')
    covs <- covs
} else {
    covs <- args$covariates
}

# solving an issue with arparser not taking optional flag and character inputs simultaneously.
# genelist: Remove blank spaces and split gene names into vector
args$genelist <- unlist(strsplit(gsub(" ", "", args$genelist), split = ",", fixed = TRUE))

# known_cancer_genes: Remove blank spaces and split gene names into vector
args$known_cancer_genes <- unlist(strsplit(gsub(" ", "", args$known_cancer_genes), split = ",", fixed = TRUE))

# Assign default value to genelist as NULL, as required by the dndscv function.
if((length(args$genelist) == 1)){
    if(toupper(args$genelist) == "NONE"){
        args$genelist <- NULL
    }
}

# read in table
mutations <- fread(args$varagg, sep='\t', header=TRUE)
# mutations <- mutations[,.(SAMPLE,CHROM,POS, REF,ALT)] # Restricting input matrix to first 5 columns
mutations <- mutations[,.(CHROM,POS, REF,ALT,SAMPLE)]
# refdb is from ensembl withouth 'chr' prefix.
setDT(mutations)[, CHROM := gsub("\\chr", "",x=CHROM)]

dndsout <- dndscv(
    args$input,  # the aggregate variantfile
    gene_list = NULL, # limit analysis to this gene list (vector)
    refdb = args$ref, # reference file to use
    sm = args$subm,  # Substitution model (precomputed models are available in the data directory) (submod_192r_3w.rda) )12r and 2r are available too.
    kc = args$kc,  # List of a-priori known cancer genes (to be excluded from the indel background model) - references the cancer gene census v81.  (cancergenes_cgc81.rda)
    cv = args$covariates, # Covariates (a matrix of covariates -columns- for each gene -rows-) [default: reference covariates] [cv=NULL runs dndscv without covariates] ()
    max_muts_per_gene_per_sample = args$maxmuts,  # If n<Inf, arbitrarily the first n mutations by chr position will be kept 
    max_coding_muts_per_sample = args$maxtmb,  # to filter out hypermutated samples.
    use_indel_sites = args$indels, # Use unique indel sites instead of the total number of indels (it tends to be more robust)
    min_indels = args$mindels, #  Minimum number of indels required to run the indel recurrence module
    maxcovs = args$maxcov, # Maximum number of covariates that will be considered (additional columns in the matrix of covariates will be excluded)
    constrain_wnon_wspl = args$constrain, # This constrains wnon==wspl (this typically leads to higher power to detect selection) nonsense=wnon, missense=wmis
    outp = args$output_type,# Output: 1 = Global dN/dS values; 2 = Global dN/dS and dNdSloc; 3 = Global dN/dS, dNdSloc and dNdScv
    numcode = args$gennum, # NCBI genetic code number (default = 1; standard genetic code). To see the list of genetic codes supported use: ? seqinr::translate. Note that the same genetic code must be used in the dndscv and buildref functions.
    outmats = FALSE, # Output the internal N and L matrices (default = F)
    mingenecovs = 500  # Minimum number of genes required to run the negative binomial regression model with covariates (default = 500), 
    dc = NULL)

write.table(
    dndsout$sel_cv, 
    file=paste0(args$outdir, '/', args$runname),
    sep='\t', 
    quote=FALSE
    )

##### TODO ###: 
# how do we reference/link the script to the rda files downloaded? their names are different:
    # i guesss we load the data and 
    # change line of code to search for hg38 in the data folder (place there when downloaded)
# how can we input other refcds files? is it the path to the .rda file?
    # (if so place into the data folder + allow specified paths.)
# are the sm, kc and cb included in the hg38 download?
# max_muts_per_gene_per_sample how does this affect the output? is it limited to non-synonymous muts only? (probably to filter false positives / repeats?)


# how can we attune this to low-frequency mutations? can we optimise the dndscv hyperparameters by including a downstream model / metric to optimise?
    # how can we link the feedback from that model back into the three tools? (would nextflow allow us to iterate on the tools and explore the searchspace?
