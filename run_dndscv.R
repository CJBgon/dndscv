#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparser"))


######### input arguments #############
p <- arg_parser("identify coding genes under positive selection")

p <- add_argument(
    p,
    "--varagg", type="character",
    help ='Path to the variant aggregate.'
)
p <- add_argument(
    p,
    "--out", type="character", default='/out/',
    help = 'Path to the output directory.'
)
p <- add_argument(
    p,
    "--ref", type="character",
    help = '"hg38, hg19: or path to the reference cds file.',
    default = "hg19"
)
p <- add_argument(
    p,
    "--subm", type="character",
    help = 'Substitution model (precomputed models are available in the data directory) (submod_192r_3w.rda,12r and 2r are available too.)',
    default ="192r_3w"
)
p <- add_argument(
    p,
    "--known_cancer_genes", type="character",
    help = 'List of a-priori known cancer genes (to be excluded from the indel background model) - references the cancer gene census v81.  ("cgc81" = cancergenes_cgc81.rda)'
)
p <- add_argument(
    p,
    "--covariates", type="character",
    help = 'Covariates (a matrix of covariates -columns- for each gene -rows-) [default: reference covariates] [cv=NULL runs dndscv without covariates] (hg19_hg38_epigenome_pcawg.rda)',
    default ="hg19"
)
p <- add_argument(
    p,
    "--maxmuts", type="integer",
    help = 'If n<Inf, arbitrarily the first n mutations by chr position will be kept',
    default = 3
)
p <- add_argument(
    p,
    "--tmbmax", type="integer",
    help = 'to filter out hypermutated samples.',
    default = 3000
)
p <- add_argument(
    p,
    "--indels", flag=TRUE,
    help = 'Use unique indel sites instead of the total number of indels (it tends to be more robust)',
)
p <- add_argument(
    p,
    "--mindels", type="integer",
    help = 'Minimum number of indels required to run the indel recurrence module.',
    default = 5
)
p <- add_argument(
    p, 
    "--amaxcov", type="integer",
    help = 'Maximum number of covariates that will be considered (additional columns in the matrix of covariates will be excluded)',
    default = 20
)
p <- add_argument(
    p,
    "--wnon_constrain", flag=TRUE,
    help = 'This constrains wnon==wspl (this typically leads to higher power to detect selection)',
)
p <- add_argument(
    p,
    "--loc_sv_type", type="integer", default=3,
    help = 'Output: 1 = Global dN/dS values; 2 = Global dN/dS and dNdSloc; 3 = Global dN/dS, dNdSloc and dNdScv',
)
p <- add_argument(
    p,
    "--gennum", type="integer", default=1,
    help = 'NCBI genetic code number (default = 1; standard genetic code). To see the list of genetic codes supported use: ? seqinr::translate. Note that the same genetic code must be used in the dndscv and buildref functions.',
)
p <- add_argument(
    p,
    "--name", type="character",
    help = 'name of the run - given to output file',
    default='dndscv'
)
p <- add_argument(
    p,
    "--genelist", type="character", nargs=Inf,
    help = 'limit the run to this vector of genes',
    default=NULL
)
p <- add_argument(
    p,
    "--mingenecovs", type="integer",
    help = 'Minimum number of genes required to run the negative binomial regression model with covariates (default = 500)',
    default=500
)

args <- parse_args(p)

#### Wrangle input #####
library("dndscv")
library(data.table)
# data("dataset_simbreast", package="dndscv") # Example dataset
# dndsout = dndscv(mutations, refdb="example_output_refcds.rda", cv=NULL)
if (args$covariates == 'covariates_hg19_hg38_epigenome_pcawg.rda') {
    load('/usr/src/app/data/covariates_hg19_hg38_epigenome_pcawg.rda')
    covs <- covs
} else {
    covs <- args$covariates
}


# read in table
mutations <- fread(args$varagg, sep='\t', header=TRUE)
# mutations <- mutations[,.(SAMPLE,CHROM,POS, REF,ALT)] # Restricting input matrix to first 5 columns
mutations <- mutations[,.(CHROM,POS, REF,ALT,SAMPLE)]
# refdb is from ensembl withouth 'chr' prefix.
setDT(mutations)[, CHROM := gsub("\\chr", "",x=CHROM)]

dndsout <- dndscv(
    mutations,  # the aggregate variantfile
    gene_list = args$genelist, # limit analysis to this gene list (vector)
    refdb = args$ref, # reference file to use
    sm = args$subm,  # Substitution model (precomputed models are available in the data directory) (submod_192r_3w.rda) )12r and 2r are available too.
    kc = args$known_cancer_genes,  # List of a-priori known cancer genes (to be excluded from the indel background model) - references the cancer gene census v81.  (cancergenes_cgc81.rda)
    cv = covs, # Covariates (a matrix of covariates -columns- for each gene -rows-) [default: reference covariates] [cv=NULL runs dndscv without covariates] ()
    max_muts_per_gene_per_sample = args$maxmuts,  # If n<Inf, arbitrarily the first n mutations by chr position will be kept 
    max_coding_muts_per_sample = args$tbmmax,  # to filter out hypermutated samples.
    use_indel_sites = args$indels, # Use unique indel sites instead of the total number of indels (it tends to be more robust)
    min_indels = args$mindels, #  Minimum number of indels required to run the indel recurrence module
    maxcovs = args$amaxcov, # Maximum number of covariates that will be considered (additional columns in the matrix of covariates will be excluded)
    constrain_wnon_wspl = args$wnon_constrain, # This constrains wnon==wspl (this typically leads to higher power to detect selection) nonsense=wnon, missense=wmis
    outp = args$loc_sv_type,# Output: 1 = Global dN/dS values; 2 = Global dN/dS and dNdSloc; 3 = Global dN/dS, dNdSloc and dNdScv
    numcode = args$gennum, # NCBI genetic code number (default = 1; standard genetic code). To see the list of genetic codes supported use: ? seqinr::translate. Note that the same genetic code must be used in the dndscv and buildref functions.
    outmats = FALSE, # Output the internal N and L matrices (default = F)
    mingenecovs = args$mingenecovs,  # Minimum number of genes required to run the negative binomial regression model with covariates (default = 500), 
    dc = NULL) # Duplex coverage per gene. Named Numeric Vector with values reflecting the mean duplex coverage per site per gene, and names corresponding to gene names. Use this argument only when running dNdScv on duplex sequencing data to use gene coverage in the offset of the regression model 

write.table(
    dndsout$sel_cv,
    file=paste0(args$out, '/', args$name),
    sep='\t',
    quote=FALSE
    )

##### TODO ###: 
# how do we reference/link the script to the rda files downloaded? their names are different:
    # i guesss we load the data and 
    # change line of code to search for hg38 in the data folder (place there when downloaded)
# how can we input other refcds files? is it the path to the .rda file?
    # (if so place into the data folder + allow specified paths.)

# max_muts_per_gene_per_sample how does this affect the output? is it limited to non-synonymous muts only? (probably to filter false positives / repeats?)


# how can we attune this to low-frequency mutations? can we optimise the dndscv hyperparameters by including a downstream model / metric to optimise?
    # how can we link the feedback from that model back into the three tools? (would nextflow allow us to iterate on the tools and explore the searchspace?
