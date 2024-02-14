library(data.table)
library(GenomicRanges)
library(seqinr)
# muts <- '/Users/christianbouwens/Documents/resources/containers/dndscv-master/synth_vcf/fake_pat_variants.bed'
muts <- '/Users/christianbouwens//Documents/resources/containers/sym_vcf/bed_files/large_fake_pat_variants.bed'
mutations <- fread(muts, sep='\t', header=TRUE)
colnames(mutations) <- c("CHROM","POS","REF","ALT",'SAMPLE')
# mutations <- mutations[,.(SAMPLE,CHROM,POS, REF,ALT)] # Restricting input matrix to first 5 columns
mutations <- mutations[,.(CHROM,POS, REF,ALT,SAMPLE)]
mutations[, c("CHROM", "REF", "ALT",'SAMPLE') := lapply(.SD, as.character), .SDcols = c("CHROM", "REF", "ALT",'SAMPLE')] # Factors to character
mutations[, POS := as.numeric(POS)]
colnames(mutations) <- c("chr","pos","ref","mut",'sampleID')
# mutations[[3]] = as.numeric(mutations[[3]]) # Chromosome position as numeric
mutations <- mutations[ref != mut] # Removing mutations with identical reference and mutant base
# colnames(mutations) = c("sampleID","chr","pos","ref","mut")
library(stringi)


setDT(mutations)[, chr := gsub("\\chr", "",x=chr)]



# TODO INPUT REFDB AS A VARIABLE.
refdb='/Users/christianbouwens/Documents/resources/containers/dndscv-master/data/RefCDS_human_GRCh38_GencodeV18_recommended.rda'
    # [Input] Reference database
refdb_class = class(refdb)
if ("character" %in% refdb_class) {
    if (refdb == "hg19") {
        data("refcds_hg19", package="dndscv")
        if (any(gene_list=="CDKN2A")) { # Replace CDKN2A in the input gene list with two isoforms
            gene_list = unique(c(setdiff(gene_list,"CDKN2A"),"CDKN2A.p14arf","CDKN2A.p16INK4a"))
        }
    } else {
        load(refdb)
    }
} else if("array" %in% refdb_class) {
    # use the user-supplied RefCDS object
    RefCDS = refdb
} else {
    stop("Expected refdb to be \"hg19\", a file path, or a RefCDS-formatted array object.")
}


# Removing NA entries from the input mutation table
indna = which(is.na(mutations),arr.ind=T)
if (nrow(indna)>0) {
    mutations = mutations[-unique(indna[,1]),] # Removing entries with an NA in any row
    warning(sprintf("%0.0f rows in the input table contained NA entries and have been removed. Please investigate.",length(unique(indna[,1]))))
}

# [Input] Gene list (The user can input a gene list as a character vector)
if (is.null(gene_list)) {
    gene_list = sapply(RefCDS, function(x) x$gene_name) # All genes [default]
} else { # Using only genes in the input gene list
    allg = sapply(RefCDS,function(x) x$gene_name)
    nonex = gene_list[!(gene_list %in% allg)]
    if (length(nonex)>0) { stop(sprintf("The following input gene names are not in the RefCDS database: %s", paste(nonex,collapse=", "))) }
    RefCDS = RefCDS[allg %in% gene_list] # Only input genes
    gr_genes = gr_genes[gr_genes$names %in% gene_list] # Only input genes
}


# [Input] Duplex coverage (dc argument)
if (!is.null(dc)) {
    if (is.vector(dc)) {
        if (any(is.na(dc)) | any(dc<=0)) {
            stop(sprintf("Invalid duplex coverage for some genes, please revisit your use of the dc argument: %s.", paste(names(dc[is.na(dc)]), collapse=", ")))
        } else {
            dc = dc[sapply(RefCDS, function(x) x$gene_name)] # Duplex coverage vector
            dc = dc / mean(dc, na.rm=T) # Normalising values relative to the mean
            Lallgenes = array(sapply(RefCDS, function(x) x$L), dim=c(192,4,length(RefCDS))) # saving original L matrices
            for (j in 1:length(RefCDS)) {
                RefCDS[[j]]$L = RefCDS[[j]]$L * dc[j] # Exposures relative to dc
            }
        }
    } else {
        stop("dc must be a Named Numeric Vector with the vector values corresponding to the mean duplex coverage per site per gene, and with the vector names corresponding to gene names")
    }
}

RefCDS <- lapply(RefCDS, function(x) {
  x$seq_cds = strsplit(as.character(x$seq_cds), split="")[[1]]
  x$seq_cds1up = strsplit(as.character(x$seq_cds1up), split="")[[1]]
  x$seq_cds1down = strsplit(as.character(x$seq_cds1down), split="")[[1]]
  
  if (!is.null(x$seq_splice)) {
    x$seq_splice = strsplit(as.character(x$seq_splice), split="")[[1]]
    x$seq_splice1up = strsplit(as.character(x$seq_splice1up), split="")[[1]]
    x$seq_splice1down = strsplit(as.character(x$seq_splice1down), split="")[[1]]
  }
  
  return(x)
})

## 2. Mutation annotation
message("[2] Annotating the mutations...")

# Subfunction: obtaining the codon positions of a coding mutation given the exon intervals
chr2cds = function(pos,cds_int,strand) {
    if (strand==1) {
        return(which(unlist(apply(cds_int, 1, function(x) x[1]:x[2])) %in% pos))
    } else if (strand==-1) {
        return(which(rev(unlist(apply(cds_int, 1, function(x) x[1]:x[2]))) %in% pos))
    }
}


annot_variant_impact_dt <- function(mutations) {
    setDT(mutations)
    mutations[, `:=` (
        impact = NA_character_,
        impind = NA_integer_,
        ref3_cod = NA_character_,
        mut3_cod = NA_character_,
        aachange = NA_character_,
        ntchange = NA_character_,
        codonsub = NA_character_,
        wrong_ref = FALSE
    )]
  
    mutations[, {
        chr = chr
        pos = pos
        ref = ref
        mut = mut
        sampleID = sampleID
        start = start
        end = end
        geneind = geneind
        gene = gene
        strand = strand
        ref_cod = ref_cod
        mut_cod = mut_cod
    # Essential Splice-site substitution
    if (pos %in% RefCDS[[geneind]]$intervals_splice) {
        pos_ind = (pos == RefCDS[[geneind]]$intervals_splice)
        cdsnt = RefCDS[[geneind]]$seq_splice[pos_ind]
        if (ref_cod != as.character(cdsnt)) {
                wrong_ref = TRUE
            } else {
                wrong_ref = FALSE
            }
        list(
            impact = "Essential_Splice",
            impind = 4,
            ref3_cod = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], RefCDS[[geneind]]$seq_splice[pos_ind], RefCDS[[geneind]]$seq_splice1down[pos_ind]),
            mut3_cod = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], mut_cod, RefCDS[[geneind]]$seq_splice1down[pos_ind]),
            aachange = ".",
            ntchange = ".",
            codonsub = ".",
            wrong_ref = wrong_ref
        )
    } else {
        # Coding substitution
        pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, strand)
        if (length(pos_ind) > 0) {
            cdsnt = RefCDS[[geneind]]$seq_cds[pos_ind]
            ref3_cod = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], RefCDS[[geneind]]$seq_cds[pos_ind], RefCDS[[geneind]]$seq_cds1down[pos_ind])
            mut3_cod = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], mut_cod, RefCDS[[geneind]]$seq_cds1down[pos_ind])
            codon_pos = c(ceiling(pos_ind/3)*3-2, ceiling(pos_ind/3)*3-1, ceiling(pos_ind/3)*3)
            old_codon = as.character(as.vector(RefCDS[[geneind]]$seq_cds[codon_pos]))
            pos_in_codon = pos_ind-(ceiling(pos_ind/3)-1)*3
            new_codon = old_codon; new_codon[pos_in_codon] = mut_cod
            old_aa = seqinr::translate(old_codon, numcode = numcode)
            new_aa = seqinr::translate(new_codon, numcode = numcode)
            aachange = sprintf('%s%0.0f%s',old_aa,ceiling(pos_ind/3),new_aa)
            ntchange = sprintf('%s%0.0f%s',ref_cod, pos_ind, mut_cod)
            codonsub = sprintf('%s>%s',paste(old_codon,collapse=""),paste(new_codon,collapse=""))
            
            # Determine the impact based on amino acid changes
            if (new_aa == old_aa){
                impact = "Synonymous"; impind = 1
            } else if (new_aa == "*"){
                impact = "Nonsense"; impind = 3
            } else if (old_aa != "*"){
                impact = "Missense"; impind = 2
            } else if (old_aa=="*") {
                impact = "Stop_loss"; impind = as.double(NA)
            }
            
            if (ref_cod != as.character(cdsnt)) {
                wrong_ref = TRUE
            } else {
                wrong_ref = FALSE
            }

            list(
            impact = impact,
            impind = impind,
            ref3_cod = ref3_cod,
            mut3_cod = mut3_cod,
            aachange = aachange,
            ntchange = ntchange,
            codonsub = codonsub,
            wrong_ref = wrong_ref
            )
        } else {
            list(
            impact = NA_character_,
            impind = NA_integer_,
            ref3_cod = NA_character_,
            mut3_cod = NA_character_,
            aachange = NA_character_,
            ntchange = NA_character_,
            codonsub = NA_character_,
            wrong_ref = NA
            )
        }
        }
    }, by = .(chr,pos, ref, mut, sampleID, start, end, geneind, gene, strand,ref_cod, mut_cod)]
}

annot_indels <- function(indels) {
    setDT(indels)
    indels[, `:=` (
        ref_cod = NA_character_,
        mut_cod = NA_character_,
        impact = NA_character_,
        ref3_cod = NA_character_,
        mut3_cod = NA_character_,
        aachange = NA_character_,
        ntchange = NA_character_,
        codonsub = NA_character_,
        pid=NA_character_
    )]
  
    indels[, {
        chr = chr
        ref = ref
        mut = mut
        sampleID = sampleID
        start = start
        end = end
        geneind = geneind
        gene = gene
        strand = strand
        ref_cod = "."
        mut_cod = "."
        ref3_cod = "."
        mut3_cod = "."
        aachange = "."
        ntchange = "."
        codonsub = "."
        impact = "no-SNV"
        pid = RefCDS[[geneind]]$protein_id
        pos = start:end
        ins = nchar(gsub("-","",ref))<nchar(gsub("-","",mut))
        del = nchar(gsub("-","",ref))>nchar(gsub("-","",mut))
        multisub = nchar(gsub("-","",ref))==nchar(gsub("-","",mut))
        l = nchar(gsub("-","",ref))-nchar(gsub("-","",mut))
        if (ins) { 
            pos = c(pos-1,pos)
        }
        pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, strand)
        if (length(pos_ind)>0) {
            inframe = (length(pos_ind) %% 3) == 0
            if (ins) { # Insertion
                ntchange = sprintf("%0.0f-%0.0f-ins%s",min(pos_ind),max(pos_ind),c("frshift","inframe")[inframe+1])
            } else if (del) { # Deletion
                ntchange = sprintf("%0.0f-%0.0f-del%s",min(pos_ind),max(pos_ind),c("frshift","inframe")[inframe+1])
            } else { # Dinucleotide and multinucleotide changes (MNVs)
                ntchange = sprintf("%0.0f-%0.0f-mnv",min(pos_ind),max(pos_ind))
            }
            list(
                mut_cod = mut_cod,
                ref_cod = ref_cod,
                impact = impact,
                ref3_cod = ref3_cod,
                mut3_cod = mut3_cod,
                aachange = aachange,
                codonsub = codonsub,
                pid = pid,
                ntchange = ntchange
            )
        } else {
            list(
                mut_cod = mut_cod,
                ref_cod = ref_cod,
                impact = impact,
                ref3_cod = ref3_cod,
                mut3_cod = mut3_cod,
                aachange = aachange,
                codonsub = codonsub,
                pid = pid,
                ntchange = ntchange
            )
        }
    }, by = .(chr,pos, ref, mut, sampleID, start, end, geneind, gene, strand)]
}

start.time <- Sys.time()
nt = c("A","C","G","T")
trinucs = paste(rep(nt,each=16,times=1),rep(nt,each=4,times=4),rep(nt,each=1,times=16), sep="")
trinucinds = setNames(1:64, trinucs)
trinucsubs = NULL
for (j in 1:length(trinucs)) {
    trinucsubs = c(trinucsubs, paste(trinucs[j], paste(substr(trinucs[j],1,1), setdiff(nt,substr(trinucs[j],2,2)), substr(trinucs[j],3,3), sep=""), sep=">"))
}
trinucsubsind = setNames(1:192, trinucsubs)

ind = setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))
gr_genes_ind = ind[gr_genes$names]

# Warning about possible unannotated dinucleotide substitutions
if (any(diff(mutations$pos)==1)) {
    warning("Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.")
}

# Warning about multiple instances of the same mutation in different sampleIDs
if (nrow(unique(mutations[,2:5])) < nrow(mutations)) {
    warning("Same mutations observed in different sampleIDs. Please verify that these are independent events and remove duplicates otherwise.")
}

# Start and end position of each mutation
mutations[, c("start", "end") := .(pos, pos)]
mutations[, end := end + nchar(ref) - 1]
# Correct 'start' for deletions based on comparison of the first character of 'ref' and 'mut'
mutations[, start := 
        start + (substr(ref, 1, 1) == substr(mut, 1, 1) 
            & nchar(ref) > nchar(mut))]

# Mapping mutations to genes
gr_muts = GenomicRanges::GRanges(mutations[,chr], IRanges::IRanges(mutations[,start],mutations[,end]))
ol = as.data.table(GenomicRanges::findOverlaps(gr_muts, gr_genes, type="any", select="all"))
mutations = mutations[ol[,queryHits]] # Duplicating subs if they hit more than one gene
mutations$geneind = gr_genes_ind[ol[,subjectHits]]
mutations$gene = sapply(RefCDS,function(x) x$gene_name)[mutations$geneind]
mutations = unique(mutations)

# Optional: Excluding samples exceeding the limit of mutations/sample [see Default parameters]
# nsampl = sort(table(mutations$sampleID))

# if (any(nsampl>max_coding_muts_per_sample)) {
#     message(sprintf('    Note: %0.0f samples excluded for exceeding the limit of mutations per sample (see the max_coding_muts_per_sample argument in dndscv). %0.0f samples left after filtering.',sum(nsampl>max_coding_muts_per_sample),sum(nsampl<=max_coding_muts_per_sample)))
#     exclsamples = names(nsampl[nsampl>max_coding_muts_per_sample])
#     mutations = mutations[!(mutations$sampleID %in% names(nsampl[nsampl>max_coding_muts_per_sample])),]
# }

# # Optional: Limiting the number of mutations per gene per sample (to minimise the impact of unannotated kataegis and other mutation clusters) [see Default parameters]
# mutrank = ave(mutations$pos, paste(mutations$sampleID,mutations$gene), FUN = function(x) rank(x))
# if (any(mutrank>max_muts_per_gene_per_sample)) {
#     message(sprintf('    Note: %0.0f mutations removed for exceeding the limit of mutations per gene per sample (see the max_muts_per_gene_per_sample argument in dndscv)',sum(mutrank>max_muts_per_gene_per_sample)))
#     exclmuts = mutations[mutrank>max_muts_per_gene_per_sample,]
#     mutations = mutations[mutrank<=max_muts_per_gene_per_sample,]
# }

# Additional annotation of substitutions
mutations$strand = sapply(RefCDS,function(x) x$strand)[mutations$geneind]
snv = (mutations$ref %in% nt & mutations$mut %in% nt)
if (!any(snv)) { stop("Zero coding substitutions found in this dataset. Unable to run dndscv. Common causes for this error are inputting only indels or using chromosome names different to those in the reference database (e.g. chr1 vs 1)") }
indels = mutations[!snv,]
mutations = mutations[snv,]
mutations$ref_cod = mutations$ref
mutations$mut_cod = mutations$mut
compnt = setNames(rev(nt), nt)
isminus = (mutations$strand==-1)
mutations$ref_cod[isminus] = compnt[mutations$ref[isminus]]
mutations$mut_cod[isminus] = compnt[mutations$mut[isminus]]

for (j in 1:length(RefCDS)) {
    RefCDS[[j]]$N = array(0, dim=c(192,4)) # Initialising the N matrices
}


# Annotating the functional impact of each substitution and populating the N matrices
mutations = annot_variant_impact_dt(mutations)
mutations$pid = sapply(RefCDS,function(x) x$protein_id)[mutations$geneind]
mutations[, trisub := trinucsubsind[paste(ref3_cod, mut3_cod, sep = ">")]]
# Filter to correct annotations only

if (any(mutations[,wrong_ref])) {
    if (mutations[,mean(wrong_ref)] < 0.1) { # If fewer than 10% of mutations have a wrong reference base, we warn the user
        warning(sprintf('%0.0f (%0.2g%%) mutations have a wrong reference base (see the affected mutations in dndsout$wrongmuts). Please identify the causes and rerun dNdScv.', mutations[,sum(wrong_ref)], 100*mutations[,mean(wrong_ref)]))
    } else { # If more than 10% of mutations have a wrong reference base, we stop the execution (likely wrong assembly or a serious problem with the data)
        stop(sprintf('%0.0f (%0.2g%%) mutations have a wrong reference base. Please confirm that you are not running data from a different assembly or species.',  mutations[,sum(wrong_ref)], 100*mutations[,mean(wrong_ref)]))
    }
    wrong_refbase = mutations[(wrong_ref), c(1:5)]
}
mutations <- mutations[!(wrong_ref) & !is.na(impind), !"wrong_ref" ]

# This is where we make another, minor, improvement in speed/resource requirements.
# instead of appending the N matrix of RefCDS mutation by mutation we first aggregate the counts and then append.
# Aggregate to count mutations by geneind, trisub, and impind
aggregated <- mutations[, .(count = .N), by = .(geneind, trisub, impind)]

# Apply updates in bulk
invisible(
    lapply(seq_len(nrow(aggregated)), function(i) {
    row <- aggregated[i]
    # <<- required as lapply won't update anything outside of the scope of the function (like global variables)
    if (
        (!is.na(row$trisub)) 
        && (row$trisub > 0) 
        && (row$geneind <= length(RefCDS))) {
        RefCDS[[row$geneind]]$N[row$trisub, row$impind] <<- RefCDS[[row$geneind]]$N[row$trisub, row$impind] + row$count
    }
    if (round(i/1e4)==(i/1e4)) { message(sprintf('    %0.3g%% ...', round(i/nrow(mutations),2)*100)) }
    })
)
# backup_indels = indels

# indels

if (any(nrow(indels))) {
    indels = annot_indels(indels)
    annot = rbind(mutations[,-c('trisub','impind')], indels)
} else {
    annot = mutations[,-c('trisub','impind')]
}
annot = annot[order(sampleID, chr,pos)]



end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken 

ncol(indels)
ncol(mutations)
