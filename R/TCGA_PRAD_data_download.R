# initial look at TCGA PRAD data
# 11/18/2026

suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(recount)))
suppressWarnings(suppressPackageStartupMessages(library(recount3)))
suppressWarnings(suppressPackageStartupMessages(library(IRanges)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicRanges)))
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicFeatures)))

# Setup recount environmental variables
recount3_cache_rm()
recount3_cache(tempdir())
options(
  "recount3_organism_human_project_homes_URL_http://duffel.rail.bio/recount3" = c(
    "data_sources/sra",
    "data_sources/gtex",
    "data_sources/tcga"
  )
)

proj_dir <- "TCGA_PRAD"
dir.create(proj_dir, showWarnings = T, recursive = T)
setwd(proj_dir)

# load in project from recount
prad_jxn <- recount3::create_rse_manual(
  project = "PRAD",
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v26",
  type = "jxn"
)

# load in recount provided exon GTF reference file
recount3_exon_gtf_file <- file_retrieve(locate_url_ann(
  type = "exon",
  organism = "human",
  annotation = metadata(prad_jxn)$annotation
))
recount3_exon_gtf <- unique(import(recount3_exon_gtf_file, format = "gtf"))
metadata(prad_jxn)$GTF_exon <- recount3_exon_gtf

dir.create("data/reference/", showWarnings = T, recursive = T)
fwrite(
  as.data.frame(metadata(prad_jxn)$GTF_exon),
  paste0("data/reference/Recount3_", metadata(prad_jxn)$annotation, "_GTF.txt"),
  sep = '\t'
)

# extract project metadata
prad_meta <- as.data.frame(colData(prad_jxn))

# subset dataset to gene of interest - AR
gtf_AR <- recount3_exon_gtf[which(recount3_exon_gtf$gene_name == "AR")]
gtf_AR_chr <- as.character(seqnames(gtf_AR)[1])
gtf_AR_str <- min(start(gtf_AR), na.rm = T)
gtf_AR_end <- max(end(gtf_AR), na.rm = T)
AR_GRanges <- GenomicRanges::GRanges(
  gtf_AR_chr,
  IRanges(gtf_AR_str, gtf_AR_end, names = "AR")
)
prad_jxn_AR_overlap <- IRanges::subsetByOverlaps(
  rowRanges(prad_jxn),
  AR_GRanges
)
prad_jxn_AR <- prad_jxn[names(prad_jxn_AR_overlap), ]
