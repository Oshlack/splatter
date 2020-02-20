#' Example GFF file
#'
#' An example GFF file containing human genes from chromosome 22 of the GRCh38
#'
#' @docType data
#' 
#' @usage data(ex_gff)
#'
#' @format A data frame with 504 rows and 9 variables:
#' \describe{
#'   \item{V1}{id}
#'   \item{V2}{source}
#'   \item{V3}{feature}
#'   \item{V4}{start}
#'   \item{V5}{end}
#'   \item{V6}{score}
#'   \item{V7}{direction}
#'   \item{V8}{frame}
#'   \item{V9}{attribute}
#'   ...
#' }
#' @source \url{ftp://ftp.ensembl.org/pub/release-98/gff3/homo_sapiens/Homo_sapiens.GRCh38.98.chromosome.22.gff3.gz
#'}
"ex_gff"


#' Example simulated genotype file (vcf)
#'
#' An example genotype file in the .vcf format, where each column is a sample
#' and each row is a SNP.
#' 
#' @docType data
#' 
#' @usage data(ex_snps)
#' 
#' @format A data frame with 141882 rows and 109 variables:
#' 
#' @source \url{ftp://ftp.ensembl.org/pub/release-98/gff3/homo_sapiens/Homo_sapiens.GRCh38.98.chromosome.22.gff3.gz
#'}
"ex_snps"


#' Example gene mean data for a whole population (from GTEx)
#'
#' An example population wide gene mean data.frame, where each column is a 
#' sample and each row is a gene. Data is from the Thyroid tissue from the GTEx
#' dataset. Example data has already been filtered to remove genes with >50% 
#' of samples have < 0.1 TPM
#' 
#' @docType data
#' 
#' @usage data(ex_means)
#' 
#' @format A data frame with 2438 rows and 850 variables:
#' 
#' @source \url{https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
#'}
"ex_means"


#' Example eQTL results from GTEx Thyroid tissue
#'
#' An example dataframe with eQTL pairs and effect sizes from the Thyroid 
#' tissue from the GTEx dataset. This example data has already been filtered to 
#' include only the top eSNP for each gene. Columns that must be included are: 
#' 'gene_id', 'pval_nominal', and 'slope'.
#' 
#' @docType data
#' 
#' @usage data(ex_pairs)
#' 
#' @format A data frame with 25,244 rows and 9 variables:
#' 
#' @source \url{https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/all_snp_gene_associations/Thyroid.allpairs.txt.gz
#'}
"ex_pairs"
