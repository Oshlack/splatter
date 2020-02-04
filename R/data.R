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