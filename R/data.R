#' Designed primers for Point Mutagenesis
#'
#' This is the result of the vignette for Point Mutagenesis with BbsI
#'
#' @format Primer object which is the result of the calculation
#' 
#' @source \url{https://github.com/ipb-halle/GoldenMutagenesis/tree/master/vignettes}
#' @name Point_Mutagenesis_BbsI_result
#' @docType data
NULL

#' Designed primers for Multiple Site Directed Mutagenesis
#'
#' This is the result of the vignette for Multiple Site Directed Mutagenesis in Level 2 with BsaI
#'
#' @format Primer object which is the result of the calculation
#' 
#' @source \url{https://github.com/ipb-halle/GoldenMutagenesis/tree/master/vignettes}
#' @name MSD_BsaI_result_lv2
#' @docType data
NULL

#' Basic setup for Point Mutagenesis
#'
#' Those are the objects of the Point Mutagenesis with BbsI vignette.
#'
#' @format The following objects are included:
#'  \describe{
#'   \item{cuf}{The Codon Usage Table which is being used to select the codon for an exchanged amino acid.}
#'   \item{input_sequence}{The sequence which should be modified. This is just an object of type character containing the sequence. }
#'   \item{mutations}{The desired mutations.}
#'   \item{recognition_site_bbsi}{The BbsI recognition site.}
#'   
#' }
#' 
#' @source \url{https://github.com/ipb-halle/GoldenMutagenesis/tree/master/vignettes}
#' @name Point_Mutagenesis_BbsI_setup
#' @docType data
NULL

#' Basic setup for Multiple Site Directed Mutagenesis
#'
#' Those are the objects of the MSD vignette.
#'
#' @format The following objects are included:
#'  \describe{
#'   \item{cuf}{The Codon Usage Table which is being used to select the codon for an exchanged amino acid.}
#'   \item{input_sequence}{The sequence which should be modified. This is just an object of type character containing the sequence. }
#'   \item{mutations}{The desired mutations.}
#'   \item{recognition_site_bsai}{The BsaI recognition site.}
#'   
#' }
#' 
#' @source \url{https://github.com/ipb-halle/GoldenMutagenesis/tree/master/vignettes}
#' @name MSD_BsaI_setup_lv2
#' @docType data
NULL