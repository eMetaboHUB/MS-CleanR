#' mscleanr: A package for cleaning and analyzing MS data
#'
#' The mscleanr package provides 2 important functions: \code{\link{clean_msdial_data}}, \code{\link{keep_top_peaks}} and \code{\link{launch_msfinder_annotation}}.
#' See the functions documentation and vignettes for more information.
#' Needs MSDial v4.00 or higher.
#'
#' @docType package
#' @name mscleanr
#'
#' @section Architecture needed for an analysis:
#' See functions documentation.
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
NULL



#' Mass differences of adducts with the original compound M, in negative mode.
#'
#' @format A data frame with 29 rows and 2 variables:
#' \describe{
#'   \item{adduct}{name of the adduct}
#'   \item{diff}{mass difference with M in negative mode}
#' }
#'
#' @seealso \code{\link{mass_adducts_pos}}, \code{\link{mass_neutral_loss_neg}}, \code{\link{mass_neutral_loss_pos}}, \code{\link{mass_isotopes}}
"mass_adducts_neg"



#' Mass differences of adducts with the original compound M, in positive mode.
#'
#' @format A data frame with 46 rows and 2 variables:
#' \describe{
#'   \item{adduct}{name of the adduct}
#'   \item{diff}{mass difference with M in positive mode}
#' }
#' @seealso \code{\link{mass_adducts_neg}}, \code{\link{mass_neutral_loss_neg}}, \code{\link{mass_neutral_loss_pos}}, \code{\link{mass_isotopes}}
"mass_adducts_pos"



#' Masses of isotopes.
#'
#' @format A data frame with 2 rows and 2 variables:
#' \describe{
#'   \item{name}{name of the isotope}
#'   \item{m}{mass}
#' }
#'
#' @seealso \code{\link{mass_adducts_neg}}, \code{\link{mass_adducts_pos}}, \code{\link{mass_neutral_loss_neg}}, \code{\link{mass_neutral_loss_pos}}
"mass_isotopes"



#' Masses of neutral losses, in negative mode.
#'
#' @format A data frame with 2 rows and 2 variables:
#' \describe{
#'   \item{name}{name of the neutral loss}
#'   \item{m}{mass in negative mode}
#' }
#'
#' @seealso \code{\link{mass_adducts_neg}}, \code{\link{mass_adducts_pos}}, \code{\link{mass_neutral_loss_pos}}, \code{\link{mass_isotopes}}
"mass_neutral_loss_neg"



#' Masses of neutral losses, in positive mode.
#'
#' @format A data frame with 2 rows and 2 variables:
#' \describe{
#'   \item{name}{name of the neutral loss}
#'   \item{m}{mass in positive mode}
#' }
#'
#' @seealso \code{\link{mass_adducts_neg}}, \code{\link{mass_adducts_pos}}, \code{\link{mass_neutral_loss_neg}}, \code{\link{mass_isotopes}}
"mass_neutral_loss_pos"
