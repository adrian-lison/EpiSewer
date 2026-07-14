#' Raw SARS-CoV-2 wastewater data from Zurich
#'
#' @description Wastewater measurements of SARS-CoV-2 RNA and case counts from
#'   the WWTP Werdhölzli in Zürich, Switzerland (January–April 2022). Data are
#'   provided by EAWAG under the CC BY 4.0 license.
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{measurements}{A `data.table` with columns `date` and
#'     `concentration` (gc/mL).}
#'   \item{flows}{A `data.table` with columns `date` and `flow` (mL/day).}
#'   \item{cases}{A `data.table` with columns `date` and `cases`
#'     (cases/day).}
#'   \item{units}{A list with concentration and flow unit strings.}
#' }
#' @source \url{https://sensors-eawag.ch/sars/}
"SARS_CoV_2_Zurich"

#' Raw Influenza A wastewater data from Zurich
#'
#' @description Wastewater measurements of Influenza A RNA from the WWTP
#'   Werdhölzli in Zürich, Switzerland (September 2022–April 2023). Data are
#'   provided by EAWAG.
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{measurements}{A `data.table` with columns `date`,
#'     `concentration` (gc/mL), `replicate_id`, `dilution`,
#'     `total_droplets`, and `positive_droplets`.}
#'   \item{flows}{A `data.table` with columns `date` and `flow` (mL/day).}
#'   \item{units}{A list with concentration and flow unit strings.}
#' }
"Influenza_A_Zurich"

#' Example sewer_data object for SARS-CoV-2 in Zurich
#'
#' @description Processed wastewater data for SARS-CoV-2 from the WWTP
#'   Werdhölzli in Zürich (January–April 2022), subsampled to Monday and
#'   Thursday measurements, ready for use with [EpiSewer()].
#'
#' @format A `sewer_data` object as returned by [sewer_data()].
#' @seealso [SARS_CoV_2_Zurich], [ww_assumptions_SARS_CoV_2_Zurich],
#'   [ww_result_SARS_CoV_2_Zurich]
"ww_data_SARS_CoV_2_Zurich"

#' Example sewer_data object for Influenza A in Zurich
#'
#' @description Processed wastewater data for Influenza A from the WWTP
#'   Werdhölzli in Zürich (September–December 2022), ready for use with
#'   [EpiSewer()].
#'
#' @format A `sewer_data` object as returned by [sewer_data()].
#' @seealso [Influenza_A_Zurich], [ww_assumptions_influenza_Zurich]
"ww_data_influenza_Zurich"

#' Example sewer_assumptions object for SARS-CoV-2 in Zurich
#'
#' @description Epidemiological assumptions for SARS-CoV-2, including
#'   generation time, incubation period, and shedding load distributions
#'   derived from the literature.
#'
#' @format A `sewer_assumptions` object as returned by [sewer_assumptions()].
#' @seealso [ww_data_SARS_CoV_2_Zurich], [ww_result_SARS_CoV_2_Zurich]
"ww_assumptions_SARS_CoV_2_Zurich"

#' Example sewer_assumptions object for Influenza A in Zurich
#'
#' @description Epidemiological assumptions for Influenza A, including
#'   generation time and shedding load distributions derived from the
#'   literature.
#'
#' @format A `sewer_assumptions` object as returned by [sewer_assumptions()].
#' @seealso [ww_data_influenza_Zurich]
"ww_assumptions_influenza_Zurich"

#' Example EpiSewerJobResult for SARS-CoV-2 in Zurich
#'
#' @description Result of running [EpiSewer()] on [ww_data_SARS_CoV_2_Zurich]
#'   with [ww_assumptions_SARS_CoV_2_Zurich]. Can be used to explore outputs
#'   and plotting functions without fitting a model.
#'
#' @format An `EpiSewerJobResult` object as returned by [EpiSewer()].
#' @seealso [ww_data_SARS_CoV_2_Zurich], [ww_assumptions_SARS_CoV_2_Zurich]
"ww_result_SARS_CoV_2_Zurich"
