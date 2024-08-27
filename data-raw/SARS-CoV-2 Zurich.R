# Load wastewater measurements for the SARS-CoV-2 N1 gene from the WWTP
# Werdhölzli in Zürich, Switzerland.
#
# Data are provided by EAWAG (Swiss Federal Institute of Aquatic Science and
# Technology) to the public domain (CC BY 4.0 license).
#
# Measurements are according to the protocol v4 (link url v2)
url <- "https://sensors-eawag.ch/sars/__data__/processed_normed_data_zurich_v2.csv"

WW_data <- data.table::setDT(readr::read_delim(url,
  delim = ";",
  skip = 1,
  col_names = c(
    "date", "load", "median_7d_sars_cov2_rna",
    "new_cases", "median_7d_new_cases", "quantification_flag", "flow"
  ),
  col_types = readr::cols(
    date = readr::col_date(), load = readr::col_double(),
    new_cases = readr::col_double()
  )
))

WW_data <- WW_data[!is.na(load), c("date", "load", "new_cases", "flow")]

# convert units
WW_data[, new_cases := new_cases * 471000/1e5] # from cases/100k pers/day to cases/day
WW_data[, load := load * 471000/1e5] # from gc/100k pers/day to gc/day
WW_data[, flow := flow * 1e+06] # from m3/day to mL/day


# interpolate flow on missing days
WW_data <- WW_data[
  data.table::CJ(date = seq.Date(min(date), max(date), by = "day"), unique = TRUE),
  on = .(date)
]
data.table::setnafill(WW_data, type = "locf", cols = "flow")
WW_data[, concentration := load / flow]
setnames(WW_data, "new_cases", "cases")

WW_data <- WW_data[
  date >= as.Date("2022-01-01") & date < as.Date("2022-05-01"),
]

SARS_CoV_2_Zurich <- list(
  measurements = WW_data[, c("date", "concentration")],
  flows = WW_data[, c("date", "flow")],
  cases = WW_data[, c("date", "cases")],
  units = list(concentration = "gc/mL", flow = "mL/day")
)

usethis::use_data(SARS_CoV_2_Zurich, overwrite = TRUE)

data_zurich <- SARS_CoV_2_Zurich
measurements_sparse <- data_zurich$measurements[,weekday := weekdays(data_zurich$measurements$date)][weekday %in% c("Monday","Thursday"),]
ww_data_SARS_CoV_2_Zurich <- sewer_data(measurements = measurements_sparse, flows = data_zurich$flows)
usethis::use_data(ww_data_SARS_CoV_2_Zurich, overwrite = TRUE)

generation_dist <- get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4, maxX = 12)
incubation_dist <- get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4, maxX = 10)
shedding_dist <- get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397, maxX = 30)
load_per_case <- 1e+11
ww_assumptions_SARS_CoV_2_Zurich <- sewer_assumptions(
  generation_dist = generation_dist,
  incubation_dist = incubation_dist,
  shedding_dist = shedding_dist,
  load_per_case = load_per_case
)
usethis::use_data(ww_assumptions_SARS_CoV_2_Zurich, overwrite = TRUE)
