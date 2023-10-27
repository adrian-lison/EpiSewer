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
