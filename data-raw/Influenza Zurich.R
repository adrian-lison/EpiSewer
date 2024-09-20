source("local/local_data.R")
library(data.table)

WW_data <- read_local_data_EAWAG()

WW_data <- WW_data[wwtp_name == "ARA Werdhoelzli" & target == "IAV-M", ]
WW_data[,concentration := gc_per_lww / 1000]
WW_data[,flow := flow * 1000 * 1000]
WW_data <- WW_data[!is.na(concentration), c("date", "concentration", "replicate_id", "dilution", "total_droplets", "flow")]

WW_data <- WW_data[
  date >= as.Date("2022-09-15") & date < as.Date("2023-05-01"),
]

# interpolate flow on missing days
WW_data <- WW_data[
  data.table::CJ(date = seq.Date(min(date), max(date), by = "day"), unique = TRUE),
  on = .(date)
]
data.table::setnafill(WW_data, type = "locf", cols = "flow")

# interpolate total droplets and dilution
WW_data[is.na(dilution) & !is.na(concentration), dilution := median(WW_data$dilution, na.rm = T)]
WW_data[is.na(total_droplets) & !is.na(concentration), total_droplets := median(WW_data$total_droplets, na.rm = T)]

Influenza_A_Zurich <- list(
  measurements = WW_data[, c("date", "concentration", "replicate_id", "dilution", "total_droplets")],
  flows = WW_data[, c("date", "flow")],
  units = list(concentration = "gc/mL", flow = "mL/day")
)

usethis::use_data(Influenza_A_Zurich, overwrite = TRUE)

measurements <- copy(Influenza_A_Zurich$measurements)
measurements[,is_outlier := FALSE]
measurements[date %in% as.Date(c("2023-01-17", "2023-02-19")), is_outlier := TRUE]
measurements <- measurements[date <= as.Date("2022-12-19")]
ww_data_influenza_Zurich <- sewer_data(measurements = measurements, flows = Influenza_A_Zurich$flows)
usethis::use_data(ww_data_influenza_Zurich, overwrite = TRUE)

generation_dist <- get_discrete_gamma_shifted(gamma_mean = 2.6, gamma_sd = 1.5)
shedding_dist <- get_discrete_gamma(gamma_mean = 2.491217, gamma_sd = 1.004283)
ww_assumptions_influenza_Zurich <- sewer_assumptions(
  generation_dist = generation_dist,
  shedding_dist = shedding_dist,
  shedding_reference = "infection"
)
usethis::use_data(ww_assumptions_influenza_Zurich, overwrite = TRUE)

