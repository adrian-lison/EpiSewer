source("local/local_data.R")

WW_data <- read_local_data_EAWAG()

WW_data <- WW_data[wwtp_name == "ARA Werdhoelzli" & target == "IAV-M", ]
WW_data[,concentration := gc_per_lww / 1000]
WW_data[,flow := flow * 1000 * 1000]
WW_data <- WW_data[!is.na(concentration), c("date", "concentration", "flow")]

WW_data <- WW_data[
  date >= as.Date("2022-09-15") & date < as.Date("2023-05-01"),
]

# interpolate flow on missing days
WW_data <- WW_data[
  data.table::CJ(date = seq.Date(min(date), max(date), by = "day"), unique = TRUE),
  on = .(date)
]
data.table::setnafill(WW_data, type = "locf", cols = "flow")

Influenza_A_Zurich <- list(
  measurements = WW_data[, c("date", "concentration")],
  flows = WW_data[, c("date", "flow")],
  units = list(concentration = "gc/mL", flow = "mL/day")
)

usethis::use_data(Influenza_A_Zurich, overwrite = TRUE)

