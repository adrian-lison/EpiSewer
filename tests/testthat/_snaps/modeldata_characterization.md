# messages for config 'calibrate_cases_direct' are stable

    Code
      invisible(modeldata_fixture_job(battery$calibrate_cases_direct()))
    Message
      i Due to non-detects at the start of the measurement time series, Rt will only be estimated from 2022-09-30 onwards. Measurements and Rt before that date are still modeled, but using an extended seeding phase. Set 'seeding = seeding_estimate_rw(extend = FALSE)' to deactivate this behaviour.

# messages for config 'default_sars' are stable

    Code
      invisible(modeldata_fixture_job(battery$default_sars()))

# messages for config 'partitions_likelihood' are stable

    Code
      invisible(modeldata_fixture_job(battery$partitions_likelihood()))

