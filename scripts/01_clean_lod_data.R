library(tidyverse)
library(here)
library(readr)

# Read in raw LOD data
lod_data <- read_csv(file = here(
  "data", "raw", "lod", "nist.lod.triplicate.averages.method.csv"
  )
  )

# Convert date to date
lod_data$analysis_date <- mdy(lod_data$analysis_date)

# Convert to long form and calculate thiamine conc from fluorescence
lod_data_long <- lod_data |>
  mutate(sample_rep = paste0(sample, "_", rep)) |> 
  pivot_longer(cols = starts_with("t"),
               names_to = "time_rep",
               values_to = "flourescence") |> 
  separate_wider_delim(time_rep, names = c("time_point", "analysis_rep"), delim = 
                         ".") |> 
  mutate(thiamine_nmol = (flourescence - run_cal_intercept) / run_cal_slope)

# Average duplicates, calculated standard deviations and RSD
summarized_lod_data_long <- lod_data_long |> 
  group_by(sample_rep, time_point, analysis_date) |> 
  summarise(
    avg_thiamine = mean(thiamine_nmol),
    sd_thiamine = sd(thiamine_nmol),
    rsd = (sd_thiamine / avg_thiamine) * 100,
    .groups = "drop")

# Convert back to wide format and calculate absolute differences
summarized_lod_data_wide <- summarized_lod_data_long |> 
  pivot_wider(names_from = time_point, 
              values_from = c(avg_thiamine, sd_thiamine, rsd)) |> 
  mutate(abs_thiamine_difference = abs(avg_thiamine_t0 - avg_thiamine_t1))

# Save cleaned datasheet for analysis
write_csv(
  summarized_lod_data_wide,
  file = here("data", "processed", "lod", "clean_lod_data.csv")
  )

# Save cleaned long datasheet for plotting
write_csv(
  summarized_lod_data_long,
  file = here("data", "processed", "lod", "clean_lod_data_long.csv")
)
  