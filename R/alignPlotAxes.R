# Alignment solution from:
# https://stackoverflow.com/questions/54153906/aligning-axes-of-r-plots-on-one-side-of-a-grid-together

# Load lung survival dataset from survival package
lung_dataset <- survival::lung |>
  # we only need the time and status (event) columns
  dplyr::select(
    time,
    status
  ) |>
  # recode the status column to 0 (censored) and 1 (event)
  dplyr::mutate(
    staus = dplyr::recode(status, `1` = 0, `2` = 1)
  )

# Fit the survival object
fit <- survival::survfit(formula = survival::Surv(time = time, event = staus) ~ 1,
                         data = lung_dataset)

# create simple KM plot and risk table
plot_and_table <- survminer::ggsurvplot(
  fit = fit,
  data = lung_dataset,
  risk.table = TRUE,
)

# Use the patchwork library for axis alignment.
# patchwork::plot_layout aligns the plot and table 
# by x-axis ticks (0 on plot directly above 0 on table, etc.)

# load the library to ensure the `+` symbols are interpreted correctly.
library(patchwork)

aligned_plot_and_table <- 
  plot_and_table$plot + 
  plot_and_table$table + 
  patchwork::plot_layout(heights = c(8, 2), nrow = 2)

# compare the two plots:
plot_and_table
aligned_plot_and_table
