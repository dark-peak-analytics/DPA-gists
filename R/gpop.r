# General population mortality and utility adjustment

#' General population mortality using a distributed method
#'
#' Computes an overall survival curve for every baseline age from 0 to 99 from
#' baseline to the time horizon in years, then uses this matrix to return
#' overall survival for males and females on the cohort level.
#'
#' @param lt data.frame. Life table including qx (probability of
#' surviving to age +1 years, given age) for each age, with columns "m"
#' and "f".
#' @param age_mean numeric of length 2. Mean age of males and females,
#' respectively.
#' @param age_sd numeric of length 2. Standard deviation of age for
#' males and females, respectively.
#' @param age_dist character of length 2. Distributional shape for
#' males and females ("normal" or "lnorm").
#' @param cycle_length numeric. Length of cycle in years.
#' @param time_horizon numeric. Length of time horizon in years.
#' @param death_at_100 logical. If TRUE, assumes people die upon
#' turning 100 years old.
#' @param out character. Output type: "os", "matrix", or "unweighted
#' matrix".
#'
#' @return Depending on `out`, returns cohort-level overall survival,
#' weighted/unweighted survival matrices.
#' @importFrom matrixStats rowCumprods
#' @examples
#' # mort_gpop_dist(lt, age_mean, age_sd, age_dist, cycle_length, time_horizon)
mort_gpop_dist <- function(
    lt,
    age_mean,
    age_sd,
    age_dist,
    cycle_length,
    time_horizon,
    death_at_100 = TRUE,
    out = "unweighted matrix") {
  # Time steps and units, then conversions:

  requireNamespace("matrixStats", quietly = TRUE)
  stopifnot(out %in% c("os", "matrix", "unweighted matrix"))
  assertthat::assert_that(
    all(
      all(names(lt) %in% c("m", "f")),
      all(names(age_mean) %in% c("m", "f")),
      all(names(age_sd) %in% c("m", "f")),
      all(names(age_dist) %in% c("m", "f"))
    ),
    msg = "lt, age_mean, age_sd, age_dist must have names 'm' and 'f'"
  )

  # There are a couple of constants involved in these computations:
  age_bl_ub <- 99
  age_ub <- age_bl_ub + 1
  age_range <- seq(0, age_ub, 1)
  age_range_no_ub <- age_range[seq_len(age_ub)]

  # making sure that the names are all ok for the inputs:
  # cycles rounds up and adds one for cycle 0
  n_cycles <- floor(time_horizon / cycle_length) + 1
  t <- seq(0, time_horizon, by = cycle_length)[seq_len(n_cycles)]
  stopifnot(t[2] == cycle_length)
  t_floor <- floor(t)

  # use the above to make a coordinate matrix for life tables in the same
  # structure as the full probability and overall survival matrix:
  age_lookup <- matrix(
    row(matrix(0, nrow = age_ub, ncol = n_cycles)),
    nrow = age_ub,
    ncol = n_cycles
  )
  age_lookup <- outer(seq_len(age_ub), t_floor, "+")
  age_lookup <- pmin(age_lookup, age_ub + 1)

  # change the time unit of the gpop mortality to match cycle_length
  lt <- 1 - exp(-1 * (-log(1 - lt[, c("m", "f")])) * cycle_length)

  # Compute a 2d overall survival matrix for each baseline age for each sex.
  # Compute the CDF at baseline of ages, then use the age reference table to
  # look up the correct mortality values from the life tables in a vectorised
  # manner
  sex_index <- setNames(c("m", "f"), nm = c("m", "f"))
  os <- lapply(sex_index, function(sex) {
    bl_dist <- switch(age_dist[sex],
      normal = dnorm(
        seq(0, 99, 1),
        mean = age_mean[sex],
        sd = age_sd[sex]
      ),
      lnorm = dlnorm(
        seq(0, 99, 1),
        meanlog = log(age_mean[sex]),
        sdlog = age_sd[sex]
      )
    )
    # equivalent of pnorm for bandwidths, but with sum == 1 is TRUE:
    bl_dist[length(bl_dist)] <- 1 - sum(bl_dist[seq_len(length(bl_dist) - 1)])

    # Use the matrix of positional indices to look up probs in the life table
    # Note that depending on death_at_100, last row is p(death)=1 or
    # p(death | 100) for all ages beyond 100th birthday
    p_matrix <- matrix(
      data = lt[age_lookup, sex],
      nrow = dim(age_lookup)[1],
      ncol = dim(age_lookup)[2]
    )
    p_matrix[age_lookup == age_ub + 1] <- 1
    s_matrix <- 1 - p_matrix
    # return the un-weighted survival matrix
    list(
      bl = bl_dist,
      os = matrixStats::rowCumprods(cbind(1, s_matrix))
    )
  })

  # depending on what the user wants to return, conduct further processing:
  if (out %in% c("os", "matrix")) {
    # OS lines weighted by baseline distribution. For utilities and multiple
    # actions, probably better to return the unweighted ones and do as a next
    # step
    weighted_os <- lapply(
      sex_index,
      function(sex) {
        diag(os[[sex]]$bl) %*% os[[sex]]$os
      }
    )
    if (out == "matrix") {
      return(weighted_os)
    }
  }

  # if simple OS lines are wanted, do that
  if (out == "os") {
    lapply(weighted_os, colSums)
  }

  # only other alternative is to produce the un-weighted survival matrix where
  # all rows start from 1
  list(
    os = os,
    index = age_lookup
  )
}

#' Uses distributed overall survival and baseline distribution of age per sex
#' to estimate health related quality of life for all still-living general
#' population sub-cohort members
qol_gpop_dist <- function(
    os_list,
    util,
    pr_fe,
    output = "simple") {
  stopifnot(output %in% c("simple", "matrix", "all"))
  stopifnot(all(names(os_list) %in% c("os", "index")))
  # contains os and index
  list2env(os_list, envir = environment())

  # The weighting to apply to utility for each baseline age cohort is the
  # proportion of their proportion of the currently surviving cohort
  rep_index <- setNames(names(os), nm = names(os))
  u_list <- lapply(rep_index, function(sex_lab) {
    sex <- os[[sex_lab]]
    # the object sex contains the un-scaled overall survival extrapolations for
    # each baseline age, as well as the baseline densities summing to 1, and a
    # reference matrix telling us which row from the utility table to apply
    #
    # With this, we can multiply the un-weighted overall survival matrix by
    # baseline weightings row-wise, and then by the column sums of the result
    # column-wise to get the weighted utility of each cohort at each time,
    # scaled such that OS sums to 1 in each column of the matrix
    surv <- sex$os
    bl <- sex$bl

    # weighted survival. `w_remsurv` has columns that sum to 1
    w_surv <- bl * surv
    w_remsurv <- w_surv %*% diag(1 / colSums(w_surv))

    # now use index to generate a utility matrix and multiply that element-wise
    # by w_remsurv
    u <- util[, sex_lab]
    u_mat <- base::matrix(
      data = u[index],
      nrow = dim(index)[1],
      ncol = dim(index)[2],
      byrow = FALSE,
      dimnames = NULL
    )

    # output matrix of utilities after applying baseline weightings and taking
    # into account the variable death rates given baseline age over time
    u_mat * w_remsurv[1:dim(u_mat)[1], 1:dim(u_mat)[2]]
  })

  switch(output,
    simple = lapply(u_list, colSums),
    matrix = u_list,
    all = list(
      u = u_list,
      os = lapply(os, function(x) x$os),
      bl = lapply(os, function(x) x$bl),
      u_cohort = lapply(u_list, colSums)
    )
  )
}

# testing area
if (FALSE) {
  library(openxlsx)
  wb <- openxlsx::loadWorkbook(
    file = "https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/datasets/nationallifetablesunitedkingdomreferencetables/current/nltuk198020213.xlsx"
  )
  sheets <- openxlsx::sheets(wb)

  # extract from the "2019-2021" sheet. Males are on the left:
  lt <- openxlsx::read.xlsx(
    wb,
    sheet = "2019-2021",
    startRow = 6,
    colNames = TRUE
  )
  lt <- lt[, c(1, which(colnames(lt) == "qx"))]
  colnames(lt) <- c("age", "m", "f")
  lt <- as.matrix(lt)

  # Set some other dummy parameters to test the function
  bl_pr_m <- 0.43
  age_mean <- setNames(c(40, 40), c("m", "f"))
  age_sd <- setNames(c(10, 10), nm = c("m", "f"))
  age_dist <- setNames(c("normal", "normal"), nm = c("m", "f"))
  cycle_length <- 28 / 365.25
  time_horizon <- 50
  death_at_100 <- TRUE
  return_matrices <- FALSE

  # Now run the function. This provides a matrix of overall survival lines
  # where each row is a baseline age, and a weighting to assign to each row
  gpop_results <- mort_gpop_dist(
    lt = lt,
    age_mean = age_mean,
    age_sd = age_sd,
    age_dist = age_dist,
    cycle_length = cycle_length,
    time_horizon = time_horizon,
    death_at_100 = death_at_100,
    out = "unweighted matrix"
  )

  # Manually, you can compute the overall survival like this:
  os_m <- colSums(diag(gpop_results$os$m$bl) %*% gpop_results$os$m$os)
  os_f <- colSums(diag(gpop_results$os$f$bl) %*% gpop_results$os$f$os)

  # Everything else required to compute general population overall survival goes
  # from there:
  pr_m <- (os_m * bl_pr_m) / ((os_m * bl_pr_m) + (os_f * (1 - bl_pr_m)))
  l_os_m <- data.table::shift(os_m, type = "lag", fill = 1)
  h_m <- 1 - (os_m / l_os_m)
  l_os_f <- data.table::shift(os_f, type = "lag", fill = 1)
  h_f <- 1 - (os_f / l_os_f)
  h_w <- (h_m * pr_m) + (h_f * (1 - pr_m))
  os_w <- cumprod(c(1, 1 - h_w[2:length(os_m)]))

  plot_dat <- data.table::as.data.table(list(
    t = cycle_length * seq_len(length(os_m)),
    os_m = os_m,
    os_f = os_f,
    os_w = os_w,
    pr_m = pr_m,
    h_m = h_m,
    h_f = h_f,
    h_w = h_w
  ))

  plot_dat <- data.table::melt(
    plot_dat,
    id.vars = "t",
    variable.name = "variable"
  )

  # OS plot:
  ggplot2::ggplot(
    data = plot_dat[grep("os", plot_dat$variable), ],
    mapping = ggplot2::aes(x = t, y = value, color = variable)
  ) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Time (years)",
      y = "Value",
      color = "Variable"
    ) +
    ggplot2::theme_minimal()

  # hazard plot: note the "birthday effect" because each time age increments it
  # affects each age (which has a different weighting!) differently in a
  # non-linear fashion. This means that those people that e.g. reach 100 years
  # old get a hazard of 1 and die immediately, meaning that the hazard jumps up
  # and back down again each birthday! Furthermore, as the mortality effect of
  # getting older is different for each age, mortality will fall throughout the
  # year as older people die faster than younger.
  ggplot2::ggplot(
    data = plot_dat[grep("h_", plot_dat$variable), ],
    mapping = ggplot2::aes(x = t, y = value, color = variable)
  ) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Time (years)",
      y = "Value",
      color = "Variable"
    ) +
    ggplot2::theme_minimal()

  # There is a separate function for computing the male OS, female OS and
  # weighted "hazard" (or estimated conditional transition probability) for
  # useful results

  # Note that in order to work out weighted average utility over time, you need
  # to use the unweighted matrix, so the first function provides the pieces you
  # need

  # if unweighted matrix, you have the most flexibility, so this is the
  # default
  plot(colSums(diag(gpop_results$os$m$bl) %*% gpop_results$os$m$os), type = "l")
  plot(gpop_results$os$m$os[50, ])

  # Something more fun:
  os_m_t <- t(gpop_results$os$m$os)
  os_m_t <- data.table::as.data.table(os_m_t)
  os_m_t$t <- cycle_length * seq_len(nrow(os_m_t))
  dt_long <- data.table::melt(
    os_m_t,
    id.vars = "t",
    variable.name = "var",
    value.name = "value"
  )
  dt_long[, baseline_age := as.integer(sub("V", "", var)) - 1]
  dt_long <- dt_long[, .(baseline_age, t, value)]
  bl_weighting <- gpop_results$os$m$bl / sum(gpop_results$os$m$bl)
  # bl_weighting[bl_weighting <= 0.95] <- 0
  # bl_weighting[bl_weighting >= 0.05] <- 0

  bl_merge <- data.table::data.table(
    baseline_age = seq(0, 99, 1),
    bl_pr = bl_weighting
  )
  dt_long <- data.table::merge.data.table(
    x = dt_long,
    y = bl_merge,
    by = "baseline_age",
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
  )
  dt_long[, dens_prop := (value / sum(value)) * bl_pr, by = c("t")]
  dt_long[, dens_rel := dens_prop / sum(dens_prop), by = c("t")]

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = dt_long,
      ggplot2::aes(
        x = t,
        y = value,
        group = baseline_age,
        colour = dens_rel,
        alpha = dens_rel
      )
    ) +
    ggplot2::geom_line(
      data = plot_dat[grep("os_m", plot_dat$variable), ],
      mapping = ggplot2::aes(x = t, y = value),
      colour = "red"
    ) +
    ggplot2::theme_classic() +
    ggplot2::scale_alpha(range = c(0.001, 1)) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.025))
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::theme(legend.position = "none")


  util <- read.csv("gpop2/util.csv")
  qol <- qol_gpop_dist(
    os_list = gpop_results,
    util = util,
    pr_fe = 1 - bl_pr_m,
    output = "all"
  )

  u_plot <- data.table::data.table(
    t = cycle_length * seq_len(length(qol$u_cohort$m)),
    data.table::as.data.table(qol$u_cohort)
  )
  u_plot <- data.table::melt.data.table(u_plot, id.vars = "t")

  p_abs <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = u_plot,
      mapping = ggplot2::aes(
        x = t,
        y = value,
        colour = variable
      )
    ) +
    ggplot2::geom_line(
      data = plot_dat[grep("os", plot_dat$variable), ],
      mapping = ggplot2::aes(x = t, y = value, color = variable)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.025))
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, max(u_plot$t)),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Absolute utility over time")

  # Relative utility plot. Note the sawtooth wave pattern. this is because
  # during each year between birthdays people at different ages are dieing at
  # non-linearly different mortality rates. This means that after a spike in
  # mortality rate among those still alive (with older people having more of a
  # jump in mortality than younger people), the attrition pushes utility up
  # again during the year (because the survivors are on average younger over
  # time within a given year).
  u_plot[, value_rel := value / value[which.min(t)], by = variable]

  p_rel <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = u_plot,
      mapping = ggplot2::aes(
        x = t,
        y = value_rel,
        colour = variable
      )
    ) +
    ggplot2::geom_line(
      data = plot_dat[grep("os", plot_dat$variable), ],
      mapping = ggplot2::aes(x = t, y = value, color = variable)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.025))
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, max(u_plot$t)),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Relative utility over time")

  (p_abs | p_rel) +
    patchwork::plot_layout(guides = "collect")
}
