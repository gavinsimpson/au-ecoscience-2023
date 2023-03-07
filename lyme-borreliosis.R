# Analysis of the Lyme borreliosis weekly incidence data from Norway

# Data are from:
# The emergence and shift in seasonality of Lyme borreliosis in Northern Europe
# Goren et al (2023) doi.org/10.1098/rspb.2022.2420

# Packages
pkgs <- c("tibble", "readr", "dplyr", "here", "tidyr", "ggplot2", "mgcv",
    "gratia")
vapply(pkgs, library, FUN.VALUE = logical(1L), character.only = TRUE,
    logical.return = TRUE)

f <- here("raw-data", "lyme-borreliosis.csv")

lyme <- read_csv(f,
    col_types = "-iiiiiiiiiii")

lyme <- lyme |>
    mutate(date = as.Date(paste(yrwk, "1"), format = "%Y%U%u"),
    fyear = factor(year),
    pop = pop / 100000) |>
    select(-yrwk)

lyme |>
    ggplot(aes(x = date, y = cases)) +
    geom_point()

lyme |>
    ggplot(aes(x = week, y = cases, group = fyear, colour = fyear)) +
    geom_line(alpha = 0.5)  +
    scale_color_viridis_c(option = "plasma")

m1 <- bam(cases ~ s(week, bs = "cc", k = 20) + s(year, k = 25) +
    offset(log(pop)),
    data = lyme, method = "fREML", knots = list(week = c(0.5, 52.5)),
    family = nb(), discrete = TRUE, nthreads = 4)

draw(m1)

appraise(m1, method = "simulate")

m2 <- update(m1, . ~ . + ti(week, year, bs = c("cc", "tp"), k = c(20, 20)))

summary(m2)

draw(m2)

appraise(m2)

m3 <- bam(cases ~ s(week, fyear, bs = "fs", xt = "cc", k = 20) +
    offset(log(pop)),
    data = lyme, method = "fREML", knots = list(week = c(0.5, 52.5)),
    family = nb(), discrete = TRUE, nthreads = 4)

m4 <- bam(cases ~ s(week, bs = "cc", k = 20) +
    s(year, k = 25) +
    s(week, fyear, bs = "fs", xt = "cc", k = 10) +
    offset(log(pop)),
    data = lyme, method = "fREML", knots = list(week = c(0.5, 52.5)),
    family = nb(), discrete = TRUE, nthreads = 4)

draw(m4, residuals = TRUE, rug = FALSE)
summary(m4)
appraise(m4, method = "simulate")

m5 <- bam(cases ~ s(week, bs = "cc", k = 20) +
    s(year, k = 25) +
    s(fyear, bs = "re") +
    ti(week, year, bs = c("cc", "tp"), k = c(20, 20)) +
    offset(log(pop)),
    data = lyme, method = "fREML", knots = list(week = c(0.5, 52.5)),
    family = nb(), discrete = TRUE, nthreads = 4)

draw(m4, residuals = FALSE, rug = FALSE)
summary(m4)
appraise(m4, method = "simulate")

draw(m5, residuals = FALSE, rug = FALSE)
summary(m5)
appraise(m5, method = "simulate")

AIC(m1, m2, m3, m4, m5)

draw(m2, residuals = FALSE, rug = FALSE)
summary(m2)
appraise(m2, method = "simulate")

# rootograms

# predict to find week of peak in each year
ds <- data_slice(m2, week = evenly(week, n = 100), year = unique(year),
    pop = log(100000))

fv <- fitted_values(m2, data = ds, scale = "response")

fv |>
    ggplot(aes(x = week, y = fitted, group = year, colour = year)) +
    geom_line(linewidth = 1) +
    scale_colour_viridis_c(option = "cividis") +
    labs(y = "Number of cases per 100,000 population", x = "Week",
        title = "Annual seasonal trend in cases of Lyme borreliosis",
        colour = NULL)

ds2 <- data_slice(m2, week = unique(week), year = unique(year),
    pop = log(100000))
ds2 <- ds2 |>
    add_column(row = seq_len(nrow(ds2)), .before = 1L)

fv2 <- fitted_values(m2, data = ds2, scale = "response")

fv2 |>
    group_by(year) |>
    slice_max(order_by = fitted) |>
    ggplot(aes(x = year, y = week)) +
    geom_point()

fs2 <- fitted_samples(m2, data = ds2, n = 10000, seed = 42, scale = "response")

fs2 <- fs2 |>
    left_join(ds2, by = join_by(row)) |>
    select(-pop)

peak_week <- fs2 |>
    group_by(year, draw) |>
    slice_max(order_by = fitted) |>
    group_by(year) |>
    summarise(peak = quantile(week, prob = 0.5),
        lower_ci = quantile(week, prob = 0.055),
        upper_ci = quantile(week, prob = 0.945))

peak_week |>
    ggplot(aes(x = year, y = peak)) +
    geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci)) +
    labs(y = "Week of peak", x = NULL,
        title = "Annual peak in cases of Lyme borreliosis",
        caption = "89% credible interval")

