## ----setup, include=FALSE, cache=FALSE----------------------------------------
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(cache = TRUE, dev = 'svg', echo = FALSE, message = FALSE, warning = FALSE,
                      fig.height=6, fig.width = 1.777777*6)

library("curl")
library('here')
library('mgcv')
library('gratia')
library('ggplot2')
library('purrr')
library('mvnfast')
library("tibble")
library('gganimate')
library('tidyr')
library("knitr")
library("viridis")
library('readr')
library('dplyr')
library('patchwork')

## plot defaults
theme_set(theme_minimal(base_size = 16, base_family = 'Fira Sans'))

## constants
anim_width <- 1000
anim_height <- anim_width / 1.77777777
anim_dev <- 'png'
anim_res <- 200


## ----smooth-fun-animation, results = FALSE, echo = FALSE----------------------
f <- function(x) {
    x^11 * (10 * (1 - x))^6 + ((10 * (10 * x)^3) * (1 - x)^10)
}

draw_beta <- function(n, k, mu = 1, sigma = 1) {
    rmvn(n = n, mu = rep(mu, k), sigma = diag(rep(sigma, k)))
}

weight_basis <- function(bf, x, n = 1, k, ...) {
    beta <- draw_beta(n = n, k = k, ...)
    out <- sweep(bf, 2L, beta, '*')
    colnames(out) <- paste0('f', seq_along(beta))
    out <- as_tibble(out)
    out <- add_column(out, x = x)
    out <- pivot_longer(out, -x, names_to = 'bf', values_to = 'y')
    out
}

random_bases <- function(bf, x, draws = 10, k, ...) {
    out <- rerun(draws, weight_basis(bf, x = x, k = k, ...))
    out <- bind_rows(out)
    out <- add_column(out, draw = rep(seq_len(draws), each = length(x) * k),
                      .before = 1L)
    class(out) <- c("random_bases", class(out))
    out
}

plot.random_bases <- function(x, facet = FALSE) {
    plt <- ggplot(x, aes(x = x, y = y, colour = bf)) +
        geom_line(lwd = 1, alpha = 0.75) +
        guides(colour = "none")
    if (facet) {
        plt + facet_wrap(~ draw)
    }
    plt
}

normalize <- function(x) {
    rx <- range(x)
    z <- (x - rx[1]) / (rx[2] - rx[1])
    z
}

set.seed(1)
N <- 500
data <- tibble(x     = runif(N),
               ytrue = f(x),
               ycent = ytrue - mean(ytrue),
               yobs  = ycent + rnorm(N, sd = 0.5))

k <- 10
knots <- with(data, list(x = seq(min(x), max(x), length = k)))
sm <- smoothCon(s(x, k = k, bs = "cr"), data = data, knots = knots)[[1]]$X
colnames(sm) <- levs <- paste0("f", seq_len(k))
basis <- pivot_longer(cbind(sm, data), -(x:yobs), names_to = 'bf')
basis

set.seed(2)
bfuns <- random_bases(sm, data$x, draws = 20, k = k)

smooth <- bfuns %>%
    group_by(draw, x) %>%
    summarise(spline = sum(y)) %>%
    ungroup()

p1 <- ggplot(smooth) +
    geom_line(data = smooth, aes(x = x, y = spline), lwd = 1.5) +
    labs(y = 'f(x)', x = 'x') +
    theme_minimal(base_size = 16, base_family = 'Fira Sans')

smooth_funs <- animate(
    p1 + transition_states(draw, transition_length = 4, state_length = 2) + 
    ease_aes('cubic-in-out'),
    nframes = 200, height = anim_height, width = anim_width, res = anim_res,
    dev = anim_dev)

anim_save('resources/spline-anim.gif', smooth_funs)


## ----basis-functions, fig.height=6, fig.width = 1.777777*6, echo = FALSE------
ggplot(basis,
       aes(x = x, y = value, colour = bf)) +
    geom_line(lwd = 2, alpha = 0.5) +
    guides(colour = "none") +
    labs(x = 'x', y = 'b(x)') +
    theme_minimal(base_size = 20, base_family = 'Fira Sans')


## ----basis-function-animation, results = 'hide', echo = FALSE-----------------
bfun_plt <- plot(bfuns) +
    geom_line(data = smooth, aes(x = x, y = spline),
              inherit.aes = FALSE, lwd = 1.5) +
    labs(x = 'x', y = 'f(x)') +
    theme_minimal(base_size = 14, base_family = 'Fira Sans')

bfun_anim <- animate(
    bfun_plt + transition_states(draw, transition_length = 4, state_length = 2) + 
    ease_aes('cubic-in-out'),
    nframes = 200, height = anim_height, width = anim_width, res = anim_res, dev = anim_dev)

anim_save('resources/basis-fun-anim.gif', bfun_anim)


## ----basis-functions-anim, results = "hide", echo = FALSE---------------------
sm2 <- smoothCon(s(x, k = k, bs = "cr"), data = data, knots = knots)[[1]]$X
beta <- coef(lm(ycent ~ sm2 - 1, data = data))
wtbasis <- sweep(sm2, 2L, beta, FUN = "*")
colnames(wtbasis) <- colnames(sm2) <- paste0("F", seq_len(k))
## create stacked unweighted and weighted basis
basis <- as_tibble(rbind(sm2, wtbasis)) %>%
    add_column(x = rep(data$x, times = 2),
               type = rep(c('unweighted', 'weighted'), each = nrow(sm2)),
               .before = 1L)
##data <- cbind(data, fitted = rowSums(scbasis))
wtbasis <- as_tibble(rbind(sm2, wtbasis)) %>%
    add_column(x      = rep(data$x, times = 2),
               fitted = rowSums(.),
               type   = rep(c('unweighted', 'weighted'), each = nrow(sm2))) %>%
    pivot_longer(-(x:type), names_to = 'bf')
basis <- pivot_longer(basis, -(x:type), names_to = 'bf')

p3 <- ggplot(data, aes(x = x, y = ycent)) +
    geom_point(aes(y = yobs), alpha = 0.2) +
    geom_line(data = basis,
              mapping = aes(x = x, y = value, colour = bf),
              lwd = 1, alpha = 0.5) +
    geom_line(data = wtbasis,
              mapping = aes(x = x, y = fitted), lwd = 1, colour = 'black', alpha = 0.75) +
    guides(colour = "none") +
    labs(y = 'f(x)', x = 'x') +
    theme_minimal(base_size = 16, base_family = 'Fira Sans')

crs_fit <- animate(p3 + transition_states(type, transition_length = 4, state_length = 2) + 
                   ease_aes('cubic-in-out'),
                   nframes = 100, height = anim_height, width = anim_width, res = anim_res,
                   dev = anim_dev)

anim_save('./resources/gam-crs-animation.gif', crs_fit)


## ----lyme-borreliosis-setup, echo = FALSE, results = "hide"-------------------
# Analysis of the Lyme borreliosis weekly incidence data from Norway

# Data are from:
# The emergence and shift in seasonality of Lyme borreliosis in Northern Europe
# Goren et al (2023) doi.org/10.1098/rspb.2022.2420

# Packages
pkgs <- c("tibble", "readr", "dplyr", "here", "tidyr", "ggplot2", "mgcv",
    "gratia")
vapply(pkgs, library, FUN.VALUE = logical(1L), character.only = TRUE,
    logical.return = TRUE)

f <- here("data", "goren-et-al", "lyme-borreliosis.csv")

lyme <- read_csv(f,
    col_types = "-iiiiiiiiiii")

lyme <- lyme |>
    mutate(date = as.Date(paste(yrwk, "1"), format = "%Y%U%u"),
    fyear = factor(year),
    pop = pop / 100000) |>
    select(-yrwk)


## ----insert-tick-img----------------------------------------------------------
knitr::include_graphics("resources/Adult_deer_tick.jpg")


## ----insert-tick-schematic, out.width = "90%"---------------------------------
knitr::include_graphics("resources/gilbert-tick-zoonoses-schematic.png")


## ----lyme-cases-plot, echo = FALSE, out.width = "95%", fig.align = "center"----
lyme |>
    ggplot(aes(x = week, y = cases, group = year, colour = year)) +
    geom_line(linewidth = 0.75)  +
    scale_color_viridis_c(option = "magma",
        guide = guide_colourbar(barheight = unit(0.5, "npc"))) +
    theme_grey(base_size = 18) +
    labs(x = "Week", y = "Number of cases per week", colour = NULL)


## ----lyme-gam-code, echo = TRUE, cache = TRUE---------------------------------
m <- bam(cases ~ s(week, bs = "cc", k = 20) +
                 s(year, k = 25) +
                 ti(week, year, bs = c("cc", "tp"), k = c(20, 20)) +
                 offset(log(pop)),
    data = lyme,
    family = nb(),
    method = "fREML",
    knots = list(week = c(0.5, 52.5)),
    discrete = TRUE,
    nthreads = 4)


## ----lyme-draw, cache = TRUE, out.width = "95%", fig.align = "center"---------
draw(m, rug = FALSE) + plot_layout(ncol = 2) & theme_grey(base_size = 14)


## ----lyme-plot-fitted-curves, cache = TRUE, out.width = "95%", fig.align = "center"----
# predict to find week of peak in each year
ds <- data_slice(m, week = evenly(week, n = 100), year = unique(year),
    pop = log(100000))

fv <- fitted_values(m, data = ds, scale = "response")

fv |>
    ggplot(aes(x = week, y = fitted, group = year, colour = year)) +
    geom_line(linewidth = 1) +
    scale_colour_viridis_c(option = "magma",
        guide = guide_colourbar(barheight = unit(0.5, "npc"))) +
    labs(y = "Number of cases per 100,000 population", x = "Week",
        title = "Annual seasonal trend in cases of Lyme borreliosis",
        colour = NULL) +
    theme_grey(base_size = 18)


## ----lyme-posterior-simulation, cache = TRUE, echo = TRUE---------------------
ds2 <- data_slice(m, week = unique(week), year = unique(year), pop = log(100000))

ds2 <- ds2 |>
    add_column(row = seq_len(nrow(ds2)), .before = 1L)

fs2 <- fitted_samples(m, data = ds2, n = 10000, seed = 42, scale = "response")

fs2 <- fs2 |>
    left_join(ds2, by = join_by(row)) |>
    select(-pop)

peak_week <- fs2 |>
    group_by(year, draw) |>
    slice_max(order_by = fitted, n = 1) |>
    group_by(year) |>
    summarise(peak = quantile(week, prob = 0.5),
        lower_ci = quantile(week, prob = 0.055),
        upper_ci = quantile(week, prob = 0.945))


## ----lyme-plot-peak-week, cache = TRUE, out.width = "95%", fig.align = "center", dependson = -1----
peak_week |>
    ggplot(aes(x = year, y = peak)) +
    geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci)) +
    labs(y = "Week of peak cases", x = NULL,
        caption = "(89% credible interval)") +
    theme_grey(base_size = 18)

