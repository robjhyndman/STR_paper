get_trend <- function(l, type=c("deterministic","stochastic")) {
  type = match.arg(type)
  if(type=="stochastic") 
    trend <- rnorm(2*l, 0, 1) %>% cumsum() %>% cumsum() %>% tail(l)
  else
    trend <- rnorm(1) * (seq(l) + rnorm(1, -1, l / 2))^2 
  trend %>% scale() %>% as.numeric()
}

get_seasonal <- function(l, sl, k=5L,  type=c("deterministic", "stochastic")) {
  type = match.arg(type)
  if(type=="stochastic") {
    season <- rnorm(sl, 0, 1) %>%
      scale() %>% rep(l %/% sl + 1) %>%
      cumsum() %>% scale() %>% cumsum() 
  } else {
    k_sin <- rnorm(k)
    k_cos <- rnorm(k)
    s <- 0
    for (i in seq(k)) {
      s <- s + k_sin[i] * sin(i * 2 * pi * seq(sl) / sl) +
               k_cos[i] * cos(i * 2 * pi * seq(sl) / sl)
    }
    season <- rep(s, rep(l %/% sl + 1))
  }
  season %>% head(l) %>% scale() %>% as.numeric()
}

get_msts_data <- function(l = 365 * 3 + 1, wl = 7, yl = 365, 
                          alpha = 1, beta = 1, gamma = 0.25,  
                          type=c("deterministic","stochastic"),
                          k = 5L) {
  type = match.arg(type)
  tibble(
      trend = get_trend(l = l, type = type),
      season_week = alpha * get_seasonal(l = l, sl = wl, type = type),
      season_year = beta * get_seasonal(l = l, sl = yl, type = type, k = k),
      remainder = gamma * rnorm(l, 0, 1)
  ) %>%
  mutate(data = trend + season_week + season_year + remainder)
}

get_predictors <- function(l=365*3+1, wl=7, yl=365, lambda1=1, lambda2=0.1) {
  list(
    Trend = list(
      name = "Trend",
      data = rep(1, l),
      times = seq(l),
      seasons = rep(1, l),
      timeKnots = seq(from = 1, to = l, length.out = yl + 1),
      seasonalStructure = list(segments = list(c(0, 1)), sKnots = list(c(1, 0))),
      lambdas = c(lambda1, 0, 0)
    ),
    WSeason = list(
      name = "Weekly seas",
      data = rep(1, l),
      times = seq(l),
      seasons = head(rep(seq(1, wl) - 1, l %/% wl + 1), l),
      timeKnots = seq(from = 1, to = l, length.out = 2),
      seasonalStructure = list(segments = list(c(0, wl)), sKnots = c(as.list(seq(wl - 1)), list(c(0, wl)))),
      lambdas = c(0, 0, -10000)
    ),
    YSeason = list(
      name = "Yearly seas",
      data = rep(1, l),
      times = seq(1, l),
      seasons = head(rep(seq(1, yl) - 1, l %/% yl + 1), l),
      timeKnots = seq(from = 1, to = l, length.out = 2),
      seasonalStructure = list(segments = list(c(0, yl)), sKnots = c(as.list(seq(yl - 1)), list(c(0, yl)))),
      lambdas = c(0, lambda2, -10000)
    )
  )
}

get_components <- function(data, method = c("STL","STR","TBATS"), l=365*3+1, wl=7, yl=365) {
  method <- match.arg(method)
  mts_data <- forecast::msts(data$data, seasonal.periods=c(wl,yl))
  if (method == "STL") {
    d <- forecast::mstl(mts_data, iterate = 5)
    trend <- as.vector(d[, 2])
    seasonal1 <- as.vector(d[, 3])
    seasonal2 <- as.vector(d[, 4])
    remainder <- as.vector(d[, 5])
  } else if (method == "STR") {
    d <- STR(
      data = as.vector(mts_data),
      predictors = get_predictors(l, wl, yl),
      gapCV = 20,
      reltol = 0.00001,
      trace = FALSE
    )
    trend <- d$output$predictors[[1]]$data
    seasonal1 <- d$output$predictors[[2]]$data
    seasonal2 <- d$output$predictors[[3]]$data
    remainder <- d$output$random$data
  } else if (method == "TBATS") {
    d <- forecast::tbats(mts_data, use.box.cox = FALSE, use.trend = TRUE)
    d_comp <- forecast::tbats.components(d)
    trend <- as.vector(d_comp[, 2])
    seasonal1 <- as.vector(d_comp[, 4])
    seasonal2 <- as.vector(d_comp[, 5])
    remainder <- as.vector(mts_data) - (trend + seasonal1 + seasonal2)
  } else {
    stop("Unknown method to run decomposition...")
  }
  return(list(trend = trend, season_week = seasonal1, 
              season_year = seasonal2, remainder = remainder))
}

compute_errors <- function(l=365*3+1, wl=7, yl=365, alpha=1, beta=1, gamma=0.25, 
                           type=c("deterministic", "stochastic")) {
  type = match.arg(type)
  d <- get_msts_data(l=l, wl=wl, yl=yl, alpha=alpha, beta=beta, gamma=gamma, type=type)
  result <- NULL
  for (method in c("STL", "TBATS", "STR")) {
    comp <- get_components(d, method=method, l=l, wl=wl, yl=yl)
    result <- rbind(result,
      data.frame(
        method = method,
        trend = d$trend - comp$trend,
        season_week = d$season_week - comp$season_week,
        season_year = d$season_year - comp$season_year,
        remainder = d$remainder - comp$remainder
      )
    )
  }
  return(result)
}

compute_n_errors <- function(n=20, l=365*3+1, wl=7, yl=365, alpha=1, beta=1, gamma=0.25, 
                                type=c("deterministic", "stochastic")) {
  result <- NULL
  for (i in seq(n)) {
    df <- compute_errors(l=l, wl=wl, yl=yl, alpha=alpha, beta=beta, gamma=gamma, type=type) %>%
      mutate(rep=i)
    result <- rbind(result, df)
  }
  as_tibble(result) %>% mutate(alpha=alpha, beta=beta, gamma=gamma, type = type)
}

get_pvalues_error <- function(errors, var, type, gamma) {
  errors <- errors %>%
    filter(type==type, gamma==gamma) %>%
    mutate(y = errors[[var]]^2)
  lm(y ~ method,  data=errors) %>%
    broom::tidy() %>%
    filter(str_detect(term,"method")) %>%
    select(term, p.value) %>%
    mutate(type=type,gamma=gamma,var=var)
}

get_pvalues <- function(errors, type, gamma) {
  bind_rows(
    get_pvalues_error(errors, "trend", type, gamma),
    get_pvalues_error(errors, "season_week", type, gamma),
    get_pvalues_error(errors, "season_year", type, gamma),
    get_pvalues_error(errors, "remainder", type, gamma)
  )
}
