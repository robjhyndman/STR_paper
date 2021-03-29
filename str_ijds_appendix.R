library(tidyverse)
library(stR)
library(fabletools)
library(distributional)
library(lubridate)
library(forecast)
library(imputeTS)
library(tsibble)

# Create predictor dist object
predictor_dist <- function(x, level) {
    dist_percentile(
      with(x, mapply(c, data, split(lower, row(lower)), split(upper, row(upper)), SIMPLIFY = FALSE)),
      list(c(50, 50 + rep(c(-1, 1), each = length(level))*rep(level, 2)/2))
    )
}


# Supermarket and grocery store log revenue in New South Wales
supermarket_nsw <- tsibbledata::aus_retail %>%
  filter(
    State == "New South Wales",
    Industry == "Supermarket and grocery stores",
    year(Month) >= 2000,
    year(Month) <= 2009
  ) %>%
  pull(Turnover) %>%
  log() %>%
  ts(start=2000, frequency=12)

# STR decomposition
d <- supermarket_nsw %>%
  AutoSTR(confidence=0.95, gapCV=12, reltol=0.00001)

# STL decomposition
d_stl = forecast::mstl(supermarket_nsw, iterate=5)

# TBATS decomposition
d_tbats = tbats(supermarket_nsw, use.box.cox=FALSE, use.trend=TRUE)

# X-13-ARIMA-SEATS decomposition
d_x13 <- seasonal::seas(supermarket_nsw, transform.function = "none")

# Make data frame with components
decomp <- bind_rows(
    tibble(
      method = "STR",
      trend = d$output$predictors[[1]]$data,
      seasonal = d$output$predictors[[2]]$data,
      remainder = d$output$random$data,
      month = yearmonth(seq(from=as.Date("2000-01-01"),by="1 month", length=120))
    ),
    as_tibble(d_stl) %>%
      transmute(
        method = "STL",
        trend=Trend,
        seasonal=Seasonal12,
        remainder=Remainder,
        month = yearmonth(seq(from=as.Date("2000-01-01"),by="1 month", length=120))
      ),
    tbats.components(d_tbats) %>%
      as_tibble() %>%
      transmute(
        method = "TBATS",
        trend = level,
        seasonal = season,
        remainder = observed - level - season,
      month = yearmonth(seq(from=as.Date("2000-01-01"),by="1 month", length=120))
      ),
    tibble(
      method = "X-13-ARIMA-SEATS",
      trend = seasonal::trend(d_x13),
      #seasonal = seasonal::series(d_x13, "seats.seasonal"),
      remainder = seasonal::irregular(d_x13),
      seasonal = supermarket_nsw - trend - remainder,
      month = yearmonth(seq(from=as.Date("2000-01-01"),by="1 month", length=120))
    )
  ) %>%
  mutate(observed = trend + seasonal + remainder)


## Plot results
decomp %>%
  select(-observed) %>%
  bind_rows(tsibble(
    method="Data",
    data = supermarket_nsw,
    month = yearmonth(seq(from=as.Date("2000-01-01"),by="1 month", length=120)),
    index=month, key=method
  )) %>%
  select(month, method, everything()) %>%
  pivot_longer(data:remainder,
               names_to = "variable", values_to="value") %>%
  mutate(
    variable = factor(variable, levels=c("data","trend","seasonal","remainder")),
    method = factor(method, levels=c("Data","STR","STL","TBATS","X-13-ARIMA-SEATS"))
  ) %>%
  ggplot(aes(x=month, y=value, group=method, col=method)) +
  geom_line() +
  facet_grid(variable ~ ., scales = "free_y") +
  scale_color_manual(values = c(Data = "black",STL="#D55E00",STR="#0072B2",TBATS="#009E73",
                                `X-13-ARIMA-SEATS` = "#CC79A7")) +
  guides(colour = guide_legend(title = "  ")) +
  labs(y="", x="Month")

