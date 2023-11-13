# Functions 

GIF.convert <- function(x, output) #Create a function to read, animate and convert the files to gif
{
  image_read(x) %>%
    image_animate(fps = 1) %>%
    image_write(output)
}

check_dots <- DHARMa:::checkDots
dispersion_fct <- DHARMa:::testDispersion
ensure_dharma <- DHARMa:::ensureDHARMa
ensure_predictor <- DHARMa:::ensurePredictor
get_p_val <- DHARMa:::getP
levene_test <- DHARMa:::leveneTest.formula
outliers_fct <- DHARMa:::testOutliers
uniformity_fct <- DHARMa:::testUniformity
test_categorical <- DHARMa:::testCategorical
test_quantiles <- DHARMa:::testQuantiles

make_dharma_res <- function(M) {
  response <- predict(M, newdata = STObj)$newdata@data$COUNT 
  preds <- simulate(M) 
  DHARMa::createDHARMa(
    simulatedResponse = preds,
    observedResponse = response,
    integerResponse = TRUE
  )
}

gg_disp_hist <- function(sim_list, alternative = c("two.sided", "greater",
                                                   "less")) {
  sim_list <- ensure_dharma(sim_list, convert = "Model")
  expected_var <- sd(sim_list$simulatedResponse) ^ 2
  spread <- function(x, sim_list, expected_var) {
    var(x - sim_list$fittedPredictedResponse) / expected_var
  }
  alternative <- match.arg(alternative)
  observed <- spread(sim_list$observedResponse, sim_list, expected_var)
  simulated <- apply(sim_list$simulatedResponse, 2, spread, sim_list,
                     expected_var)
  p_val_ <- get_p_val(simulated = simulated, observed = observed,
                      alternative = alternative)
  x_lab_ <- paste("Simulated values, red line = fitted model. p-value (",
                  alternative, ") = ", p_val_, sep = "")
  data.frame(simulated = simulated) %>%
    ggplot(data = .) +
    geom_histogram(mapping = aes(x = simulated), colour = "grey60",
                   bins = max(round(sim_list$nSim / 5), 20)) +
    geom_vline(xintercept = observed, colour = "tomato") +
    labs(x = x_lab_, y = "",
         title = "DHARMa nonparametric dispersion test via sd of residuals",
         subtitle = "Fitted vs. Simulated") +
    theme_bw() +
    theme(plot.title = element_text(size = 9),
          plot.subtitle = element_text(size = 9))
}


plot_qq_unif <- function(sim_listrix, test_uniformity = TRUE,
                         test_outliers = TRUE, test_dispersion = TRUE, ...) {
  sim_listrix <- ensure_dharma(sim_listrix, convert = "Model")
  res_ <- data.frame(y = sim_listrix$scaledResiduals) %>%
    dplyr::mutate(x = seq_len(dplyr::n()) / (dplyr::n() + 1))
  res_ <- qqplot(res_$x, res_$y, plot.it = FALSE) %>%
    as.data.frame
  p_ <- ggplot(data = res_) +
    geom_point(mapping = aes(x = x, y = y), shape = 16,
               colour = "skyblue", alpha = 0.8) +
    geom_abline(slope = 1, linetype = 2) +
    labs(x = "Expected", y = "Observed", subtitle = "QQ plot residuals") +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    theme_bw() +
    theme(plot.title = element_text(size = 9))
  if (test_uniformity) {
    tmp_ks <- uniformity_fct(sim_listrix, plot = FALSE)
    k_lab_a <- paste("KS test: p =", round(tmp_ks$p.value, digits = 5))
    k_lab_b <- paste("Deviation:", ifelse(tmp_ks$p.value < 0.05, "Significant",
                                          "N.S."))
    p_ <- p_ +
      annotate("text", x = 0, y = 0.98, hjust = 0, vjust = 0.5, label = k_lab_a,
               colour = ifelse(tmp_ks$p.value < 0.05, "tomato", "black"),
               size = 3) +
      annotate("text", x = 0, y = 0.90, hjust = 0, vjust = 0.5, label = k_lab_b,
               colour = ifelse(tmp_ks$p.value < 0.05, "tomato", "black"),
               size = 3)
  }
  if (test_outliers) {
    tmp_ot <- outliers_fct(sim_listrix, plot = FALSE)
    ot_lab_a <- paste("Outlier test: p =", round(tmp_ot$p.value, digits = 5))
    ot_lab_b <- paste("Deviation:", ifelse(tmp_ot$p.value < 0.05, "Significant",
                                           "N.S."))
    p_ <- p_ +
      annotate("text", x = 0, y = 0.78, hjust = 0, vjust = 0.5, size = 3,
               label = ot_lab_a, colour = ifelse(tmp_ot$p.value < 0.05,
                                                 "tomato", "black")) +
      annotate("text", x = 0, y = 0.70, hjust = 0, vjust = 0.5, size = 3,
               label = ot_lab_b, colour = ifelse(tmp_ot$p.value < 0.05,
                                                 "tomato", "black"))
  }
  p_
}

plot_residuals <- function(sim_list, form = NULL, quantreg = NULL, rank = TRUE,
                           as_factor = NULL, smooth_scatter = NULL,
                           quantiles = c(0.25, 0.5, 0.75), ...) {
  ##### Checks #####
  if (is.null(form)) {
    xlab_ <- check_dots("xlab", ifelse(rank, "Model predictions (rank transformed)", "Model predictions"), ...)
  } else {
    xlab_ <- "Predictor"
  }
  sim_list <- ensure_dharma(sim_list, convert = TRUE)
  res <- sim_list$scaledResiduals
  if (inherits(form, "DHARMa")) {
    stop("DHARMa::plot_residuals > argument form cannot be of class DHARMa. ",
         "Note that the syntax of plot_residuals has changed since DHARMa ",
         "0.3.0. See ?plot_residuals.")
  }
  pred <- ensure_predictor(sim_list, form)
  
  ##### Rank transform and factor conversion #####
  if (!is.factor(pred)) {
    if (rank) {
      pred <- rank(pred, ties.method = "average")
      pred <- pred / max(pred)
    }
    nuniq <- length(unique(pred))
    ndata <- length(pred)
    if (is.null(as_factor)) {
      as_factor <- (nuniq == 1) | (nuniq < 10 & ndata / nuniq > 10)
    }
    if (as_factor) {
      pred <- factor(pred)
    }
  }
  
  ##### Residual scatter plots #####
  if (is.null(quantreg)) {
    if (length(res) > 2000) {
      quantreg <- FALSE
    } else {
      quantreg <- TRUE
    }
  }
  switch_scatter <- 1e4
  if (is.null(smooth_scatter)) {
    if (length(res) > switch_scatter) {
      smooth_scatter <- TRUE
    } else {
      smooth_scatter <- FALSE
    }
  }
  blackcol <- rgb(0, 0, 0,
                  alpha = max(0.1, 1 - 3 * length(res) / switch_scatter))
  if (is.factor(pred)) {
    out <- cat_gg_res(sim_list = sim_list, cat_pred = pred,
                      quantiles = quantiles)
  } else if (smooth_scatter) {
    def_col_ <- ifelse(res == 0 | res == 1, 2, blackcol)
    df_ <- data.frame(x = pred, y = res, col = def_col_ == 2)
    out <- ggplot(data = df_, mapping = aes(x = x, y = y)) +
      stat_density2d(mapping = aes(fill = ..density..^0.25),
                     geom = "tile", contour = FALSE, n = 200) +
      geom_point(mapping = aes(x = x, y = y), size = 0.2, shape = 16,
                 colour = "black") +
      geom_point(data = df_ %>% dplyr::filter(col),
                 mapping = aes(x = x, y = y), size = 0.2, shape = 16,
                 colour = "tomato") +
      scale_fill_continuous(low = "white", high = "dodgerblue4") +
      ylim(c(0, 1)) +
      theme_bw()
  } else {
    symbol <- ifelse(res == 0 | res == 1, "Yes", "No")
    df_ <- data.frame(x = pred, y = res, extreme = symbol)
    out <- ggplot() +
      geom_point(data = df_ %>% dplyr::filter(extreme == "No"),
                 mapping = aes(x = x, y = y), colour = "black", shape = 1,
                 size = 1) +
      geom_point(data = df_ %>% dplyr::filter(extreme == "Yes"),
                 mapping = aes(x = x, y = y), colour = "tomato", shape = 8,
                 size = 1) +
      ylim(c(0, 1)) +
      labs(y = "Standardized residual", x = xlab_) +
      theme_bw()
  }
  
  # ##### Quantile regressions #####
  tit_ <- check_dots("main", "Residual vs. predicted", ...)
  tmp_ <- NULL
  if (is.numeric(pred)) {
    if (!quantreg) {
      out <- out +
        labs(title = tit_) +
        geom_hline(yintercept = quantiles, linetype = 2, colour = "lightgrey")
      sspl_ <- try({
        smooth.spline(pred, res, df = 10)
      }, silent = TRUE)
      if (inherits(sspl_, "try-error")) {
        out <- out +
          geom_line(data = data.frame(x = sspl_$x, y = sspl_$y),
                    mapping = aes(x = x, y = y), linetype = 2,
                    colour = "tomato") +
          geom_hline(yintercept = 0.5, linetype = 1, colour = "tomato")
      }
    } else {
      tmp_ <- test_quantiles(sim_list, pred, quantiles = quantiles,
                             plot = FALSE)
      if (any(tmp_$pvals < 0.05, na.rm = TRUE)) {
        tit_ <- paste(tit_, "Quantile deviations detected (red curves)",
                      sep = "\n")
        if (tmp_$p.value <= 0.05) {
          tit_ <- paste(tit_, "Combined adjusted quantile test: Significant",
                        sep = "\n")
        } else {
          tit_ <- paste(tit_, "Combined adjusted quantile test: N.S.",
                        sep = "\n")
        }
        maincol <- "tomato"
      } else {
        tit_ <- paste(tit_, "No significant problems detected", sep = "\n")
        maincol <- "black"
      }
      out <- out +
        labs(title = tit_) +
        theme(plot.title = element_text(size = 9, colour = maincol))
      for (i in seq_along(quantiles)) {
        line_col <- ifelse(tmp_$pvals[i] <= 0.05 & !(is.na(tmp_$pvals[i])),
                           "tomato", "black")
        pol_df_ <- data.frame(
          x = c(tmp_$predictions$pred, rev(tmp_$predictions$pred)),
          y = c(tmp_$predictions[, 2 * i] - tmp_$predictions[, 2 * i + 1],
                rev(tmp_$predictions[, 2 * i] + tmp_$predictions[, 2 * i + 1]))
        )
        pred_df_ <- data.frame(x = tmp_$predictions$pred,
                               y = tmp_$predictions[, 2 * i])
        out <- out +
          geom_hline(yintercept = quantiles[i], colour = line_col, lwd = 0.5,
                     linetype = 2) +
          geom_polygon(data = pol_df_, mapping = aes(x = x, y = y),
                       fill = "#00000020") +
          geom_line(data = pred_df_, mapping = aes(x = x, y = y),
                    colour = line_col, lwd = 0.5)        
      }
      subt_ <- c(paste("Quantile test: p =", round(tmp_$p.value, digits = 5)),
                 paste("Deviation:", ifelse(tmp_$p.value < 0.05, "Significant",
                                            "N.S.")))
      subt_col_ <- ifelse(tmp_$p.value < 0.05, "tomato", "black")
      out <- out +
        labs(subtitle = subt_) +
        theme(plot.subtitle = element_text(size = 9, colour = subt_col_))
    }
  }
  out
}

gg_dharma <- function(x, ...) {
  a_ <- plot_qq_unif(x)
  b_ <- plot_residuals(x, ...)
  c_ <- gg_disp_hist(x, alternative = "greater")
  (a_ + b_) / c_ + patchwork::plot_annotation(title = "DHARMa residual diagnostics") + plot_layout(nrow = 2, widths = c(2, 1))
}

#### Predictive measures 
## 95% coverage
coverage95 <- function(z, lower, upper) {
  sum((z < upper) & (z > lower)) / length(z)
}

## 95% interval score
IS95 <- function(true, lower, upper) {
  alpha = 0.1
  pred95l <- lower 
  pred95u <- upper
  ISs <- (pred95u - pred95l) + 2/alpha * (pred95l - true) * (true < pred95l) +
    2/alpha * (true - pred95u) * (true > pred95u)
  mean(ISs)
}

## Root-mean-squared prediction error
RMSPE <- function(z,pred) {
  Y <- (z - pred)^2
  sqrt(mean(Y))
}
