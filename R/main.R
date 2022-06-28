#' segen
#'
#' @param df A data frame with time features on columns
#' @param seq_len Positive integer. Time-step number of the forecasting sequence. Default: NULL (automatic selection between 2 and max limit).
#' @param similarity Positive numeric. Degree of similarity between two sequences, based on quantile conversion of distance. Default: NULL (automatic selection between 0.01, maximal difference, and 0.99, minimal difference).
#' @param dist_method String. Method for calculating distance among sequences. Available options are: "euclidean", "manhattan", "maximum", "minkowski". Default: NULL (random search).
#' @param rescale Logical. Flag to TRUE for min-max scaling of distances. Default: NULL (random search).
#' @param smoother Logical. Flag to TRUE for loess smoothing. Default: FALSE.
#' @param ci Confidence interval for prediction. Default: 0.8
#' @param error_scale String. Scale for the scaled error metrics (for continuous variables). Two options: "naive" (average of naive one-step absolute error for the historical series) or "deviation" (standard error of the historical series). Default: "naive".
#' @param error_benchmark String. Benchmark for the relative error metrics (for continuous variables). Two options: "naive" (sequential extension of last value) or "average" (mean value of true sequence). Default: "naive".
#' @param n_windows Positive integer. Number of validation windows to test prediction error. Default: 10.
#' @param n_samp Positive integer. Number of samples for random search. Default: 30.
#' @param dates Date. Vector with dates for time features.
#' @param seed Positive integer. Random seed. Default: 42.
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
#' @return This function returns a list including:
#' \itemize{
#' \item exploration: list of all not-null models, complete with predictions and error metrics
#' \item history: a table with the sampled models, hyper-parameters, validation errors
#' \item best_model: results for the best selected model according to the weighted average rank, including:
#' \itemize{
#' \item predictions: for continuous variables, min, max, q25, q50, q75, quantiles at selected ci, mean, sd, mode, skewness, kurtosis, IQR to range, risk ratio, upside probability and divergence for each point fo predicted sequences; for factor variables, min, max, q25, q50, q75, quantiles at selected ci, proportions, difformity (deviation of proportions normalized over the maximum possible deviation), entropy, upgrade probability and divergence for each point fo predicted sequences
#' \item testing_errors: testing errors for each time feature for the best selected model (for continuous variables: me, mae, mse, rmsse, mpe, mape, rmae, rrmse, rame, mase, smse, sce, gmrae; for factor variables: czekanowski, tanimoto, cosine, hassebrook, jaccard, dice, canberra, gower, lorentzian, clark)
#' \item plots: standard plots with confidence interval for each time feature
#' }
#' \item time_log
#' }
#'
#' @export
#'
#' @import purrr
#' @import tictoc
#' @importFrom readr parse_number
#' @importFrom lubridate seconds_to_period is.Date as.duration
#' @importFrom stats weighted.mean ecdf na.omit
#' @importFrom imputeTS na_kalman
#' @importFrom fANCOVA loess.as
#' @importFrom modeest mlv1
#' @import ggplot2
#' @importFrom moments skewness kurtosis
#' @importFrom stats quantile sd lm rnorm pnorm fft runif
#' @importFrom scales number
#' @importFrom utils tail head
#' @import greybox
#' @importFrom philentropy distance
#' @importFrom entropy entropy
#' @importFrom Rfast Dist
#' @importFrom narray split


#'@examples
#'segen(time_features, seq_len = 30, similarity = 0.7, n_windows = 3, n_samp = 1)
#'
#'

segen <- function(df, seq_len = NULL, similarity = NULL, dist_method = NULL, rescale = NULL, smoother = FALSE, ci = 0.8,
error_scale = "naive", error_benchmark = "naive", n_windows = 10, n_samp = 30, dates = NULL, seed = 42)
{
  tic.clearlog()
  tic("time")

  set.seed(seed)
  n_length <- nrow(df)
  n_feats <- ncol(df)
  feat_names <- colnames(df)

  feat_n_class <- NULL
  feat_level_names <- NULL

  all_classes <- all(map_lgl(df, ~ is.factor(.x) | is.character(.x)))
  all_numbers <- all(map_lgl(df, ~ is.numeric(.x) | is.integer(.x)))

  if(!all_classes & !all_numbers){stop("only numbers or factors, but not both")}

  event_class <- all_classes & !all_numbers
  if(anyNA(df) & event_class == FALSE){df <- as.data.frame(na_kalman(df)); message("kalman imputation on time features\n")}
  if(smoother == TRUE & event_class == FALSE){df <- as.data.frame(map(df, ~ suppressWarnings(loess.as(x=1:n_length, y=.x)$fitted))); message("performing optimal smoothing\n")}

  if(event_class == TRUE)
  {
    df <- as.data.frame(map(df, ~ as.factor(.x)))
    feat_level_names <- map(df, ~ levels(.x))
    feat_n_class <- map_dbl(feat_level_names, ~ length(.x))
    df <- as.data.frame(data.matrix(df)-1)
  }

  if(length(seq_len) == 1 & length(similarity) == 1 & length(rescale) == 1){n_samp <- 1}

  deriv <- map_dbl(df, ~ best_deriv(.x))
  max_limit <- floor((n_length/(n_windows + 1) - max(deriv))/2)###AT LEAST TWO ROW IN THE SEGMENTED FRAME
  sqln_set <- sampler(seq_len, n_samp, range = c(2, max_limit), integer = TRUE)
  sml_set <- sampler(similarity, n_samp, range = c(0.01, 0.99), integer = FALSE)
  dst_set <- sampler(dist_method, n_samp, range = c( "euclidean", "manhattan", "maximum", "minkowski"), integer = FALSE, multi = n_feats)
  rscl_set <- as.logical(sampler(rescale, n_samp, range = c(0, 1), integer = TRUE))

  if(any(sqln_set < 2)){sqln_set[sqln_set < 2] <- 2; message("setting min seq_len to 2")}
  if(any(sqln_set > max_limit)){sqln_set[sqln_set > max_limit] <- max_limit; message(paste0("setting max seq_len to ", max_limit))}

  exploration <- mapply(function(ftn) pmap(list(sqln_set, sml_set, dst_set, rscl_set), ~ windower(df[, ftn], seq_len = ..1, similarity = ..2, dist_method = ..3[ftn], rescale = ..4, n_windows, ci, error_scale, error_benchmark, feat_n_class[ftn], dates, feat_name = feat_names[ftn], feat_level_names[[ftn]], seed)), ftn = 1:n_feats, SIMPLIFY = FALSE)
  exploration <- transpose(exploration)

  models <- map_depth(exploration, 2, ~.x$quant_pred)
  errors <- map_depth(exploration, 2, ~.x$errors)
  if(n_feats > 1){aggr_errors <- map(errors, ~ colMeans(Reduce(rbind, .x)))} else {aggr_errors <- errors}
  aggr_errors <- t(as.data.frame(aggr_errors))

  history <- data.frame(seq_len = unlist(sqln_set))
  history$similarity <- sml_set
  history$dist_method <- dst_set
  history$rescale <- rscl_set
  history <- data.frame(history, round(aggr_errors, 4))
  rownames(history) <- NULL

  if(event_class == FALSE){history <- ranker(history, focus = -c(1, 2, 3, 4), inverse = NULL, absolute = c("me", "mpe", "sce"), reverse = FALSE)}
  if(event_class == TRUE){history <- ranker(history, focus = -c(1, 2, 3, 4), inverse = NULL, absolute = NULL, reverse = FALSE)}

  best_index <- as.numeric(rownames(history[1,]))
  predictions <- models[[best_index]]

  testing_errors <- t(as.data.frame(errors[[best_index]]))
  if(event_class == FALSE){plots <- pmap(list(predictions, feat_names), ~ plotter(..1, ci, df[, ..2], n_class = NULL, level_names = NULL, dates, ..2))}
  if(event_class == TRUE){plots <- pmap(list(predictions, feat_names, feat_n_class, feat_level_names), ~ plotter(..1, ci, df[, ..2], n_class = ..3, level_names = ..4, dates, ..2))}

  names(predictions) <- feat_names
  rownames(testing_errors) <- feat_names
  names(plots) <- feat_names

  best_model <- list(exploration = exploration, predictions = predictions, testing_errors = testing_errors, plots = plots)

  toc(log = TRUE)
  time_log <- seconds_to_period(round(parse_number(unlist(tic.log())), 0))

  outcome <- list(exploration = exploration, history = history, best_model = best_model, time_log = time_log)

  return(outcome)
}


windower <- function(ts, seq_len, similarity, dist_method, rescale = FALSE, n_windows = 10, ci = 0.8, error_scale = "naive", error_benchmark = "naive", n_class = NULL, dates = NULL, feat_name = NULL, level_names = NULL, seed = 42)
{
  n_length <- length(ts)
  idx <- c(rep(1, n_length%%(n_windows + 1)), rep(1:(n_windows + 1), each = n_length/(n_windows + 1)))

  window_results <- map(1:n_windows, ~ engine(ts[idx <= .x], seq_len, similarity, dist_method, rescale, ci, holdout = head(ts[idx == (.x + 1)], seq_len), error_scale, error_benchmark, n_class, dates, feat_name, level_names, seed))
  errors <- colMeans(Reduce(rbind, map(window_results, ~ .x$testing_error)))
  pred_scores <- rowMeans(Reduce(cbind, map(window_results, ~ .x$quant_pred$pred_scores)))
  model <- engine(ts, seq_len, similarity, dist_method, rescale, ci, holdout = NULL, error_scale, error_benchmark, n_class, dates, feat_name, level_names, seed)
  quant_pred <- model$quant_pred
  quant_pred <- cbind(quant_pred, pred_scores = pred_scores)
  plot <- model$plot

  outcome <- list(quant_pred = quant_pred, errors = errors, plot = plot)

  return(outcome)
}



engine <- function(ts, seq_len, similarity, dist_method, rescale = FALSE, ci = 0.8, holdout = NULL, error_scale = "naive", error_benchmark = "naive", n_class = NULL, dates = NULL, feat_name = NULL, level_names = NULL, seed = 42)
{
  diffmodel <- recursive_diff(ts, best_deriv(ts))
  dts <- diffmodel$vector
  reframed <- smart_reframer(dts, seq_len, seq_len)
  n_row <- nrow(reframed)
  dmat <- as.matrix(Dist(reframed, method = dist_method, p = 3))
  dmat[upper.tri(dmat)]<- NA
  smat <- Reduce(rbind, map(split(dmat, along = 1), ~ .x < quantile(.x, probs = 1 - similarity, na.rm = TRUE)))
  rownames(smat) <- NULL
  smat[1, 1] <- TRUE
  location <- apply(smat, 2, mean, na.rm = TRUE)
  deviation <- apply(smat, 2, sd, na.rm = TRUE)
  deviation[n_row] <- 0
  weights <- replicate(1000, rnorm(n_row, location, deviation), simplify = FALSE)
  if(rescale == FALSE){weights <- map(weights, ~ {.x[.x < 0] <- 0; return(.x)})}
  if(rescale == TRUE){weights <- map(weights, ~ minmax(.x, 0, 1))}
  raw_pred <- t(as.data.frame(map(weights, ~ colSums((reframed * .x))/sum(.x))))
  rownames(raw_pred) <- NULL
  raw_pred <- t(as.data.frame(map(split(raw_pred, along = 1), ~ invdiff(.x, diffmodel$tail_value))))

  quant_pred <- qpred_fun(raw_pred, holdout, ts, ci, error_scale, error_benchmark, n_class, level_names, seed)

  return(quant_pred)
}


minmax <- function(x, min_v, max_v)
{
  span <- (x - min(x))/(diff(range(x)))
  rescaled <- span * (max_v - min_v) + min_v
  return(rescaled)
}

qpred_fun <- function(raw_pred, holdout = NULL, ts, ci, error_scale = "naive", error_benchmark = "naive", n_class = NULL, level_names = NULL, seed = 42)
{
  set.seed(seed)

  if(!is.null(raw_pred))
  {
    raw_pred <- doxa_filter(ts, raw_pred, n_class)
    quants <- sort(unique(c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)))

    if(is.null(n_class))
    {
      p_stats <- function(x){c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), mode = suppressWarnings(mlv1(x[is.finite(x)], method = "shorth")), kurtosis = suppressWarnings(kurtosis(x[is.finite(x)], na.rm = TRUE)), skewness = suppressWarnings(skewness(x[is.finite(x)], na.rm = TRUE)))}
      quant_pred <- as.data.frame(t(as.data.frame(apply(raw_pred, 2, p_stats))))
      p_value <- apply(raw_pred, 2, function(x) ecdf(x)(seq(min(raw_pred), max(raw_pred), length.out = 1000)))
      divergence <- c(NA, apply(p_value[,-1, drop = FALSE] - p_value[,-ncol(p_value), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE))))
      upside_prob <- c(NA, colMeans(apply(raw_pred[,-1, drop = FALSE]/raw_pred[,-ncol(raw_pred), drop = FALSE], 2, function(x) x > 1)))
      iqr_to_range <- (quant_pred[, "75%"] - quant_pred[, "25%"])/(quant_pred[, "max"] - quant_pred[, "min"])
      median_range_ratio <- (quant_pred[, "max"] - quant_pred[, "50%"])/(quant_pred[, "50%"] - quant_pred[, "min"])
      quant_pred <- round(cbind(quant_pred, iqr_to_range, median_range_ratio, upside_prob, divergence), 4)
    }

    if(is.numeric(n_class) & !is.null(level_names))
    {
      freq <- function(x) {sapply(0:(n_class-1), function(n) sum(x == n))}
      p_stats <- function(x){c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), props = freq(x)/length(x), difformity = sd(freq(x)/length(x)), entropy = entropy(freq(x)))}
      quant_pred <- as.data.frame(t(as.data.frame(apply(raw_pred, 2, p_stats))))
      p_value <- apply(raw_pred, 2, function(x) ecdf(x)(0:(n_class - 1)))
      divergence <- c(NA, apply(p_value[,-1, drop = FALSE] - p_value[,-ncol(p_value), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE))))
      upgrade_prob <- c(NA, colMeans(apply((raw_pred[,-1, drop = FALSE] + 1)/(raw_pred[,-ncol(raw_pred), drop = FALSE] + 1), 2, function(x) x > 1)))###ADDING ONE TO SOLVE ISSUE WITH CLASS ZERO
      iqr_to_range <- (quant_pred[, "75%"] - quant_pred[, "25%"])/(quant_pred[, "max"] - quant_pred[, "min"])
      median_range_ratio <- (quant_pred[, "max"] - quant_pred[, "50%"])/(quant_pred[, "50%"] - quant_pred[, "min"])
      quant_pred <- round(cbind(quant_pred, iqr_to_range, median_range_ratio, upgrade_prob, divergence), 4)

      ###FIXING FACTOR LABELS
      n_quants <- length(quants) + 2###QUANTILES + MIN & MAX
      m_index <- as.matrix(quant_pred[, 1:n_quants] + 1)###INDEXING FROM O TO 1
      releveled <- as.data.frame(matrix(level_names[m_index], nrow(m_index), ncol(m_index)))
      colnames(releveled) <- c("min", paste0(quants * 100, "%"), "max")
      quant_pred <- as.data.frame(cbind(releveled, quant_pred[, c(paste0("props", 1:n_class), "difformity", "entropy", "upgrade_prob", "divergence")]))
      colnames(quant_pred) <- c(colnames(releveled), paste0("prop_", level_names), "difformity", "entropy", "upgrade_prob", "divergence")
    }
  }

  testing_error <- NULL
  if(!is.null(holdout))
  {
    if(is.null(n_class)){
      mean_pred <- colMeans(raw_pred)
      testing_error <- my_metrics(holdout, mean_pred, ts, error_scale, error_benchmark, n_class)}

    if(is.numeric(n_class) & !is.null(level_names)){
      testing_error <- apply(apply(raw_pred, 1, function(x) my_metrics(holdout, x, n_class = n_class)), 1, function(x) mean(x, na.rm = TRUE))}

    pred_scores <- round(prediction_score(raw_pred, holdout), 4)
    quant_pred <- cbind(quant_pred, pred_scores = pred_scores)
  }

  outcome <- list(quant_pred = quant_pred, testing_error = testing_error)
  return(outcome)
}

###
doxa_filter <- function(ts, mat, n_class = NULL)
{
  discrete_check <- all(ts%%1 == 0)
  all_positive_check <- all(ts >= 0)
  all_negative_check <- all(ts <= 0)
  monotonic_increase_check <- all(diff(ts) >= 0)
  monotonic_decrease_check <- all(diff(ts) <= 0)

  monotonic_fixer <- function(x, mode)
  {
    model <- recursive_diff(x, 1)
    vect <- model$vector
    if(mode == 0){vect[vect < 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    if(mode == 1){vect[vect > 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    return(vect)
  }

  if(is.null(n_class))
  {
    if(all_positive_check){mat[mat < 0] <- 0}
    if(all_negative_check){mat[mat > 0] <- 0}
    if(discrete_check){mat <- round(mat)}
    if(monotonic_increase_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 0)))}
    if(monotonic_decrease_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 1)))}
  }

  if(is.numeric(n_class)){mat <- round(mat); mat[mat > (n_class - 1)] <- (n_class - 1); mat[mat < 1] <- 0}
  mat <- na.omit(mat)

  return(mat)
}

###
recursive_diff <- function(vector, deriv)
{
  vector <- unlist(vector)
  head_value <- vector("numeric", deriv)
  tail_value <- vector("numeric", deriv)
  if(deriv==0){head_value = NULL; tail_value = NULL}
  if(deriv > 0){for(i in 1:deriv){head_value[i] <- head(vector, 1); tail_value[i] <- tail(vector, 1); vector <- diff(vector)}}
  outcome <- list(vector = vector, head_value = head_value, tail_value = tail_value)
  return(outcome)
}

###
invdiff <- function(vector, heads, add = FALSE)
{
  vector <- unlist(vector)
  if(is.null(heads)){return(vector)}
  for(d in length(heads):1){vector <- cumsum(c(heads[d], vector))}
  if(add == FALSE){return(vector[-c(1:length(heads))])} else {return(vector)}
}

###
my_metrics <- function(holdout, forecast, actuals, error_scale = "naive", error_benchmark = "naive", n_class = NULL)
{

  if(is.null(n_class))
  {
    scale <- switch(error_scale, "deviation" = sd(actuals), "naive" = mean(abs(diff(actuals))))
    benchmark <- switch(error_benchmark, "average" = rep(mean(forecast), length(forecast)), "naive" = rep(tail(actuals, 1), length(forecast)))
    me <- ME(holdout, forecast, na.rm = TRUE)
    mae <- MAE(holdout, forecast, na.rm = TRUE)
    mse <- MSE(holdout, forecast, na.rm = TRUE)
    rmsse <- RMSSE(holdout, forecast, scale, na.rm = TRUE)
    mre <- MRE(holdout, forecast, na.rm = TRUE)
    mpe <- MPE(holdout, forecast, na.rm = TRUE)
    mape <- MAPE(holdout, forecast, na.rm = TRUE)
    rmae <- rMAE(holdout, forecast, benchmark, na.rm = TRUE)
    rrmse <- rRMSE(holdout, forecast, benchmark, na.rm = TRUE)
    rame <- rAME(holdout, forecast, benchmark, na.rm = TRUE)
    mase <- MASE(holdout, forecast, scale, na.rm = TRUE)
    smse <- sMSE(holdout, forecast, scale, na.rm = TRUE)
    sce <- sCE(holdout, forecast, scale, na.rm = TRUE)
    gmrae <- GMRAE(holdout, forecast, benchmark, na.rm = TRUE)
    out <- round(c(me = me, mae = mae, mse = mse, rmsse = rmsse, mpe = mpe, mape = mape, rmae = rmae, rrmse = rrmse, rame = rame, mase = mase, smse = smse, sce = sce, gmrae = gmrae), 3)
  }

  if(is.numeric(n_class))
  {
    avg <- suppressMessages(distance(rbind(holdout, forecast), method = "avg"))
    tanimoto <- suppressMessages(distance(rbind(holdout, forecast), method = "tanimoto"))
    hassebrook <- 1 - suppressMessages(distance(rbind(holdout, forecast), method = "hassebrook"))
    jaccard <- suppressMessages(distance(rbind(holdout, forecast), method = "jaccard"))
    taneja <- suppressMessages(distance(rbind(holdout, forecast), method = "taneja"))
    canberra <- suppressMessages(distance(rbind(holdout, forecast), method = "canberra"))
    gower <- suppressMessages(distance(rbind(holdout, forecast), method = "gower"))
    lorentzian <- suppressMessages(distance(rbind(holdout, forecast), method = "lorentzian"))
    clark <- suppressMessages(distance(rbind(holdout, forecast), method = "clark"))

    out <- round(c(avg, tanimoto, hassebrook, jaccard, taneja, canberra, gower, lorentzian, clark), 4)
  }

  return(out)
}

###
prediction_score <- function(integrated_preds, ground_truth)
{
  pfuns <- apply(integrated_preds, 2, ecdf)
  pvalues <- map2_dbl(pfuns, ground_truth, ~ .x(.y))
  scores <- 1 - 2 * abs(pvalues - 0.5)
  return(scores)
}

###

best_deriv <- function(ts, max_diff = 3, thresh = 0.001)
{
  pvalues <- vector(mode = "double", length = as.integer(max_diff))

  for(d in 1:(max_diff + 1))
  {
    model <- lm(ts ~ t, data.frame(ts, t = 1:length(ts)))
    pvalues[d] <- with(summary(model), pf(fstatistic[1], fstatistic[2], fstatistic[3],lower.tail=FALSE))
    ts <- diff(ts)
  }

  best <- tail(cumsum(pvalues < thresh), 1)

  return(best)
}

###

ranker <- function(df, focus, inverse = NULL, absolute = NULL, reverse = FALSE)
{
  rank_set <- df[, focus, drop = FALSE]
  if(!is.null(inverse)){rank_set[, inverse] <- - rank_set[, inverse]}###INVERSION BY COL NAMES
  if(!is.null(absolute)){rank_set[, absolute] <- abs(rank_set[, absolute])}###ABS BY COL NAMES
  index <- apply(scale(rank_set), 1, mean, na.rm = TRUE)
  if(reverse == FALSE){df <- df[order(index),]}
  if(reverse == TRUE){df <- df[order(-index),]}
  return(df)
}

###
ts_graph <- function(x_hist, y_hist, x_forcat, y_forcat, lower = NULL, upper = NULL, line_size = 1.3, label_size = 11,
                     forcat_band = "seagreen2", forcat_line = "seagreen4", hist_line = "gray43", label_x = "Horizon", label_y= "Forecasted Var", dbreak = NULL, date_format = "%b-%d-%Y")
{
  if(is.character(y_hist)){y_hist <- as.factor(y_hist)}
  if(is.character(y_forcat)){y_forcat <- factor(y_forcat, levels = levels(y_hist))}
  if(is.character(lower)){lower <- factor(lower, levels = levels(y_hist))}
  if(is.character(upper)){upper <- factor(upper, levels = levels(y_hist))}

  n_class <- NULL
  if(is.factor(y_hist)){class_levels <- levels(y_hist); n_class <- length(class_levels)}

  all_data <- data.frame(x_all = c(x_hist, x_forcat), y_all = as.numeric(c(y_hist, y_forcat)))
  forcat_data <- data.frame(x_forcat = x_forcat, y_forcat = as.numeric(y_forcat))

  if(!is.null(lower) & !is.null(upper)){forcat_data$lower <- as.numeric(lower); forcat_data$upper <- as.numeric(upper)}

  plot <- ggplot()+ geom_line(data = all_data, aes_string(x = "x_all", y = "y_all"), color = hist_line, size = line_size)
  if(!is.null(lower) & !is.null(upper)){plot <- plot + geom_ribbon(data = forcat_data, aes_string(x = "x_forcat", ymin = "lower", ymax = "upper"), alpha = 0.3, fill = forcat_band)}
  plot <- plot + geom_line(data = forcat_data, aes_string(x = "x_forcat", y = "y_forcat"), color = forcat_line, size = line_size)
  if(!is.null(dbreak)){plot <- plot + scale_x_date(name = paste0("\n", label_x), date_breaks = dbreak, date_labels = date_format)}
  if(is.null(dbreak)){plot <- plot + xlab(label_x)}
  if(is.null(n_class)){plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), labels = number)}
  if(is.numeric(n_class)){plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), breaks = 1:n_class, labels = class_levels)}
  plot <- plot + ylab(label_y)  + theme_bw()
  plot <- plot + theme(axis.text=element_text(size=label_size), axis.title=element_text(size=label_size + 2))

  return(plot)
}

###
sampler <- function(vect, n_samp, range = NULL, integer = FALSE, fun = NULL, multi = NULL)
{
  if(is.null(vect) & is.null(fun))
  {
    if(!is.character(range)){if(integer){set <- min(range):max(range)} else {set <- seq(min(range), max(range), length.out = 1000)}} else {set <- range}
    if(is.null(multi)){samp <- sample(set, n_samp, replace = TRUE)}
    if(is.numeric(multi)){samp <- replicate(n_samp, sample(set, multi, replace = TRUE), simplify = FALSE)}
  }

  if(is.null(vect) & !is.null(fun)){samp <- fun}

  if(is.null(multi)){
    if(length(vect)==1){samp <- rep(vect, n_samp)}
    if(length(vect) > 1){samp <- sample(vect, n_samp, replace = TRUE)}}

  if(is.numeric(multi)){
    if(length(vect)==1){samp <- replicate(n_samp, rep(vect, multi), simplify = FALSE)}
    if(length(vect) > 1){samp <- replicate(n_samp, sample(vect, multi, replace = TRUE), simplify = FALSE)}}

  return(samp)
}

###

plotter <- function(quant_pred, ci, ts, n_class = NULL, level_names = NULL, dates = NULL, feat_name)
{

  seq_len <- nrow(quant_pred)
  n_ts <- length(ts)

  if(is.Date(dates))
  {
    new_dates<- seq.Date(tail(dates, 1), tail(dates, 1) + seq_len * mean(diff(dates)), length.out = seq_len)
    x_hist <- dates
    x_forcat <- new_dates
    rownames(quant_pred) <- as.character(new_dates)
  }
  else
  {
    x_hist <- 1:n_ts
    x_forcat <- (n_ts + 1):(n_ts + seq_len)
    rownames(quant_pred) <- paste0("t", 1:seq_len)
  }

  quant_pred <- as.data.frame(quant_pred)
  x_lab <- paste0("Forecasting Horizon for sequence n = ", seq_len)
  y_lab <- paste0("Forecasting Values for ", feat_name)

  if(is.numeric(n_class) & !is.null(level_names)){ts <- level_names[ts + 1]}
  lower_b <- paste0((1-ci)/2 * 100, "%")
  upper_b <- paste0((ci+(1-ci)/2) * 100, "%")

  plot <- ts_graph(x_hist = x_hist, y_hist = ts, x_forcat = x_forcat, y_forcat = quant_pred[, "50%"], lower = quant_pred[, lower_b], upper = quant_pred[, upper_b], label_x = x_lab, label_y = y_lab)
  return(plot)
}

###
smart_reframer <- function(ts, seq_len, stride)
{
  n_length <- length(ts)
  if(seq_len > n_length | stride > n_length){stop("vector too short for sequence length or stride")}
  if(n_length%%seq_len > 0){ts <- tail(ts, - (n_length%%seq_len))}
  n_length <- length(ts)
  idx <- base::seq(from = 1, to = (n_length - seq_len + 1), by = 1)
  reframed <- t(sapply(idx, function(x) ts[x:(x+seq_len-1)]))
  if(seq_len == 1){reframed <- t(reframed)}
  idx <- rev(base::seq(nrow(reframed), 1, - stride))
  reframed <- reframed[idx,,drop = FALSE]
  colnames(reframed) <- paste0("t", 1:seq_len)
  return(reframed)
}

