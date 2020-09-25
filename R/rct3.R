#' Run a calibrated regression to predict recruitment
#'
#' Function to run a calibrated regression to predict recruitment using the
#' method decribed by Shepherd (1997)
#'
#' @param formula a formula to define which surveys to use in the recruitment
#'        estimation.
#' @param data a dataframe with one column named 'yearclass' and other columns
#'        with the recruitment and the survey index relavent for that recruitment
#'        value
#' @param predictions which yearclasses to make recruitment predictions for
#' @param shrink shrink predictions to the VPA mean?
#' @param power the power to use 0 - no weighting, 2 - bisquare, 3 - tricubic
#' @param range the year range to use in the time tapered weighting
#' @param min.se the minimum standard error used in the weighting of predictions
#' @param old default TRUE, defines how to treat zero values.  In the origional
#'        implmentation al values were transformed using log(x + 1), old=TRUE
#'        maintains this.
#'
#'
#' @return
#' Object of class \verb{rct3}.
#'
#' @note
#'
#' This function was written based on the publication by Shepherd (1997)
#'  with additional reverse engeneering by comparing results to
#' previous examples run using the RCT3 ver3.1 dos program
#'
#' @seealso
#'
#' \code{\link{rct3-package}} gives an overview of the package.
#'
#' @examples
#' # load recruitment data
#' data(recdata)
#'
#' formula <- recruitment ~ NT1 + NT2 + NT3 +
#'                          NAK1 + NAK2 + NAK3 +
#'                          RT1 + RT2 + RT3 +
#'                          EC01 + ECO2 + ECO3
#'
#' my_rct3 <- rct3(formula, recdata, predictions = 2012:2017, shrink = TRUE)
#'
#' # see a short summary
#' my_rct3
#'
#' # for a full summary do:
#' summary(my_rct3)
#'
#' # the components are here:
#' my_rct3$rct3
#' my_rct3$rct3.summary
#'
#' # predicted recruitment
#' t(my_rct3$rct3.summary["WAP"])
#'
#' @references
#'
#' J. G. Shepherd, Prediction of year-class strength by calibration
#' regression analysis of multiple recruit index series, ICES
#' Journal of Marine Science, Volume 54, Issue 5, October 1997,
#' Pages 741–752, \url{https://doi.org/10.1006/jmsc.1997.0222}
#'
#' @importFrom stats coef lm model.matrix
#'
#' @export

rct3 <- function(formula, data, predictions = NULL, shrink = FALSE,
                 power = 3, range = 20, min.se = 0.2,
                 old = TRUE)
{
  form <- formula[[3]]
  bits <- list()
  while(length(form)>1) {
    bits <- c(bits, form[[3]])
    form <- form[[2]]
  }
  bits <- rev(c(bits, form))
  formulas <- lapply(bits, function(x) {tmp <- formula; tmp[[3]] <- tmp[[2]]; tmp[[2]] <- x; tmp})
  formulas2 <- lapply(bits, function(x) {tmp <- formula; tmp[[3]] <- x; tmp})

  weight <- function(y, y0, D, p) pmax(0, (1 - ((y0 - y)/D)^p)^p)

  log.data <- data
  if (old) {
    log.data[names(data) != "yearclass"] <- log(data[names(data) != "yearclass"] + 1)
  } else # think of something to do with zeros
  {
    log.data[names(data) != "yearclass"] <- log(data[names(data) != "yearclass"])
  }

  # fit one model at a time
  do.one.prediction <- function(i, predict.yr) {
    wk.data <- log.data[log.data$yearclass < predict.yr,]
    yr_diff <- max(wk.data$yearclass) - wk.data$yearclass
    wk.data$wts <- (1 - (pmin(range, yr_diff)/range)^power)^power
    if (nrow(wk.data) < 3) stop("too few data points in one survey!")
    m <- lm(formulas[[i]], data = wk.data, weights = wts)
    b <- {function(x) c(-x[1], 1)/x[2] }(unname(coef(m)))
    wts <- wk.data[names(m$residuals),"wts"]
    rss <- sum( wts * m$residuals^2 )
    mss <- sum(wts * (m$fitted.values - mean(m$fitted.values))^2)
    sigma <- b[2] * sqrt(rss / (sum(wts) - 2))
    rsqr <- mss / (rss + mss)
    n <- m$df.residual + 2

    Xp <- unname(model.matrix(formulas2[[i]][c(1,3)], log.data[log.data$yearclass == predict.yr,]))
    if (nrow(Xp)) {
      X <- unname(model.matrix(formulas2[[i]], wk.data))
      XXinv <- solve(t(X) %*% diag(wts) %*% X)
      pred <- drop(Xp %*% b)
      se.pred <- sqrt(sum(wts) / (sum(wts)-2)) * sigma * sqrt(1 + drop(Xp %*% XXinv %*% t(Xp)))
      index <- Xp[,2]
    } else {
      index <- pred <- se.pred <- NA
    }

    data.frame(index = as.character(formulas[[i]][[2]]),
               slope = b[2], intercept = b[1],
               se = sigma, rsquare = rsqr, n = n,
               indices = index, prediction = pred,
               se.pred = se.pred)
  }

  if (is.null(predictions)) {
    y <- eval(formula[[2]], log.data)
    predictions <- log.data $ yearclass[is.na(y)]
  }

  out <-
    lapply(predictions,
      function(yr)
      {
        out <- do.call( rbind, lapply(1:length(formulas), do.one.prediction, predict.yr = yr))
        # drop years with no indices
        out$slope[is.na(out$indices)] <- NA
        out$intercept[is.na(out$indices)] <- NA
        out$se[is.na(out$indices)] <- NA
        out$rsquare[is.na(out$indices)] <- NA
        out$n[is.na(out$indices)] <- NA
        # calculate vpa contribution
        vpa <- eval(formula[[2]], log.data)[log.data $ yearclass < yr]
        yrs <- log.data$yearclass[log.data$yearclass < yr]
        notNA <- !is.na(vpa)
        wts <- (1 - (pmin(range, max(yrs) - yrs)/range)^power)^power
        vpa <- vpa[notNA]
        yrs <- yrs[notNA]
        wts <- wts[notNA]
        out <- rbind(out, data.frame(index = "VPA Mean",
                                     slope = NA, intercept = NA,
                                     se = NA, rsquare = NA, n = length(vpa),
                                     indices = NA,
                                     prediction = sum(wts * vpa) / sum(wts),
                                     se.pred = sqrt(sum(wts * (vpa - mean(vpa))^2) / (sum(wts)-1))
                                     ))
        if (shrink)
        {
          se.pred <- pmax(out$se.pred, min.se)
          out$WAP.weights <- (1/se.pred^2) / sum(1/se.pred^2, na.rm = TRUE)
        }
        else
        {
          se.pred <- pmax(out[1:(nrow(out)-1),]$se.pred, min.se)
          out $ WAP.weights <- c((1/se.pred^2) / sum(1/se.pred^2, na.rm = TRUE), 0)
        }

        out
      })
  names(out) <- paste("yearclass", predictions, sep=":")

  summarise.rct3 <- function(tmp)
  {
    pred <- with(tmp, sum(prediction * WAP.weights, na.rm = TRUE))
    int.se <- 1/sqrt(sum(1/tmp $ se.pred^2, na.rm = TRUE))

    data.frame(WAP = exp(pred), logWAP = pred, int.se = int.se)
  }

  out <- list(stock = attr(data, "stock"),
              info = c(length(bits), nrow(data), range(log.data $ yearclass)),
              rct3 = out,
              rct3.summary = do.call(rbind, lapply(out, summarise.rct3)),
              shrink = shrink,
              power = power,
              range = range,
              min.se = min.se)

  class(out) <- "rct3"
  out
}
