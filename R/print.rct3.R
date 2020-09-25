#' Print an rct3 fit
#'
#' Print an rct3 fit showing the model settings and predicted recruitments
#'
#' @param x an object of class rct3 - an output from the rct3 function.
#' @param digits optional integer for how much to round the values in the
#'        output tables.
#' @param ... additional arguments to print.data.frame
#'
#' @return
#' invisibly returns a summary data frame.
#'
#' @seealso
#'
#' \code{\link{rct3}} run a calibrated regression to predict rectruitment.
#'
#' \code{\link{summary.rct3}} summarise a rct3 object
#'
#' \code{\link{rct3-package}} gives an overview of the package.
#'
#' @examples
#' # load recruitment data
#' data(recdata)
#'
#' formula <-
#'  recruitment ~ NT1 + NT2 + NT3 + NAK1 + NAK2 + NAK3 +
#'                RT1 + RT2 + RT3 + EC01 + ECO2 + ECO3
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
#' @export

print.rct3 <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  hdr <- with(x,
  c("Analysis by RCT3 ver4.0\n",
    stock,
    paste("\nData for ", info[1]," surveys over ", info[2]," year classes : ", info[3]," - ", info[4], sep=""),
    "Regression type = C",
    if (x $ power > 0) {
    paste("Tapered time weighting applied\npower =", power, "over", range, "years")
    } else {
    "Tapered time weighting not applied"
    },
    "Survey weighting not applied",
    if (shrink) {
      "Final estimates ARE shrunk towards mean"
    } else {
      "Final estimates not shrunk towards mean"
    },
    "Estimates with S.E.'S greater than that of mean included",
    paste("Minimum S.E. for any survey taken as   ", min.se),
    "Minimum of   3 points used for regression\n",
    "Forecast/Hindcast variance correction used.\n"
  ))

  cat(paste(hdr, collapse = "\n"), "\n")

  print.data.frame(x $ rct3.summary, digits = digits)
}
