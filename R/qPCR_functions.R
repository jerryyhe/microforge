# qPCR_functions.R

# FUNCTIONS

#' Generate a Standard Curve
#'
#'
#'
#' @param data An input dataframe containing at least a column called "Cq" and a column with copies called "Copies"
#' @return A list containing key parameters of the standard curve and a plot of the curve.
#' @export
calculate_std_curve <- function(data) {
  # Take input dataframe containing at least a column called "Cq" and a column with copies called "Copies"

  df <- data %>%
    mutate(log10_copies = log10(Copies))

  # Create visual plot
  p1 <- df %>%
    ggplot(aes(x = Cq, y = log10_copies)) +
    geom_point(colour = "steelblue") +
    geom_smooth(method = "lm", formula = y ~ x)

  # Calculate linear regression
  fit <- lm(Cq ~ log10_copies, data = df)
  slope <- coef(fit)[["log10_copies"]]
  intercept <- coef(fit)[["(Intercept)"]]
  efficiency <- 10^(-1/slope) - 1
  r2 <- summary(fit)$r.squared

  # Package into a list
  result <- list(
    plot = ggplotly(p1),
    model = fit,
    slope = slope,
    intercept = intercept,
    efficiency = efficiency,
    r2 = r2
  )

  return(result)
}

#' Generate a Standard Curve
#'
#' Take an input dataframe that has unknown samples and Cq values and uses std_list input (direct output of `calculate_std_curve`) to convert Cq values to Copy Number
#'
#' @param data An input dataframe that has unknown samples and Cq values
#' @param std_list A standard curve list, direct output of `calculate_std_curve`, to convert Cq values to Copy Number
#' @return A dataframe with a new column called `log10_copies` and `copies_rxn` (copies per reaction)
#' @export
calculate_copies <- function(data, std_list) {
  # Take an input dataframe that has unknown samples and Cq values and uses std_list input (direct output of `calculate_std_curve`) to convert Cq values to Copy Number
  intercept <- std_list$intercept
  slope <- std_list$slope

  df <- data %>%
    mutate(
      log10_copies = (Cq - intercept)/slope,
      copies_rxn = 10^log10_copies
    )

  return(df)

}
