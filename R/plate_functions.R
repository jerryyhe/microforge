# plate_functions.R


# FUNCTIONS

#' Turn a 384-Well Plate Into a Tidy Dataframe
#'
#' Turns a 384-well plate where replicates are in quadrants (e.g. A1, A2, B1, B2) into a tidy dataframe where each quadrant has a unique ID
#'
#' @param data An input dataframe. Columns are numbers (1-25). The first column are the plate row names (A-P).
#' @param plate_id A unique identifier for a plate.
#' @return A tidy dataframe with a unique plate and quadrant identifier. Can be used to join metadata onto.
#' @export
tidy_plate384 <- function(data, plate_id = "P1") {

  df <- data %>%
    select(-`...26`) %>%
    dplyr::rename(row = 1) %>%
    mutate(across(
      where(is.character) & !any_of("row"),
      ~as.double(if_else(.x == "OVRFLW", "4", .x)))
    ) %>%
    pivot_longer(
      cols = 2:25,
      names_to = "column",
      values_to = "AU"
    ) %>%
    mutate(
      plate_id = plate_id,
      column = as.double(column),
      well = sprintf("%s%02d", row, column),
      row_idx = match(row, LETTERS),                                     # maps a number to each letter. LETTERS is an environmental object

      quad_row = ceiling(row_idx / 2),                                   # Rounds the row index up to nearest integer
      quad_col = ceiling(column / 2),                                    # Rounds the column up to nearest integer
      row96_letter = LETTERS[quad_row],                                  # Essentially, what wells this quadrant maps to on a 96-well plate

      quadrant_id = sprintf("%s_Q%02d_%02d", plate_id, quad_row, quad_col),
      single_well_id = sprintf("%s_%s%02d", plate_id, row96_letter, quad_col)
    )

  return(df)

}

#' Turn a 96-Well Plate Into a Tidy Dataframe
#'
#'
#'
#' @param data An input dataframe. Columns are numbers (1-24). The first column are the plate row names (A-P).
#' @param plate_id A unique identifier for a plate.
#' @return A tidy dataframe with a unique plate and quadrant identifier. Can be used to join metadata onto.
#' @export
tidy_plate96 <- function(data, plate_id = "P1") {

  df <- data %>%
    select(-`...14`) %>%
    dplyr::rename(row = 1) %>%
    pivot_longer(
      -row,
      names_to = "column",
      values_to = "AU"
    ) %>%
    mutate(
      plate_id = plate_id,
      column = as.double(column),
      well = sprintf("%s%02d", row, column),
      row_idx = match(row, LETTERS),
      single_well_id = sprintf("%s_%s", plate_id, well)
    )
}


