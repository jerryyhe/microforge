# deseq2_functions

#' Generate a list of files to load
#'
#'
#'
#' @param directory A directory containing all files to load.
#' @param files A string vector of file names to load.
#' @return A list of files to load.
#' @export
load_files <- function(directory, files) {
  # For the following chunk of code we are leaving two options for the user to use: (1) The user can input a directory where all the files are
  # found and there will be a function that will load them into a string vector, (2) The user can input a string vector of the the relative paths
  # where the files they want to load are found

  if(missing(directory)) { # missing() checks if the directory parameter is given by the use or not, if it is not given it assumes the 'files=' argument was given
    strain_files <- files
  } else {
    strain_files <- list.files(path = directory, pattern = "*.tsv", full.names = TRUE, recursive = FALSE)
    # list.files load the relative paths of files in a directory into a string vector!
  }
  return(strain_files)
}

#' Load files into R
#'
#'
#'
#' @param strain_files Output from `load_files()`. A string vector of relative paths to files.
#' @return A dataframe.
#' @export
files_to_array <- function(strain_files) {
  # input: a vector of strings that contains the relative paths of the .tsv files
  # output: a dataframe that contains categorical variables (System, subsystem, role, etc) of RAST genes and then columns
  #         with the name of each strain (Which is extracted from the basename of the file) that will have '1' if they contain
  #         the genes in a row or '0' if they do not contain that gene

  df_exists <- FALSE # This will allow us to initiate the dataframe with the first file and add the following files to the existing df

  # Dont take into consideration the ABSOLUTE feature count -> some genes have more than 1 version inside a genome
  for(f in strain_files){
    if(df_exists == FALSE){ # This conditions is met for the first file as there is no existing data frame (df)
      df <- read_tsv(file = f, col_names = c("x", "count", "Function")) # load f which is the object from the iteration
      df_formated <- format_df(df, f) # See helper function below
      df_exists <- TRUE # change to true so the condition is not met in second iteration
    }
    else{
      df <- read_tsv(file = f, col_names = c("x", "count", "Function"))
      df_to_add <- format_df(df, f)
      df_formated <- full_join(df_formated, df_to_add) # adds new dataframe to old dataframe
    }
  }

  df_formated <- as.data.frame(mutate_all(df_formated, ~replace(., is.na(.), 0))) # Change NAs to 0 and set tibble as data frame
  return(df_formated)
}

#' Format dataframes to be joined
#'
#'
#'
#' @param df A dataframe output from `files_to_array()`.
#' @param file An input file?
#' @return A dataframe.
#' @export
format_df <- function(df, file){
  # purpose: This is a helper function that formats the dataframe from the file loaded using read_tsv()
  #          such that other df can be joined to it to create an array
  # Input: The input is a dataframe from the loaded file and the file path as a string
  # output: a dataframe with the 'Features' column removed and a column with the basename of the file (without extension) with the number 1
  #         per row, which signifies that the specific gene is found for that strain

  file_name <- tools::file_path_sans_ext(basename(file)) # This gives the file name without extension or file path
  df[(file_name)] <- df %>% select(count) %>% pull # add the new column with the extracted 'file_name' as its tittle
  df <- df %>% select(-x, -count)
  return(df)
}

#' Add labels to the DESeq2 results dataframe
#'
#'
#'
#' @param dat A DESeq2 results dataframe.
#' @return A dataframe with labels.
#' @export
add_labels <- function(dat){

  dat_l2fc <- dat %>%
    rownames_to_column("Function") %>%
    arrange(desc(log2FoldChange)) %>%
    dplyr::select(Function, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
    mutate(abs_l2fc = abs(log2FoldChange), up_down = if_else(log2FoldChange <= 0, "down", "up")) %>%
    mutate(l2fc = case_when(abs_l2fc >=4 ~ ">4",
                            abs_l2fc >=3 ~ "3",
                            abs_l2fc >=2 ~ "2",
                            abs_l2fc >=1 ~ "1",
                            TRUE ~ "<1"
    )) %>%
    mutate(l2fc = factor(l2fc))

  return(dat_l2fc)

}

#' Run DESeq2 and Export to and RDS File or Load DESeqDataSet Object from an Existing RDS File
#'
#'
#'
#' @param rdsName The desired name of your RDS or existing RDS file.
#' @param directory The directory to write or load your RDS file.
#' @return A dataframe.
#' @export
run_deseq2 <- function(rdsName, directory) {
  # Define name of DESeqDataSet RDS file
  rdsName <- rdsName
  ddsName <- paste0(directory, rdsName)

  # Check if an existing dds file exists, if so then load it. Otherwise, perform DESeq on the constructed DDS object
  if (file.exists(ddsName)) {

    # Load an existing DDS object
    dds <- readRDS(ddsName)
    res <- results(dds)
    print("Existing DESeqDataSet object loaded.")

  } else {

    # Perform DESeq2 on the constructed DDS object
    dds <- DESeq(dds)
    res <- results(dds)

    # Save DDS object to an RDS file
    saveRDS(dds, file = ddsName)
    print("Saving DDS object to an RDS file...")

  }

  ddsResult <- list(dds = dds, res = res)
  return(ddsResult)
}


