# Check and Call directory argument
suppressPackageStartupMessages(library(optparse))

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package \"optparse\" not installed. Please install it to run this script.", call. = FALSE)
}

option_list <- list(
  make_option(c("-d", "--directory"), type="character", action = "store", default=NULL, 
              help="path of the working directory",metavar="PATH"),
  make_option(c("-m", "--mode"), type="character", action = "store", default=NULL,
              help="Mode of Overrepresentation analysis [accepted mode: Outliers, Enrichment]", metavar="string"),
  make_option(c("-a", "--approximation"), type="character", action = "store", default=NULL,
              help="Define the IQR approximation that is going to be used to defined outliers in 'Outliers' mode [accepted approximation: Standard, Skew].", metavar="string"),
  make_option(c("-n", "--coefficient"), type="float", action = "store", default=1.5,
              help="Define the coefficient for fence definition in 'Outliers' mode [default %default] [range = 1.5 - 3].", metavar="number"),
  make_option(c("-c", "--column"), type="character", action = "store", default=NULL,
              help="column name of the variable in the metadata file to be used.", metavar="string"),
  make_option(c("-v", "--value"), type="character", action = "store", default=NULL,
              help="value from the variable to be compare against the rest.", metavar="string")
)

# Parse the command-line arguments
arguments <- parse_args(OptionParser(option_list = option_list))

# Check for required arguments and flags
if (is.null(arguments$directory) | is.null(arguments$directory)) {
  stop("Error: a mandatory argument was not provided", call.=FALSE)
}

# Check for argument mode and their sequential mandatory variables
if(arguments$mode != "Outliers" & arguments$mode != "Enrichment") {
  stop("Error: Provided mode is not accepted", call.=FALSE)
} else {
  if(arguments$mode == "Outliers") {
    if (is.null(arguments$approximation)) {
      stop("Error: the mandatory argument 'approximation' was not provided for 'Outliers' mode.", call.=FALSE)
    } else {
      if(arguments$approximation != "Standard" & arguments$approximation != "Skew") {
        stop("Error: Provided approximation is not accepted", call.=FALSE)
      }
    }
    if(arguments$coefficient < 1.5 | arguments$coefficient > 3){
      stop("Error: Provided coefficient for 'Otliers' mode is out of accepted range", call.=FALSE)
    }
  } else if(arguments$mode == "Enrichment"){
    if (is.null(arguments$column) | is.null(arguments$value)) {
      stop("Error: a mandatory argument was not provided for 'Enrichment' mode.", call.=FALSE)
    }
  }
}

# Check software installation
suppressPackageStartupMessages(library(syntenet))

if (!requireNamespace("syntenet", quietly = TRUE)) {
   stop("Package \"syntenet\" not installed. Please install it to run this script.", call. = FALSE)
}
