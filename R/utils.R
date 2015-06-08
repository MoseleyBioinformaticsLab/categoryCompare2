#' choose apply
#'
#' based on whether the parallel package is loaded and options(mc.cores) is not 1,
#' passes data on to the base lapply or the mclapply function, and returns
#' the results.
#'
#' @param x the list object to work on
#' @param fun the function to be called
#' @param ... other parameters
#' @export
#' @return list results
have_parallel <- function(){
  has_parallel <- FALSE
  
  loaded_packages <- loadedNamespaces()
  num_cores <- getOption("mc.cores")

  if (is.null(num_cores)){
    num_cores <- 1
  }

  if (("parallel" %in% loaded_packages) && (num_cores > 1)){
    has_parallel <- TRUE
  } 

  has_parallel
}

