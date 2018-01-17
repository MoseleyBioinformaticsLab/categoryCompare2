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

#' install executables
#' 
#' move executables to user location, default is ~/bin and changes their
#' permissions to make them executable.
#' 
#' @param path the path to put the executable scripts
#' 
#' @export
#' 
#' @return the listing of the files.
#' 
install_executables <- function(path = "~/bin"){
  stopifnot(dir.exists(path))
  exec_locs <- system.file("executables", package = "categoryCompare2")
  
  if (!is.null(exec_locs)) {
    exec_files <- dir(path = exec_locs, pattern = "*.R")
    
    file.copy(from = file.path(exec_locs, exec_files), to = file.path(path, exec_files), overwrite = TRUE)
    
    Sys.chmod(file.path(path, exec_files), mode = "0750")
    message(cat(file.path(path, exec_files), sep = "\n"))
  } else {
    stop("No scripts to move!")
  }
  
}

#' list executables
#' 
#' Show the path to the executables, so the user can add them to whatever they
#' want.
#' 
#' @export
#' @return NULL
#' 
list_executables <- function(){
  exec_locs <- system.file("executables", package = "categoryCompare2")
  exec_files <- dir(exec_locs, pattern = "*.R", full.names = TRUE)
  
  message(cat(exec_files, sep = "\n"))
}
