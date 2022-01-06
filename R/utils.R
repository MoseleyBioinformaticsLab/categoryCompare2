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
  path_exists = dir.exists(path)
  
  if (!path_exists) {
    stop(paste0("Provided path ", path, " does not exist!"))
  }
  
  exec_locs <- system.file("exec", package = "categoryCompare2")
  
  if (dir.exists(exec_locs)) {
    exec_files <- dir(path = exec_locs, pattern = "*.R")
    
    file.copy(from = file.path(exec_locs, exec_files), to = file.path(path, exec_files), overwrite = TRUE)
    
    message(cat(file.path(path, exec_files), sep = "\n"))
  } else {
    stop("Exec directory does not exist, and therefore there are no scripts to move!")
  }
  
}

#' executable path
#' 
#' Show the path to the executables, so the user can add them to whatever they
#' want.
#' 
#' @export
#' @return NULL
#' 
executable_path <- function(){
  exec_path <- system.file("exec", package = "categoryCompare2")
  
  message(paste0("To access executables, add this path to your $PATH variable: \n",
                 exec_path))
}
