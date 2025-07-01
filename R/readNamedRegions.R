#' Read all named ranges from a given Excel file
#' 
#' This function reads all named ranges from a given Excel file and returns them
#' as a list of data-frames or vectors length 1, each named after the named range in the Excel file.
#' If a named range is in the Excel file but is empty or missing "Named Range Missing" is returned 
#' for that named range. This is done to avoid error-ing out on irrelevant/legacy named ranges
#' while still giving the user useful information.
#'
#' @param file_path The path to the Excel file
#'
#' @return A list of data-frames or vectors length 1, each named after the named range in the Excel file
# Requires openxlsx R package
readNamedRegions <- function(file_path) {
  
  # check file exists
  if (!file.exists(file_path)) {
    stop("File does not exist")
  }
  
  # Identify all named ranges
  v_named_ranges <-
    openxlsx::getNamedRegions(file_path) |> as.vector()
  
  # Ensure there are named ranges
  if (length(v_named_ranges) == 0) {
    stop("No named ranges found - check Excel file")
  }
  
  # Read all named ranges into a list
  l_input_data <- sapply(
    X = v_named_ranges,
    USE.NAMES = TRUE,
    simplify = FALSE,
    FUN = function(named_range) {
      tryCatch({
        # the named range is loaded as a data-frame, regardless of dimensions
        df_range <- openxlsx::read.xlsx(xlsxFile = file_path, namedRegion = named_range)
        # if the named range is empty (one of the dimensions is zero), return the name only
        if(nrow(df_range) * ncol(df_range) == 0) return(names(df_range))
        # if the named range is not empty, return the data-frame
        return(df_range)
      }, error = function(e) {
        # If the named range can't be found, or an error occurs, return a warning in that named range's place
        return("Named Range Missing")
      })
    }
  )
  
  # create a vector of the missing/empty ones
  v_empty_named_range <- names(which(l_input_data == "Named Range Missing"))
  
  # if there are missing/empty named ranges, print a warning
  if (length(v_empty_named_range) > 0) {
    warning(paste("The following named ranges are empty or missing: ", paste(v_empty_named_range, collapse = ", ")))
  }
  
  return(l_input_data)
  
}
