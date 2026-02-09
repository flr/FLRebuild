# Script to fix FLCandy Collate field issue
# This removes or fixes the Collate field in FLCandy's DESCRIPTION file

fix_flcandy_collate <- function(flcandy_path = NULL) {
  # Try to find FLCandy
  if (is.null(flcandy_path)) {
    # Common locations
    possible_paths <- c(
      "../FLCandy",
      "../../FLCandy",
      "FLCandy",
      file.path(Sys.getenv("HOME"), "FLCandy"),
      file.path(Sys.getenv("USERPROFILE"), "FLCandy")
    )
    
    for (path in possible_paths) {
      if (file.exists(file.path(path, "DESCRIPTION"))) {
        flcandy_path <- path
        break
      }
    }
  }
  
  if (is.null(flcandy_path) || !file.exists(file.path(flcandy_path, "DESCRIPTION"))) {
    stop("Could not find FLCandy package. Please specify the path to FLCandy.")
  }
  
  desc_file <- file.path(flcandy_path, "DESCRIPTION")
  desc_lines <- readLines(desc_file)
  
  # Find Collate field
  collate_start <- grep("^Collate:", desc_lines)
  
  if (length(collate_start) == 0) {
    cat("No Collate field found in FLCandy DESCRIPTION. Nothing to fix.\n")
    return(invisible(NULL))
  }
  
  # Find where Collate field ends (next field or end of file)
  collate_end <- collate_start
  for (i in (collate_start + 1):length(desc_lines)) {
    if (grepl("^[A-Za-z]+:", desc_lines[i]) || desc_lines[i] == "") {
      collate_end <- i - 1
      break
    }
    collate_end <- i
  }
  
  cat("Found Collate field at lines", collate_start, "to", collate_end, "\n")
  cat("Removing Collate field...\n")
  
  # Remove Collate field
  new_desc <- desc_lines[-seq(collate_start, collate_end)]
  
  # Write back
  writeLines(new_desc, desc_file)
  
  cat("✓ Fixed FLCandy DESCRIPTION file\n")
  cat("  Removed Collate field referencing missing abi.R\n")
  cat("  You can now install FLCandy or FLRebuild\n")
  
  return(invisible(NULL))
}

# Run if called directly
if (!interactive()) {
  fix_flcandy_collate()
}
