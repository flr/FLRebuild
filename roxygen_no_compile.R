# Helper script to run roxygen2 without compiling C++ code
# Usage: Rscript roxygen_no_compile.R

# Temporarily move src directory to avoid compilation
if (dir.exists("src")) {
  if (dir.exists("src_backup")) {
    unlink("src_backup", recursive = TRUE)
  }
  file.rename("src", "src_backup")
  on.exit({
    if (dir.exists("src_backup")) {
      if (dir.exists("src")) {
        unlink("src", recursive = TRUE)
      }
      file.rename("src_backup", "src")
    }
  }, add = TRUE)
}

# Run roxygen2
roxygen2::roxygenise()

cat("\n✓ Documentation updated successfully!\n")
cat("✓ NAMESPACE file regenerated\n")
cat("✓ .Rd files generated in man/ directory\n")
