a <- commandArgs(trailingOnly = TRUE)
before <- readRDS(a[1]); after <- readRDS(a[2])
stopifnot(identical(names(before), names(after)))
bad <- character(0)
for (nm in names(before)) {
  if (identical(before[[nm]], after[[nm]])) next
  bad <- c(bad, nm)
  cat("### DIFFERS:", nm, "\n")
  b <- before[[nm]]; x <- after[[nm]]
  for (f in union(names(b), names(x))) {
    if (!identical(b[[f]], x[[f]])) {
      cat("   field:", f, "\n")
      cat("     before:", paste(utils::head(as.character(unlist(b[[f]])), 12), collapse=" | "), "\n")
      cat("     after :", paste(utils::head(as.character(unlist(x[[f]])), 12), collapse=" | "), "\n")
    }
  }
}
cat("\n", length(before) - length(bad), "of", length(before), "scenarios identical\n")
if (length(bad)) cat("differing:", paste(bad, collapse=", "), "\n") else cat("GOLDEN MASTER IDENTICAL\n")
