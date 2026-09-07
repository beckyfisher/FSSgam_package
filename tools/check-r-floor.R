# Checks that DESCRIPTION's declared "Depends: R (>= x)" is reachable given what
# this package's hard dependencies require, transitively.
#
# This is the direction R CMD check cannot see. Installing a package rejects a
# floor set ABOVE the running R and accepts any floor BELOW it, so no matrix of
# R versions detects a floor that is too low -- which is the whole of
# FSSgam_package#31, where "R (>= 3.5)" was declared while mgcv and MuMIn both
# required R (>= 4.4.0). Reverting to 3.5 leaves every R CMD check green.
#
# The closure is walked transitively, not just over direct Imports. Matrix
# declares R (>= 4.4) and is not a direct dependency of this package: it arrives
# through mgcv and MuMIn. A check of the direct Imports alone therefore counts
# one fewer package at the binding floor, and would report "reachable" for a
# floor that is not, for a package whose only route to Matrix is transitive.
#
# An uninstalled dependency is an error, not a skip. Silently contributing no
# floor is the same failure this script exists to catch, one level down, and
# CLAUDE.md section 3 records that this repository's environments do sometimes
# lack a declared dependency.

# Returns list(status, floor, raw), where status is one of:
#   "none"     the field declares no R constraint at all
#   "ok"       an R constraint that was parsed, floor is a package_version
#   "unparsed" an R constraint in a form this does not read
#
# The three are distinguished deliberately. Returning NULL for both "no R
# constraint" and "an R constraint I could not read" is how a check comes to
# report success without having checked, which is the defect this script exists
# to catch: an unread constraint would print as "-", the same glyph as a package
# that genuinely declares nothing, and contribute no floor.
parse_r_constraint <- function(depends) {
  none <- list(status = "none", floor = NULL, strict = NA, raw = NA_character_)
  if (is.null(depends) || length(depends) == 0 || all(is.na(depends))) return(none)
  fields <- trimws(strsplit(depends, ",")[[1]])
  # Anchored on R immediately followed by "(", so a package named Rcpp or R6 is
  # not mistaken for the R entry. The space before "(" is optional and is
  # written both ways in practice.
  r <- grep("^R[[:space:]]*\\(", fields, value = TRUE)
  if (length(r) == 0) return(none)
  raw <- r[1]

  # ">" as well as ">=", and the two are kept apart rather than conflated.
  # "R (> 4.4.0)" is declared by real packages and does not mean the same as
  # "R (>= 4.4.0)": a declared floor of exactly 4.4.0 satisfies the second and
  # not the first. Reading one as the other would understate the requirement,
  # so the operator is recorded and each constraint is tested as written. Any
  # other operator is reported rather than guessed at.
  m <- regmatches(raw, regexec("^R[[:space:]]*\\([[:space:]]*(>=?)[[:space:]]*([^)]*)\\)", raw))[[1]]
  if (length(m) != 3) return(list(status = "unparsed", floor = NULL, strict = NA, raw = raw))
  v <- trimws(m[3])
  if (!grepl("^[0-9]+([.-][0-9]+)*$", v)) return(list(status = "unparsed", floor = NULL, strict = NA, raw = raw))
  list(status = "ok", floor = package_version(v), strict = identical(m[2], ">"), raw = raw)
}

# Package names from a Depends/Imports/LinkingTo field, version constraints and
# the R entry removed.
dep_names <- function(field) {
  if (is.null(field) || length(field) == 0 || all(is.na(field))) return(character(0))
  x <- trimws(strsplit(gsub("[[:space:]]+", " ", field), ",")[[1]])
  x <- trimws(sub("\\(.*$", "", x))
  x <- x[nzchar(x)]
  setdiff(x, "R")
}

# Base-priority packages ship with R itself and impose no floor of their own.
# Recommended packages (Matrix, MASS, nlme, mgcv, nnet) are NOT excluded: they
# are versioned independently and are where the binding floor is found here.
base_pkgs <- rownames(utils::installed.packages(priority = "base"))

# read.dcf() returns only the fields present, so index by name through a helper
# rather than directly: a package with no LinkingTo has no such column.
desc_raw <- read.dcf("DESCRIPTION")[1, ]
field <- function(name) if (name %in% names(desc_raw)) desc_raw[[name]] else NULL

declared_c <- parse_r_constraint(field("Depends"))
if (declared_c$status == "unparsed") {
  stop(sprintf("DESCRIPTION declares \"%s\", which this script does not read. Write it as R (>= x).",
               declared_c$raw), call. = FALSE)
}
if (declared_c$status == "none") stop("DESCRIPTION declares no Depends: R (>= x) floor.", call. = FALSE)
declared <- declared_c$floor

# Hard dependencies are Depends, Imports and LinkingTo. Suggests is excluded:
# a Suggests floor does not stop anyone installing the package, so it does not
# bound this field. Depends is included -- moving a package from Imports to
# Depends must not disable the check.
seed <- unique(c(dep_names(field("Depends")), dep_names(field("Imports")),
                 dep_names(field("LinkingTo"))))
seed <- setdiff(seed, base_pkgs)

seen <- character(0)
rows <- list()
queue <- seed
missing <- character(0)

while (length(queue) > 0) {
  p <- queue[1]
  queue <- queue[-1]
  if (p %in% seen) next
  seen <- c(seen, p)

  d <- tryCatch(utils::packageDescription(p), error = function(e) NA)
  if (length(d) == 1 && is.na(d)) { missing <- c(missing, p); next }

  con <- parse_r_constraint(d$Depends)
  rows[[p]] <- list(package = p,
                    version = as.character(d$Version),
                    direct = p %in% seed,
                    status = con$status,
                    raw = con$raw,
                    strict = con$strict,
                    floor = con$floor)
  nxt <- unique(c(dep_names(d$Depends), dep_names(d$Imports), dep_names(d$LinkingTo)))
  queue <- c(queue, setdiff(nxt, c(seen, base_pkgs)))
}

cat(sprintf("DESCRIPTION declares R (>= %s)\n\n", declared))

# Every guard runs before any exit path. Ordering matters and has been wrong:
# an empty-closure exit placed above the uninstalled-dependency check accepted
# the declared floor, having read nothing, for a package all of whose
# dependencies were missing -- reporting "no hard dependencies outside base R"
# when it declared two. That is the failure this script exists to catch,
# occurring inside it. Add new guards above the reporting, never below it.
if (length(missing) > 0) {
  stop(sprintf(paste0("Cannot check the declared floor: %s not installed, so the R version ",
                      "they declare could not be read. Install the package's hard ",
                      "dependencies and re-run."),
               paste(sort(unique(missing)), collapse = ", ")), call. = FALSE)
}

unread <- Filter(function(x) x$status == "unparsed", rows)
if (length(unread) > 0) {
  stop(sprintf(paste0("Cannot check the declared floor: the R constraint of %s is in a form ",
                      "this script does not read (%s). Extend parse_r_constraint() rather ",
                      "than letting the constraint go uncounted."),
               paste(vapply(unread, `[[`, character(1), "package"), collapse = ", "),
               paste(vapply(unread, `[[`, character(1), "raw"), collapse = "; ")), call. = FALSE)
}

# Every field is coerced to a length-one string before printing. sprintf()
# returns character(0) if any argument is zero-length, so a package whose
# DESCRIPTION omits Version would otherwise drop its whole row from the table
# rather than printing a blank.
one <- function(x, blank = "?") if (length(x) == 0 || is.na(x[1])) blank else as.character(x[1])

if (length(rows) > 0) {
  cat(sprintf("%-14s %-10s %-8s %s\n", "package", "version", "direct", "declared R floor"))
  for (x in rows[order(names(rows))]) {
    shown <- switch(x$status,
                    ok = paste0(if (isTRUE(x$strict)) "> " else ">= ", one(as.character(x$floor))),
                    none = "-",
                    unparsed = paste0("UNREAD: ", one(x$raw)))
    cat(sprintf("%-14s %-10s %-8s %s\n", one(x$package), one(x$version),
                if (isTRUE(x$direct)) "yes" else "-", shown))
  }
} else {
  cat("This package declares no hard dependencies outside base R.\n")
}

floors <- Filter(Negate(is.null), lapply(rows, `[[`, "floor"))

if (length(floors) == 0) {
  cat("\nNo hard dependency declares an R floor, so nothing bounds the declared one\n")
  cat("from below. Something other than a dependency is imposing it. Record what,\n")
  cat("so that the next person does not lower it.\n")
  quit(status = 0)
}

# package_version compares componentwise and pads, so "4.4" and "4.4.0" compare
# equal and "4.10.0" is above "4.9.0". A numeric coercion would order the second
# pair wrongly, and "4.4.0" is not a number at all.
needed <- floors[[1]]
for (v in floors) if (v > needed) needed <- v

binding <- names(floors)[vapply(floors, function(v) v == needed, logical(1))]
cat(sprintf("\n%d hard dependencies, %d of them direct. Highest floor: R (>= %s), from %s\n",
            length(rows), sum(vapply(rows, function(x) isTRUE(x$direct), logical(1))),
            needed, paste(sort(binding), collapse = ", ")))

# Each constraint is tested as written rather than against a single computed
# maximum, so that "> x" and ">= x" are not conflated: a declared floor of
# exactly x satisfies the second and fails the first.
satisfied <- vapply(rows, function(x) {
  if (x$status != "ok") return(TRUE)
  if (isTRUE(x$strict)) declared > x$floor else declared >= x$floor
}, logical(1))

if (any(!satisfied)) {
  unmet <- rows[!satisfied]
  stop(sprintf(paste0("DESCRIPTION declares R (>= %s), which does not satisfy %s. The ",
                      "declared floor is not reachable: raise it, or drop the dependency ",
                      "that imposes the higher one."),
               declared,
               paste(vapply(unmet, function(x)
                 sprintf("%s's %s", x$package, x$raw), character(1)), collapse = ", ")),
       call. = FALSE)
}
if (declared > needed) {
  cat("\nThe declared floor is above what the dependencies require. That is allowed,\n")
  cat("but something other than a dependency is imposing it. Record what, so that the\n")
  cat("next person does not lower it.\n")
}
cat("\nOK: the declared floor is reachable.\n")
