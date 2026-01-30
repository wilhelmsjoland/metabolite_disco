# 0) Fresh R session recommended

# 1) Bioconductor repos
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
options(repos = BiocManager::repositories())

# 2) Confirm build toolchain (Rtools)
if (!requireNamespace("pkgbuild", quietly = TRUE)) install.packages("pkgbuild")
if (!pkgbuild::has_build_tools(debug = TRUE)) stop("Rtools toolchain not found/working.")

# 3) OpenBabel bundle root (you renamed it to this)
OB_ROOT <- "C:/openbabel3"

# 4) Helper path converters
msys_path <- function(p) {
    p <- normalizePath(p, winslash = "/", mustWork = TRUE)
    if (grepl("^[A-Za-z]:", p)) p <- paste0(tolower(substr(p, 1, 1)), substr(p, 2, nchar(p)))
    p
}
win_path <- function(p) normalizePath(p, winslash = "\\", mustWork = TRUE)

# 5) Bundle layout (matches what your compiler output is using)
OB_INC <- file.path(OB_ROOT, "src/x64/include")     # has openbabel/obutil.h
OB_LIB <- file.path(OB_ROOT, "bin/x64/bin")         # has libopenbabel.a
OB_EIG <- file.path(OB_ROOT, "deps/eigen-3.4.0")    # has Eigen/Core

stopifnot(file.exists(file.path(OB_INC, "openbabel/obutil.h")))
stopifnot(file.exists(file.path(OB_LIB, "libopenbabel.a")))
stopifnot(file.exists(file.path(OB_EIG, "Eigen/Core")))

# 6) Install zlibbioc (needed because your OpenBabel static lib expects bioc_* zlib symbols)
#    Use the R-4.5 Windows binary from the Bioconductor 3.21 archive. :contentReference[oaicite:1]{index=1}
ZLIB_ZIP <- "https://bioconductor.statistik.tu-dortmund.de/packages/3.21/bioc/bin/windows/contrib/4.5/zlibbioc_1.54.0.zip"
install.packages(ZLIB_ZIP, repos = NULL, type = "win.binary")
stopifnot(requireNamespace("zlibbioc", quietly = TRUE))

# 7) Export env vars for ChemmineOB configure/Makevars
Sys.setenv(
    OPEN_BABEL_INCDIR = msys_path(OB_INC),
    OPEN_BABEL_LIBDIR = msys_path(OB_LIB),
    EIGEN_DIR         = msys_path(OB_EIG),
    PATH = paste(win_path(OB_LIB), Sys.getenv("PATH"), sep = ";")
)

# 8) Download ChemmineOB source and patch Makevars to link against zlibbioc (libzbioc) instead of -lz
BIOC_VER <- as.character(BiocManager::version())
CHEM_VER <- "1.48.0"  # match your install log; change only if your log shows a different version
CHEM_URL <- sprintf("https://bioconductor.org/packages/%s/bioc/src/contrib/ChemmineOB_%s.tar.gz", BIOC_VER, CHEM_VER)

td  <- file.path(tempdir(), "chemmineob_build")
dir.create(td, showWarnings = FALSE, recursive = TRUE)
tgz <- file.path(td, sprintf("ChemmineOB_%s.tar.gz", CHEM_VER))
download.file(CHEM_URL, tgz, mode = "wb", quiet = FALSE)
untar(tgz, exdir = td)
PKGDIR <- file.path(td, "ChemmineOB")
stopifnot(dir.exists(PKGDIR))

# 9) Patch all Makevars* in src:
#    - remove '-lz' (wrong zlib)
#    - add zlibbioc pkgconfig flags (recommended on Windows) :contentReference[oaicite:2]{index=2}
makevars_files <- list.files(file.path(PKGDIR, "src"), pattern = "^Makevars", full.names = TRUE)
stopifnot(length(makevars_files) > 0)

for (f in makevars_files) {
    x <- readLines(f, warn = FALSE)
    
    # remove system zlib
    x <- gsub("\\s+-lz(\\s|$)", " ", x)
    
    # add zlibbioc include + libs
    add_cpp <- 'PKG_CPPFLAGS += $(shell echo \'zlibbioc::pkgconfig("PKG_CFLAGS")\' | "${R_HOME}/bin/R" --vanilla --slave)'
    add_lib <- 'PKG_LIBS     += $(shell echo \'zlibbioc::pkgconfig("PKG_LIBS_shared")\' | "${R_HOME}/bin/R" --vanilla --slave)'
    
    if (!any(grepl('zlibbioc::pkgconfig\\("PKG_CFLAGS"\\)', x))) x <- c(x, add_cpp)
    if (!any(grepl('zlibbioc::pkgconfig\\("PKG_LIBS_shared"\\)', x))) x <- c(x, add_lib)
    
    writeLines(x, f, useBytes = TRUE)
}

# 10) Install from patched source
RCMD <- file.path(R.home("bin"), "Rcmd.exe")
system2(RCMD, c("INSTALL", "--preclean", shQuote(PKGDIR)))

# 11) Smoke test
library(ChemmineOB)
packageVersion("ChemmineOB")
