install.packages(
    "https://bioconductor.statistik.tu-dortmund.de/packages/3.21/bioc/bin/windows/contrib/4.5/zlibbioc_1.54.0.zip",
    repos = NULL,
    type = "win.binary"
)

Sys.setenv(
    OPEN_BABEL_SRC = "C:/Users/wilhelm/Documents/openbabel3r42/openbabel3r42/src/",
    OPEN_BABEL_BIN = "C:/Users/wilhelm/Documents/openbabel3r42/openbabel3r42/bin/"
)

file.exists(file.path(Sys.getenv("OPEN_BABEL_SRC"), "x64/include/openbabel/obutil.h"))
file.exists(file.path(Sys.getenv("OPEN_BABEL_BIN"), "x64/bin/libopenbabel.a"))

# download + unpack
download.file(
    "https://bioconductor.org/packages/3.22/bioc/src/contrib/ChemmineOB_1.48.0.tar.gz",
    "ChemmineOB_1.48.0.tar.gz", mode = "wb"
)
untar("ChemmineOB_1.48.0.tar.gz")

# edit src/Makevars.win:
# - comment out the line: PKG_LIBS += -lz
# - add ONE of the following lines:

# (A) safest: static zlibbioc (no runtime DLL hunt)
# PKG_LIBS += $(shell "${R_HOME}/bin${R_ARCH_BIN}/Rscript" -e "zlibbioc::pkgconfig('PKG_LIBS_static')")

# (B) shared zlibbioc (Windows-recommended in zlibbioc doc)
# PKG_LIBS += $(shell "${R_HOME}/bin${R_ARCH_BIN}/Rscript" -e "zlibbioc::pkgconfig('PKG_LIBS_shared')")

system("R CMD build ChemmineOB")
install.packages("ChemmineOB_1.48.0.tar.gz", repos = NULL, type = "source")

options(repos = BiocManager::repositories())
BiocManager::install("ChemmineOB", type="source", force=TRUE, ask=FALSE, update=FALSE)
# 