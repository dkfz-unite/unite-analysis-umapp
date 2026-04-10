install.packages("pak", repos = "https://r-lib.github.io/p/pak/stable/")

if (!file.exists("packages.txt")) stop("packages.txt not found")
packages <- readLines("packages.txt")
packages <- packages[nchar(trimws(packages)) > 0]

pak::pkg_install(packages, ask = FALSE)