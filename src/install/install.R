install.packages("pak", repos = "https://r-lib.github.io/p/pak/stable/")

packages <- readLines("packages.txt")
packages <- packages[nchar(trimws(packages)) > 0]

pak::pkg_install(packages, ask = FALSE)