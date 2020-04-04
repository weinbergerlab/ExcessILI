# ExcessILI

ExcessILI facilitates formatting line list data from syndromic surveillance
datasets into time series, and enables analysis of these data to detect
increases in reporting above the seasonal baseline. For US data, there is an
option to automatically adjust the data for state-specific flu activity (using
data from [NREVSS](https://www.cdc.gov/surveillance/nrevss/index.html) and/or
state-specific RSV activity (based on Google search volume). The user can
either start with line list data, or formatted time series data. An
[rshiny](https://shiny.rstudio.com/) app is provided to examine data products.

# Installation

```r
# Currently, ExcessILI is not availble on CRAN
# install.packages("ExcessILI")

# Install development version from GitHub. This requires that the 'devtools'
# package be installed.
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("weinbergerlab/ExcessILI")
```

# Usage

- To see how to go from line list data to an interactive app, see
  `vignette("ESSENCE")`

- For an example of how to use the analysis functions with CDC ILINet
  data, see `vignette("ILINet")`.

# Resources

- [Ask a question](mailto:daniel.weinberger@yale.edu?subject=ExcessILI)
- [Open an issue](https://github.com/weinbergerlab/ExcessILI/issues) (GitHub
  issues for bug reports, feature requests)
