# ExcessILI

The goal for this package is to facilitate the formatting of line list data
from syndromic surveillance datasets into time series and then the analysis of
these data to detect increases above the seasonal baseline. For US data, there
is an option to automatically adjust the data for state-specific flu activity
(using data from [NREVSS](https://www.cdc.gov/surveillance/nrevss/index.html)
and/or state-specific RSV activity (based on Google search volume). The user
can either start with line list data or formatted time series data. An rshiny
app is provided to examine data products.

# Installation

```r
# Currently, ExcessILI is not availble on CRAN
# install.packages("ExcessILI")

# Install development version from GitHub
devtools::install_github("weinbergerlab/ExcessILI")
```

ExcessILI's dependencies require the packages `gdal` and `udunits`, which
cannot be installed using R's package manager. If ExcessILI installation fails
for this reason, you'll need to install these by hand. For installation of
these two packages on OS X, we recommend using a package manager such as
[brew](https://brew.sh).

```bash
brew install gdal udunits
```

# Usage

Try some examples in the vignettes `vignette("YourData")` and 
`vignette("ILINet")`

