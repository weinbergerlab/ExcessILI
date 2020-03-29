library(cdcfluview)

test_that("ILINet data is fetched and plotted correctly", {

  ili.data <- ilinet(region = c("state"))
  ili.data$state <- state.abb[match(ili.data$region, state.name)]

  ili.data <- ili.data[, c("state", "week_start", "ilitotal", "total_patients")]
  ili.data <- ili.data[ili.data$state %in% c("CA", "NY", "WA", "NJ", "CT"), ]

  excess_cases1 <-
    excessCases(ds = ili.data,
                datevar = "week_start",
                agevar = "none",
                statevar = "state",
                denom.var = "total_patients",
                use.syndromes = c("ilitotal"),
                time.res = "week")

  # Use package 'shinytest' to test this. Not implemented yet!
  # dashboardPlot(excess_cases1)

  unexplained.cases <-
    excessExtract(ds = excess_cases1,
                  syndrome = "ilitotal",
                  extract.quantity = "unexplained.cases")

  excess.rr <-
    excessExtract(ds = excess_cases1,
                  syndrome = "ilitotal",
                  extract.quantity = "resid1")

  expect_true(!is.null(excess.rr))
  expect_true(!is.null(unexplained.cases))

  par(mfrow = c(1, 1))

  expect_silent(matplot(exp(excess.rr[, , 1]), type = "l"))
})
