############################################################
## Some tests on the ordering/sorting of the results.
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

## Compare the performance of doing the sorting within R or
## directly in the SQL query.
dontrun_test_ordering_performance <- function() {
}
