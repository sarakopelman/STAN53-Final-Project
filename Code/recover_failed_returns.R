library(quantmod)

# List of tickers to recover
failed_tickers <- c("VSURE.ST", "PLSVO.ST", "ENVIPO.ST", "NORSEO.ST", "CAPSLO.ST", "CAVENO.ST")

# Safe download function with fallback start date
safe_get <- function(tkr) {
  # Try full history from 2010
  x <- tryCatch(getSymbols(tkr, src = "yahoo", from = "2010-01-01", auto.assign = FALSE),
                error = function(e) NULL)
  if(!is.null(x)) return(list(data = x, start_used = "2010-01-01"))
  
  # If failed, try last 2 years
  recent_start <- Sys.Date() - 365*2
  x <- tryCatch(getSymbols(tkr, src = "yahoo", from = recent_start, auto.assign = FALSE),
                error = function(e) NULL)
  if(!is.null(x)) return(list(data = x, start_used = as.character(recent_start)))
  
  # Still failed
  return(list(data = NULL, start_used = NA))
}

# Download the tickers
recovered_data <- lapply(failed_tickers, safe_get)
names(recovered_data) <- failed_tickers

# Separate recovered vs still missing
recovered <- names(recovered_data)[sapply(recovered_data, function(x) !is.null(x$data))]
still_missing <- names(recovered_data)[sapply(recovered_data, function(x) is.null(x$data))]

cat("Recovered tickers:\n")
print(recovered)
cat("Still missing tickers:\n")
print(still_missing)

# Optional: compute daily log returns for recovered tickers
if(length(recovered) > 0){
  returns_list <- lapply(recovered_data[recovered], function(x) dailyReturn(Ad(x$data), type = "log"))
  returns <- Reduce(function(x,y) merge(x,y, all = TRUE), returns_list)
  colnames(returns) <- recovered
  
  # Save results separately
  write.csv(as.data.frame(returns), "recovered_failed_tickers_returns.csv", row.names = TRUE)
  saveRDS(returns, "recovered_failed_tickers_returns.rds")
  
  # Also print first available date for each recovered ticker
  first_dates <- sapply(recovered_data[recovered], function(x) index(x$data)[1])
  print(first_dates)
  
  cat("Done — recovered returns saved.\n")
}
