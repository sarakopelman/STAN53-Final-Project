##########################################################
# R script: Read ticker Excel file, clean, and download prices
# - Handles tickers that IPO after 2010 automatically
# - Saves under new filenames
##########################################################

library(readxl)
library(dplyr)
library(stringr)
library(quantmod)

# 1. Read the Excel file
tickers_df <- read_excel("Data/tickers.xlsx")  # adjust path if needed
head(tickers_df)

# Assuming the tickers are in a column named "symbol"
tickers <- tickers_df$symbol

# 2. Clean tickers for Yahoo Finance
yahoo_tickers <- tickers %>%
  str_replace_all("\\.", "-") %>%  # SEB.B -> SEB-B
  paste0(".ST")                     # add Yahoo .ST suffix

# 3. Batch download historical prices
batch_size <- 20
start_date <- "2010-01-01"

safe_get <- function(tkr) {
  tryCatch({
    # Yahoo will return all available data from earliest date if later IPO
    getSymbols(tkr, src = "yahoo", from = start_date, auto.assign = FALSE)
  }, error = function(e) {
    message("FAILED: ", tkr)
    return(NULL)
  })
}

all_data <- list()
batches <- split(yahoo_tickers, ceiling(seq_along(yahoo_tickers)/batch_size))

for (i in seq_along(batches)) {
  batch <- batches[[i]]
  cat("Downloading batch", i, "of", length(batches), "\n")
  
  res <- lapply(batch, safe_get)
  names(res) <- batch
  res <- res[!sapply(res, is.null)]
  all_data <- c(all_data, res)
  
  Sys.sleep(1)  # avoid throttling Yahoo
}

if (length(all_data) == 0) stop("No tickers downloaded successfully!")

# 4. Track failed tickers
failed_tickers <- names(all_data)[!names(all_data) %in% names(res)]
if(length(failed_tickers) > 0){
  cat("Tickers failed to download:\n")
  print(failed_tickers)
  write.csv(failed_tickers, "failed_tickers_new.csv", row.names = FALSE)
}

# 5. Compute daily log returns using quantmod
returns_list <- lapply(all_data, function(x) {
  dailyReturn(Ad(x), type = "log")
})

# 6. Merge all returns into one xts object
returns <- Reduce(function(x, y) merge(x, y, all = TRUE), returns_list)

# 7. Rename columns to ticker names
colnames(returns) <- names(returns_list)

cat("Merged daily log returns: ", nrow(returns), "rows ×", ncol(returns), "columns\n")

# 8. Save returns only under new filenames
write.csv(as.data.frame(returns), "sweden_returns_v2.csv", row.names = TRUE)
saveRDS(returns, "sweden_returns_v2.rds")

cat("Done — returns saved under new files (v2).\n")
