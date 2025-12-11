## Loading data

library(dplyr)

df <- read.csv("Data/sweden_returns_merged.csv")

colnames(df)[colnames(df) == "X2CUREX.ST"] <- "2CUREX.ST"
## Choosing time period: 2019 - 2021

df_2019_2021 <- df %>%
  mutate(X = as.Date(X)) %>%   # if this gives NA, see note below
  filter(X >= as.Date("2019-01-01"),
         X <= as.Date("2022-12-31"))


## Checking for missing values

na_counts <- sort(colSums(is.na(df_2019_2021)), decreasing = TRUE)
na_counts[na_counts > 0]

## Removing all assets that have number of NA > 1 for time period

na_counts <- colSums(is.na(df_2019_2021))
keep_cols <- names(na_counts)[na_counts <= 1]
keep_cols <- union("X", setdiff(keep_cols, "X"))  # ensure X stays first/kept
df_clean <- df_2019_2021 %>% select(all_of(keep_cols))

## Interpolating/Extrapolating the two missing values for the two assets with only one value missing

df_clean$INCOAX.ST[1] <- 0 # assign 0 as if stock did not move
df_clean$RO.ST[219] <- 0 # assign 0 as if stock did not move
df_clean <- df_clean[, setdiff(names(df_clean), "BOTX.ST")] # removing this stock since it was only zeros
df_clean <- df_clean[, setdiff(names(df_clean), "MTG.A.ST")] # removing this stock since it was only zeros
df_clean <- df_clean[, setdiff(names(df_clean), "KOGO.ST")] # removing this stock since it was only zeros
df_clean <- df_clean[, setdiff(names(df_clean), "G2M.ST")] # removing this stock since it was only zeros
df_clean <- df_clean[, setdiff(names(df_clean), "HAKI.A.ST")] # removing this stock since it was only zeros
df_clean <- df_clean[, setdiff(names(df_clean), "SBOK.ST")] # removing this stock since it was only zeros
df_clean <- df_clean[, setdiff(names(df_clean), "KAKEL.ST")] # removing this stock since it was only zeros


