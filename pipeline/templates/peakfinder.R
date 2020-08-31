

peakfinder <- function(data_umi) {
  dd <- density(data_umi)
  #define window size
  smallBins <-
    ifelse(round(length(dd$x) * 0.01, 0) %% 2,
           round(length(dd$x) * 0.01, 0),
           round(length(dd$x) * 0.01, 0) + 1)
  #add tails to begining and end of the density plot
  zero_beg <-
    seq(min(dd$x) - 0.01, min(dd$x), length.out = ((1 + smallBins) / 2))
  zero_end <-
    seq(max(dd$x) + 0.01, max(dd$x), length.out = ((1 + smallBins) / 2))
  dd$x <- c(zero_beg, dd$x, zero_end)
  dd$y <-
    c(rep(0, (1 + smallBins) / 2), dd$y, rep(0, (1 + smallBins) / 2))
  #define range for local maximums searching
  isRange <-
    ((1 + smallBins) / 2):(length(dd$y) - (smallBins - 1) / 2)
  #find local maximums inside defined window
  isLocalMax <-
    sapply(isRange, function(i)
      which.max(dd$y[(i - ((smallBins / 2) - 0.5)):(i + ((smallBins / 2) - 0.5))]) == (smallBins / 2) + 0.5)
  #filter local maximums with too small values
  peaks_x <- dd$x[which(isLocalMax)]
  peaks_y <- dd$y[which(isLocalMax)]
  peaks <- as.data.frame(cbind(peaks_x, peaks_y))
  peaks <- peaks %>% filter(peaks_y > max(dd$y) * 0.1)
  peaks <- peaks$peaks_x
  peaks
}
