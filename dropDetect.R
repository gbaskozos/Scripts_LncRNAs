###################################
# Detect sudden drops in coverage #
###################################
# Function adapted from Jean-Paul van Brakel, https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa #

dropDetect <- function(coverage, start, seqnames, strand, lag,threshold, length, influence, intron_identification) {
require(IRanges)
require(GenomicRanges)

if (length(coverage) > lag) {
  
  y <- as.vector(coverage)
  x <- rev(y)

  signals <- rep(0,length(y))
  signals_x <- rep(0,length(x))


  filteredY <- y[0:lag]
  avgFilterY <- NULL
  stdFilterY <- NULL
  avgFilterY[lag] <- median(y[0:lag])
  stdFilterY[lag] <- mad(y[0:lag])

  filteredX <- x[0:lag]
  avgFilterX <- NULL
  stdFilterX <- NULL
  avgFilterX[lag] <- median(x[0:lag])
  stdFilterX[lag] <- mad(x[0:lag])

  for (i in (lag+1):length(y)){
    if (abs(y[i]-avgFilterY[i-1]) > threshold*stdFilterY[i-1]) {
      if (y[i] > avgFilterY[i-1]) {
        signals[i] <- 1;
      } else {
        signals[i] <- -1;
      }
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0
      filteredY[i] <- y[i]
    }
    avgFilterY[i] <- median(filteredY[(i-lag):i])
    stdFilterY[i] <- mad(filteredY[(i-lag):i])
  }

  for (i in (lag+1):length(x)){
    if (abs(x[i]-avgFilterX[i-1]) > threshold*stdFilterX[i-1]) {
      if (x[i] > avgFilterX[i-1]) {
        signals_x[i] <- 1;
      } else {
        signals_x[i] <- -1;
      }
      filteredX[i] <- influence*x[i]+(1-influence)*filteredX[i-1]
    } else {
      signals_x[i] <- 0
      filteredX[i] <- x[i]
    }
    avgFilterX[i] <- median(filteredX[(i-lag):i])
    stdFilterX[i] <- mad(filteredX[(i-lag):i])
  }

comb <- Rle(signals + rev(signals_x))

start_drop <- Views(comb, start = NULL, end = NULL)
end_drop <- Views(comb, start = NULL, end = NULL)

if (runValue(comb)[1] <= -1 ) {
start_drop <- Views(comb, start = 1, width = runLength(comb)[1])
}

if (runValue(comb)[length(runValue(comb))] <= -1) {
end_drop <- Views(comb, start =  length(comb) - runLength(comb)[length(runValue(comb))]+1, width = runLength(comb)[length(runValue(comb))])
}

if (intron_identification) {
drops <- slice(comb, upper = -2)
drops <- drops[width(drops) > length]
drops <- reduce(c(drops, start_drop, end_drop))
} else {
drops <- reduce(c(start_drop, end_drop))
 }

if (length(drops) > 0) {
return( GRanges(seqnames, IRanges(start(drops)+start, width = width(drops)), strand=strand) )
} else {
GRanges(seqnames, 0, strand=strand)
} 

} else {
GRanges(seqnames, 0, strand=strand)
	}

}
