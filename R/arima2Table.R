#'
#'Construct table of AIC for all combinations 0<=p<=P and 0<=q<=Q
#'@param data a time series object
#'@param P a positive integer value
#'@param Q a positive integer value
aicTable <- function(data,P,Q, D=0){
  table <- matrix(NA,(P+1),(Q+1))
  for(p in 0:P) {
    for(q in 0:Q) {
      table[p+1,q+1] <- arima2(data,order=c(p,D,q))$aic
    }
  }
  dimnames(table) <- list(paste("AR",0:P, sep=""),paste("MA",0:Q,sep=""))
  table
}
