#' Data Driven Smooth Test for Stochastic Dominance in Two Samples
#'
#' Performs the data driven smooth test for detection of
#' the stochastic ordering, as described
#' in detail in Ledwina and Wyłupek (2012).
#' Suppose that we have random samples from two distributions F and G.
#' The null hypothesis is that F(x) >= G(x) for all x while the alternative is that at
#' F(x) < G(x) for some x.
#' Detailed description of the test statistic is provided in Ledwina and Wylupek (2012).
#'
#' @param x a (non-empty) numeric vector of data
#' @param y a (non-empty) numeric vector of data
#' @param K.N an integer specifying a level of complexity of the grid considered, only for advanced users
#' @param alpha a significance level
#' @param t an alpha-dependent tunning parameter in the penalty in the model selection rule
#' @param nr an integer specifying the number of runs for a p-value and a critical value computation if any
#' @param compute.p a logical value indicating whether to compute a p-value or not
#' @param compute.cv a logical value indicating whether to compute a critical value corresponding to the significance level alpha or not
#'
#' @references Nonparametric tests for stochastic ordering. Ledwina and Wyłupek (2012) \url{https://doi.org/10.1007/s11749-011-0278-7}
#' @export
#' @examples
#' set.seed(7)
#' library("rmutil", warn.conflicts = FALSE)
#' # 1. Pareto(1)/Pareto(1.5)
#' # H0 is false
#' x <- rpareto(50, 2, 2)
#' y <- rpareto(50, 1.5, 1.5)
#' t <- ddst.forstochdom.test(x, y, t = 2.2, K.N = 4)
#' t
#' plot(t)
#'
#' # 2. Laplace(0,1)/Laplace(1,25)
#' # H0 is false
#' x <- rlaplace(50, 0, 1)
#' y <- rlaplace(50, 1, 25)
#' t <- ddst.forstochdom.test(x, y, t = 2.2, K.N = 4)
#' t
#' plot(t)
#'
#' # 3. LN(0.85,0.6)/LN(1.2,0.2)
#' # H0 is true
#' x <- rlnorm(50, 0.85, 0.6)
#' y <- rlnorm(50, 1.2, 0.2)
#' t <- ddst.forstochdom.test(x, y, t = 2.2, K.N = 4)
#' t
#' plot(t)
#'
#' \dontrun{
#' # Generate distribution of test statistic
#' N <- 1000
#' samp <- replicate(N, {
#'    x <- runif(30)
#'    y <- runif(30)
#'    # statistics with Schwartz penalty
#'    ddst.forstochdom.test(x, y)$statistic
#' })
#' quantile(samp, 0.95)
#' plot(ecdf(samp))
#' }
#'
#' @keywords htest
`ddst.forstochdom.test` <-
  function(x,
           y,
           K.N = floor(log(length(x)+length(y),2))-1, # d
           alpha = 0.05,
           t, #  = 2.2
           nr = 100000,
           compute.p = TRUE,
           compute.cv = TRUE) {
    m = length(x)
    n = length(y)

### INTERNAL FUNCTIONS
    l.j = function(k, i, z) {
      a.j = (2*i-1)/2^(k+1)
      left = -sqrt((1 - a.j)/a.j)
      right = sqrt(a.j/(1 - a.j))
      score = ifelse(z < a.j, left, right)
      return(score)
    }

    L.j = function(k, i, x, y, m, n) {
      N = m + n
      H = ecdf(c(x,y))
      rx = H(x) - 1/(2*N)
      ry = H(y) - 1/(2*N)
      cx = sum( l.j(k,i,rx) ) / m
      cy = sum( l.j(k,i,ry) ) / n
      score = (cy - cx)*sqrt(m*n/N)
      return(score)
    }

    test.M.d = function(x,y,m,n,k.N) {
      N = m+n
      d.N = 2^(k.N+1)-1
      L.vec = matrix(0,1,d.N)
      for(k in 0:k.N) {
        for(i in 1:(2^k)){
          j = 2^k - 1 + i
          L.vec[1,j] = L.j(k,i,x,y,m,n)
        }
      }
      M.d = min( L.vec[1,] )
      result = c(M.d,L.vec)
      return(result)
    }

    test.Q = function(x,y,m,n,k.N,t){
      N = m+n
      d.N = 2^(k.N+1)-1
      L.vec = matrix(0,1,d.N)
      L.vec.trun = matrix(0,1,d.N)
      for(k in 0:k.N) {
        for(i in 1:(2^k)) {
          j = 2^k - 1 + i
          L.vec[1,j] = L.j(k,i,x,y,m,n)
          L.vec.trun[1,j] = max( - L.vec[1,j], 0 )
        }
      }
      Q.d.vec = matrix(0,1,k.N+1)
      D.vec = 2^(0:k.N+1) - 1
      for(k in 1:(k.N+1)){
        Q.d.vec[1,k] = L.vec.trun[1,1:D.vec[k]] %*% L.vec.trun[1,1:D.vec[k]]
      }
      S = which.max( Q.d.vec[1,] - D.vec*log(N) )
      M = which.max( Q.d.vec[1,] )
      if( max( -L.vec[1,] ) <= sqrt(t*log(N)) ){
        T = S
      }else{
        T = M
      }
      Q.S = Q.d.vec[1,S]
      Q.T = Q.d.vec[1,T]
      result = list(Q.T, Q.S, T, S, L = L.vec.trun[1,])
      return(result)
    }

## END OF INTERNAL FUNCTIONS

    # T, S
    statistics <- test.Q(x, y, m, n, K.N, t)
    names(statistics) <- c("stat.T","stat.S","ncoord.T","ncoord.S", "L")
    # coords <- test.M.d(x,y,m,n,k.N)

    l = statistics$stat.T
    attr(l, "names") = "QT"
    qt = statistics$ncoord.T
    attr(qt, "names") = "T"
    result = list(statistic = l,
                  parameter = qt,
                  coordinates = statistics$L,
                  method = "Data Driven Stochastic Ordering Test")
    result$data.name = paste(paste(as.character(substitute(x)),
                             as.character(substitute(y)), collapse = " "),
                             ", t: ",
                             t,
                             ", K.N: ",
                             K.N,
                             sep = "")
    class(result) = c("htest", "ddst.test", "ddst.stochasticorder.test")

    result
}
