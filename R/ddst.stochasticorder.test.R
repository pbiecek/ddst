#' Data Driven Nonparametric Test for Stochastic Ordering
#'
#' Performs the data driven smooth test for detection of
#' the stochastic ordering, as described
#' in detail in Ledwina and Wyłupek (2012).
#' Suppose that we have random samples from two distributions F and G.
#' The null hypothesis is that F(x) >= G(x) for all x while the alternative is that at
#' F(x) < G(x) for some x.
#' Detailed description of the test statistic is provided in Ledwina and Wylupek (2012).
#'
#' @param x a (non-empty) numeric vector of data values
#' @param y a (non-empty) numeric vector of data values
#' @param t a positive number, penalty for model selection rule, see package description
#' @param d an integer, number of coordinates that measure potential deviation from null hypothesis
#' @param ... further arguments
#'
#' @references Nonparametric tests for stochastic ordering. Ledwina and Wyłupek (2012) \url{https://doi.org/10.1007/s11749-011-0278-7}
#' @export
#' @examples
#' set.seed(13)
#' library("rmutil", warn.conflicts = FALSE)
#' # 1. Pareto(1)/Pareto(1.5)
#' #uzyc parametrow z tabeli 3, p. 742, np dla m = n = 50,
#' x <- rpareto(50, 2, 2)
#' y <- rpareto(50, 1.5, 1.5)
#' t <- ddst.stochasticorder.test(x, y)
#' t
#' plot(t)
#'
#' # 2. Laplace(0,1)/Laplace(1,25)
#' x <- rlaplace(50, 0, 1)
#' y <- rlaplace(50, 1, 25)
#' t <- ddst.stochasticorder.test(x, y)
#' t
#' plot(t)
#'
#' # 3. LN(0.85,0.6)/LN(1.2,0.2)
#' x <- rlnorm(50, 0.85, 0.6)
#' y <- rlnorm(50, 1.2, 0.2)
#' t <- ddst.stochasticorder.test(x, y)
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
#'    ddst.stochasticorder.test(x, y)$statistic
#' })
#' quantile(samp, 0.95)
#' plot(ecdf(samp))
#' }
#'
#' @keywords htest
`ddst.stochasticorder.test` <-
  function(x,
           y,
           t = 2.2,
           d = 4) {
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
    statistics <- test.Q(x, y, m, n, d, t)
    names(statistics) <- c("stat.T","stat.S","ncoord.T","ncoord.S", "L")
    # coords <- test.M.d(x,y,m,n,k.N)

    l = statistics$stat.T
    attr(l, "names") = "Q.T"
    t = statistics$ncoord.T
    attr(t, "names") = "T"
    result = list(statistic = l,
                  parameter = t,
                  coordinates = statistics$L,
                  method = "Data Driven Stochastic Ordering Test")
    class(result) = c("htest", "ddst.test", "ddst.stochasticorder.test")

    result
}
